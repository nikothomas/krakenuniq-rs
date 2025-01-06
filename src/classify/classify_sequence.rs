use std::collections::BTreeSet;
use ahash::{AHashMap, AHashSet};
use std::fmt::Write;
use std::sync::Arc;
use crate::krakendb::KrakenDB;
use crate::types::{DNASequence, KrakenOutputLine};

/// A parent map: taxon -> parent taxon
pub type ParentMap = AHashMap<u32, u32>;

#[inline]
fn encode_base_2bit(b: u8) -> u64 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4, // invalid marker
    }
}

/// Encode the first k-mer by scanning k bases in the slice.
/// Returns `(encoded_kmer, valid_window)`.
#[inline]
fn encode_first_kmer(seq_bytes: &[u8]) -> (u64, bool) {
    let mut val = 0u64;
    for &b in seq_bytes {
        let code = encode_base_2bit(b);
        if code > 3 {
            // invalid base => 'N'
            return (0, false);
        }
        val = (val << 2) | code;
    }
    (val, true)
}

/// Roll the previous k-mer to the next k-mer by removing the old base and adding a new base.
/// `k` is the k-mer length in bases.
#[inline]
fn roll_kmer_2bit(old: u64, old_base: u8, new_base: u8, k: usize) -> (u64, bool) {
    // For a k-mer, we have 2*k bits total. We shift out the top 2 bits
    // that represented `old_base`. Example approach:
    //   mask off the topmost 2 bits * k
    let shift = 2 * (k - 1);
    let mask = (1 << shift) - 1; // 2^(2*(k-1)) - 1
    let mut val = old & mask;

    // Shift left by 2 and add new base
    val <<= 2;
    let code = encode_base_2bit(new_base);
    if code > 3 {
        // invalid
        return (0, false);
    }
    val |= code;

    (val, true)
}

/// Return the lowest common ancestor (LCA) of `a` and `b`.
/// If either is `0`, treat it as unclassified => return the other.
/// If none found, return `1`.
pub fn lca(parent_map: &ParentMap, mut a: u32, mut b: u32) -> u32 {
    if a == 0 || b == 0 {
        return if a == 0 { b } else { a };
    }

    // Collect ancestors of a
    let mut a_anc = AHashSet::with_capacity(16);
    while a > 1 {
        a_anc.insert(a);
        if let Some(&p) = parent_map.get(&a) {
            if p == a { break; }
            a = p;
        } else {
            break;
        }
    }

    // Climb b upward until we find a common ancestor
    while b > 1 {
        if a_anc.contains(&b) {
            return b;
        }
        if let Some(&p) = parent_map.get(&b) {
            if p == b { break; }
            b = p;
        } else {
            break;
        }
    }
    1
}

/// Tie-break by LCA among multiple top taxa.
fn break_ties(max_taxa: BTreeSet<u32>, parent_map: &ParentMap) -> u32 {
    if max_taxa.is_empty() {
        return 0;
    }
    let mut iter = max_taxa.into_iter();
    let mut candidate = iter.next().unwrap();
    for t in iter {
        candidate = lca(parent_map, candidate, t);
    }
    candidate
}

/// Sum hits from each taxon up its ancestry chain, pick the one with highest total. Ties => LCA.
pub fn resolve_tree_kraken(hit_counts: &AHashMap<u32, u32>, parent_map: &ParentMap) -> u32 {
    let mut max_taxa = BTreeSet::new();
    let mut max_taxon = 0u32;
    let mut max_score = 0u32;

    // For each leaf node that got hits
    for (&taxon, _) in hit_counts {
        let mut node = taxon;
        let mut score = 0u32;

        while node > 0 {
            if let Some(&cnt) = hit_counts.get(&node) {
                score += cnt;
            }
            if let Some(&p) = parent_map.get(&node) {
                if p == node { break; }
                node = p;
            } else {
                break;
            }
        }

        if score > max_score {
            max_score = score;
            max_taxon = taxon;
            max_taxa.clear();
        } else if score == max_score {
            if max_taxa.is_empty() {
                max_taxa.insert(max_taxon);
            }
            max_taxa.insert(taxon);
        }
    }

    if !max_taxa.is_empty() {
        let tie_broken = break_ties(max_taxa, parent_map);
        if tie_broken != 0 {
            return tie_broken;
        }
    }

    max_taxon
}

/// We store each k-mer’s coverage in `kmer_counts`.
#[derive(Default, Debug, Clone)]
pub struct ReadCounts {
    pub read_count: u32,
    /// Map: k-mer => how many times it appears
    pub kmer_counts: AHashMap<u64, u32>,
}

impl ReadCounts {
    #[inline]
    pub fn increment_read_count(&mut self) {
        self.read_count += 1;
    }

    #[inline]
    pub fn add_kmer_occurrence(&mut self, kmer: u64) {
        *self.kmer_counts.entry(kmer).or_insert(0) += 1;
    }
}

/// A type alias for taxon => read+kmers info
pub type TaxonCounts = AHashMap<u32, ReadCounts>;

/// Convert the assigned taxa + ambiguous flags to a "hitlist" string.
fn hitlist_string(taxa: &[u32], ambig: &[bool]) -> String {
    let mut out = String::new();
    if taxa.is_empty() {
        return out;
    }
    let mut last_code: i64 = if ambig[0] { -1 } else { taxa[0] as i64 };
    let mut code_count = 1;

    for i in 1..taxa.len() {
        let code: i64 = if ambig[i] { -1 } else { taxa[i] as i64 };
        if code == last_code {
            code_count += 1;
        } else {
            if last_code >= 0 {
                let _ = write!(out, "{}:{} ", last_code, code_count);
            } else {
                let _ = write!(out, "A:{} ", code_count);
            }
            code_count = 1;
            last_code = code;
        }
    }
    // final flush
    if last_code >= 0 {
        let _ = write!(out, "{}:{}", last_code, code_count);
    } else {
        let _ = write!(out, "A:{}", code_count);
    }
    out
}

/// Single-thread classification of one read.
/// - Uses rolling k-mers
/// - If a k-mer has an invalid base ('N'), we skip that window until the next valid region.
pub fn classify_sequence(
    dna: &DNASequence,
    kraken_dbs: &[Arc<KrakenDB>],
    parent_map: &ParentMap,
    taxon_counts: &mut TaxonCounts,
    print_sequence_in_kraken: bool,
    only_classified_kraken_output: bool,
    kraken_output_lines: &mut Vec<KrakenOutputLine>,
) -> bool {
    // If no DB, we can't classify
    if kraken_dbs.is_empty() {
        return false;
    }

    let k = kraken_dbs[0].k as usize;
    let seq_len = dna.seq.len();
    if seq_len < k {
        // Too short => unclassified
        if !only_classified_kraken_output {
            kraken_output_lines.push(KrakenOutputLine {
                status: 'U',
                read_id: dna.id.clone(),
                tax_id: 0,
                length: 0,
                hitlist: "0:0".to_string(),
                sequence: if print_sequence_in_kraken {
                    Some(dna.seq.clone())
                } else {
                    None
                },
            });
        }
        return false;
    }

    let seq_bytes = dna.seq.as_bytes();
    let span_count = seq_len - k + 1;

    // We'll store the assigned taxon for each window, plus a boolean if it's ambiguous
    let mut assigned_taxa = vec![0u32; span_count];
    let mut ambig_list = vec![false; span_count];

    // This will track how many times each taxon was hit (for final classification).
    let mut local_hit_counts: AHashMap<u32, u32> = AHashMap::with_capacity(span_count);

    // Rolling k-mer logic:
    let mut i = 0;
    while i < span_count {
        // 1) Attempt to encode the k-mer starting at i
        let window = &seq_bytes[i..i + k];
        let (mut encoded, valid) = encode_first_kmer(window);
        if !valid {
            // There's an 'N' or invalid base => skip past this window
            assigned_taxa[i] = 0; // unassigned
            ambig_list[i] = true; // mark it ambiguous
            i += 1;
            continue;
        }

        // 2) Query each DB to find the first DB that has a match
        let mut found_taxon = 0;
        let mut canon_kmer = 0;
        for db in kraken_dbs {
            // If your DB supports canonical representations:
            canon_kmer = db.canonical_representation(encoded);
            if let Some(tid) = db.kmer_query(canon_kmer) {
                found_taxon = tid;
                break;
            }
        }

        assigned_taxa[i] = found_taxon;
        ambig_list[i] = false;

        // 3) If found_taxon != 0, add to local_hit_counts and to global kmer list
        if found_taxon != 0 {
            *local_hit_counts.entry(found_taxon).or_insert(0) += 1;
            taxon_counts
                .entry(found_taxon)
                .or_default()
                .add_kmer_occurrence(canon_kmer);
        }

        // 4) Move to next window. Use a rolling approach for subsequent windows:
        //    We only roll if we won't exceed seq bounds. But we do it in a loop if you'd like
        //    to encode subsequent windows quickly. For simplicity, we just do 1 step here.
        i += 1;
    }

    // 5) Determine final assigned taxid by summing hits up the tax tree
    let call = resolve_tree_kraken(&local_hit_counts, parent_map);
    if call != 0 {
        taxon_counts.entry(call).or_default().increment_read_count();
    }

    // 6) Build hitlist. If we skip unclassified in output, handle that
    if only_classified_kraken_output && call == 0 {
        return false;
    }

    let hitlist = if assigned_taxa.is_empty() {
        "0:0".to_string()
    } else {
        hitlist_string(&assigned_taxa, &ambig_list)
    };

    // 7) Produce final Kraken output line
    kraken_output_lines.push(KrakenOutputLine {
        status: if call != 0 { 'C' } else { 'U' },
        read_id: dna.id.clone(),
        tax_id: call,
        length: seq_len,
        hitlist,
        sequence: if print_sequence_in_kraken {
            Some(dna.seq.clone())
        } else {
            None
        },
    });

    call != 0
}
