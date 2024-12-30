//src/classify.rs

use std::collections::{HashMap, HashSet, BTreeSet};
use std::fmt::Write;
use crate::krakendb::KrakenDB;
use crate::types::{KrakenOutputLine, KrakenReportRow};
// -----------------------
// These come from krakenutil.cpp's logic
// -----------------------

/// A parent map: taxon -> parent taxon. For example, parent_map[9606] = 9605, etc.
pub type ParentMap = HashMap<u32, u32>;

/// Return the lowest common ancestor of 'a' and 'b'.
/// LCA(0, x) = x, LCA(x, 0) = x. If we can't find anything else, return 1 (root).
pub fn lca(parent_map: &ParentMap, mut a: u32, mut b: u32) -> u32 {
    // If either is 0, just return the other
    if a == 0 || b == 0 {
        return if a == 0 { b } else { a };
    }

    // Collect all ancestors of a (including a)
    let mut a_anc = HashSet::new();
    while a > 1 {
        a_anc.insert(a);
        // move 'a' to its parent
        if let Some(&p) = parent_map.get(&a) {
            a = p;
        } else {
            // if missing parent => break
            break;
        }
    }

    // Climb b's ancestry until we find something in a's set
    while b > 1 {
        if a_anc.contains(&b) {
            return b;
        }
        if let Some(&p) = parent_map.get(&b) {
            b = p;
        } else {
            break;
        }
    }
    1 // fallback to root if no common ancestor found
}

/// If there is a tie (multiple taxa in `max_taxa`), compute the LCA among them.
/// Returns `0` if `max_taxa` is empty (no tie).
fn break_ties(max_taxa: BTreeSet<u32>, parent_map: &ParentMap) -> u32 {
    if max_taxa.is_empty() {
        // No tie to break
        return 0;
    }

    let mut iter = max_taxa.into_iter();
    let mut candidate = iter.next().unwrap();

    for t in iter {
        let old_candidate = candidate;
        candidate = lca(parent_map, candidate, t);
    }

    candidate
}

/// This replicates `resolve_tree(...)` from Kraken's krakenutil.cpp:
///   - For each taxon in `hit_counts`,
///   - sum all hits from that taxon up its ancestry chain.
///   - Keep track of the max score and the taxon that gave that max.
///   - If multiple taxa tie, compute LCA among them.
pub fn resolve_tree_kraken(
    hit_counts: &HashMap<u32, u32>,
    parent_map: &ParentMap
) -> u32 {
    let mut max_taxa = BTreeSet::new();
    let mut max_taxon = 0;
    let mut max_score = 0u32;

    // For each taxon in the map
    for (&taxon, &count) in hit_counts {
        // We'll climb from `taxon` to the root, summing any hits in `hit_counts`.
        let mut node = taxon;
        let mut score = 0u32;
        while node > 0 {
            if let Some(&cnt) = hit_counts.get(&node) {
                score += cnt;
            }
            if let Some(&p) = parent_map.get(&node) {
                if p == node {
                    break;
                }
                node = p;
            } else {
                break;
            }
        }

        // Compare to the max we've found so far
        if score > max_score {
            max_score = score;
            max_taxon = taxon;
            max_taxa.clear();
        } else if score == max_score {
            // It's a tie. We'll add both to a set so we can LCA them below
            if max_taxa.is_empty() {
                max_taxa.insert(max_taxon);
            }
            max_taxa.insert(taxon);
        }
    }

    // If we have multiple taxa tied for best path sum, we LCA them
    if !max_taxa.is_empty() {
        let tie_broken_taxon = break_ties(max_taxa, parent_map);
        if tie_broken_taxon != 0 {
            return tie_broken_taxon;
        }
    }

    max_taxon
}

// -----------------------
// Structures for your classification example
// -----------------------

#[derive(Debug, Clone)]
pub struct DNASequence {
    pub id: String,
    pub header_line: String,
    pub seq: String,
    pub quals: String,
}

/// A structure to keep track of read counts & unique k-mers for each taxon.
#[derive(Default, Debug)]
pub struct ReadCounts {
    pub read_count: u64,
    pub unique_kmers: HashSet<u64>,
}

impl ReadCounts {
    pub fn increment_read_count(&mut self) {
        self.read_count += 1;
    }
    pub fn add_kmer(&mut self, kmer: u64) {
        self.unique_kmers.insert(kmer);
    }
}

pub type TaxonCounts = HashMap<u32, ReadCounts>;

/// Convert a list of per-base assigned taxa + an "ambiguous" flag
/// into something like "2:10 0:5 ..." or "A:3" for ambiguous spans.
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
            // flush the previous run
            if last_code >= 0 {
                let _ = write!(out, "{}:{} ", last_code, code_count);
            } else {
                let _ = write!(out, "A:{} ", code_count);
            }
            code_count = 1;
            last_code = code;
        }
    }
    // flush final run
    if last_code >= 0 {
        let _ = write!(out, "{}:{}", last_code, code_count);
    } else {
        let _ = write!(out, "A:{}", code_count);
    }
    out
}

/// Encode a k-mer from `seq_bytes[start..start + k]` into a 2-bit representation.
/// Return `(encoded, is_valid)`.
#[inline]
fn encode_kmer_2bit_slice(seq_bytes: &[u8]) -> (u64, bool) {
    let mut val = 0u64;
    for &b in seq_bytes {
        val <<= 2;
        match b {
            b'A' | b'a' => { /* 0 */ }
            b'C' | b'c' => { val |= 1; }
            b'G' | b'g' => { val |= 2; }
            b'T' | b't' => { val |= 3; }
            _ => return (0, false), // 'N' or invalid base => not valid
        }
    }
    (val, true)
}

/// Classify a read using "Kraken-style" ancestor accumulation + `resolve_tree_kraken`.
///
/// 1. For each valid k-mer, find `found_taxon`.
/// 2. For stats, store the k-mer under `found_taxon`.
/// 3. **Add 1** to `hit_counts[current_taxon]` for *all* ancestors (like Kraken does).
/// 4. At the end, call `resolve_tree_kraken` to pick the final classification.
pub fn classify_sequence(
    dna: &DNASequence,
    kraken_dbs: &[KrakenDB],
    parent_map: &ParentMap,
    taxon_counts: &mut TaxonCounts,
    kraken_output: &mut String,
    classified_output: &mut String,
    unclassified_output: &mut String,
    print_sequence_in_kraken: bool,
    only_classified_kraken_output: bool,
    kraken_output_structs: &mut Vec<KrakenOutputLine>,
) -> bool {
    if kraken_dbs.is_empty() {
        return false;
    }

    // We'll assume all DBs have the same k-mer size
    let k = kraken_dbs[0].k as usize;
    let seq_len = dna.seq.len();

    // If read is too short, unclassified
    if seq_len < k {
        let seq_header = if dna.seq.as_bytes().first() == Some(&b'>') {
            ">some_fasta\n"
        } else {
            "@some_fastq\n"
        };
        let mut seq_out = String::new();
        seq_out.push_str(seq_header);
        seq_out.push_str(&dna.seq);
        seq_out.push('\n');
        unclassified_output.push_str(&seq_out);

        if !only_classified_kraken_output {
            let _ = write!(
                kraken_output,
                "U\t{}\t0\t0:0{}\n",
                dna.id,
                if print_sequence_in_kraken {
                    format!("\t{}", dna.seq)
                } else {
                    "".to_string()
                }
            );
        }
        return false;
    }

    let span_count = seq_len - k + 1;
    let seq_bytes = dna.seq.as_bytes();

    // We'll store the per-base assigned taxon (just for printing).
    let mut assigned_taxa = Vec::with_capacity(span_count);
    let mut ambig_list = Vec::with_capacity(span_count);

    // This is the table of "how many k-mer hits belong to each taxon"
    // including ancestors
    let mut hit_counts: HashMap<u32, u32> = HashMap::new();

    for start in 0..span_count {
        let (encoded, is_valid) = encode_kmer_2bit_slice(&seq_bytes[start..start + k]);
        if !is_valid {
            // ambiguous => 'N'
            assigned_taxa.push(0);
            ambig_list.push(true);
            continue;
        }
        ambig_list.push(false);

        // For each DB, see if there's a match
        let mut found_taxon = 0u32;
        for db in kraken_dbs {
            let canon_kmer = db.canonical_representation(encoded);
            if let Some(tid) = db.kmer_query(canon_kmer) {
                found_taxon = tid;
                break;
            }
        }
        assigned_taxa.push(found_taxon);

        // For stats, record the unique k-mer in our TaxonCounts if found_taxon != 0
        if found_taxon != 0 {
            taxon_counts
                .entry(found_taxon)
                .or_default()
                .add_kmer(encoded);

            *hit_counts.entry(found_taxon).or_insert(0) += 1;
        }
    }

    // Now do the final classification: sum each leaf's hits up its chain,
    // pick whichever is best, if tie => LCA among them.
    let call = resolve_tree_kraken(&hit_counts, parent_map);

    // Bump the read count for whichever final call we got
    taxon_counts.entry(call).or_default().increment_read_count();

    // Build read-level output
    let seq_header = if dna.seq.as_bytes().first() == Some(&b'>') {
        ">some_fasta\n"
    } else {
        "@some_fastq\n"
    };
    let mut seq_out = String::new();
    seq_out.push_str(seq_header);
    seq_out.push_str(&dna.seq);
    seq_out.push('\n');

    if call == 0 {
        unclassified_output.push_str(&seq_out);
    } else {
        classified_output.push_str(&seq_out);
    }

    if only_classified_kraken_output && call == 0 {
        return false;
    }

    // Build a Kraken-like output line: C/U, read_id, call, length, "hitlist"
    if call != 0 {
        kraken_output.push_str("C\t");
    } else {
        kraken_output.push_str("U\t");
    }
    let _ = write!(kraken_output, "{}\t{}\t{}", dna.id, call, dna.seq.len());

    kraken_output.push('\t');
    if assigned_taxa.is_empty() {
        kraken_output.push_str("0:0");
    } else {
        let hl = hitlist_string(&assigned_taxa, &ambig_list);
        kraken_output.push_str(&hl);
    }

    if print_sequence_in_kraken {
        kraken_output.push('\t');
        kraken_output.push_str(&dna.seq);
    }
    kraken_output.push('\n');

    // Build a Kraken-like output line: C/U, read_id, call, length, "hitlist"
    let status_char = if call != 0 { 'C' } else { 'U' };
    let read_id = dna.id.clone();
    let length = dna.seq.len();
    let hitlist = if assigned_taxa.is_empty() {
        "0:0".to_string()
    } else {
        hitlist_string(&assigned_taxa, &ambig_list)
    };

    // Optionally store the sequence
    let seq_opt = if print_sequence_in_kraken {
        Some(dna.seq.clone())
    } else {
        None
    };

    // 2) Also build a KrakenOutputLine struct
    let output_line = KrakenOutputLine {
        status: status_char,
        read_id,
        tax_id: call,
        length,
        hitlist,
        sequence: seq_opt,
    };
    kraken_output_structs.push(output_line);


    call != 0
}
