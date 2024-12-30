// src/classify_sequence.rs

use std::collections::BTreeSet;
use ahash::{AHashMap, AHashSet};
use std::fmt::Write;
use crate::krakendb::KrakenDB;
use crate::types::KrakenOutputLine;

/// A parent map: taxon -> parent taxon (e.g. parent_map[9606] = 9605).
pub type ParentMap = AHashMap<u32, u32>;

/// Return the lowest common ancestor (LCA) of `a` and `b`.
/// If either is `0`, we treat it as “unclassified,” so LCA(0,x)=x, etc.
/// If we can’t find another ancestor, return `1` (root).
pub fn lca(parent_map: &ParentMap, mut a: u32, mut b: u32) -> u32 {
    if a == 0 || b == 0 {
        return if a == 0 { b } else { a };
    }

    // Collect all ancestors of `a`
    let mut a_anc = AHashSet::new();
    while a > 1 {
        a_anc.insert(a);
        if let Some(&p) = parent_map.get(&a) {
            a = p;
        } else {
            break;
        }
    }

    // Climb `b` upward until we find a common ancestor
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
    1
}

/// If multiple taxa tie for best “score,” break ties by LCA.
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

/// Sum hits from each taxon up its ancestry chain (Kraken’s approach).
/// Then pick the taxon with the highest total. If ties => LCA.
pub fn resolve_tree_kraken(
    hit_counts: &AHashMap<u32, u32>,
    parent_map: &ParentMap
) -> u32 {
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
                if p == node {
                    break;
                }
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

/// A minimal representation of a read.
#[derive(Debug, Clone)]
pub struct DNASequence {
    pub id: String,
    pub header_line: String,
    pub seq: String,
    pub quals: String,
}

/// We now store each k-mer’s coverage in `kmer_counts`.
#[derive(Default, Debug, Clone)]
pub struct ReadCounts {
    pub read_count: u32,

    /// Map: k-mer => how many times it appears
    pub kmer_counts: AHashMap<u64, u32>,
}

impl ReadCounts {
    pub fn increment_read_count(&mut self) {
        self.read_count += 1;
    }

    /// Add one occurrence for this `kmer`.
    pub fn add_kmer_occurrence(&mut self, kmer: u64) {
        *self.kmer_counts.entry(kmer).or_insert(0) += 1;
    }

    /// How many *distinct* k-mers were found in this taxon?
    pub fn distinct_kmers(&self) -> usize {
        self.kmer_counts.len()
    }

    /// How many total k-mer occurrences (coverage) in this taxon?
    pub fn total_coverage(&self) -> u64 {
        self.kmer_counts.values().map(|&cnt| cnt as u64).sum()
    }
}

/// A type alias for “taxon => read+kmers info.”
pub type TaxonCounts = AHashMap<u32, ReadCounts>;

/// Build a string like “2:10 0:5 ...” from the per-base assigned taxa + ambig flags.
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
    // Flush final run
    if last_code >= 0 {
        let _ = write!(out, "{}:{}", last_code, code_count);
    } else {
        let _ = write!(out, "A:{}", code_count);
    }
    out
}

/// Convert a slice of A/C/G/T (ASCII) to a 2-bit encoded k-mer.
/// Return `(encoded_kmer, is_valid)`.
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
            _ => return (0, false), // “N” or invalid base => invalid
        }
    }
    (val, true)
}

/// Modified to return classification status and avoid string building
pub fn classify_sequence(
    dna: &DNASequence,
    kraken_dbs: &[KrakenDB],
    parent_map: &ParentMap,
    taxon_counts: &mut TaxonCounts,
    print_sequence_in_kraken: bool,
    only_classified_kraken_output: bool,
    kraken_output_lines: &mut Vec<KrakenOutputLine>,
) -> bool {
    // If no DB provided, nothing to do
    if kraken_dbs.is_empty() {
        return false;
    }

    let k = kraken_dbs[0].k as usize;
    let seq_len = dna.seq.len();

    // If read is too short => unclassified
    if seq_len < k {
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

    let span_count = seq_len.saturating_sub(k) + 1;
    let seq_bytes = dna.seq.as_bytes();
    let mut assigned_taxa = Vec::with_capacity(span_count);
    let mut ambig_list = Vec::with_capacity(span_count);
    let mut hit_counts: AHashMap<u32, u32> = AHashMap::new();

    // Process each k-mer
    for start in 0..span_count {
        let (encoded, is_valid) = encode_kmer_2bit_slice(&seq_bytes[start..start + k]);
        if !is_valid {
            assigned_taxa.push(0);
            ambig_list.push(true);
            continue;
        }
        ambig_list.push(false);

        // Check each DB
        let mut found_taxon = 0;
        let mut canon_kmer = 0;
        for db in kraken_dbs {
            canon_kmer = db.canonical_representation(encoded);
            if let Some(tid) = db.kmer_query(canon_kmer) {
                found_taxon = tid;
                break;
            }
        }

        assigned_taxa.push(found_taxon);

        if found_taxon != 0 {
            *hit_counts.entry(found_taxon).or_insert(0) += 1;
            let entry = taxon_counts.entry(found_taxon).or_default();
            entry.add_kmer_occurrence(canon_kmer);
        }
    }

    // Classify using the hit counts
    let call = resolve_tree_kraken(&hit_counts, parent_map);

    // Update taxon counts for the final call
    if call != 0 {
        taxon_counts.entry(call).or_default().increment_read_count();
    }

    // Skip unclassified output if requested
    if only_classified_kraken_output && call == 0 {
        return false;
    }

    // Create output line
    let hitlist = if assigned_taxa.is_empty() {
        "0:0".to_string()
    } else {
        hitlist_string(&assigned_taxa, &ambig_list)
    };

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
