// src/classify_sequence.rs

use std::collections::BTreeSet;
use ahash::{AHashMap, AHashSet};
use std::fmt::Write;
use std::sync::Arc;
use crate::krakendb::KrakenDB;
use crate::types::{DNASequence, KrakenOutputLine};
use rayon::prelude::*; // <-- Add Rayon

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

/// If multiple taxa tie for best "score," break ties by LCA.
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
            // We have a tie => track both in max_taxa
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

/// A type alias for "taxon => read+kmers info."
pub type TaxonCounts = AHashMap<u32, ReadCounts>;

/// Build a string like "2:10 0:5 ..." from the per-base assigned taxa + ambig flags.
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
    // Flush the final run
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
            _ => return (0, false), // "N" or invalid base => invalid
        }
    }
    (val, true)
}

/// Parallelized version of classify_sequence, using Rayon at the k-mer level.
/// If you'd rather parallelize over reads, move the parallel logic outside this function.
pub fn classify_sequence(
    dna: &DNASequence,
    kraken_dbs: &[Arc<KrakenDB>],
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

    // ----------------------------------------------------------------
    // 1) Parallel: process each k-mer
    //    We'll store partial results in (assigned_taxon, is_ambig, local_hits).
    // ----------------------------------------------------------------
    let (assigned_info, local_maps): (Vec<_>, Vec<_>) =
        (0..span_count)
            .into_par_iter()
            .map(|start| {
                let (encoded, is_valid) = encode_kmer_2bit_slice(&seq_bytes[start..start + k]);
                if !is_valid {
                    // Return "unassigned" plus no hits
                    return (
                        (start, 0u32, true), // assigned_taxon=0, is_ambig=true
                        AHashMap::new(),     // local AHashMap<taxid, Vec<kmer>>
                    );
                }

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

                if found_taxon == 0 {
                    // No match => assigned_taxon=0
                    return ((start, 0u32, false), AHashMap::new());
                }

                // Found a taxon => store kmer in a local map
                let mut local_hits: AHashMap<u32, Vec<u64>> = AHashMap::default();
                local_hits.insert(found_taxon, vec![canon_kmer]);

                ((start, found_taxon, false), local_hits)
            })
            // fold partial results
            .fold(
                || (Vec::new(), Vec::new()),
                |mut acc, (info, local_map)| {
                    acc.0.push(info);
                    acc.1.push(local_map);
                    acc
                },
            )
            // reduce across threads
            .reduce(
                || (Vec::new(), Vec::new()),
                |(mut assigned_a, mut local_a), (mut assigned_b, mut local_b)| {
                    assigned_a.append(&mut assigned_b);
                    local_a.append(&mut local_b);
                    (assigned_a, local_a)
                },
            );

    // ----------------------------------------------------------------
    // 2) Reconstruct assigned_taxa + ambig_list in correct order
    // ----------------------------------------------------------------
    let mut assigned_taxa = vec![0u32; span_count];
    let mut ambig_list = vec![false; span_count];
    for (idx, taxon, is_ambig) in assigned_info {
        assigned_taxa[idx] = taxon;
        ambig_list[idx] = is_ambig;
    }

    // ----------------------------------------------------------------
    // 3) Merge local_maps into a single hit_counts and also into global taxon_counts
    // ----------------------------------------------------------------
    let mut hit_counts: AHashMap<u32, u32> = AHashMap::new();
    for lm in local_maps {
        for (taxid, kmer_list) in lm {
            // Bump the raw count
            *hit_counts.entry(taxid).or_insert(0) += 1;

            // Also store k-mers in the shared taxon_counts
            let entry = taxon_counts.entry(taxid).or_default();
            for k in kmer_list {
                entry.add_kmer_occurrence(k);
            }
        }
    }

    // ----------------------------------------------------------------
    // 4) Classify using the hit_counts => final assigned taxid
    // ----------------------------------------------------------------
    let call = resolve_tree_kraken(&hit_counts, parent_map);

    // Update taxon counts for the final call
    if call != 0 {
        taxon_counts.entry(call).or_default().increment_read_count();
    }

    // Skip unclassified output if requested
    if only_classified_kraken_output && call == 0 {
        return false;
    }

    // ----------------------------------------------------------------
    // 5) Build the final "hitlist" string for output
    // ----------------------------------------------------------------
    let hitlist = if assigned_taxa.is_empty() {
        "0:0".to_string()
    } else {
        hitlist_string(&assigned_taxa, &ambig_list)
    };

    // ----------------------------------------------------------------
    // 6) Append to kraken_output_lines
    // ----------------------------------------------------------------
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
