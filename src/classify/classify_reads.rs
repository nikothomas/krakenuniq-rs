// src/classify/classify_reads.rs

use rayon::prelude::*;
use std::sync::Arc;

use super::{classify_sequence, TaxonCounts};
use super::KrakenDB;
use super::ParentMap;
use crate::types::{DNASequence, KrakenOutputLine};

/// Parallel classification of multiple reads at once.
///
/// - `db` and `parent_map` are wrapped in `Arc` so they can be shared by threads.
/// - `all_reads` is a slice of reads to classify.
/// - Returns a tuple with:
///   1) `Vec<DNASequence>` of classified reads,
///   2) `Vec<DNASequence>` of unclassified reads,
///   3) `Vec<KrakenOutputLine>` of classification lines,
///   4) `TaxonCounts` (global map of taxon => read & k-mer counts).
pub fn classify_reads_parallel(
    db: Arc<KrakenDB>,
    parent_map: Arc<ParentMap>,
    all_reads: &[DNASequence],
    print_sequence_in_kraken: bool,
    only_classified_kraken_output: bool,
) -> (Vec<DNASequence>, Vec<DNASequence>, Vec<KrakenOutputLine>, TaxonCounts) {
    // Use Rayon to process reads in parallel
    all_reads
        .par_iter()
        .map(|dna| {
            // Each thread holds local partial results
            let mut local_taxon_counts = TaxonCounts::default();
            let mut local_kraken_lines = Vec::new();

            // Perform the existing single-read classification
            let was_classified = classify_sequence(
                dna,
                &[db.clone()],    // pass slice of 1 DB or more if needed
                &parent_map,
                &mut local_taxon_counts,
                print_sequence_in_kraken,
                only_classified_kraken_output,
                &mut local_kraken_lines,
            );

            // Return everything needed for merging
            (dna.clone(), was_classified, local_taxon_counts, local_kraken_lines)
        })
        // Now fold partial results within each worker thread
        .fold(
            // Identity for each thread
            || (
                Vec::new(),                // classified reads
                Vec::new(),                // unclassified reads
                Vec::new(),                // kraken output lines
                TaxonCounts::default(),    // taxon counts
            ),
            // Fold function merges partials within a thread
            |mut acc, (read, was_classified, local_map, local_lines)| {
                // Track read
                if was_classified {
                    acc.0.push(read);
                } else {
                    acc.1.push(read);
                }

                // Merge local taxon counts
                for (taxid, rc2) in local_map {
                    let rc1 = acc.3.entry(taxid).or_default();
                    rc1.read_count += rc2.read_count;
                    // Merge k-mer coverage
                    for (kmer, count) in rc2.kmer_counts {
                        *rc1.kmer_counts.entry(kmer).or_insert(0) += count;
                    }
                }

                // Merge classification lines
                acc.2.extend(local_lines);

                acc
            },
        )
        // Finally reduce across threads
        .reduce(
            // Another identity
            || (
                Vec::new(),                 // classified
                Vec::new(),                 // unclassified
                Vec::new(),                 // lines
                TaxonCounts::default(),     // taxon counts
            ),
            // Merge function merges two thread partials
            |(mut c1, mut u1, mut lines1, mut t1),
             (mut c2, mut u2, mut lines2, mut t2)| {
                c1.append(&mut c2);
                u1.append(&mut u2);
                lines1.append(&mut lines2);

                // Merge taxon maps
                for (taxid, rc2) in t2 {
                    let rc1 = t1.entry(taxid).or_default();
                    rc1.read_count += rc2.read_count;
                    for (kmer, count) in rc2.kmer_counts {
                        *rc1.kmer_counts.entry(kmer).or_insert(0) += count;
                    }
                }
                (c1, u1, lines1, t1)
            },
        )
}
