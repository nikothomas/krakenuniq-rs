use rayon::prelude::*;
use std::sync::Arc;

use super::{classify_sequence, TaxonCounts};
use super::KrakenDB;
use super::ParentMap;
use crate::types::{DNASequence, KrakenOutputLine};

/// Parallel classification of multiple reads at once.
///
/// This version avoids cloning the entire `DNASequence`, instead returning
/// references to `DNASequence` in the result. We also pre-allocate space in
/// vectors and parallel fold over reads (only one level of parallelism).
pub fn classify_reads_parallel<'a>(
    db: Arc<KrakenDB>,
    parent_map: Arc<ParentMap>,
    all_reads: &'a [DNASequence],
    print_sequence_in_kraken: bool,
    only_classified_kraken_output: bool,
) -> (
    Vec<&'a DNASequence>,   // classified
    Vec<&'a DNASequence>,   // unclassified
    Vec<KrakenOutputLine>,  // kraken output lines
    TaxonCounts,            // global taxon counts
) {
    // Optional: chunk the reads if extremely large and each is small.
    // Uncomment and adjust if beneficial for your workload:
    // return all_reads
    //     .par_chunks(1024)
    //     .map(|chunk| {
    //         classify_reads_chunk(
    //             db.clone(), parent_map.clone(), chunk,
    //             print_sequence_in_kraken, only_classified_kraken_output
    //         )
    //     })
    //     .reduce(
    //         || (
    //             Vec::with_capacity(1024),
    //             Vec::with_capacity(1024),
    //             Vec::with_capacity(1024),
    //             TaxonCounts::default(),
    //         ),
    //         merge_partial_results
    //     );

    // If not chunking, we just do par_iter over all reads:
    all_reads
        .par_iter()
        .fold(
            // Thread-local identity
            || (
                Vec::with_capacity(256),     // classified
                Vec::with_capacity(256),     // unclassified
                Vec::with_capacity(256),     // local output lines
                TaxonCounts::default(),       // local taxon counts
            ),
            // Per-read classification
            |mut acc, dna| {
                let mut local_taxon_counts = TaxonCounts::default();
                let mut local_kraken_lines = Vec::with_capacity(1);

                let was_classified = classify_sequence(
                    dna,
                    &[db.clone()],
                    &parent_map,
                    &mut local_taxon_counts,
                    print_sequence_in_kraken,
                    only_classified_kraken_output,
                    &mut local_kraken_lines,
                );

                if was_classified {
                    acc.0.push(dna);
                } else {
                    acc.1.push(dna);
                }

                // Merge local taxon counts
                for (taxid, rc2) in local_taxon_counts {
                    let rc1 = acc.3.entry(taxid).or_default();
                    rc1.read_count += rc2.read_count;
                    for (kmer, count) in rc2.kmer_counts {
                        *rc1.kmer_counts.entry(kmer).or_insert(0) += count;
                    }
                }

                // Merge local classification lines
                acc.2.extend(local_kraken_lines);
                acc
            },
        )
        .reduce(
            // Another identity
            || (
                Vec::with_capacity(256),
                Vec::with_capacity(256),
                Vec::with_capacity(256),
                TaxonCounts::default(),
            ),
            // Merge partial results from different threads
            merge_partial_results,
        )
}

/// Merges two partial classification results from different threads.
fn merge_partial_results<'a>(
    mut a: (Vec<&'a DNASequence>, Vec<&'a DNASequence>, Vec<KrakenOutputLine>, TaxonCounts),
    mut b: (Vec<&'a DNASequence>, Vec<&'a DNASequence>, Vec<KrakenOutputLine>, TaxonCounts),
) -> (Vec<&'a DNASequence>, Vec<&'a DNASequence>, Vec<KrakenOutputLine>, TaxonCounts) {
    // Merge vectors
    a.0.append(&mut b.0);
    a.1.append(&mut b.1);
    a.2.append(&mut b.2);

    // Merge taxon counts
    let t1 = &mut a.3;
    t1.reserve(b.3.len());
    for (taxid, rc2) in b.3 {
        let rc1 = t1.entry(taxid).or_default();
        rc1.read_count += rc2.read_count;
        for (kmer, count) in rc2.kmer_counts {
            *rc1.kmer_counts.entry(kmer).or_insert(0) += count;
        }
    }

    a
}

/// Optional helper if chunking reads above.
#[allow(dead_code)]
fn classify_reads_chunk<'a>(
    db: Arc<KrakenDB>,
    parent_map: Arc<ParentMap>,
    chunk: &'a [DNASequence],
    print_sequence_in_kraken: bool,
    only_classified_kraken_output: bool,
) -> (
    Vec<&'a DNASequence>,
    Vec<&'a DNASequence>,
    Vec<KrakenOutputLine>,
    TaxonCounts
) {
    let mut classified = Vec::with_capacity(chunk.len());
    let mut unclassified = Vec::with_capacity(chunk.len());
    let mut lines = Vec::with_capacity(chunk.len());
    let mut tax_counts = TaxonCounts::default();

    for dna in chunk {
        let mut local_taxon_counts = TaxonCounts::default();
        let mut local_lines = Vec::with_capacity(1);

        let was_classified = classify_sequence(
            dna,
            &[db.clone()],
            &parent_map,
            &mut local_taxon_counts,
            print_sequence_in_kraken,
            only_classified_kraken_output,
            &mut local_lines,
        );

        if was_classified {
            classified.push(dna);
        } else {
            unclassified.push(dna);
        }

        for (taxid, rc2) in local_taxon_counts {
            let rc1 = tax_counts.entry(taxid).or_default();
            rc1.read_count += rc2.read_count;
            for (kmer, count) in rc2.kmer_counts {
                *rc1.kmer_counts.entry(kmer).or_insert(0) += count;
            }
        }
        lines.extend(local_lines);
    }

    (classified, unclassified, lines, tax_counts)
}
