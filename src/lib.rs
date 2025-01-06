// src/lib.rs

pub mod types;
pub mod fastq;
pub mod krakendb;
pub mod taxdb;
pub mod classify;

use ahash::AHashMap;
use rayon::prelude::*;
use std::fmt::Write as FmtWrite;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use crate::classify::classify_reads::classify_reads_parallel;
use crate::classify::classify_sequence::TaxonCounts;
use crate::classify::classify_stats::build_kraken_report;
use crate::fastq::read_fastq_records;
use crate::krakendb::{KrakenDB, KrakenDBIndex};
use crate::taxdb::{read_taxdb_and_counts, NameMap, RankMap, TaxonomyData};
use crate::types::{DNASequence, KrakenOutputLine, KrakenReportRow};

/// A struct to hold classification results, with minimal duplication.
pub struct ClassificationResults {
    pub kraken_output_lines: Vec<KrakenOutputLine>,
    pub classified_reads: Vec<DNASequence>,
    pub unclassified_reads: Vec<DNASequence>,
    pub taxon_counts: TaxonCounts,
    pub kraken_report_rows: Option<Vec<KrakenReportRow>>,
    pub name_map: Option<Arc<NameMap>>,
    pub rank_map: Option<Arc<RankMap>>,
}

impl ClassificationResults {
    /// Generate Kraken output text on demand
    pub fn get_kraken_output(&self) -> String {
        // Since we could have as many lines as we have reads, reserve the approximate capacity
        // for the output string. This is an estimate: say ~100 bytes per line on average.
        // If you have extremely long sequences, adjust as needed.
        let capacity_estimate = self.kraken_output_lines.len().saturating_mul(100);
        let mut output = String::with_capacity(capacity_estimate);

        for line in &self.kraken_output_lines {
            // If the sequence was requested, append it
            let seq_part = match &line.sequence {
                Some(s) => format!("\t{}", s),
                None => String::new(),
            };

            // Using `writeln!` for convenience; you could manually push_str if you want
            // finer control or performance.
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}{}",
                line.status,
                line.read_id,
                line.tax_id,
                line.length,
                line.hitlist,
                seq_part
            ).unwrap();
        }

        output
    }

    /// Generate classified FASTQ output
    pub fn get_classified_reads_text(&self) -> String {
        // Reserve ~ (read length + header length + 10) per read
        let cap = self.classified_reads.len().saturating_mul(150);
        let mut output = String::with_capacity(cap);

        for read in &self.classified_reads {
            writeln!(output, "@{}\n{}\n+\n{}", read.header_line, read.seq, read.quals).unwrap();
        }
        output
    }

    /// Generate unclassified FASTQ output
    pub fn get_unclassified_reads_text(&self) -> String {
        let cap = self.unclassified_reads.len().saturating_mul(150);
        let mut output = String::with_capacity(cap);

        for read in &self.unclassified_reads {
            writeln!(output, "@{}\n{}\n+\n{}", read.header_line, read.seq, read.quals).unwrap();
        }
        output
    }

    /// Generate Kraken report text on demand
    pub fn get_kraken_report(&self) -> Option<String> {
        if let Some(rows) = &self.kraken_report_rows {
            // Rough estimate for lines, say ~80 bytes per row
            let cap = rows.len().saturating_mul(80);
            let mut output = String::with_capacity(cap);

            output.push_str("%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName\n");

            for row in rows {
                let mut indented_name = String::new();
                // If you expect big depth, you could even do `String::with_capacity(...)`.
                for _ in 0..row.depth {
                    indented_name.push('\t');
                }
                indented_name.push_str(&row.tax_name);

                writeln!(
                    output,
                    "{:.4}\t{}\t{}\t{}\t{:.4}\t{:.6}\t{}\t{}\t{}",
                    row.pct,
                    row.reads,
                    row.tax_reads,
                    row.kmers,
                    row.dup,
                    row.cov,
                    row.tax_id,
                    row.rank,
                    indented_name
                ).unwrap();
            }
            Some(output)
        } else {
            None
        }
    }
}

/// Unified function to classify reads from one or multiple files.
pub fn classify_reads(
    db_path: &str,
    db_index_path: &str,
    db_counts_path: &str,
    taxdb_path: &str,
    reads_paths: Vec<PathBuf>,
    print_sequence_in_kraken: bool,
    only_classified_kraken_output: bool,
    generate_report: bool,
) -> Result<ClassificationResults, Box<dyn std::error::Error>> {
    // STEP 1. Load & setup DB
    let t0 = Instant::now();
    let mut db = KrakenDB::new().open_file(db_path)?;
    let idx_data = std::fs::read(db_index_path)?;
    let db_index = KrakenDBIndex::from_slice(Arc::from(idx_data))?;
    db.index_ptr = Some(Arc::new(db_index));
    eprintln!("Step 1 (Load & setup DB) took: {:?}", t0.elapsed());

    // STEP 2. Parse taxonomy DB & counts
    let t1 = Instant::now();
    let taxonomy_data: TaxonomyData = read_taxdb_and_counts(taxdb_path, db_counts_path)?;
    let parent_map = Arc::new(taxonomy_data.parent_map);
    let name_map = Arc::new(taxonomy_data.name_map);
    let rank_map = Arc::new(taxonomy_data.rank_map);
    eprintln!("Step 2 (Parse taxonomy DB & counts) took: {:?}", t1.elapsed());

    // STEP 3. Read all FASTQ data in parallel.
    // - If you have many files, consider chunking them if needed.
    let t2 = Instant::now();
    // We can guess an upper bound if we know about how many reads each file might have.
    // For instance, if you typically expect 10k reads per file, with 10 files => 100k total.
    // That helps reduce reallocation. Adjust as appropriate.
    let files_count = reads_paths.len();
    let mut all_reads = Vec::with_capacity(files_count.saturating_mul(10_000));

    // We read each path in parallel, each returning a `Vec<DNASequence>`.
    // Then flatten. This approach is mostly the same as your original code,
    // but we do an extra local `collect` into a temporary so we can `extend` into `all_reads` (to pre-allocate).
    let all_reads_vecs = reads_paths
        .into_par_iter()
        .map(|path| read_fastq_records(path))
        .collect::<Result<Vec<Vec<DNASequence>>, _>>()?;

    // Flatten them without too many re-allocations:
    for read_vec in all_reads_vecs {
        all_reads.reserve(read_vec.len());
        all_reads.extend(read_vec);
    }
    eprintln!("Step 3 (Read all FASTQ data) took: {:?}", t2.elapsed());

    // STEP 4. Run parallel classification
    let t3 = Instant::now();
    let db_arc = Arc::new(db); // share DB across threads
    let (classified_reads, unclassified_reads, kraken_output_lines, taxon_counts) =
        classify_reads_parallel(
            db_arc,
            parent_map.clone(),
            &all_reads,
            print_sequence_in_kraken,
            only_classified_kraken_output,
        );
    eprintln!("Step 4 (Parallel classification) took: {:?}", t3.elapsed());

    // STEP 5. Build the Kraken report (if requested)
    let t4 = Instant::now();
    let kraken_report_rows = if generate_report {
        let total_reads = all_reads.len() as u32;

        let (rows, _text) = build_kraken_report(
            &taxon_counts,
            &parent_map,
            Some(&name_map),
            Some(&rank_map),
            1, // root taxid
            total_reads,
            &taxonomy_data.total_counts,
            &taxonomy_data.direct_counts,
        );
        Some(rows)
    } else {
        None
    };
    eprintln!("Step 5 (Build report) took: {:?}", t4.elapsed());
    // Now convert to Vec<DNASequence> by cloning each read
    let classified_return: Vec<DNASequence> =
        classified_reads.iter().map(|r| (*r).clone()).collect();

    let unclassified_return: Vec<DNASequence> =
        unclassified_reads.iter().map(|r| (*r).clone()).collect();

    // STEP 6. Return final results
    Ok(ClassificationResults {
        kraken_output_lines,
        classified_reads: classified_return,
        unclassified_reads: unclassified_return,
        taxon_counts,
        kraken_report_rows,
        name_map: Some(name_map),
        rank_map: Some(rank_map),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn test_classify_reads_api() {
        // 1. Gather all *.fastq or *.fastq.gz files under "barcode03"
        let read_files: Vec<PathBuf> = fs::read_dir("./barcode03")
            .expect("Cannot read directory 'barcode03'")
            .filter_map(|entry| {
                let path = entry.ok()?.path();
                let filename = path.file_name()?.to_string_lossy().to_lowercase();
                if filename.ends_with(".fastq") || filename.ends_with(".fastq.gz") {
                    Some(path)
                } else {
                    None
                }
            })
            .collect();

        // 2. Call the classify_reads function
        let results = classify_reads(
            "./pr2/database.kdb",
            "./pr2/database.idx",
            "./pr2/database.kdb.counts",
            "./pr2/taxDB",
            read_files,            // Pass all FASTQ paths
            false,
            false,
            true
        ).expect("Classification failed");

        // 3. Write outputs using the new getters
        fs::write("kraken_output.txt", results.get_kraken_output())
            .expect("Could not write kraken_output.txt");

        fs::write("classified_reads.fq", results.get_classified_reads_text())
            .expect("Could not write classified_reads.fq");

        fs::write("unclassified_reads.fq", results.get_unclassified_reads_text())
            .expect("Could not write unclassified_reads.fq");

        // 4. Write report if generated
        if let Some(report_text) = results.get_kraken_report() {
            fs::write("kraken_report.txt", report_text)
                .expect("Could not write kraken_report.txt");
        }

        // 5. Basic checks
        assert!(!results.kraken_output_lines.is_empty(), "Expected some classification output");

        let total_reads = results.classified_reads.len() + results.unclassified_reads.len();
        assert!(total_reads > 0, "Expected some reads to be processed");

        // Basic checks on classification lines
        for line in &results.kraken_output_lines {
            match line.status {
                'C' => assert!(
                    results.classified_reads.iter().any(|r| r.id == line.read_id),
                    "Classified read not found in classified_reads"
                ),
                'U' => assert!(
                    results.unclassified_reads.iter().any(|r| r.id == line.read_id),
                    "Unclassified read not found in unclassified_reads"
                ),
                _ => panic!("Invalid classification status: {}", line.status),
            }
        }

        // 6. Check taxonomy report
        if let Some(rows) = &results.kraken_report_rows {
            assert!(!rows.is_empty(), "Expected some taxonomy report rows");
        }

        // 7. Compare the output to a known reference, ignoring line order
        let kraken_output = fs::read_to_string("kraken_output.txt").expect("Could not read kraken_output.txt");
        let real_report = fs::read_to_string("realoutput.txt").expect("Could not read realoutput.txt");

        let mut lines_kraken: Vec<_> = kraken_output.lines().collect();
        lines_kraken.sort_unstable();

        let mut lines_real: Vec<_> = real_report.lines().collect();
        lines_real.sort_unstable();

        assert_eq!(
            lines_kraken, lines_real,
            "kraken_output.txt differs from realoutput.txt (ignoring order)."
        );
    }
}
