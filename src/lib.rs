// src/lib.rs
pub mod types;
pub mod classifications_stats;
pub mod classify_sequence;
pub mod fastq;
pub mod krakendb;
pub mod taxdb;
pub mod db_counts;

use ahash::AHashMap;
use std::fmt::Write as FmtWrite;
use std::sync::Arc;
use std::path::PathBuf;

use crate::types::{KrakenOutputLine, KrakenReportRow};
use crate::classifications_stats::build_kraken_report;
use crate::classify_sequence::{classify_sequence, DNASequence, TaxonCounts};
use crate::fastq::read_fastq_records;
use crate::krakendb::{KrakenDB, KrakenDBIndex};
use crate::taxdb::{parse_taxdb, NameMap, RankMap};

/// A struct to hold classification results with minimal duplication.
/// Now we only store structured data and generate text on demand.
pub struct ClassificationResults {
    /// Structured version of the Kraken output lines (one per read)
    pub kraken_output_lines: Vec<KrakenOutputLine>,

    /// Classified and unclassified reads stored as DNASequence to save memory
    pub classified_reads: Vec<DNASequence>,
    pub unclassified_reads: Vec<DNASequence>,

    /// A map from `taxid` -> read counts & unique k-mers (`ReadCounts`)
    pub taxon_counts: TaxonCounts,

    /// Structured version of the taxonomy report rows (if generated)
    pub kraken_report_rows: Option<Vec<KrakenReportRow>>,

    /// Optional lookups (weak references to avoid duplication)
    pub name_map: Option<Arc<NameMap>>,
    pub rank_map: Option<Arc<RankMap>>,
}

impl ClassificationResults {
    /// Generate Kraken output text on demand
    pub fn get_kraken_output(&self) -> String {
        let mut output = String::new();
        for line in &self.kraken_output_lines {
            let seq_part = line.sequence.as_ref().map(|s| format!("\t{}", s)).unwrap_or_default();
            writeln!(output, "{}\t{}\t{}\t{}\t{}{}",
                     line.status, line.read_id, line.tax_id, line.length, line.hitlist, seq_part
            ).unwrap();
        }
        output
    }

    /// Generate classified reads text on demand
    pub fn get_classified_reads_text(&self) -> String {
        let mut output = String::new();
        for read in &self.classified_reads {
            writeln!(output, "@{}\n{}\n+\n{}", read.header_line, read.seq, read.quals).unwrap();
        }
        output
    }

    /// Generate unclassified reads text on demand
    pub fn get_unclassified_reads_text(&self) -> String {
        let mut output = String::new();
        for read in &self.unclassified_reads {
            writeln!(output, "@{}\n{}\n+\n{}", read.header_line, read.seq, read.quals).unwrap();
        }
        output
    }

    /// Generate Kraken report text on demand
    pub fn get_kraken_report(&self) -> Option<String> {
        if let Some(rows) = &self.kraken_report_rows {
            let mut output = String::new();
            output.push_str("%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName\n");

            for row in rows {
                let mut indented_name = String::new();
                for _ in 0..row.depth {
                    indented_name.push('\t');
                }
                indented_name.push_str(&row.tax_name);

                writeln!(output,
                         "{:.4}\t{}\t{}\t{}\t{:.4}\t{:.6}\t{}\t{}\t{}",
                         row.pct, row.reads, row.tax_reads, row.kmers,
                         row.dup, row.cov, row.tax_id, row.rank, indented_name
                ).unwrap();
            }
            Some(output)
        } else {
            None
        }
    }
}

/// Unified function to classify reads from one or multiple files
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
    // 1. Load and setup database
    let mut db = KrakenDB::new().open_file(db_path)?;
    let idx_data = std::fs::read(db_index_path)?;
    let db_index = KrakenDBIndex::from_slice(Arc::from(idx_data))?;
    db.index_ptr = Some(Arc::new(db_index));

    // 2. Parse taxonomy database
    let (parent_map, name_map, rank_map) = parse_taxdb(taxdb_path)?;
    let name_map = Arc::new(name_map);
    let rank_map = Arc::new(rank_map);

    // 3. Process all reads
    let mut all_reads = Vec::new();
    for path in reads_paths {
        let reads = read_fastq_records(path)?;
        all_reads.extend(reads);
    }

    // 4. Prepare outputs
    let mut classified_reads = Vec::new();
    let mut unclassified_reads = Vec::new();
    let mut kraken_output_lines = Vec::new();
    let mut taxon_counts: TaxonCounts = AHashMap::new();

    // 5. Classify each read
    for dna in &all_reads {
        let was_classified = classify_sequence(
            dna,
            &[db.clone()],
            &parent_map,
            &mut taxon_counts,
            print_sequence_in_kraken,
            only_classified_kraken_output,
            &mut kraken_output_lines,
        );

        if was_classified {
            classified_reads.push(dna.clone());
        } else {
            unclassified_reads.push(dna.clone());
        }
    }

    // 6. Generate report if requested
    let kraken_report_rows = if generate_report {
        let total_reads = all_reads.len() as u32;
        let (rows, _) = build_kraken_report(
            &taxon_counts,
            &parent_map,
            Some(&name_map),
            Some(&rank_map),
            1,
            total_reads,
            taxdb_path,
            db_counts_path,
        );
        Some(rows)
    } else {
        None
    };

    Ok(ClassificationResults {
        kraken_output_lines,
        classified_reads,
        unclassified_reads,
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
        // Example usage with modern PathBuf handling
        let results = classify_reads(
            "./pr2/database.kdb",
            "./pr2/database.idx",
            "./pr2/database.kdb.counts",
            "./pr2/taxDB",
            vec![PathBuf::from("reads.fastq")],
            false,
            false,
            true
        ).expect("Classification failed");

        // Write outputs using the new getters
        fs::write("kraken_output.txt", results.get_kraken_output())
            .expect("Could not write kraken_output.txt");

        fs::write("classified_reads.fq", results.get_classified_reads_text())
            .expect("Could not write classified_reads.fq");

        fs::write("unclassified_reads.fq", results.get_unclassified_reads_text())
            .expect("Could not write unclassified_reads.fq");

        // Write report if generated
        if let Some(report_text) = results.get_kraken_report() {
            fs::write("kraken_report.txt", report_text)
                .expect("Could not write kraken_report.txt");
        }

        // Test the structured data
        assert!(!results.kraken_output_lines.is_empty(), "Expected some classification output");

        // Check classified/unclassified reads
        let total_reads = results.classified_reads.len() + results.unclassified_reads.len();
        assert!(total_reads > 0, "Expected some reads to be processed");

        // Inspect classification results
        for line in &results.kraken_output_lines {
            eprintln!("Read: {}, status: {}, taxID={}", line.read_id, line.status, line.tax_id);

            // Validate classification status matches the reads vectors
            match line.status {
                'C' => assert!(results.classified_reads.iter().any(|r| r.id == line.read_id)),
                'U' => assert!(results.unclassified_reads.iter().any(|r| r.id == line.read_id)),
                _ => panic!("Invalid classification status"),
            }
        }

        // Check taxonomy report
        if let Some(rows) = &results.kraken_report_rows {
            assert!(!rows.is_empty(), "Expected some taxonomy report rows");

            // Validate first row
            if let Some(first_row) = rows.first() {
                eprintln!("First row => tax {} => pct={:.2}", first_row.tax_id, first_row.pct);

                // Basic sanity checks
                assert!(first_row.pct >= 0.0 && first_row.pct <= 100.0, "Invalid percentage");
                assert!(first_row.reads > 0, "Expected non-zero reads in first row");
                assert!(!first_row.rank.is_empty(), "Expected non-empty rank");
                assert!(!first_row.tax_name.is_empty(), "Expected non-empty taxonomy name");
            }
        }

        // Verify taxon counts
        assert!(!results.taxon_counts.is_empty(), "Expected some taxon counts");

        // Verify name and rank maps are present if report was generated
        if results.kraken_report_rows.is_some() {
            assert!(results.name_map.is_some(), "Expected name map with report");
            assert!(results.rank_map.is_some(), "Expected rank map with report");
        }
    }
}