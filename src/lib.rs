pub mod classify_sequence;
pub mod krakendb;
pub mod taxdb;
pub mod fastq;
pub mod classifications_stats;
pub mod types;

use ahash::AHashMap;
use std::fmt::Write as FmtWrite;
use std::sync::Arc;

// We'll import our newly created structured types:
use crate::types::{KrakenOutputLine, KrakenReportRow};

use crate::classifications_stats::build_kraken_report;
use crate::classify_sequence::{classify_sequence, DNASequence, TaxonCounts};
use crate::fastq::read_fastq_records;
use crate::krakendb::{KrakenDB, KrakenDBIndex};
use crate::taxdb::{parse_taxdb, NameMap, RankMap};

/// A simple struct to hold the final outputs in memory.
pub struct ClassificationResults {
    /// Kraken-style output lines as text (C/U <tab> read_id <tab> taxid <tab> length <tab> "hitlist").
    pub kraken_output: String,

    /// FASTQ-like text of all classified reads.
    pub classified_reads: String,

    /// FASTQ-like text of all unclassified reads.
    pub unclassified_reads: String,

    /// A map from `taxid` -> read counts & unique k-mers (`ReadCounts`).
    pub taxon_counts: TaxonCounts,

    /// CSV-style lines summarizing `(taxid, total_reads, unique_kmer_count)`.
    pub taxon_counts_csv: String,

    /// Optionally, the final Kraken-style taxonomy report as text (if generated).
    pub kraken_report: Option<String>,

    /// New: the structured version of the Kraken output lines (one per read).
    pub kraken_output_lines: Vec<KrakenOutputLine>,

    /// New: the structured version of the taxonomy report rows (if generated).
    pub kraken_report_rows: Option<Vec<KrakenReportRow>>,

    /// Optionally, name and rank maps if you need them later.
    pub name_map: Option<NameMap>,
    pub rank_map: Option<RankMap>,
}

/// Classify reads in one go, returning in-memory results.
///
/// # Arguments
/// - `db_path`: Path to a `.kdb` file (Kraken DB).
/// - `db_index_path`: Path to a matching `.idx` file.
/// - `taxdb_path`: Path to a tab-delimited "taxDB" file (format: `<taxid>\t<parent>\t<name>\t<rank>`).
/// - `reads_path`: Path to a FASTQ file containing reads.
/// - `print_sequence_in_kraken`: If `true`, append the read sequence to each `kraken_output` line.
/// - `only_classified_kraken_output`: If `true`, omit lines for unclassified reads from `kraken_output`.
/// - `generate_report`: If `true`, also compute a Kraken-style taxonomic report.
///
/// # Returns
/// On success, returns a `ClassificationResults` struct containing:
///   - `kraken_output` + `kraken_output_lines`
///   - `classified_reads`
///   - `unclassified_reads`
///   - `taxon_counts` (HashMap of taxid -> ReadCounts)
///   - `taxon_counts_csv` (a CSV summarizing the taxon counts)
///   - `kraken_report` + `kraken_report_rows` (if `generate_report=true`)
///   - `name_map` and `rank_map` (optional lookups)
///
/// On failure, returns `Err` with a boxed error.
pub fn classify_in_one_call(
    db_path: &str,
    db_index_path: &str,
    taxdb_path: &str,
    reads_path: &str,
    print_sequence_in_kraken: bool,
    only_classified_kraken_output: bool,
    generate_report: bool,
) -> Result<ClassificationResults, Box<dyn std::error::Error>> {
    // 1. Load Kraken DB
    let mut db = KrakenDB::new().open_file(db_path)?;

    // 2. Load the .idx file into memory, parse it, attach to DB
    let idx_data = std::fs::read(db_index_path)?;
    let db_index = KrakenDBIndex::from_slice(Arc::from(idx_data))?;
    db.index_ptr = Some(Arc::new(db_index));

    // 3. Parse taxDB => (parent_map, name_map, rank_map)
    let (parent_map, name_map, rank_map) = parse_taxdb(taxdb_path)?;

    // 4. Read FASTQ
    let reads = read_fastq_records(reads_path)?;

    // Prepare classification-time outputs
    let mut kraken_output = String::new();
    let mut classified_reads = String::new();
    let mut unclassified_reads = String::new();
    let mut taxon_counts: TaxonCounts = AHashMap::new();
    // We'll store the structured lines here:
    let mut kraken_output_structs = Vec::new(); // Vec<KrakenOutputLine>

    // 5. Classify each read
    for dna in &reads {
        classify_sequence(
            dna,
            &[db.clone()],
            &parent_map,
            &mut taxon_counts,
            &mut kraken_output,
            &mut classified_reads,
            &mut unclassified_reads,
            print_sequence_in_kraken,
            only_classified_kraken_output,
            &mut kraken_output_structs, // <-- collecting structured lines
        );
    }

    // 6. Build a CSV summary of taxon counts
    let mut counts_csv = String::new();
    for (taxid, counts) in &taxon_counts {
        let _ = writeln!(
            counts_csv,
            "{},{},{}",
            taxid,
            counts.read_count,
            counts.unique_kmers.len()
        );
    }

    // 7. Optionally build a Kraken-style taxonomy report (both text & structured rows)
    let (kraken_report_rows, kraken_report_text) = if generate_report {
        let total_reads = reads.len() as u64;
        let (rows, text) = build_kraken_report(
            &taxon_counts,
            &parent_map,
            Some(&name_map),
            Some(&rank_map),
            1, // root taxid
            total_reads,
        );
        (Some(rows), Some(text))
    } else {
        (None, None)
    };

    Ok(ClassificationResults {
        kraken_output,
        classified_reads,
        unclassified_reads,
        taxon_counts,
        taxon_counts_csv: counts_csv,

        // We'll store the text from the taxonomy report (if generated)
        kraken_report: kraken_report_text,

        // Our newly added structured results:
        kraken_output_lines: kraken_output_structs,
        kraken_report_rows,

        name_map: Some(name_map),
        rank_map: Some(rank_map),
    })
}

/// Classify reads from multiple FASTQ files in one go, returning in-memory results.
///
/// # Arguments
/// - `db_path`: Path to a `.kdb` file (Kraken DB).
/// - `db_index_path`: Path to a matching `.idx` file.
/// - `taxdb_path`: Path to a tab-delimited "taxDB" file (`<taxid>\t<parent>\t<name>\t<rank>`).
/// - `reads_paths`: A list of paths to FASTQ or FASTQ.GZ files.
/// - `print_sequence_in_kraken`: If `true`, append the read sequence to each `kraken_output` line.
/// - `only_classified_kraken_output`: If `true`, omit lines for unclassified reads from `kraken_output`.
/// - `generate_report`: If `true`, compute a Kraken-style taxonomic report.
///
/// # Returns
/// A `ClassificationResults` struct on success, or `Err` on failure.
pub fn classify_multiple_in_one_call(
    db_path: &str,
    db_index_path: &str,
    taxdb_path: &str,
    reads_paths: Vec<String>,
    print_sequence_in_kraken: bool,
    only_classified_kraken_output: bool,
    generate_report: bool,
) -> Result<ClassificationResults, Box<dyn std::error::Error>> {
    // 1. Load Kraken DB
    let mut db = KrakenDB::new().open_file(db_path)?;

    // 2. Load the .idx file into memory
    let idx_data = std::fs::read(db_index_path)?;
    let db_index = KrakenDBIndex::from_slice(Arc::from(idx_data))?;
    db.index_ptr = Some(Arc::new(db_index));

    // 3. Parse taxDB => (parent_map, name_map, rank_map)
    let (parent_map, name_map, rank_map) = parse_taxdb(taxdb_path)?;

    // 4. Read all FASTQ/FASTQ.GZ files
    let mut all_reads: Vec<DNASequence> = Vec::new();
    for path in reads_paths {
        let reads = read_fastq_records(&*path.to_string())?;
        all_reads.extend(reads);
    }

    // Prepare classification-time outputs
    let mut kraken_output = String::new();
    let mut classified_reads = String::new();
    let mut unclassified_reads = String::new();
    let mut taxon_counts: TaxonCounts = AHashMap::new();
    let mut kraken_output_structs = Vec::new(); // collect structured lines

    // 5. Classify each read
    for dna in &all_reads {
        classify_sequence(
            dna,
            &[db.clone()],
            &parent_map,
            &mut taxon_counts,
            &mut kraken_output,
            &mut classified_reads,
            &mut unclassified_reads,
            print_sequence_in_kraken,
            only_classified_kraken_output,
            &mut kraken_output_structs,
        );
    }

    // 6. Build a CSV summary of taxon counts
    let mut counts_csv = String::new();
    for (taxid, counts) in &taxon_counts {
        let _ = writeln!(
            counts_csv,
            "{},{},{}",
            taxid,
            counts.read_count,
            counts.unique_kmers.len()
        );
    }

    // 7. Optionally build a Kraken-style taxonomy report
    let (kraken_report_rows, kraken_report_text) = if generate_report {
        let total_reads = all_reads.len() as u64;
        let (rows, text) = build_kraken_report(
            &taxon_counts,
            &parent_map,
            Some(&name_map),
            Some(&rank_map),
            1, // root taxid
            total_reads,
        );
        (Some(rows), Some(text))
    } else {
        (None, None)
    };

    Ok(ClassificationResults {
        kraken_output,
        classified_reads,
        unclassified_reads,
        taxon_counts,
        taxon_counts_csv: counts_csv,
        kraken_report: kraken_report_text,
        kraken_output_lines: kraken_output_structs,
        kraken_report_rows,
        name_map: Some(name_map),
        rank_map: Some(rank_map),
    })
}

// -------------------
// Unit Test
// -------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_classify_in_one_call_api() {
        // Example usage:
        let results = classify_in_one_call(
            "./pr2/database.kdb",
            "./pr2/database.idx",
            "./pr2/taxDB",
            "reads.fastq",
            /* print_sequence_in_kraken = */ false,
            /* only_classified_kraken_output = */ false,
            /* generate_report = */ true,
        ).expect("Classification failed");

        // (Optional) Write the outputs
        fs::write("kraken_output.txt", &results.kraken_output)
            .expect("Could not write kraken_output.txt");
        fs::write("classified_reads.fq", &results.classified_reads)
            .expect("Could not write classified_reads.fq");
        fs::write("unclassified_reads.fq", &results.unclassified_reads)
            .expect("Could not write unclassified_reads.fq");
        fs::write("taxon_counts.csv", &results.taxon_counts_csv)
            .expect("Could not write taxon_counts.csv");

        if let Some(report_text) = &results.kraken_report {
            fs::write("kraken_report.txt", report_text)
                .expect("Could not write kraken_report.txt");
        }

        // Basic checks
        assert!(!results.kraken_output.is_empty(), "Expected some classification output");
        assert!(results.taxon_counts_csv.contains(','), "Expected CSV data for taxon counts");

        // Inspect the structured lines in memory:
        assert!(!results.kraken_output_lines.is_empty(), "Expected some structured Kraken lines");
        for line in &results.kraken_output_lines {
            eprintln!("Read: {}, status: {}, taxID={}", line.read_id, line.status, line.tax_id);
        }

        if let Some(rows) = &results.kraken_report_rows {
            // We have a structured taxonomy report
            assert!(!rows.is_empty(), "Expected some tax report rows");
            // Example: show top row
            if let Some(first_row) = rows.first() {
                eprintln!("First row => tax {} => pct={:.2}", first_row.tax_id, first_row.pct);
            }
        }
    }
}
