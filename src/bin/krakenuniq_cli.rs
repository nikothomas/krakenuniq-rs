use indicatif::{ProgressBar, ProgressStyle};
use std::fs;
use std::path::PathBuf;

// Replace with your actual krakenuniq-rs import
use krakenuniq_rs::classify_reads;

fn main() {
    // 1. Spinner for gathering *.fastq or *.fastq.gz files
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&[
                "⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏",
            ])
            .template("{spinner:.blue} {msg}")
            .expect("Invalid spinner template"),
    );
    spinner.set_message("Gathering FASTQ files under 'barcode03'...");

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

    spinner.finish_with_message(format!("Found {} FASTQ file(s).", read_files.len()));

    // 2. Spinner for classifying reads
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&[
                "⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏",
            ])
            .template("{spinner:.green} {msg}")
            .expect("Invalid spinner template"),
    );
    spinner.set_message("Classifying reads...");

    let results = classify_reads(
        "./pr2/database.kdb",
        "./pr2/database.idx",
        "./pr2/database.kdb.counts",
        "./pr2/taxDB",
        read_files,
        false,
        false,
        true,
    )
        .expect("Classification failed");

    spinner.finish_with_message("Classification finished.");

    // 3. Spinner for writing outputs
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&[
                "⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏",
            ])
            .template("{spinner:.yellow} {msg}")
            .expect("Invalid spinner template"),
    );
    spinner.set_message("Writing output files...");

    fs::write("kraken_output.txt", results.get_kraken_output())
        .expect("Could not write kraken_output.txt");

    fs::write("classified_reads.fq", results.get_classified_reads_text())
        .expect("Could not write classified_reads.fq");

    fs::write("unclassified_reads.fq", results.get_unclassified_reads_text())
        .expect("Could not write unclassified_reads.fq");

    // Optionally write the Kraken report if generated
    if let Some(report_text) = results.get_kraken_report() {
        fs::write("kraken_report.txt", report_text)
            .expect("Could not write kraken_report.txt");
    }

    spinner.finish_with_message("Output files created.");

    // 4. Final message
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&[
                "⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏",
            ])
            .template("{spinner:.cyan} {msg}")
            .expect("Invalid spinner template"),
    );
    spinner.set_message("All done!");
    spinner.finish_with_message("All done!");
}
