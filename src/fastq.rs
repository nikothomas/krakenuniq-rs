use std::fs::File;
use std::io::{BufRead, BufReader};
use crate::classify_sequence::DNASequence;

// Add these new imports:
use flate2::read::MultiGzDecoder;
use std::path::Path;

/// Minimal FASTQ read function that also supports .gz
pub fn read_fastq_records(path: &str) -> std::io::Result<Vec<DNASequence>> {
    // Open the file:
    let f = File::open(path)?;

    // If the file ends with ".gz", wrap it in a MultiGzDecoder
    let is_gz = Path::new(path)
        .extension()
        .map(|ext| ext == "gz")
        .unwrap_or(false);

    let mut reader: Box<dyn BufRead> = if is_gz {
        Box::new(BufReader::new(MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    };

    let mut sequences = Vec::new();
    let mut line = String::new();

    // We'll reuse the same loop logic, but reading from `reader` instead of `br`
    loop {
        line.clear();
        // 1) read header
        if reader.read_line(&mut line)? == 0 {
            break; // EOF
        }
        let header_line = line.trim_end().to_string();
        if !header_line.starts_with('@') {
            // Not a valid FASTQ header; skip or handle error
            continue;
        }
        // remove '@'
        let header_str = &header_line[1..];

        // 2) read sequence
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            // malformed?
            break;
        }
        let seq_str = line.trim_end().to_string();

        // 3) read plus line
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        // 4) read quality
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let qual_str = line.trim_end().to_string();

        // Create DNASequence
        let dna = DNASequence {
            id: header_str.split(' ').next().unwrap().to_string(),
            header_line: header_str.to_string(),
            seq: seq_str,
            quals: qual_str,
        };
        sequences.push(dna);
    }

    Ok(sequences)
}
