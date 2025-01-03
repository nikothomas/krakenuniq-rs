//src/types.rs

/// A structured representation of one row in the Kraken taxonomy report.
/// For example:
///  %  reads  taxReads  kmers  taxKmers  kmersDB  taxKmersDB  dup  cov  taxID  rank  taxName
#[derive(Debug, Clone)]
pub struct KrakenReportRow {
    pub pct: f32,
    pub reads: u32,
    pub tax_reads: u32,
    pub kmers: u32,
    pub dup: f32,
    pub cov: f32,
    pub tax_id: u32,
    pub rank: String,
    pub tax_name: String,
    pub depth: usize,                // Still needed for indentation
    pub parent_tax_id: Option<u32>,  // Keep for tree structure
    pub children_tax_ids: Vec<u32>,  // Keep for tree structure
}

/// A structured representation of one Kraken output line.
#[derive(Debug, Clone)]
pub struct KrakenOutputLine {
    pub status: char,          // 'C' or 'U'
    pub read_id: String,
    pub tax_id: u32,
    pub length: usize,
    pub hitlist: String,
    pub sequence: Option<String>, // if you store the read's sequence
}

/// A minimal representation of a read.
#[derive(Debug, Clone)]
pub struct DNASequence {
    pub id: String,
    pub header_line: String,
    pub seq: String,
    pub quals: String,
}
