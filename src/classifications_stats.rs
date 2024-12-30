// src/classification_stats.rs

use ahash::AHashMap;
use crate::types::KrakenReportRow;
use super::classify_sequence::TaxonCounts;
use super::taxdb::{ParentMap, NameMap, RankMap};

/// Per-node stats that store:
///   - read counts (self & clade)
///   - distinct k-mer counts & total coverage (self & clade)
///   - optional placeholders for DB-based metrics
#[derive(Default, Debug, Clone)]
pub struct NodeStats {
    /// Reads directly assigned to this node
    pub self_reads: u64,
    /// Reads in node + descendants
    pub clade_reads: u64,

    /// How many **distinct** k-mers in this node alone
    pub self_num_kmers_distinct: usize,
    /// Total coverage for this node's k-mers (sum of counts)
    pub self_num_kmers_coverage: u64,

    /// Distinct k-mers in clade (node + descendants)
    pub clade_num_kmers_distinct: usize,
    /// Coverage in clade (node + descendants)
    pub clade_num_kmers_coverage: u64,

    // Optional placeholders for DB-based metrics
    pub self_kmers_db: u64,
    pub clade_kmers_db: u64,
    pub self_tax_kmers_db: u64,
    pub clade_tax_kmers_db: u64,
}

impl NodeStats {
    /// Return duplication factor for this clade = coverage / distinct.
    pub fn clade_duplication(&self) -> f64 {
        if self.clade_num_kmers_distinct == 0 {
            0.0
        } else {
            self.clade_num_kmers_coverage as f64
                / self.clade_num_kmers_distinct as f64
        }
    }
}

/// Build a map of `parent -> Vec<child>` for traversing the taxonomy.
pub fn build_children_map(parent_map: &ParentMap) -> AHashMap<u32, Vec<u32>> {
    let mut children_map: AHashMap<u32, Vec<u32>> = AHashMap::new();

    // Initialize every known taxid as a key to avoid missing entries
    for &taxid in parent_map.keys() {
        children_map.entry(taxid).or_default();
    }

    // Populate children
    for (&child, &parent) in parent_map {
        if parent > 0 && child != parent {
            children_map.entry(parent).or_default().push(child);
        }
    }
    children_map
}

/// Initialize each NodeStats from `TaxonCounts`, using coverage-based approach:
/// - `self_reads` from `ReadCounts.read_count`
/// - `self_num_kmers_distinct` = `read_counts.kmer_counts.len()`
/// - `self_num_kmers_coverage` = sum of values in `kmer_counts`
pub fn init_node_stats(taxon_counts: &TaxonCounts) -> AHashMap<u32, NodeStats> {
    let mut stats_map = AHashMap::new();

    for (&taxid, read_counts) in taxon_counts {
        let mut node_stats = NodeStats::default();

        node_stats.self_reads = read_counts.read_count;
        node_stats.clade_reads = read_counts.read_count;

        // Distinct k-mers for this taxon
        let distinct = read_counts.kmer_counts.len();
        // Sum coverage
        let coverage = read_counts
            .kmer_counts
            .values()
            .fold(0u64, |acc, &v| acc + v as u64);

        node_stats.self_num_kmers_distinct = distinct;
        node_stats.self_num_kmers_coverage = coverage;

        // Initially, clade == self
        node_stats.clade_num_kmers_distinct = distinct;
        node_stats.clade_num_kmers_coverage = coverage;

        // DB placeholders
        node_stats.self_kmers_db = 0;
        node_stats.clade_kmers_db = 0;
        node_stats.self_tax_kmers_db = 0;
        node_stats.clade_tax_kmers_db = 0;

        stats_map.insert(taxid, node_stats);
    }

    stats_map
}

/// Recursively sum children’s stats into the parent (reads, distinct kmers, coverage, etc.)
pub fn accumulate_clade_stats(
    taxid: u32,
    children_map: &AHashMap<u32, Vec<u32>>,
    stats_map: &mut AHashMap<u32, NodeStats>,
) -> (u64, usize, u64, u64, u64) {
    // Ensure an entry exists
    if !stats_map.contains_key(&taxid) {
        stats_map.insert(taxid, NodeStats::default());
    }
    let cur = stats_map[&taxid].clone();

    // Start with self
    let mut total_reads = cur.self_reads;
    let mut total_distinct = cur.self_num_kmers_distinct;
    let mut total_coverage = cur.self_num_kmers_coverage;
    let mut total_db = cur.self_kmers_db;
    let mut total_tax_db = cur.self_tax_kmers_db;

    // Recurse children
    if let Some(kids) = children_map.get(&taxid) {
        for &child in kids {
            let (c_reads, c_distinct, c_cov, c_db, c_tax_db) =
                accumulate_clade_stats(child, children_map, stats_map);

            total_reads += c_reads;
            total_distinct += c_distinct;
            total_coverage += c_cov;
            total_db += c_db;
            total_tax_db += c_tax_db;
        }
    }

    // Update parent's clade fields
    if let Some(node) = stats_map.get_mut(&taxid) {
        node.clade_reads = total_reads;
        node.clade_num_kmers_distinct = total_distinct;
        node.clade_num_kmers_coverage = total_coverage;
        node.clade_kmers_db = total_db;
        node.clade_tax_kmers_db = total_tax_db;
    }

    (total_reads, total_distinct, total_coverage, total_db, total_tax_db)
}

/// Generate a Kraken-style report (both structured rows and text).
pub fn generate_kraken_style_report(
    root_taxid: u32,
    stats_map: &AHashMap<u32, NodeStats>,
    children_map: &AHashMap<u32, Vec<u32>>,
    name_map: Option<&NameMap>,
    rank_map: Option<&RankMap>,
    total_reads: u64,
) -> (Vec<KrakenReportRow>, String) {
    let mut report_lines = Vec::new();
    let mut report_text = String::new();

    // Minimal header
    report_text.push_str("%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName\n");

    fn dfs(
        taxid: u32,
        parent_taxid: Option<u32>,
        depth: usize,
        stats_map: &AHashMap<u32, NodeStats>,
        children_map: &AHashMap<u32, Vec<u32>>,
        name_map: Option<&NameMap>,
        rank_map: Option<&RankMap>,
        total_reads: u64,
        report_lines: &mut Vec<KrakenReportRow>,
        report_text: &mut String,
    ) {
        let default_stats = NodeStats::default();
        let stats = stats_map.get(&taxid).unwrap_or(&default_stats);

        if stats.clade_reads == 0 || total_reads == 0 {
            return;
        }

        let pct = 100.0 * (stats.clade_reads as f64) / (total_reads as f64);
        // duplication = coverage / distinct
        let dup = stats.clade_duplication();

        // coverage = fraction of DB k-mers (if used)
        let cov = if stats.clade_kmers_db > 0 {
            stats.clade_tax_kmers_db as f64 / stats.clade_kmers_db as f64
        } else {
            0.0
        };

        // Possibly retrieve rank & name
        let rank_str = rank_map
            .and_then(|rm| rm.get(&taxid))
            .cloned()
            .unwrap_or_default();
        let raw_name = name_map
            .and_then(|nm| nm.get(&taxid))
            .cloned()
            .unwrap_or_default();

        // Indent name based on depth
        let mut indented_name = String::new();
        for _ in 0..depth {
            indented_name.push('\t');
        }
        indented_name.push_str(&raw_name);

        // Sort children by clade_reads desc
        let mut kids = children_map.get(&taxid).cloned().unwrap_or_default();
        kids.sort_by_key(|child_id| {
            let cstats = stats_map.get(child_id).unwrap_or(&default_stats);
            std::cmp::Reverse(cstats.clade_reads)
        });

        // We’ll store “kmers” in the final report as “clade_num_kmers_distinct”
        let clade_kmers = stats.clade_num_kmers_distinct as u64;

        // Build structured row
        let row = KrakenReportRow {
            pct,
            reads: stats.clade_reads,
            tax_reads: stats.self_reads,
            kmers: clade_kmers,
            dup,
            cov,
            tax_id: taxid,
            rank: rank_str.clone(),
            tax_name: raw_name.clone(),
            depth,
            parent_tax_id: parent_taxid,
            children_tax_ids: kids.clone(),
        };
        report_lines.push(row);

        // Write text line
        use std::fmt::Write as _;
        let _ = write!(
            report_text,
            "{:.4}\t{}\t{}\t{}\t{:.4}\t{:.5}\t{}\t{}\t{}\n",
            pct,
            stats.clade_reads,
            stats.self_reads,
            clade_kmers,
            dup,
            cov,
            taxid,
            rank_str,
            indented_name
        );

        // Recurse on children
        for child in kids {
            dfs(
                child,
                Some(taxid),
                depth + 1,
                stats_map,
                children_map,
                name_map,
                rank_map,
                total_reads,
                report_lines,
                report_text,
            );
        }
    }

    dfs(
        root_taxid,
        None,
        0,
        stats_map,
        children_map,
        name_map,
        rank_map,
        total_reads,
        &mut report_lines,
        &mut report_text,
    );

    (report_lines, report_text)
}

/// The main pipeline to build a Kraken-style report:
///  1) Build children map
///  2) init_node_stats
///  3) accumulate_clade_stats
///  4) generate_kraken_style_report
pub fn build_kraken_report(
    taxon_counts: &TaxonCounts,
    parent_map: &ParentMap,
    name_map: Option<&NameMap>,
    rank_map: Option<&RankMap>,
    root_taxid: u32,
    total_reads: u64,
) -> (Vec<KrakenReportRow>, String) {
    let children_map = build_children_map(parent_map);
    let mut stats_map = init_node_stats(taxon_counts);
    accumulate_clade_stats(root_taxid, &children_map, &mut stats_map);
    generate_kraken_style_report(
        root_taxid,
        &stats_map,
        &children_map,
        name_map,
        rank_map,
        total_reads,
    )
}
