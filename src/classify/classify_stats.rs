// src/classify/classify_stats.rs

use ahash::AHashMap;
use crate::types::KrakenReportRow;
use super::classify_sequence::TaxonCounts;
use super::{ParentMap, NameMap, RankMap};

/// Per-node stats that store:
///   - read counts (self & clade)
///   - distinct k-mer counts & total coverage (self & clade)
///   - placeholders for DB-based metrics:
///       * `clade_kmers_db`      = total k-mers in the DB for this clade
///       * `clade_tax_kmers_db`  = how many of those DB k-mers we actually observed in the reads
///
/// Therefore, coverage = clade_tax_kmers_db / clade_kmers_db
#[derive(Default, Debug, Clone)]
pub struct NodeStats {
    /// Reads directly assigned to this node
    pub self_reads: u32,
    /// Reads in node + descendants
    pub clade_reads: u32,

    /// How many **distinct** k-mers in this node alone
    pub self_num_kmers_distinct: usize,
    /// Total coverage for this node's k-mers (sum of counts)
    pub self_num_kmers_coverage: u32,

    /// Distinct k-mers in clade (node + descendants)
    pub clade_num_kmers_distinct: usize,
    /// Coverage in clade (node + descendants)
    pub clade_num_kmers_coverage: u32,

    // -------------------------
    // DB-based metrics
    // -------------------------
    /// Number of DB k-mers assigned specifically to this node
    pub self_kmers_db: u32,
    /// Number of DB k-mers assigned to the entire clade (node + descendants)
    pub clade_kmers_db: u32,
    /// Number of node-level k-mers (in DB) we actually found in the reads
    pub self_tax_kmers_db: u32,
    /// Number of clade-level k-mers (in DB) we actually found in the reads
    pub clade_tax_kmers_db: u32,
}


impl NodeStats {
    /// Return duplication factor for this clade = coverage / distinct in the reads, i.e.
    /// `(clade_num_kmers_coverage) / (clade_num_kmers_distinct)`.
    pub fn clade_duplication(&self) -> f32 {
        if self.clade_num_kmers_distinct == 0 {
            0.0
        } else {
            self.clade_num_kmers_coverage as f32
                / self.clade_num_kmers_distinct as f32
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

pub fn init_node_stats(
    taxon_counts: &TaxonCounts,
    total_db_counts: &AHashMap<u32, u32>,
    direct_db_counts: &AHashMap<u32, u32>,
    parent_map: &ParentMap,
) -> AHashMap<u32, NodeStats> {
    let mut stats_map = AHashMap::new();

    // First, initialize all nodes in the taxonomy with their DB counts
    for (&taxid, _) in parent_map {
        let mut node_stats = NodeStats::default();

        // Set DB counts even if no reads
        if let Some(&direct_count) = direct_db_counts.get(&taxid) {
            node_stats.self_kmers_db = direct_count;
        }
        if let Some(&total_count) = total_db_counts.get(&taxid) {
            node_stats.clade_kmers_db = total_count;
        }

        stats_map.insert(taxid, node_stats);
    }

    // Then update nodes that have read counts
    for (&taxid, read_counts) in taxon_counts {
        let node_stats = stats_map.entry(taxid).or_default();

        // Read-based stats
        let distinct = read_counts.kmer_counts.len();
        let coverage = read_counts.kmer_counts.values().map(|&v| v as u32).sum();

        node_stats.self_reads = read_counts.read_count;
        node_stats.clade_reads = read_counts.read_count;
        node_stats.self_num_kmers_distinct = distinct;
        node_stats.self_num_kmers_coverage = coverage;
        node_stats.clade_num_kmers_distinct = distinct;
        node_stats.clade_num_kmers_coverage = coverage;

        // Observed DB k-mers
        node_stats.self_tax_kmers_db = distinct as u32;
        node_stats.clade_tax_kmers_db = distinct as u32;
    }

    stats_map
}

/// Recursively sum only the read-based statistics (not DB counts which are pre-calculated)
pub fn accumulate_clade_stats(
    taxid: u32,
    children_map: &AHashMap<u32, Vec<u32>>,
    stats_map: &mut AHashMap<u32, NodeStats>,
) -> (u32, usize, u32, u32) {
    // Ensure an entry exists
    if !stats_map.contains_key(&taxid) {
        stats_map.insert(taxid, NodeStats::default());
    }
    let cur = stats_map[&taxid].clone();

    // Start with self
    let mut total_reads = cur.self_reads;
    let mut total_distinct = cur.self_num_kmers_distinct;
    let mut total_coverage = cur.self_num_kmers_coverage;
    let mut total_tax_db = cur.self_tax_kmers_db;

    // Recurse children
    if let Some(kids) = children_map.get(&taxid) {
        for &child in kids {
            let (c_reads, c_distinct, c_cov, c_tax_db) =
                accumulate_clade_stats(child, children_map, stats_map);

            total_reads += c_reads;
            total_distinct += c_distinct;
            total_coverage += c_cov;
            total_tax_db += c_tax_db;
        }
    }

    // Update parent's clade fields
    if let Some(node) = stats_map.get_mut(&taxid) {
        node.clade_reads = total_reads;
        node.clade_num_kmers_distinct = total_distinct;
        node.clade_num_kmers_coverage = total_coverage;
        node.clade_tax_kmers_db = total_tax_db;
    }

    (total_reads, total_distinct, total_coverage, total_tax_db)
}

/// Build a Kraken-style report, but we no longer re-read the DB counts:
/// instead, we take `&total_db_counts` and `&direct_db_counts` from the caller.
pub fn build_kraken_report(
    taxon_counts: &TaxonCounts,
    parent_map: &ParentMap,
    name_map: Option<&NameMap>,
    rank_map: Option<&RankMap>,
    root_taxid: u32,
    total_reads: u32,
    // <-- Added these two:
    total_db_counts: &AHashMap<u32, u32>,
    direct_db_counts: &AHashMap<u32, u32>,
) -> (Vec<KrakenReportRow>, String) {
    // Build children map once
    let children_map = build_children_map(parent_map);

    // Use total_db_counts, direct_db_counts from the caller
    let mut stats_map = init_node_stats(
        taxon_counts,
        total_db_counts,
        direct_db_counts,
        parent_map,
    );

    // Accumulate child data
    accumulate_clade_stats(root_taxid, &children_map, &mut stats_map);

    // Generate final structured rows & text
    generate_kraken_style_report(
        root_taxid,
        &stats_map,
        &children_map,
        name_map,
        rank_map,
        total_reads,
    )
}

/// Generate a Kraken-style report (both structured rows and text).
pub fn generate_kraken_style_report(
    root_taxid: u32,
    stats_map: &AHashMap<u32, NodeStats>,
    children_map: &AHashMap<u32, Vec<u32>>,
    name_map: Option<&NameMap>,
    rank_map: Option<&RankMap>,
    total_reads: u32,
) -> (Vec<KrakenReportRow>, String) {
    let mut report_lines = Vec::new();
    let mut report_text = String::new();

    // Minimal header
    // Note that `cov` = coverage of k-mers in this clade (fraction in DB actually observed).
    report_text.push_str("%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName\n");

    fn dfs(
        taxid: u32,
        parent_taxid: Option<u32>,
        depth: usize,
        stats_map: &AHashMap<u32, NodeStats>,
        children_map: &AHashMap<u32, Vec<u32>>,
        name_map: Option<&NameMap>,
        rank_map: Option<&RankMap>,
        total_reads: u32,
        report_lines: &mut Vec<KrakenReportRow>,
        report_text: &mut String,
    ) {
        let default_stats = NodeStats::default();
        let stats = stats_map.get(&taxid).unwrap_or(&default_stats);

        // If no reads in clade => skip
        if stats.clade_reads == 0 && stats.self_reads == 0 {
            return;
        }

        // % of reads in entire dataset
        let pct = 100.0 * (stats.clade_reads as f32) / (total_reads as f32);
        let dup = stats.clade_duplication();

        // -- coverage = fraction of DB k-mers for this clade that were found in reads --
        // i.e. clade_tax_kmers_db / clade_kmers_db
        let cov = if stats.clade_kmers_db > 0 {
            stats.clade_tax_kmers_db as f32 / stats.clade_kmers_db as f32
        } else {
            0.0 as f32
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

        // Indent name based on tree depth in the final text
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

        let clade_kmers = stats.clade_num_kmers_distinct as u32;

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
            "{:.4}\t{}\t{}\t{}\t{:.4}\t{:.6}\t{}\t{}\t{}\n",
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

    // DFS starting at root
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
