// src/classification_stats.rs

use std::collections::HashMap;
use crate::types::KrakenReportRow;
use super::classify_sequence::{TaxonCounts};
use super::taxdb::{ParentMap, NameMap, RankMap};

/// Holds per-node stats:
///   - `self_reads`: how many reads directly assigned to this node
///   - `clade_reads`: how many reads in total (node + descendants)
///   - `self_kmers`: how many *unique* k-mers directly assigned to this node
///   - `clade_kmers`: how many *unique* k-mers in node + descendants
///
/// Additionally, we show placeholders for "database" metrics:
///   - `self_kmers_db` / `clade_kmers_db`: for how many k-mers matched in some external DB
///   - `self_tax_kmers_db` / `clade_tax_kmers_db`: for how many *taxonomic* k-mers matched in DB
#[derive(Default, Debug)]
pub struct NodeStats {
    pub self_reads: u64,
    pub clade_reads: u64,
    pub self_kmers: usize,
    pub clade_kmers: usize,

    // Placeholders for DB metrics (not actually computed by default):
    pub self_kmers_db: u64,
    pub clade_kmers_db: u64,
    pub self_tax_kmers_db: u64,
    pub clade_tax_kmers_db: u64,
}

/// Build a map of parent -> Vec<child> from the existing child->parent `parent_map`.
/// This helps us traverse the taxonomy tree from the root(s) downward.
pub fn build_children_map(parent_map: &ParentMap) -> HashMap<u32, Vec<u32>> {
    let mut children_map: HashMap<u32, Vec<u32>> = HashMap::new();

    // Initialize every known taxid as a key to avoid missing entries
    for &taxid in parent_map.keys() {
        children_map.entry(taxid).or_default();
    }

    // Populate children
    for (&child, &parent) in parent_map {
        // skip if parent=0 or child=parent (in case of weird data)
        if parent > 0 && child != parent {
            children_map.entry(parent).or_default().push(child);
        }
    }
    children_map
}

/// Initialize `NodeStats` for each taxid in `taxon_counts`.
/// The `self_reads` and `self_kmers` come from classification results.
/// `clade_reads` and `clade_kmers` start off equal to the self values
/// but will be updated once we accumulate child stats.
///
/// For the “DB” fields (`_kmers_db`), this example initializes them to 0;
/// if you have logic to populate them, add it here.
pub fn init_node_stats(taxon_counts: &TaxonCounts) -> HashMap<u32, NodeStats> {
    let mut stats_map: HashMap<u32, NodeStats> = HashMap::new();

    for (&taxid, rc) in taxon_counts {
        let mut node_stats = NodeStats::default();
        node_stats.self_reads = rc.read_count;
        node_stats.clade_reads = rc.read_count;
        node_stats.self_kmers = rc.unique_kmers.len();
        node_stats.clade_kmers = node_stats.self_kmers;

        // If you actually track “kmersDB” or “taxKmersDB”, populate them here.
        node_stats.self_kmers_db = 0;
        node_stats.clade_kmers_db = 0;
        node_stats.self_tax_kmers_db = 0;
        node_stats.clade_tax_kmers_db = 0;

        stats_map.insert(taxid, node_stats);
    }

    stats_map
}

/// Recursively sum children's stats into the parent's `clade_reads`, `clade_kmers`,
/// and the DB fields (`clade_kmers_db`, etc.) if used.
pub fn accumulate_clade_stats(
    taxid: u32,
    children_map: &HashMap<u32, Vec<u32>>,
    stats_map: &mut HashMap<u32, NodeStats>,
) -> (u64, usize, u64, u64) {
    if !stats_map.contains_key(&taxid) {
        stats_map.insert(taxid, NodeStats::default());
    }

    let mut total_reads = stats_map[&taxid].clade_reads;
    let mut total_kmers = stats_map[&taxid].clade_kmers;
    let mut total_kmers_db = stats_map[&taxid].clade_kmers_db;
    let mut total_tax_kmers_db = stats_map[&taxid].clade_tax_kmers_db;

    // Recurse into children
    if let Some(children) = children_map.get(&taxid) {
        for &child_taxid in children {
            let (child_clade_reads, child_clade_kmers, child_kmers_db, child_tax_kmers_db) =
                accumulate_clade_stats(child_taxid, children_map, stats_map);
            total_reads += child_clade_reads;
            total_kmers += child_clade_kmers;
            total_kmers_db += child_kmers_db;
            total_tax_kmers_db += child_tax_kmers_db;
        }
    }

    // Update our NodeStats with the accumulated totals
    if let Some(parent_stats) = stats_map.get_mut(&taxid) {
        parent_stats.clade_reads = total_reads;
        parent_stats.clade_kmers = total_kmers;
        parent_stats.clade_kmers_db = total_kmers_db;
        parent_stats.clade_tax_kmers_db = total_tax_kmers_db;
    }

    (total_reads, total_kmers, total_kmers_db, total_tax_kmers_db)
}

pub fn generate_kraken_style_report(
    root_taxid: u32,
    stats_map: &HashMap<u32, NodeStats>,
    children_map: &HashMap<u32, Vec<u32>>,
    name_map: Option<&NameMap>,
    rank_map: Option<&RankMap>,
    total_reads: u64,
) -> (Vec<KrakenReportRow>, String) {
    let mut report_lines = Vec::new();
    let mut report_text = String::new();

    // Header row in text
    report_text.push_str(
        "%\treads\ttaxReads\tkmers\ttaxKmers\tkmersDB\ttaxKmersDB\tdup\tcov\ttaxID\trank\ttaxName\n"
    );

    fn dfs(
        taxid: u32,
        depth: usize,
        stats_map: &HashMap<u32, NodeStats>,
        children_map: &HashMap<u32, Vec<u32>>,
        name_map: Option<&NameMap>,
        rank_map: Option<&RankMap>,
        total_reads: u64,
        report_lines: &mut Vec<KrakenReportRow>,
        report_text: &mut String,
    ) {
        let node_stats_default = &NodeStats::default();
        let stats = stats_map.get(&taxid).unwrap_or(&node_stats_default);
        if stats.clade_reads == 0 || total_reads == 0 {
            return;
        }

        let pct = 100.0 * (stats.clade_reads as f64) / (total_reads as f64);
        let dup = if stats.self_kmers == 0 {
            0.0
        } else {
            stats.clade_kmers as f64 / stats.self_kmers as f64
        };
        let cov = (stats.clade_kmers as f64) / 10_000.0;

        let rank_str = rank_map
            .and_then(|rm| rm.get(&taxid))
            .cloned()
            .unwrap_or_default();
        let raw_name = name_map
            .and_then(|nm| nm.get(&taxid))
            .cloned()
            .unwrap_or_default();

        // Indent with tabs
        let mut indented_name = String::new();
        for _ in 0..depth {
            indented_name.push('\t');
        }
        indented_name.push_str(&raw_name);

        // ----- Build structured row
        let row = KrakenReportRow {
            pct,
            reads: total_reads,
            tax_reads: stats.self_reads,
            kmers: stats.clade_kmers as u64,
            tax_kmers: stats.self_kmers as u64,
            kmers_db: stats.clade_kmers_db,
            tax_kmers_db: stats.self_kmers_db,
            dup,
            cov,
            tax_id: taxid,
            rank: rank_str.clone(),
            tax_name: raw_name.clone(),
            depth,
        };
        report_lines.push(row);

        // ----- Build the text line
        use std::fmt::Write as FmtWrite;
        let _ = write!(
            report_text,
            "{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.5}\t{}\t{}\t{}\n",
            pct,
            total_reads,
            stats.self_reads,
            stats.clade_kmers,
            stats.self_kmers,
            stats.clade_kmers_db,
            stats.self_kmers_db,
            dup,
            cov,
            taxid,
            rank_str,
            indented_name
        );

        // Descend into children
        let mut kids = children_map
            .get(&taxid)
            .cloned()
            .unwrap_or_default();
        kids.sort_by_key(|child_id| std::cmp::Reverse(stats_map.get(child_id).map_or(0, |c| c.clade_reads)));

        for child_id in kids {
            dfs(child_id, depth + 1, stats_map, children_map, name_map, rank_map,
                total_reads, report_lines, report_text);
        }
    }

    // Start DFS
    dfs(root_taxid, 0, stats_map, children_map, name_map, rank_map, total_reads, &mut report_lines, &mut report_text);

    (report_lines, report_text)
}

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
