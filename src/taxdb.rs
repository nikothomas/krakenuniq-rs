// src/taxdb.rs

use ahash::AHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind};
use std::path::Path;

/// A type alias for your parent, name, and rank maps
pub type ParentMap = AHashMap<u32, u32>;
pub type NameMap = AHashMap<u32, String>;
pub type RankMap = AHashMap<u32, String>;

/// This struct holds all data after reading the merged tax DB + counts
pub struct TaxonomyData {
    pub parent_map: ParentMap,
    pub name_map: NameMap,
    pub rank_map: RankMap,

    /// The total counts after propagating up the tree
    pub total_counts: AHashMap<u32, u32>,
    /// Direct (per-node only) counts from the counts file
    pub direct_counts: AHashMap<u32, u32>,
}

#[derive(Debug)]
struct TaxNode {
    parent_id: u32,
    taxname: String,
    rank: String,
    count: u32, // direct count from counts file (if any)
}

/// Reads a combined "taxonomy DB + rank + name" file plus a "counts" file,
/// merges them into a single `TaxonomyData`.
///
/// - `taxdb_path` lines each have: `taxid\tparentid\ttaxname\trank`
/// - `counts_path` lines each have: `taxid\tcount`
pub fn read_taxdb_and_counts<P: AsRef<Path>, Q: AsRef<Path>>(
    taxdb_path: P,
    counts_path: Q,
) -> Result<TaxonomyData, Box<dyn std::error::Error>> {
    // 1) Read the taxdb
    let file = File::open(&taxdb_path)?;
    let reader = BufReader::new(file);

    let mut tax_tree: AHashMap<u32, TaxNode> = AHashMap::new();
    let mut ordered_tax_ids = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 {
            continue;
        }

        let taxid_str = parts[0].trim();
        let parentid_str = parts[1].trim();
        let name_str = parts[2].trim();
        let rank_str = parts[3].trim();

        let tax_id: u32 = match taxid_str.parse() {
            Ok(n) => n,
            Err(_) => continue,
        };
        let parent_id: u32 = parentid_str.parse().unwrap_or(0);

        // Insert into our structure
        ordered_tax_ids.push(tax_id);
        tax_tree.insert(
            tax_id,
            TaxNode {
                parent_id,
                taxname: name_str.to_string(),
                rank: rank_str.to_string(),
                count: 0, // fill from the counts file next
            },
        );
    }

    // 2) Read the counts file
    let file_counts = File::open(&counts_path)?;
    let reader_counts = BufReader::new(file_counts);
    let mut direct_counts = AHashMap::new();

    for line_res in reader_counts.lines() {
        let line = line_res?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            continue;
        }

        let tax_id: u32 = fields[0].parse()?;
        let count: u32 = fields[1].parse()?;

        // Update the node if present
        if let Some(node) = tax_tree.get_mut(&tax_id) {
            node.count = count;
        }
        direct_counts.insert(tax_id, count);
    }

    // 3) Compute total_counts by propagating up the tree
    let mut total_counts = AHashMap::new();
    // init with direct count
    for &tax_id in &ordered_tax_ids {
        let c = tax_tree.get(&tax_id).map(|n| n.count).unwrap_or(0);
        total_counts.insert(tax_id, c);
    }

    // propagate
    for &tax_id in &ordered_tax_ids {
        if let Some(node) = tax_tree.get(&tax_id) {
            let mut parent = node.parent_id;
            let count = total_counts.get(&tax_id).copied().unwrap_or(0);
            if parent != tax_id {
                while let Some(parent_node) = tax_tree.get(&parent) {
                    *total_counts.entry(parent).or_insert(0) += count;
                    if parent == parent_node.parent_id {
                        // root or cycle
                        break;
                    }
                    parent = parent_node.parent_id;
                }
            }
        }
    }

    // 4) Build parent_map, name_map, rank_map
    let mut parent_map = AHashMap::default();
    let mut name_map = AHashMap::default();
    let mut rank_map = AHashMap::default();

    for (&tax_id, node) in &tax_tree {
        parent_map.insert(tax_id, node.parent_id);
        name_map.insert(tax_id, node.taxname.clone());
        rank_map.insert(tax_id, node.rank.clone());
    }

    // 5) Return everything
    Ok(TaxonomyData {
        parent_map,
        name_map,
        rank_map,
        total_counts,
        direct_counts,
    })
}
