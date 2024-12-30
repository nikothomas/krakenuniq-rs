use ahash::AHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::error::Error;
use std::path::Path;

#[derive(Debug)]
struct TaxNode {
    parent_id: u32,
    count: u32,
}

pub fn process_taxonomy<P: AsRef<Path>, Q: AsRef<Path>>(
    taxdb_path: P,
    counts_path: Q,
) -> Result<(AHashMap<u32, u32>, AHashMap<u32, u32>), Box<dyn Error>> {
    // Read taxonomy database
    let tax_file = File::open(&taxdb_path)?;
    let reader = BufReader::new(tax_file);
    let mut tax_tree: AHashMap<u32, TaxNode> = AHashMap::new();
    let mut ordered_tax_ids: Vec<u32> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            continue;
        }

        let tax_id: u32 = fields[0].parse()?;
        let parent_id: u32 = fields[1].parse()?;

        ordered_tax_ids.push(tax_id);
        tax_tree.insert(
            tax_id,
            TaxNode {
                parent_id,
                count: 0,
            },
        );
    }

    // Read count file and initialize counts
    let count_file = File::open(&counts_path)?;
    let reader = BufReader::new(count_file);
    let mut direct_counts: AHashMap<u32, u32> = AHashMap::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 {
            continue;
        }

        let tax_id: u32 = fields[0].parse()?;
        let count: u32 = fields[1].parse()?;

        if let Some(node) = tax_tree.get_mut(&tax_id) {
            node.count = count;
            direct_counts.insert(tax_id, count);
        }
    }

    // Calculate total counts
    let mut total_counts: AHashMap<u32, u32> = AHashMap::new();

    // First pass: Initialize total_counts with direct counts
    for &tax_id in &ordered_tax_ids {
        if let Some(node) = tax_tree.get(&tax_id) {
            total_counts.insert(tax_id, node.count);
        }
    }

    // Second pass: Propagate counts up the tree
    for &tax_id in &ordered_tax_ids {
        if let Some(node) = tax_tree.get(&tax_id) {
            let mut parent_id = node.parent_id;
            let count = total_counts.get(&tax_id).copied().unwrap_or(0);

            // Skip self-referential nodes
            if parent_id != tax_id {
                while let Some(parent_node) = tax_tree.get(&parent_id) {
                    let current_count = total_counts.entry(parent_id).or_insert(0);
                    *current_count += count;

                    if parent_id == parent_node.parent_id {
                        break; // Reached root or cycle
                    }
                    parent_id = parent_node.parent_id;
                }
            }
        }
    }

    Ok((total_counts, direct_counts))
}
