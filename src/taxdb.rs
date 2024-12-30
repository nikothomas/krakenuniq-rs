//src/taxdb.rs

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

pub type ParentMap = HashMap<u32, u32>;
pub type NameMap = HashMap<u32, String>;
pub type RankMap = HashMap<u32, String>;

/// Parses a taxDB file in the format:
/// ```text
/// <taxid>\t<parentid>\t<taxname>\t<rank>
/// ```
/// Returns:
/// - a `ParentMap` mapping child_taxid -> parent_taxid
/// - a `NameMap` mapping taxid -> taxname
/// - a `RankMap` mapping taxid -> rank
pub fn parse_taxdb<P: AsRef<Path>>(
    filepath: P,
) -> io::Result<(ParentMap, NameMap, RankMap)> {
    let file = File::open(filepath)?;
    let reader = BufReader::new(file);

    let mut parent_map: ParentMap = HashMap::new();
    let mut name_map: NameMap = HashMap::new();
    let mut rank_map: RankMap = HashMap::new();

    for line_result in reader.lines() {
        let line = line_result?;
        // Expecting 4 tab-separated fields: taxid, parentid, taxname, rank
        // e.g. "2   1   Eukaryota   domain"
        let parts: Vec<&str> = line.split('\t').collect();

        // Skip malformed lines
        if parts.len() < 4 {
            continue;
        }

        let taxid_str = parts[0].trim();
        let parentid_str = parts[1].trim();
        let taxname_str = parts[2].trim();
        let rank_str = parts[3].trim();

        // Attempt to parse numeric fields
        let taxid: u32 = taxid_str.parse().unwrap_or(0);
        let parentid: u32 = parentid_str.parse().unwrap_or(0);

        // Store in maps
        if taxid != 0 {
            parent_map.insert(taxid, parentid);
            name_map.insert(taxid, taxname_str.to_string());
            rank_map.insert(taxid, rank_str.to_string());
        }
    }
    Ok((parent_map, name_map, rank_map))
}