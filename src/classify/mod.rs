pub mod classify_reads;
pub mod classify_sequence;
pub mod classify_stats;

use classify_sequence::*;
use super::krakendb::KrakenDB;
use super::taxdb::{ParentMap, RankMap, NameMap};