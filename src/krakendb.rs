//src/krakendb.rs

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Error as IoError, ErrorKind, Read, Seek, SeekFrom};
use std::ops::Range;
use std::str;
use std::sync::Arc;
use ahash::AHashMap;

/// File type code for Jellyfish/Kraken DBs.
pub const DATABASE_FILE_TYPE: &str = "JFLISTDN";

/// File type code on Kraken DB index.
pub const KRAKEN_INDEX_STRING: &str = "KRAKIDX";

/// File type code for Kraken DB index (v2).
pub const KRAKEN_INDEX2_STRING: &str = "KRAKIX2";

/// XOR mask for minimizer bin keys (allows for better distribution).
pub const INDEX2_XOR_MASK: u64 = 0xe37e28c4271b5a2d;

/// A fully safe, slice-based representation of a Kraken DB index.
pub struct KrakenDBIndex {
    /// The raw index data.
    /// In C++ this was `char *fptr`; here it’s a safe slice.
    pub index_data: Arc<[u8]>,
    /// Index version (1 => `KRAKIDX`, 2 => `KRAKIX2`).
    pub idx_type: u8,
    /// Number of indexed nucleotides in the minimizer.
    pub nt: u8,

    /// A parsed array of offsets (the “bin starts”).
    /// Equivalent to `get_array()` in the original C++ code.
    offset_array: Arc<[u64]>,
}

/// A fully safe, slice-based representation of the Kraken DB file.
#[derive(Clone)]
pub struct KrakenDB {
    /// Entire DB file loaded in memory.
    db_data: Arc<[u8]>,
    /// Number of k-mer/taxon pairs in the DB.
    pub key_ct: u64,
    /// Byte length of each taxon value (e.g. 4).
    pub val_len: u64,
    /// Byte length of each k-mer representation.
    pub key_len: u64,
    /// Number of bits used to represent a k-mer.
    pub key_bits: u64,
    /// The “k” value (k-mer length in nucleotides).
    pub k: u8,

    /// The total file size (for reference).
    pub filesize: u64,

    /// Optional index pointer (similar to `set_index` in C++).
    pub index_ptr: Option<Arc<KrakenDBIndex>>,

    // ---------- CHUNKING FIELDS ----------
    /// The number of chunks, if chunked loading is used.
    pub chunks: u32,
    /// Minimizer boundaries in the index space.
    pub idx_chunk_bounds: Vec<u64>,
    /// K-mer/taxon pair boundaries in the DB space.
    pub dbx_chunk_bounds: Vec<u64>,
    /// The offset inside the DB slice used during chunk queries.
    pub data_offset: u64,
}

impl KrakenDB {
    /// Creates an empty `KrakenDB`.
    pub fn new() -> Self {
        Self {
            db_data: Arc::new([]),
            key_ct: 0,
            val_len: 0,
            key_len: 0,
            key_bits: 0,
            k: 0,
            filesize: 0,
            index_ptr: None,

            chunks: 0,
            idx_chunk_bounds: Vec::new(),
            dbx_chunk_bounds: Vec::new(),
            data_offset: 0,
        }
    }

    /// Opens a DB file, reads it fully into memory, validates, and parses header fields.
    pub fn open_file(mut self, path: &str) -> Result<Self, IoError> {
        let mut file = File::open(path)?;
        let filesize = file.metadata()?.len();

        if filesize < 8 {
            return Err(IoError::new(
                ErrorKind::InvalidData,
                "File too small to be a valid Kraken DB",
            ));
        }

        // Read the entire file into a buffer
        let mut buf = vec![0u8; filesize as usize];
        file.seek(SeekFrom::Start(0))?;
        file.read_exact(&mut buf)?;

        // Wrap the buffer in Arc<[u8]>
        let db_data = Arc::<[u8]>::from(buf);

        // Check DB signature
        let sig_len = DATABASE_FILE_TYPE.len();
        if filesize < sig_len as u64 {
            return Err(IoError::new(
                ErrorKind::InvalidData,
                "File truncated, cannot read DB signature",
            ));
        }
        let file_type_bytes = &db_data[..sig_len];
        let file_type_str = match std::str::from_utf8(file_type_bytes) {
            Ok(s) => s,
            Err(_) => {
                return Err(IoError::new(
                    ErrorKind::InvalidData,
                    "DB signature not valid UTF-8",
                ))
            }
        };
        if file_type_str != DATABASE_FILE_TYPE {
            return Err(IoError::new(
                ErrorKind::InvalidData,
                format!("DB in improper format. Found {file_type_str}"),
            ));
        }

        // Parse metadata from known offsets
        let key_bits = read_u64_le(&db_data[8..16]);
        let val_len = read_u64_le(&db_data[16..24]);
        let key_ct = read_u64_le(&db_data[48..56]);

        if val_len != 4 {
            return Err(IoError::new(
                ErrorKind::InvalidData,
                "Only 4-byte DB values supported in this example",
            ));
        }
        let k = (key_bits / 2) as u8;
        let key_len = (key_bits / 8) + u64::from(key_bits % 8 != 0);

        log::info!(
            "Loaded database with {} keys, k={}, val_len={}, key_len={}",
            key_ct,
            k,
            val_len,
            key_len
        );

        self.db_data = db_data;
        self.key_bits = key_bits;
        self.val_len = val_len;
        self.key_ct = key_ct;
        self.k = k;
        self.key_len = key_len;
        self.filesize = filesize;

        Ok(self)
    }

    /// Returns the size of the DB header, from which `(kmer, taxon)` pairs start.
    /// Matches the original C++ logic: `72 + 2 * (4 + 8 * key_bits)`
    pub fn header_size(&self) -> usize {
        (72 + 2 * (4 + 8 * self.key_bits)) as usize
    }

    /// Returns the “pair size” = `key_len + val_len` in bytes.
    pub fn pair_size(&self) -> usize {
        (self.key_len + self.val_len) as usize
    }

    /// Returns a slice for all `(kmer, taxon)` pairs (i.e., from header end to the file end).
    fn pairs_slice(&self) -> &[u8] {
        let start = self.header_size();
        if start >= self.db_data.len() {
            return &[];
        }
        &self.db_data[start..]
    }

    // -----------------------------------------------------------------------
    //  Counting Taxons
    // -----------------------------------------------------------------------

    /// Counts how many times each taxon ID occurs across all k-mers,
    /// akin to `count_taxons()` in the original code.
    pub fn count_taxons(&self) -> Result<BTreeMap<u32, u64>, IoError> {
        let mut taxon_counts = BTreeMap::new();

        let pairs = self.pairs_slice();
        let pair_sz = self.pair_size();
        let key_ct = self.key_ct as usize;
        let total_bytes_needed = key_ct.checked_mul(pair_sz)
            .ok_or_else(|| IoError::new(ErrorKind::InvalidData, "Overflow in key_ct * pair_sz"))?;

        if total_bytes_needed > pairs.len() {
            return Err(IoError::new(
                ErrorKind::InvalidData,
                format!("File truncated? Need {total_bytes_needed} bytes in pair data."),
            ));
        }

        for i in 0..key_ct {
            if i % 10_000_000 == 1 {
                eprint!("\r{:.2}%", 100.0 * (i as f64 / key_ct as f64));
            }
            // The taxon offset is the second half of each pair
            let offset = i * pair_sz + self.key_len as usize;
            let taxon_slice = &pairs[offset..(offset + self.val_len as usize)];
            let taxon_val = u32::from_le_bytes(
                taxon_slice
                    .try_into()
                    .map_err(|_| IoError::new(ErrorKind::InvalidData, "Bad taxon slice"))?,
            );
            *taxon_counts.entry(taxon_val).or_insert(0) += 1;
        }
        eprint!("\r");
        Ok(taxon_counts)
    }

    // -----------------------------------------------------------------------
    //  K-mer / Minimizer Manipulations
    // -----------------------------------------------------------------------

    /// Reverse complement of a k-mer, for `n` nucleotides.
    /// Mirrors the Jellyfish-based logic from the original code.
    pub fn reverse_complement_n(&self, mut kmer: u64, n: u8) -> u64 {
        // The bit manipulations from Jellyfish (C++).
        // This is a pure function, safe to do in Rust.
        kmer = ((kmer >> 2) & 0x3333333333333333) | ((kmer & 0x3333333333333333) << 2);
        kmer = ((kmer >> 4) & 0x0F0F0F0F0F0F0F0F) | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);
        kmer = ((kmer >> 8) & 0x00FF00FF00FF00FF) | ((kmer & 0x00FF00FF00FF00FF) << 8);
        kmer = ((kmer >> 16) & 0x0000FFFF0000FFFF) | ((kmer & 0x0000FFFF0000FFFF) << 16);
        kmer = (kmer >> 32) | (kmer << 32);
        (((u64::MAX) - kmer) >> (8 * std::mem::size_of::<u64>() - ((n as usize) << 1))) as u64
    }

    /// Reverse complement using the DB’s k.
    pub fn reverse_complement(&self, kmer: u64) -> u64 {
        self.reverse_complement_n(kmer, self.k)
    }

    /// Lexicographically smaller of `(kmer, revcom(kmer))`, for `n`.
    pub fn canonical_representation_n(&self, kmer: u64, n: u8) -> u64 {
        let rc = self.reverse_complement_n(kmer, n);
        if kmer < rc {
            kmer
        } else {
            rc
        }
    }

    /// Lexicographically smaller of `(kmer, revcom(kmer))`, using the DB’s k.
    pub fn canonical_representation(&self, kmer: u64) -> u64 {
        let rc = self.reverse_complement(kmer);
        if kmer < rc {
            kmer
        } else {
            rc
        }
    }

    /// Compute the bin key for the given `kmer` with length = `idx_nt`.
    /// Mirrors `bin_key(uint64_t kmer, uint64_t idx_nt)` in C++.
    pub fn bin_key_n(&self, mut kmer: u64, idx_nt: u64) -> u64 {
        let use_idx_type = if let Some(ref idx) = self.index_ptr {
            idx.idx_type
        } else {
            1
        };
        let xor_mask = if use_idx_type == 2 { INDEX2_XOR_MASK } else { 0 };
        let mask = (1u64 << (idx_nt * 2)) - 1;
        let xor_mask = xor_mask & mask;

        let mut min_bin_key = u64::MAX;
        let iterations = self.key_bits / 2 - idx_nt + 1;
        for _ in 0..iterations {
            let rep = self.canonical_representation_n(kmer & mask, idx_nt as u8);
            let tmp = xor_mask ^ rep;
            if tmp < min_bin_key {
                min_bin_key = tmp;
            }
            kmer >>= 2;
        }
        min_bin_key
    }

    /// Shortcut bin key that picks `nt` from the DB index or falls back to `self.k`.
    pub fn bin_key(&self, kmer: u64) -> u64 {
        if let Some(ref index) = self.index_ptr {
            self.bin_key_n(kmer, index.nt as u64)
        } else {
            self.bin_key_n(kmer, self.k as u64)
        }
    }

    // -----------------------------------------------------------------------
    //  Query logic using the index
    // -----------------------------------------------------------------------

    /// Simplified version of `kmer_query`:
    ///  1) compute `b_key`
    ///  2) get bin range from index
    ///  3) do a binary search in `[start, end)` among pairs
    pub fn kmer_query(&self, kmer: u64) -> Option<u32> {
        let b_key = self.bin_key(kmer);
        let idx = self.index_ptr.as_ref()?;

        let range = idx.bin_range(b_key)?;
        self.binary_search_in_range(kmer, range)
    }

    /// Binary search for `kmer` in `[start, end)` among `(kmer, taxon)` pairs.
    fn binary_search_in_range(&self, kmer: u64, range: Range<u64>) -> Option<u32> {
        let pairs = self.pairs_slice();
        let pair_sz = self.pair_size();

        let start = range.start as usize;
        let end = range.end as usize;
        let total_bytes = end.checked_mul(pair_sz)?;

        // Sanity check
        if total_bytes > pairs.len() {
            // The file or index might be corrupt/truncated
            return None;
        }

        // We'll do a standard binary search
        let mut left = start as i64;
        let mut right = end as i64 - 1;

        while left <= right {
            let mid = (left + right) / 2;
            let offset = mid as usize * pair_sz;

            // Read kmer from the first `key_len` bytes
            let comp_kmer = read_kmer_bits(&pairs[offset..offset + self.key_len as usize], self.key_bits);
            if kmer > comp_kmer {
                left = mid + 1;
            } else if kmer < comp_kmer {
                right = mid - 1;
            } else {
                // Found match => read the taxon
                let taxon_offset = offset + self.key_len as usize;
                let taxon_slice = &pairs[taxon_offset..taxon_offset + self.val_len as usize];
                let taxon = u32::from_le_bytes(taxon_slice.try_into().ok()?);
                return Some(taxon);
            }
        }

        None
    }

    // -----------------------------------------------------------------------
    //  Chunking logic (highly simplified)
    // -----------------------------------------------------------------------

    /// Demonstrates a naive approach to chunking: we simply specify some boundaries.
    /// In real code, you’d do a more sophisticated approach to limit memory usage, etc.
    pub fn prepare_chunking(&mut self, max_bytes_for_db: usize) {
        // In a real system, you'd chunk the index into segments.
        // For demonstration, we only store chunk boundaries.
        // We do not allocate or copy partial data to `self.data`.
        self.chunks = 2;
        self.idx_chunk_bounds = vec![0, 42, 84];
        self.dbx_chunk_bounds = vec![0, 100, 200];

        log::info!(
            "prepare_chunking with max bytes {} => stub with 2 chunks",
            max_bytes_for_db
        );
    }

    /// Load a chunk, if you needed to do partial I/O or partial mapping.
    /// Here we just stub it out.
    pub fn load_chunk(&mut self, db_chunk_id: u32) {
        log::info!("Stub: load_chunk {db_chunk_id}");
        // Real code would copy partial data from `pairs_slice()` into a smaller buffer.
    }

    /// Check if `minimizer` is within chunk `db_chunk_id` (based on `idx_chunk_bounds`).
    pub fn is_minimizer_in_chunk(&self, minimizer: u64, db_chunk_id: u32) -> bool {
        if (db_chunk_id as usize) >= self.idx_chunk_bounds.len() - 1 {
            return false;
        }
        let start = self.idx_chunk_bounds[db_chunk_id as usize];
        let end = self.idx_chunk_bounds[db_chunk_id as usize + 1];
        (start..end).contains(&minimizer)
    }
}

impl KrakenDBIndex {
    /// Build a new index from an Arc slice.
    /// In the original code, we’d parse pointers in `from_mmap`.
    pub fn from_slice(index_data: Arc<[u8]>) -> Result<Self, IoError> {
        // We must check if this is KRAKEN_INDEX_STRING or KRAKEN_INDEX2_STRING, etc.
        let s1_len = KRAKEN_INDEX_STRING.len();
        let s2_len = KRAKEN_INDEX2_STRING.len();
        let max_len = s1_len.max(s2_len);

        if index_data.len() < max_len {
            return Err(IoError::new(
                ErrorKind::InvalidData,
                "Index file too short",
            ));
        }

        let prefix1 = &index_data[..s1_len];
        let prefix2 = &index_data[..s2_len];

        let idx_type = if prefix1 == KRAKEN_INDEX_STRING.as_bytes() {
            1
        } else if prefix2 == KRAKEN_INDEX2_STRING.as_bytes() {
            2
        } else {
            return Err(IoError::new(
                ErrorKind::InvalidData,
                "Illegal Kraken DB index format",
            ));
        };

        // The next byte is `nt`
        let nt_offset = if idx_type == 1 { s1_len } else { s2_len };
        if nt_offset + 1 > index_data.len() {
            return Err(IoError::new(
                ErrorKind::InvalidData,
                "Index truncated; cannot read nt",
            ));
        }
        let nt = index_data[nt_offset];

        // The rest of the index is a u64 array of offsets
        let arr_offset = nt_offset + 1;
        if arr_offset > index_data.len() {
            return Err(IoError::new(
                ErrorKind::InvalidData,
                "Index truncated; cannot read offset array",
            ));
        }
        // Convert the remainder to a vector of u64
        if (index_data.len() - arr_offset) % 8 != 0 {
            return Err(IoError::new(
                ErrorKind::InvalidData,
                "Index array must be a multiple of 8 bytes",
            ));
        }

        let mut tmp_vec = Vec::new();
        let mut i = arr_offset;
        while i < index_data.len() {
            let next64 = read_u64_le(&index_data[i..i + 8]);
            tmp_vec.push(next64);
            i += 8;
        }

        let offset_array = Arc::<[u64]>::from(tmp_vec);

        Ok(Self {
            index_data,
            idx_type,
            nt,
            offset_array,
        })
    }

    /// Returns how many nucleotides the index uses for bin keys.
    pub fn indexed_nt(&self) -> u8 {
        self.nt
    }

    /// Returns `[start, end)` for a given bin key (`b_key`).
    /// Equivalent to `at(b_key)` and `at(b_key+1)` in the C++ code.
    pub fn bin_range(&self, b_key: u64) -> Option<Range<u64>> {
        let start = self.at(b_key)?;
        let end = self.at(b_key + 1)?;
        Some(start..end)
    }

    /// Safe accessor to the offset array (like `KrakenDBIndex::at(idx)`).
    /// This is equivalent to `fptr + offset` in the original code.
    pub fn at(&self, idx: u64) -> Option<u64> {
        self.offset_array.get(idx as usize).copied()
    }
}

// ---------------------------------------------------------------------------
//  Helper functions
// ---------------------------------------------------------------------------

/// Reads a `u64` from a little-endian byte slice.
fn read_u64_le(bytes: &[u8]) -> u64 {
    let mut arr = [0u8; 8];
    arr.copy_from_slice(bytes);
    u64::from_le_bytes(arr)
}

/// Reads a kmer from the front of `bytes` up to `key_len` bytes,
/// then mask it by `(1 << key_bits) - 1` to handle any leftover bits.
fn read_kmer_bits(bytes: &[u8], key_bits: u64) -> u64 {
    // Convert up to 8 bytes
    let mut arr = [0u8; 8];
    for (i, b) in bytes.iter().enumerate().take(8) {
        arr[i] = *b;
    }
    let val = u64::from_le_bytes(arr);

    let mask = (1u64 << key_bits) - 1;
    val & mask
}

pub fn read_db_tax_kmer_counts(path: &str) -> Result<AHashMap<u32, u64>, Box<dyn std::error::Error>> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);

    let mut result = AHashMap::new();
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<_> = line.trim().split('\t').collect();
        if parts.len() != 2 {
            continue; // skip malformed lines
        }
        let taxid: u32 = parts[0].parse().unwrap_or(0);
        let count: u64 = parts[1].parse().unwrap_or(0);
        if taxid != 0 && count > 0 {
            result.insert(taxid, count);
        }
    }
    Ok(result)
}

/// Return `true` if the given raw `kmer` is assigned exactly to `taxid` in the Kraken DB.
/// For example, if `db.kmer_query(...)` returns `Some(t) == taxid`, we say it "belongs."
///
/// If you need to consider ancestor/descendant relationships, you would do more complex logic
/// (e.g., climb parents to see if `taxid` is an ancestor).
pub fn kmer_belongs_to_taxon(db: &KrakenDB, raw_kmer: u64, taxid: u32) -> bool {
    // 1) Convert to canonical form
    let canon_kmer = db.canonical_representation(raw_kmer);

    // 2) Query the DB
    if let Some(hit_taxid) = db.kmer_query(canon_kmer) {
        // 3) Check if exact match
        return hit_taxid == taxid;
    }
    false
}