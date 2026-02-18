use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use memmap2::Mmap;

use super::InputError;

const HEADER_SIZE_V1: usize = 256;
const MAGIC: &[u8; 4] = b"KORG";
const VERSION_MAJOR: u16 = 1;
const VERSION_MINOR: u16 = 0;
const ENDIAN_TAG: u32 = 0x1234_5678;
const CRC64_ECMA_POLY: u64 = 0x42F0_E1EB_A9EA_3693;

#[derive(Debug, Clone)]
pub struct SharedCacheData {
    pub path: std::path::PathBuf,
    pub n_genes: u64,
    pub n_cells: u64,
    pub nnz: u64,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
    pub col_ptr: Vec<u64>,
    pub row_idx: Vec<u32>,
    pub values_u32: Vec<u32>,
    pub mmap: Arc<Mmap>,
}

#[derive(Debug, Clone, Copy)]
struct HeaderV1 {
    n_genes: u64,
    n_cells: u64,
    nnz: u64,
    genes_table_offset: u64,
    genes_table_bytes: u64,
    barcodes_table_offset: u64,
    barcodes_table_bytes: u64,
    col_ptr_offset: u64,
    row_idx_offset: u64,
    values_u32_offset: u64,
    n_blocks: u64,
    blocks_offset: u64,
    file_bytes: u64,
    header_crc64: u64,
    data_crc64: u64,
}

pub fn read_shared_cache(path: &Path) -> Result<SharedCacheData, InputError> {
    let file = File::open(path)?;
    let mmap = {
        // SAFETY: The file is opened read-only and the mapping is used immutably.
        unsafe { Mmap::map(&file) }
            .map_err(|e| InputError::Parse(format!("failed to mmap cache file: {e}")))?
    };

    if mmap.len() < HEADER_SIZE_V1 {
        return Err(InputError::Parse(
            "cache file is smaller than header".to_string(),
        ));
    }

    let header = parse_header(&mmap[..HEADER_SIZE_V1])?;

    if header.file_bytes != mmap.len() as u64 {
        return Err(InputError::Parse(format!(
            "file_bytes mismatch: header={}, actual={}",
            header.file_bytes,
            mmap.len()
        )));
    }

    validate_header_crc(&mmap[..HEADER_SIZE_V1], header.header_crc64)?;
    validate_section_bounds(&header, mmap.len() as u64)?;

    let genes = parse_string_table(
        &mmap,
        header.genes_table_offset,
        header.genes_table_bytes,
        header.n_genes,
        "genes",
    )?;
    let barcodes = parse_string_table(
        &mmap,
        header.barcodes_table_offset,
        header.barcodes_table_bytes,
        header.n_cells,
        "barcodes",
    )?;

    let col_ptr = parse_u64_array(
        &mmap,
        header.col_ptr_offset,
        (header.n_cells + 1)
            .checked_mul(8)
            .ok_or_else(|| InputError::Parse("col_ptr byte length overflow".to_string()))?,
    )?;
    let row_idx = parse_u32_array(
        &mmap,
        header.row_idx_offset,
        header
            .nnz
            .checked_mul(4)
            .ok_or_else(|| InputError::Parse("row_idx byte length overflow".to_string()))?,
    )?;
    let values_u32 = parse_u32_array(
        &mmap,
        header.values_u32_offset,
        header
            .nnz
            .checked_mul(4)
            .ok_or_else(|| InputError::Parse("values_u32 byte length overflow".to_string()))?,
    )?;

    validate_csc(
        &col_ptr,
        &row_idx,
        header.n_genes,
        header.n_cells,
        header.nnz,
    )?;

    Ok(SharedCacheData {
        path: path.to_path_buf(),
        n_genes: header.n_genes,
        n_cells: header.n_cells,
        nnz: header.nnz,
        genes,
        barcodes,
        col_ptr,
        row_idx,
        values_u32,
        mmap: Arc::new(mmap),
    })
}

fn parse_header(header: &[u8]) -> Result<HeaderV1, InputError> {
    if &header[0..4] != MAGIC {
        return Err(InputError::Parse("cache magic is invalid".to_string()));
    }

    let version_major = read_u16_le(header, 4)?;
    let version_minor = read_u16_le(header, 6)?;
    if version_major != VERSION_MAJOR {
        return Err(InputError::Parse(format!(
            "unsupported cache major version: {version_major}"
        )));
    }
    if version_minor != VERSION_MINOR {
        return Err(InputError::Parse(format!(
            "unsupported cache minor version: {version_minor}"
        )));
    }

    let endian_tag = read_u32_le(header, 8)?;
    if endian_tag != ENDIAN_TAG {
        return Err(InputError::Parse("cache endian tag is invalid".to_string()));
    }

    let header_size = read_u32_le(header, 12)?;
    if header_size as usize != HEADER_SIZE_V1 {
        return Err(InputError::Parse(format!(
            "cache header_size invalid: {header_size}"
        )));
    }

    let parsed = HeaderV1 {
        n_genes: read_u64_le(header, 16)?,
        n_cells: read_u64_le(header, 24)?,
        nnz: read_u64_le(header, 32)?,
        genes_table_offset: read_u64_le(header, 40)?,
        genes_table_bytes: read_u64_le(header, 48)?,
        barcodes_table_offset: read_u64_le(header, 56)?,
        barcodes_table_bytes: read_u64_le(header, 64)?,
        col_ptr_offset: read_u64_le(header, 72)?,
        row_idx_offset: read_u64_le(header, 80)?,
        values_u32_offset: read_u64_le(header, 88)?,
        n_blocks: read_u64_le(header, 96)?,
        blocks_offset: read_u64_le(header, 104)?,
        file_bytes: read_u64_le(header, 112)?,
        header_crc64: read_u64_le(header, 120)?,
        data_crc64: read_u64_le(header, 128)?,
    };

    if parsed.n_blocks != 0 {
        return Err(InputError::Parse("n_blocks must be 0 in v1".to_string()));
    }
    if parsed.blocks_offset != 0 {
        return Err(InputError::Parse(
            "blocks_offset must be 0 in v1".to_string(),
        ));
    }
    if parsed.data_crc64 != 0 {
        return Err(InputError::Parse("data_crc64 must be 0 in v1".to_string()));
    }

    Ok(parsed)
}

fn validate_header_crc(header: &[u8], expected_crc: u64) -> Result<(), InputError> {
    let mut header_for_crc = header.to_vec();
    header_for_crc[120..128].fill(0);
    let actual_crc = crc64_ecma(&header_for_crc);
    if actual_crc != expected_crc {
        return Err(InputError::Parse(format!(
            "cache header CRC mismatch: expected={expected_crc}, actual={actual_crc}"
        )));
    }
    Ok(())
}

fn validate_section_bounds(header: &HeaderV1, file_bytes: u64) -> Result<(), InputError> {
    let sections = [
        (
            "genes_table",
            header.genes_table_offset,
            header.genes_table_bytes,
        ),
        (
            "barcodes_table",
            header.barcodes_table_offset,
            header.barcodes_table_bytes,
        ),
        (
            "col_ptr",
            header.col_ptr_offset,
            (header.n_cells + 1)
                .checked_mul(8)
                .ok_or_else(|| InputError::Parse("col_ptr bytes overflow".to_string()))?,
        ),
        (
            "row_idx",
            header.row_idx_offset,
            header
                .nnz
                .checked_mul(4)
                .ok_or_else(|| InputError::Parse("row_idx bytes overflow".to_string()))?,
        ),
        (
            "values_u32",
            header.values_u32_offset,
            header
                .nnz
                .checked_mul(4)
                .ok_or_else(|| InputError::Parse("values_u32 bytes overflow".to_string()))?,
        ),
    ];

    for (name, offset, bytes) in sections {
        let end = offset
            .checked_add(bytes)
            .ok_or_else(|| InputError::Parse(format!("{name} section overflow")))?;
        if end > file_bytes {
            return Err(InputError::Parse(format!(
                "{name} section out of bounds: end={end}, file={file_bytes}"
            )));
        }
    }

    Ok(())
}

fn parse_string_table(
    file_bytes: &[u8],
    offset: u64,
    bytes_len: u64,
    expected_count: u64,
    name: &str,
) -> Result<Vec<String>, InputError> {
    let table = slice_at(file_bytes, offset, bytes_len)?;

    if table.len() < 4 {
        return Err(InputError::Parse(format!("{name} table too small")));
    }

    let count = read_u32_le(table, 0)? as usize;
    if count as u64 != expected_count {
        return Err(InputError::Parse(format!(
            "{name} table count mismatch: table={count}, expected={expected_count}"
        )));
    }

    let offsets_bytes =
        4usize
            .checked_add((count + 1).checked_mul(4).ok_or_else(|| {
                InputError::Parse(format!("{name} table offsets length overflow"))
            })?)
            .ok_or_else(|| InputError::Parse(format!("{name} table layout overflow")))?;

    if table.len() < offsets_bytes {
        return Err(InputError::Parse(format!(
            "{name} table too small for offsets"
        )));
    }

    let mut offsets: Vec<u32> = Vec::with_capacity(count + 1);
    for i in 0..=count {
        offsets.push(read_u32_le(table, 4 + i * 4)?);
    }

    for i in 1..offsets.len() {
        if offsets[i] < offsets[i - 1] {
            return Err(InputError::Parse(format!(
                "{name} table offsets are not monotonic"
            )));
        }
    }

    let blob = &table[offsets_bytes..];
    let blob_len = blob.len() as u32;
    if offsets[count] != blob_len {
        return Err(InputError::Parse(format!(
            "{name} table terminal offset mismatch"
        )));
    }

    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let start = offsets[i] as usize;
        let end = offsets[i + 1] as usize;
        let s = std::str::from_utf8(&blob[start..end])
            .map_err(|_| InputError::Parse(format!("{name} table contains invalid UTF-8")))?;
        out.push(s.to_string());
    }

    Ok(out)
}

fn validate_csc(
    col_ptr: &[u64],
    row_idx: &[u32],
    n_genes: u64,
    n_cells: u64,
    nnz: u64,
) -> Result<(), InputError> {
    if col_ptr.len() != n_cells as usize + 1 {
        return Err(InputError::Parse("col_ptr length is invalid".to_string()));
    }
    if row_idx.len() != nnz as usize {
        return Err(InputError::Parse("row_idx length is invalid".to_string()));
    }
    if col_ptr.first().copied().unwrap_or(1) != 0 {
        return Err(InputError::Parse("col_ptr[0] must be 0".to_string()));
    }
    if col_ptr.last().copied().unwrap_or(u64::MAX) != nnz {
        return Err(InputError::Parse(
            "col_ptr[n_cells] must equal nnz".to_string(),
        ));
    }

    for i in 1..col_ptr.len() {
        if col_ptr[i] < col_ptr[i - 1] {
            return Err(InputError::Parse("col_ptr must be monotonic".to_string()));
        }
    }

    for c in 0..n_cells as usize {
        let start = col_ptr[c] as usize;
        let end = col_ptr[c + 1] as usize;
        let mut prev: Option<u32> = None;
        for &row in &row_idx[start..end] {
            if row as u64 >= n_genes {
                return Err(InputError::Parse("row_idx out of bounds".to_string()));
            }
            if let Some(p) = prev {
                if row <= p {
                    return Err(InputError::Parse(
                        "row_idx must be strictly increasing per column".to_string(),
                    ));
                }
            }
            prev = Some(row);
        }
    }

    Ok(())
}

fn parse_u64_array(file_bytes: &[u8], offset: u64, bytes_len: u64) -> Result<Vec<u64>, InputError> {
    let raw = slice_at(file_bytes, offset, bytes_len)?;
    if raw.len() % 8 != 0 {
        return Err(InputError::Parse(
            "u64 array has invalid byte length".to_string(),
        ));
    }
    let mut out = Vec::with_capacity(raw.len() / 8);
    for i in 0..out.capacity() {
        out.push(read_u64_le(raw, i * 8)?);
    }
    Ok(out)
}

fn parse_u32_array(file_bytes: &[u8], offset: u64, bytes_len: u64) -> Result<Vec<u32>, InputError> {
    let raw = slice_at(file_bytes, offset, bytes_len)?;
    if raw.len() % 4 != 0 {
        return Err(InputError::Parse(
            "u32 array has invalid byte length".to_string(),
        ));
    }
    let mut out = Vec::with_capacity(raw.len() / 4);
    for i in 0..out.capacity() {
        out.push(read_u32_le(raw, i * 4)?);
    }
    Ok(out)
}

fn slice_at(file_bytes: &[u8], offset: u64, bytes_len: u64) -> Result<&[u8], InputError> {
    let start = usize::try_from(offset)
        .map_err(|_| InputError::Parse("offset does not fit usize".to_string()))?;
    let len = usize::try_from(bytes_len)
        .map_err(|_| InputError::Parse("bytes length does not fit usize".to_string()))?;
    let end = start
        .checked_add(len)
        .ok_or_else(|| InputError::Parse("slice bounds overflow".to_string()))?;
    file_bytes
        .get(start..end)
        .ok_or_else(|| InputError::Parse("slice out of bounds".to_string()))
}

fn read_u16_le(bytes: &[u8], offset: usize) -> Result<u16, InputError> {
    let arr: [u8; 2] = bytes
        .get(offset..offset + 2)
        .ok_or_else(|| InputError::Parse("read_u16 out of bounds".to_string()))?
        .try_into()
        .map_err(|_| InputError::Parse("read_u16 failed".to_string()))?;
    Ok(u16::from_le_bytes(arr))
}

fn read_u32_le(bytes: &[u8], offset: usize) -> Result<u32, InputError> {
    let arr: [u8; 4] = bytes
        .get(offset..offset + 4)
        .ok_or_else(|| InputError::Parse("read_u32 out of bounds".to_string()))?
        .try_into()
        .map_err(|_| InputError::Parse("read_u32 failed".to_string()))?;
    Ok(u32::from_le_bytes(arr))
}

fn read_u64_le(bytes: &[u8], offset: usize) -> Result<u64, InputError> {
    let arr: [u8; 8] = bytes
        .get(offset..offset + 8)
        .ok_or_else(|| InputError::Parse("read_u64 out of bounds".to_string()))?
        .try_into()
        .map_err(|_| InputError::Parse("read_u64 failed".to_string()))?;
    Ok(u64::from_le_bytes(arr))
}

fn crc64_ecma(bytes: &[u8]) -> u64 {
    let mut crc = 0u64;
    for &b in bytes {
        crc ^= (b as u64) << 56;
        for _ in 0..8 {
            if (crc & 0x8000_0000_0000_0000) != 0 {
                crc = (crc << 1) ^ CRC64_ECMA_POLY;
            } else {
                crc <<= 1;
            }
        }
    }
    crc
}
