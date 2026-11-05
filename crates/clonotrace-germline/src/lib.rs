//! Experimental CIGAR-driven reference projection. No alignment or D inference.
pub mod v2;
use std::collections::BTreeMap;
/// Only exact unique FASTA IDs are resolved. Source sequences remain external.
pub fn fasta(text: &str) -> Result<BTreeMap<String, Vec<u8>>, String> {
    let mut out = BTreeMap::new();
    let mut current: Option<String> = None;
    for line in text.lines() {
        if let Some(id) = line.strip_prefix('>') {
            if id.trim().is_empty()
                || id.bytes().any(|b| b.is_ascii_whitespace())
                || out.contains_key(id)
            {
                return Err("blank/duplicate/noncanonical reference ID".into());
            }
            current = Some(id.into());
            out.insert(id.to_owned(), Vec::new());
        } else {
            let id = current.as_ref().ok_or("sequence before FASTA header")?;
            if line.is_empty() || !line.bytes().all(|b| b"ACGT".contains(&b)) {
                return Err("reference requires uppercase ACGT".into());
            }
            out.get_mut(id).unwrap().extend_from_slice(line.as_bytes());
        }
    }
    if out.is_empty() || out.values().any(Vec::is_empty) {
        return Err("empty reference".into());
    }
    Ok(out)
}
fn operations(cigar: &str) -> Result<Vec<(usize, u8)>, &'static str> {
    let bytes = cigar.as_bytes();
    let mut ops = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        let start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if start == i || bytes[start] == b'0' || i == bytes.len() || !b"SMNX=".contains(&bytes[i]) {
            return Err("unsupported_cigar");
        }
        let n: usize = cigar[start..i].parse().map_err(|_| "unsupported_cigar")?;
        ops.push((n, bytes[i]));
        i += 1;
    }
    if ops.is_empty() {
        return Err("unsupported_cigar");
    }
    // AIRR clipping/skipping belongs at alignment ends, around one contiguous aligned block.
    let first = ops
        .iter()
        .position(|&(_, op)| b"M=X".contains(&op))
        .ok_or("no_aligned_segment")?;
    let last = ops
        .iter()
        .rposition(|&(_, op)| b"M=X".contains(&op))
        .unwrap();
    if ops[first..=last]
        .iter()
        .any(|&(_, op)| !b"M=X".contains(&op))
    {
        return Err("unsupported_cigar");
    }
    for end in [&ops[..first], &ops[last + 1..]] {
        let codes: Vec<u8> = end.iter().map(|&(_, op)| op).collect();
        if !matches!(codes.as_slice(), [] | [b'S'] | [b'N'] | [b'S', b'N']) {
            return Err("unsupported_cigar");
        }
    }
    Ok(ops)
}
pub fn project(
    sequence: &[u8],
    v_call: &str,
    j_call: &str,
    v_cigar: &str,
    j_cigar: &str,
    refs: &BTreeMap<String, Vec<u8>>,
) -> Result<Vec<u8>, &'static str> {
    if sequence.is_empty() || !sequence.iter().all(|b| b"ACGT".contains(b)) {
        return Err("invalid_sequence");
    }
    let v = refs.get(v_call).ok_or("missing_reference")?;
    let j = refs.get(j_call).ok_or("missing_reference")?;
    let mut out = vec![b'N'; sequence.len()];
    let mut used = vec![false; sequence.len()];
    for (reference, cigar) in [(v, v_cigar), (j, j_cigar)] {
        let mut query = 0usize;
        let mut target = 0usize;
        let mut aligned = 0usize;
        for (n, op) in operations(cigar)? {
            match op {
                b'S' => {
                    query = query.checked_add(n).ok_or("query_bounds")?;
                    if query > out.len() {
                        return Err("query_bounds");
                    }
                }
                b'N' => {
                    target = target.checked_add(n).ok_or("reference_bounds")?;
                    if target > reference.len() {
                        return Err("reference_bounds");
                    }
                }
                _ => {
                    let end = query.checked_add(n).ok_or("query_bounds")?;
                    let rend = target.checked_add(n).ok_or("reference_bounds")?;
                    if end > out.len() {
                        return Err("query_bounds");
                    }
                    if rend > reference.len() {
                        return Err("reference_bounds");
                    }
                    if used[query..end].iter().any(|x| *x) {
                        return Err("overlapping_segments");
                    }
                    out[query..end].copy_from_slice(&reference[target..rend]);
                    used[query..end].fill(true);
                    query = end;
                    target = rend;
                    aligned += n;
                }
            }
        }
        if aligned == 0 {
            return Err("no_aligned_segment");
        }
    }
    Ok(out)
}
