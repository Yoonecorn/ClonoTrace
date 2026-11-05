//! Exact reference-derived V(D)J alignment reconstruction from supplied AIRR CIGARs.
//! The output is conditional on the supplied annotation/reference, not an inferred ancestor.
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DPolicy {
    Mask,
    Include,
}
impl DPolicy {
    pub fn name(self) -> &'static str {
        match self {
            Self::Mask => "mask-d",
            Self::Include => "include-d",
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct SegmentInput<'a> {
    pub call: &'a str,
    pub cigar: &'a str,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SegmentDiagnostic {
    pub segment: &'static str,
    pub call: String,
    pub cigar: String,
    /// Zero-based half-open query coordinates, including insertions.
    pub query_start: usize,
    pub query_end: usize,
    /// Zero-based half-open reference coordinates, including deletions.
    pub reference_start: usize,
    pub reference_end: usize,
    pub aligned_bases: usize,
    pub inserted_bases: usize,
    pub deleted_bases: usize,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Reconstruction {
    /// Exactly input query length. Insertions and unassigned positions are N.
    pub query_projection: Vec<u8>,
    /// Original query with '-' inserted for reference bases deleted from query.
    pub query_alignment: Vec<u8>,
    /// Reference bases, '-' for query insertions, N for unassigned query bases.
    pub germline_alignment: Vec<u8>,
    pub segments: Vec<SegmentDiagnostic>,
}
fn ops(cigar: &str) -> Result<Vec<(usize, u8)>, &'static str> {
    let bytes = cigar.as_bytes();
    let mut at = 0;
    let mut result = Vec::new();
    while at < bytes.len() {
        let start = at;
        while at < bytes.len() && bytes[at].is_ascii_digit() {
            at += 1;
        }
        if start == at || bytes[start] == b'0' || at == bytes.len() {
            return Err("unsupported_cigar");
        }
        let n: usize = cigar[start..at].parse().map_err(|_| "unsupported_cigar")?;
        let op = bytes[at];
        if !b"SMNX=ID".contains(&op) {
            return Err("unsupported_cigar");
        }
        result.push((n, op));
        at += 1;
    }
    let first = result
        .iter()
        .position(|&(_, o)| b"M=XID".contains(&o))
        .ok_or("unsupported_cigar")?;
    let last = result
        .iter()
        .rposition(|&(_, o)| b"M=XID".contains(&o))
        .ok_or("unsupported_cigar")?;
    for end in [&result[..first], &result[last + 1..]] {
        let codes: Vec<u8> = end.iter().map(|&(_, o)| o).collect();
        if !matches!(codes.as_slice(), [] | [b'S'] | [b'N'] | [b'S', b'N']) {
            return Err("unsupported_cigar");
        }
    }
    let body = &result[first..=last];
    if body.iter().any(|&(_, o)| !b"M=XID".contains(&o)) {
        return Err("unsupported_cigar");
    }
    if body.iter().any(|&(_, o)| o == b'M') && body.iter().any(|&(_, o)| b"=X".contains(&o)) {
        return Err("mixed_match_encoding");
    }
    Ok(result)
}
/// Reconstruct an alignment without inventing alleles, alignment positions or junction bases.
/// Mask deliberately ignores D annotation. Include requires D on IGH and forbids it on light chains.
/// Reference IDs are exact. Calls with commas or whitespace are not resolved by this API.
pub fn reconstruct(
    sequence: &[u8],
    locus: &str,
    v: SegmentInput<'_>,
    j: SegmentInput<'_>,
    d: Option<SegmentInput<'_>>,
    policy: DPolicy,
    references: &BTreeMap<String, Vec<u8>>,
) -> Result<Reconstruction, &'static str> {
    reconstruct_bounded(sequence, locus, v, j, d, policy, references, 20_000)
}
/// The explicit alignment bound is checked before allocating output alignment arrays.
#[allow(clippy::too_many_arguments)]
pub fn reconstruct_bounded(
    sequence: &[u8],
    locus: &str,
    v: SegmentInput<'_>,
    j: SegmentInput<'_>,
    d: Option<SegmentInput<'_>>,
    policy: DPolicy,
    references: &BTreeMap<String, Vec<u8>>,
    max_alignment_length: usize,
) -> Result<Reconstruction, &'static str> {
    if sequence.len() > max_alignment_length {
        return Err("max_alignment_length_exceeded");
    }
    if sequence.is_empty() || !sequence.iter().all(|b| b"ACGT".contains(b)) {
        return Err("invalid_sequence");
    }
    if !matches!(locus, "IGH" | "IGK" | "IGL") {
        return Err("outside_locus_scope");
    }
    if policy == DPolicy::Include {
        if locus == "IGH" && d.is_none() {
            return Err("missing_D_segment");
        }
        if locus != "IGH" && d.is_some() {
            return Err("D_not_supported_for_light_chain");
        }
    }
    let mut segments = vec![("V", v)];
    if policy == DPolicy::Include {
        if let Some(d) = d {
            segments.push(("D", d));
        }
    }
    segments.push(("J", j));
    let mut projected = vec![b'N'; sequence.len()];
    let mut aligned_reference = vec![None; sequence.len()];
    // Reference deletion events are anchored immediately before the query position.
    let mut deleted: BTreeMap<usize, (usize, Vec<u8>)> = BTreeMap::new();
    let mut diagnostics = Vec::new();
    let mut previous_end = 0;
    let mut projected_deletions = 0usize;
    for (segment_index, (segment, input)) in segments.into_iter().enumerate() {
        if input.call.is_empty()
            || input
                .call
                .bytes()
                .any(|b| b == b',' || b.is_ascii_whitespace())
        {
            return Err("invalid_reference_call");
        }
        let reference = references.get(input.call).ok_or("missing_reference")?;
        if reference.is_empty() || !reference.iter().all(|b| b"ACGT".contains(b)) {
            return Err("invalid_reference");
        }
        let operations = ops(input.cigar)?;
        let mut query = 0usize;
        let mut target = 0usize;
        let mut start = None;
        let mut reference_start = None;
        let mut end = 0;
        let mut reference_end = 0;
        let mut matches = 0usize;
        let mut insertions = 0usize;
        let mut deletions = 0usize;
        for (n, operation) in operations {
            match operation {
                b'S' => query = query.checked_add(n).ok_or("query_bounds")?,
                b'N' => target = target.checked_add(n).ok_or("reference_bounds")?,
                b'D' => {
                    let target_end = target.checked_add(n).ok_or("reference_bounds")?;
                    if target_end > reference.len() {
                        return Err("reference_bounds");
                    }
                    if query > sequence.len() {
                        return Err("query_bounds");
                    }
                    projected_deletions = projected_deletions
                        .checked_add(n)
                        .ok_or("alignment_size_overflow")?;
                    if sequence
                        .len()
                        .checked_add(projected_deletions)
                        .ok_or("alignment_size_overflow")?
                        > max_alignment_length
                    {
                        return Err("max_alignment_length_exceeded");
                    }
                    let event = deleted
                        .entry(query)
                        .or_insert_with(|| (segment_index, Vec::new()));
                    if event.0 != segment_index {
                        return Err("overlapping_deletion_anchors");
                    }
                    event.1.extend_from_slice(&reference[target..target_end]);
                    start.get_or_insert(query);
                    reference_start.get_or_insert(target);
                    target = target_end;
                    end = query;
                    reference_end = target;
                    deletions = deletions.checked_add(n).ok_or("reference_bounds")?;
                }
                b'I' => {
                    let query_end = query.checked_add(n).ok_or("query_bounds")?;
                    if query_end > sequence.len() {
                        return Err("query_bounds");
                    }
                    if aligned_reference[query..query_end]
                        .iter()
                        .any(Option::is_some)
                    {
                        return Err("overlapping_segments");
                    }
                    aligned_reference[query..query_end].fill(Some(b'-'));
                    start.get_or_insert(query);
                    reference_start.get_or_insert(target);
                    query = query_end;
                    end = query;
                    reference_end = target;
                    insertions = insertions.checked_add(n).ok_or("query_bounds")?;
                }
                _ => {
                    let query_end = query.checked_add(n).ok_or("query_bounds")?;
                    let target_end = target.checked_add(n).ok_or("reference_bounds")?;
                    if query_end > sequence.len() {
                        return Err("query_bounds");
                    }
                    if target_end > reference.len() {
                        return Err("reference_bounds");
                    }
                    if aligned_reference[query..query_end]
                        .iter()
                        .any(Option::is_some)
                    {
                        return Err("overlapping_segments");
                    }
                    for (q, r) in (query..query_end).zip(target..target_end) {
                        if operation == b'=' && sequence[q] != reference[r] {
                            return Err("cigar_match_conflict");
                        }
                        if operation == b'X' && sequence[q] == reference[r] {
                            return Err("cigar_mismatch_conflict");
                        }
                        projected[q] = reference[r];
                        aligned_reference[q] = Some(reference[r]);
                    }
                    start.get_or_insert(query);
                    reference_start.get_or_insert(target);
                    query = query_end;
                    target = target_end;
                    end = query;
                    reference_end = target;
                    matches = matches.checked_add(n).ok_or("query_bounds")?;
                }
            }
            if query > sequence.len() {
                return Err("query_bounds");
            }
            if target > reference.len() {
                return Err("reference_bounds");
            }
        }
        if matches == 0 {
            return Err("no_aligned_segment");
        }
        let start = start.ok_or("no_aligned_segment")?;
        if start < previous_end {
            return Err("segment_order");
        }
        previous_end = end;
        diagnostics.push(SegmentDiagnostic {
            segment,
            call: input.call.to_owned(),
            cigar: input.cigar.to_owned(),
            query_start: start,
            query_end: end,
            reference_start: reference_start.ok_or("no_aligned_segment")?,
            reference_end,
            aligned_bases: matches,
            inserted_bases: insertions,
            deleted_bases: deletions,
        });
    }
    let total_deleted = deleted
        .values()
        .try_fold(0usize, |a, (_, b)| a.checked_add(b.len()))
        .ok_or("alignment_size_overflow")?;
    let alignment_len = sequence
        .len()
        .checked_add(total_deleted)
        .ok_or("alignment_size_overflow")?;
    if alignment_len > max_alignment_length {
        return Err("max_alignment_length_exceeded");
    }
    let mut query_alignment = Vec::with_capacity(alignment_len);
    let mut germline_alignment = Vec::with_capacity(alignment_len);
    for q in 0..=sequence.len() {
        if let Some((_, bases)) = deleted.get(&q) {
            query_alignment.resize(query_alignment.len() + bases.len(), b'-');
            germline_alignment.extend_from_slice(bases);
        }
        if q < sequence.len() {
            query_alignment.push(sequence[q]);
            germline_alignment.push(aligned_reference[q].unwrap_or(b'N'));
        }
    }
    Ok(Reconstruction {
        query_projection: projected,
        query_alignment,
        germline_alignment,
        segments: diagnostics,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn refs() -> BTreeMap<String, Vec<u8>> {
        [
            ("IGHV1", "AACCGGTT"),
            ("IGHD1", "TGGA"),
            ("IGHJ1", "TTGCAA"),
        ]
        .into_iter()
        .map(|(a, b)| (a.into(), b.as_bytes().to_vec()))
        .collect()
    }
    fn run(sequence: &str, cigar: &str) -> Reconstruction {
        reconstruct(
            sequence.as_bytes(),
            "IGH",
            SegmentInput {
                call: "IGHV1",
                cigar,
            },
            SegmentInput {
                call: "IGHJ1",
                cigar: "10S2N3M1S",
            },
            None,
            DPolicy::Mask,
            &refs(),
        )
        .unwrap()
    }
    #[test]
    fn reference_indels_preserve_both_alignment_sides() {
        let ins = run("TTACACGTTTGCAT", "2S1N2M1I2M7S");
        assert_eq!(ins.query_projection, b"NNACNCGNNNGCAN");
        assert_eq!(ins.query_alignment, b"TTACACGTTTGCAT");
        assert_eq!(ins.germline_alignment, b"NNAC-CGNNNGCAN");
        let del = run("TTACGGTTTTGCAT", "2S1N2M1D2M8S");
        assert_eq!(del.query_projection, b"NNACGGNNNNGCAN");
        assert_eq!(del.query_alignment, b"TTAC-GGTTTTGCAT");
        assert_eq!(del.germline_alignment, b"NNACCGGNNNNGCAN");
    }
    #[test]
    fn explicit_d_and_mask_have_distinct_declared_meanings() {
        let v = SegmentInput {
            call: "IGHV1",
            cigar: "2S1N4M8S",
        };
        let j = SegmentInput {
            call: "IGHJ1",
            cigar: "10S2N3M1S",
        };
        let d = SegmentInput {
            call: "IGHD1",
            cigar: "7S1N2M5S",
        };
        let x = reconstruct(
            b"TTACCGTTTTGCAT",
            "IGH",
            v,
            j,
            Some(d),
            DPolicy::Include,
            &refs(),
        )
        .unwrap();
        assert_eq!(x.query_projection, b"NNACCGNGGNGCAN");
        assert_eq!(
            reconstruct(
                b"TTACCGTTTTGCAT",
                "IGH",
                v,
                j,
                None,
                DPolicy::Include,
                &refs()
            )
            .unwrap_err(),
            "missing_D_segment"
        );
        assert_eq!(
            reconstruct(
                b"TTACCGTTTTGCAT",
                "IGH",
                v,
                j,
                Some(SegmentInput {
                    call: "missing",
                    cigar: "bad"
                }),
                DPolicy::Mask,
                &refs()
            )
            .unwrap()
            .query_projection,
            b"NNACCGNNNNGCAN"
        );
    }
    #[test]
    fn malformed_public_inputs_return_errors() {
        let j = SegmentInput {
            call: "IGHJ1",
            cigar: "10S2N3M1S",
        };
        for c in [
            "",
            "0M",
            "1M1S1M",
            "1M1N1M",
            "1N1S1M",
            "999999999999999999999999999M",
            "1H3M",
            "1M1P1M",
        ] {
            assert!(
                reconstruct(
                    b"TTACCGTTTTGCAT",
                    "IGH",
                    SegmentInput {
                        call: "IGHV1",
                        cigar: c
                    },
                    j,
                    None,
                    DPolicy::Mask,
                    &refs()
                )
                .is_err(),
                "{c}"
            );
        }
        assert_eq!(
            reconstruct(
                b"AAAAAAAAAAAAAA",
                "IGH",
                SegmentInput {
                    call: "IGHV1",
                    cigar: "2S1N4=8S"
                },
                j,
                None,
                DPolicy::Mask,
                &refs()
            )
            .unwrap_err(),
            "cigar_match_conflict"
        );
    }
    #[test]
    fn alignment_limit_and_invalid_map_are_explicit_errors() {
        let v = SegmentInput {
            call: "IGHV1",
            cigar: "2S1N2M1D2M8S",
        };
        let j = SegmentInput {
            call: "IGHJ1",
            cigar: "10S2N3M1S",
        };
        assert_eq!(
            reconstruct_bounded(
                b"TTACGGTTTTGCAT",
                "IGH",
                v,
                j,
                None,
                DPolicy::Mask,
                &refs(),
                14
            )
            .unwrap_err(),
            "max_alignment_length_exceeded"
        );
        let mut bad = refs();
        bad.insert("IGHV1".into(), b"AANCCGTT".to_vec());
        assert_eq!(
            reconstruct(b"TTACGGTTTTGCAT", "IGH", v, j, None, DPolicy::Mask, &bad).unwrap_err(),
            "invalid_reference"
        );
        let mut huge = refs();
        huge.insert("IGHV1".into(), vec![b'A'; 100000]);
        let v = SegmentInput {
            call: "IGHV1",
            cigar: "1M99999D",
        };
        assert_eq!(
            reconstruct_bounded(
                b"AAAAAAAAAAAAAA",
                "IGH",
                v,
                j,
                None,
                DPolicy::Mask,
                &huge,
                20000
            )
            .unwrap_err(),
            "max_alignment_length_exceeded"
        );
    }
}
