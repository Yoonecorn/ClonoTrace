//! Checked canonical-amino-acid and standard-code translation APIs.
use std::collections::BTreeMap;
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AaSequence(Vec<u8>);
impl AaSequence {
    pub fn new(sequence: &[u8]) -> Result<Self, &'static str> {
        if sequence.is_empty() || !sequence.iter().all(|b| b"ACDEFGHIKLMNPQRSTVWY".contains(b)) {
            return Err("invalid_amino_acid_sequence");
        }
        Ok(Self(sequence.to_vec()))
    }
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }
    pub fn into_bytes(self) -> Vec<u8> {
        self.0
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AmbiguityPolicy {
    Reject,
    ResolveUnanimous,
}
#[derive(Clone, Copy, Debug)]
pub struct TranslationConfig {
    pub frame: usize,
    pub ambiguity: AmbiguityPolicy,
}
fn bases(b: u8, policy: AmbiguityPolicy) -> Option<&'static [u8]> {
    let exact: Option<&'static [u8]> = match b {
        b'A' => Some(b"A"),
        b'C' => Some(b"C"),
        b'G' => Some(b"G"),
        b'T' => Some(b"T"),
        _ => None,
    };
    if exact.is_some() || policy == AmbiguityPolicy::Reject {
        return exact;
    }
    match b {
        b'R' => Some(b"AG"),
        b'Y' => Some(b"CT"),
        b'S' => Some(b"CG"),
        b'W' => Some(b"AT"),
        b'K' => Some(b"GT"),
        b'M' => Some(b"AC"),
        b'B' => Some(b"CGT"),
        b'D' => Some(b"AGT"),
        b'H' => Some(b"ACT"),
        b'V' => Some(b"ACG"),
        b'N' => Some(b"ACGT"),
        _ => None,
    }
}
fn codon(a: u8, b: u8, c: u8) -> u8 {
    let base = |b| match b {
        b'T' => 0,
        b'C' => 1,
        b'A' => 2,
        b'G' => 3,
        _ => unreachable!(),
    };
    b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
        [base(a) * 16 + base(b) * 4 + base(c)]
}
pub fn translate(sequence: &[u8], config: TranslationConfig) -> Result<AaSequence, &'static str> {
    if sequence.is_empty()
        || sequence
            .iter()
            .any(|&b| bases(b, config.ambiguity).is_none())
    {
        return Err("invalid_nucleotide_sequence");
    }
    if config.frame > 2
        || sequence.len() <= config.frame
        || !(sequence.len() - config.frame).is_multiple_of(3)
    {
        return Err("empty_or_partial_codon");
    }
    let mut out = Vec::with_capacity((sequence.len() - config.frame) / 3);
    for c in sequence[config.frame..].as_chunks::<3>().0 {
        let mut result = None;
        for &a in bases(c[0], config.ambiguity).unwrap() {
            for &b in bases(c[1], config.ambiguity).unwrap() {
                for &z in bases(c[2], config.ambiguity).unwrap() {
                    let x = codon(a, b, z);
                    if result.is_some_and(|r| r != x) {
                        return Err("ambiguous_codon");
                    }
                    result = Some(x);
                }
            }
        }
        let x = result.unwrap();
        if x == b'*' {
            return Err("stop_codon");
        }
        out.push(x);
    }
    Ok(AaSequence(out))
}
pub fn normalized_hamming(a: &AaSequence, b: &AaSequence) -> Result<f64, &'static str> {
    if a.0.len() != b.0.len() {
        return Err("unequal_amino_acid_lengths");
    }
    Ok(a.0.iter().zip(&b.0).filter(|(x, y)| x != y).count() as f64 / a.0.len() as f64)
}
/// Exact identity; indices remain independent observations, including duplicates.
pub fn identity(sequences: &[AaSequence]) -> Vec<Vec<usize>> {
    let mut by_sequence = BTreeMap::new();
    for (i, s) in sequences.iter().enumerate() {
        by_sequence.entry(&s.0).or_insert_with(Vec::new).push(i);
    }
    let mut groups: Vec<Vec<usize>> = by_sequence.into_values().collect();
    groups.sort_by_key(|g| g[0]);
    groups
}
/// Separate length buckets, ordinary normalized-Hamming single linkage.
/// Empty batches yield an empty partition. Threshold validation never panics.
pub fn single_linkage(
    sequences: &[AaSequence],
    threshold: f64,
) -> Result<Vec<Vec<usize>>, &'static str> {
    if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
        return Err("threshold must be finite decimal in [0,1]");
    }
    let mut parent: Vec<usize> = (0..sequences.len()).collect();
    fn root(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            p[i] = p[p[i]];
            i = p[i];
        }
        i
    }
    for i in 0..sequences.len() {
        for j in i + 1..sequences.len() {
            if sequences[i].0.len() == sequences[j].0.len()
                && normalized_hamming(&sequences[i], &sequences[j])? <= threshold
            {
                let a = root(&mut parent, i);
                let b = root(&mut parent, j);
                parent[b] = a;
            }
        }
    }
    let mut groups = BTreeMap::new();
    for i in 0..sequences.len() {
        groups
            .entry(root(&mut parent, i))
            .or_insert_with(Vec::new)
            .push(i);
    }
    let mut groups: Vec<Vec<usize>> = groups.into_values().collect();
    groups.sort_by_key(|g| g[0]);
    Ok(groups)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn table_and_policies() {
        for (nt, aa) in [
            (b"GCTTGTAAA".as_slice(), b"ACK".as_slice()),
            (b"GCCTGCAAG", b"ACK"),
        ] {
            assert_eq!(
                translate(
                    nt,
                    TranslationConfig {
                        frame: 0,
                        ambiguity: AmbiguityPolicy::Reject
                    }
                )
                .unwrap()
                .as_bytes(),
                aa
            );
        }
        for (nt, aa) in [
            (b"GCN".as_slice(), b"A".as_slice()),
            (b"TGY", b"C"),
            (b"AAR", b"K"),
        ] {
            assert_eq!(
                translate(
                    nt,
                    TranslationConfig {
                        frame: 0,
                        ambiguity: AmbiguityPolicy::ResolveUnanimous
                    }
                )
                .unwrap()
                .as_bytes(),
                aa
            );
        }
        for nt in [b"TAN".as_slice(), b"NNN", b"TGR"] {
            assert_eq!(
                translate(
                    nt,
                    TranslationConfig {
                        frame: 0,
                        ambiguity: AmbiguityPolicy::ResolveUnanimous
                    }
                ),
                Err("ambiguous_codon")
            );
        }
        assert_eq!(
            translate(
                b"TAR",
                TranslationConfig {
                    frame: 0,
                    ambiguity: AmbiguityPolicy::ResolveUnanimous
                }
            ),
            Err("stop_codon")
        );
        assert!(translate(
            b"AGCTTGTAAA",
            TranslationConfig {
                frame: 1,
                ambiguity: AmbiguityPolicy::Reject
            }
        )
        .is_ok());
        assert_eq!(
            translate(
                b"NGCTTGTAAA",
                TranslationConfig {
                    frame: 1,
                    ambiguity: AmbiguityPolicy::Reject
                }
            ),
            Err("invalid_nucleotide_sequence")
        );
    }
    #[test]
    fn hamming_boundaries() {
        let a = AaSequence::new(b"ACK").unwrap();
        let b = AaSequence::new(b"ACE").unwrap();
        let t = 1.0 / 3.0;
        assert_eq!(
            single_linkage(&[a.clone(), b.clone()], t).unwrap(),
            vec![vec![0, 1]]
        );
        assert_eq!(
            single_linkage(&[a, b], f64::from_bits(t.to_bits() - 1)).unwrap(),
            vec![vec![0], vec![1]]
        );
        assert_eq!(identity(&[]), Vec::<Vec<usize>>::new());
    }
}
