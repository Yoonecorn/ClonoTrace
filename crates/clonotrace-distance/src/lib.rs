//! Binary64 normalized Hamming threshold behavior for qualified uppercase ACGT.
pub fn linked(a: &[u8], b: &[u8], threshold: f64) -> bool {
    assert!(!a.is_empty() && a.len() == b.len());
    let n = a.iter().zip(b).filter(|(x, y)| x != y).count();
    (n as f64) / (a.len() as f64) <= threshold
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn binary64_boundaries() {
        let t = 1f64 / 3.;
        assert!(linked(b"AAA", b"CAA", t));
        assert!(!linked(b"AAA", b"CAA", f64::from_bits(t.to_bits() - 1)));
        assert!(linked(b"AAAA", b"TTTT", 1.));
    }
    #[test]
    fn exhaustive_three_bases() {
        let seq: Vec<Vec<u8>> = (0..64)
            .map(|n| (0..3).map(|k| b"ACGT"[(n >> (2 * k)) & 3]).collect())
            .collect();
        for a in &seq {
            for b in &seq {
                let mismatches = (0..3).filter(|&i| a[i] != b[i]).count();
                for n in 0..=3 {
                    assert_eq!(linked(a, b, n as f64 / 3.), mismatches <= n);
                }
            }
        }
    }
}
