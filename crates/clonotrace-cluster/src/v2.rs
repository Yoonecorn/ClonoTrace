//! CT hierarchical-v2: deterministic dense agglomeration over original row weights.
//! Callers establish canonical byte-order record IDs before passing matrix rows.
pub const MAX_ROWS: usize = 512;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Linkage {
    Single,
    Average,
    Complete,
}
impl Linkage {
    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "single" => Ok(Self::Single),
            "average" => Ok(Self::Average),
            "complete" => Ok(Self::Complete),
            _ => Err("unsupported hierarchical linkage".into()),
        }
    }
}
#[derive(Clone, Debug)]
pub struct Merge {
    pub left: usize,
    pub right: usize,
    pub height: f64,
    pub size: usize,
}
#[derive(Clone, Debug)]
pub struct Dendrogram {
    pub rows: usize,
    pub merges: Vec<Merge>,
}
fn validate(d: &[Vec<f64>], limit: usize) -> Result<(), String> {
    let n = d.len();
    if limit == 0 || limit > MAX_ROWS || n == 0 || n > limit {
        return Err("hierarchical dense row bound exceeded".into());
    }
    for (i, row) in d.iter().enumerate() {
        if row.len() != n {
            return Err("distance matrix must be square".into());
        }
        for (j, &value) in row.iter().enumerate() {
            if !value.is_finite() || value < 0.0 || (i == j && value != 0.0) {
                return Err(
                    "distance matrix requires finite nonnegative values and zero diagonal".into(),
                );
            }
        }
    }
    for (i, row) in d.iter().enumerate() {
        for (j, other) in d.iter().enumerate().skip(i + 1) {
            if row[j] != other[i] {
                return Err("distance matrix must be exactly symmetric".into());
            }
        }
    }
    Ok(())
}
pub fn hamming_matrix(seqs: &[&[u8]], limit: usize) -> Result<Vec<Vec<f64>>, String> {
    if seqs.is_empty()
        || seqs.len() > limit
        || limit == 0
        || limit > MAX_ROWS
        || seqs[0].is_empty()
        || seqs.iter().any(|s| s.len() != seqs[0].len())
    {
        return Err("Hamming requires bounded nonempty equal-length sequences".into());
    }
    let n = seqs.len();
    let mut d = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let v = seqs[i].iter().zip(seqs[j]).filter(|(a, b)| a != b).count() as f64
                / seqs[i].len() as f64;
            d[i][j] = v;
            d[j][i] = v;
        }
    }
    Ok(d)
}
pub fn build(distances: &[Vec<f64>], linkage: Linkage, limit: usize) -> Result<Dendrogram, String> {
    validate(distances, limit)?;
    let n = distances.len();
    let mut groups: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut nodes: Vec<usize> = (0..n).collect();
    let mut merges = Vec::with_capacity(n - 1);
    while groups.len() > 1 {
        let mut best = None;
        for i in 0..groups.len() {
            for j in i + 1..groups.len() {
                let mut value = if linkage == Linkage::Single {
                    f64::INFINITY
                } else {
                    0.0
                };
                for &a in &groups[i] {
                    for &b in &groups[j] {
                        let d = distances[a][b];
                        value = match linkage {
                            Linkage::Single => value.min(d),
                            Linkage::Complete => value.max(d),
                            Linkage::Average => value + d,
                        };
                    }
                }
                if linkage == Linkage::Average {
                    value /= (groups[i].len() * groups[j].len()) as f64;
                }
                if !value.is_finite() {
                    return Err("hierarchical distance arithmetic overflow".into());
                }
                if best.is_none_or(|(d, _, _)| value < d) {
                    best = Some((value, i, j));
                }
            }
        }
        let (height, i, j) = best.ok_or("hierarchical empty merge set")?;
        let tail = groups.remove(j);
        groups[i].extend(tail);
        groups[i].sort_unstable();
        merges.push(Merge {
            left: nodes[i],
            right: nodes[j],
            height,
            size: groups[i].len(),
        });
        nodes[i] = n + merges.len() - 1;
        nodes.remove(j);
    }
    Ok(Dendrogram { rows: n, merges })
}
impl Dendrogram {
    fn cut(&self, threshold: Option<f64>, count: Option<usize>) -> Result<Vec<Vec<usize>>, String> {
        if self.rows == 0 || self.rows > MAX_ROWS || self.merges.len() != self.rows - 1 {
            return Err("invalid dendrogram shape".into());
        }
        // Public dendrogram fields may be externally constructed; validate the complete tree
        // before honoring an early cut, including nodes beyond the requested threshold.
        let mut active = vec![Some(1usize); self.rows];
        for (step, m) in self.merges.iter().enumerate() {
            if !m.height.is_finite()
                || m.height < 0.0
                || m.left >= self.rows + step
                || m.right >= self.rows + step
                || m.left == m.right
            {
                return Err("invalid dendrogram merge".into());
            }
            let a = active[m.left].take().ok_or("reused dendrogram node")?;
            let b = active[m.right].take().ok_or("reused dendrogram node")?;
            if a + b != m.size {
                return Err("invalid dendrogram merge size".into());
            }
            active.push(Some(m.size));
        }
        let mut nodes: Vec<Option<Vec<usize>>> = (0..self.rows).map(|i| Some(vec![i])).collect();
        let mut remaining = self.rows;
        for (step, m) in self.merges.iter().enumerate() {
            if threshold.is_some_and(|t| m.height > t)
                || count.is_some_and(|k| remaining <= k.min(self.rows))
            {
                break;
            }
            if !m.height.is_finite()
                || m.height < 0.0
                || m.left >= self.rows + step
                || m.right >= self.rows + step
                || m.left == m.right
            {
                return Err("invalid dendrogram merge".into());
            }
            let mut left = nodes[m.left].take().ok_or("reused dendrogram node")?;
            let right = nodes[m.right].take().ok_or("reused dendrogram node")?;
            left.extend(right);
            left.sort_unstable();
            if left.len() != m.size {
                return Err("invalid dendrogram merge size".into());
            }
            nodes.push(Some(left));
            remaining -= 1;
        }
        let mut result: Vec<Vec<usize>> = nodes.into_iter().flatten().collect();
        result.sort_by_key(|g| g[0]);
        Ok(result)
    }
    pub fn cut_threshold(&self, threshold: f64) -> Result<Vec<Vec<usize>>, String> {
        if !threshold.is_finite() || threshold < 0.0 {
            return Err("hierarchical threshold must be finite and nonnegative".into());
        }
        self.cut(Some(threshold), None)
    }
    /// The requested count saturates at component row count, including singleton work units.
    pub fn cut_count(&self, count: usize) -> Result<Vec<Vec<usize>>, String> {
        if count == 0 {
            return Err("hierarchical cluster count must be positive".into());
        }
        self.cut(None, Some(count))
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn linkage_weighting_and_cuts() {
        let ss: Vec<&[u8]> = vec![b"AAAA", b"AAAA", b"AAAC", b"AACC"];
        let d = hamming_matrix(&ss, MAX_ROWS).unwrap();
        let tree = build(&d, Linkage::Average, MAX_ROWS).unwrap();
        assert_eq!(
            tree.cut_threshold(0.375).unwrap(),
            vec![vec![0, 1, 2], vec![3]]
        );
        assert_eq!(tree.cut_count(1).unwrap(), vec![vec![0, 1, 2, 3]]);
        assert_eq!(
            tree.cut_count(8).unwrap(),
            vec![vec![0], vec![1], vec![2], vec![3]]
        );
        assert!(tree.cut_count(0).is_err());
    }
    #[test]
    fn forged_tree_is_rejected_before_early_cut() {
        let tree = Dendrogram {
            rows: 3,
            merges: vec![
                Merge {
                    left: 0,
                    right: 1,
                    height: 0.5,
                    size: 2,
                },
                Merge {
                    left: 3,
                    right: 99,
                    height: 0.8,
                    size: 3,
                },
            ],
        };
        assert!(tree.cut_count(3).is_err());
        assert!(tree.cut_threshold(0.0).is_err());
        let tree = Dendrogram {
            rows: 2,
            merges: vec![Merge {
                left: 0,
                right: 1,
                height: f64::NAN,
                size: 2,
            }],
        };
        assert!(tree.cut_count(2).is_err());
    }
    #[test]
    fn matrix_validation() {
        assert!(build(&[vec![0.0, 1.0], vec![0.5, 0.0]], Linkage::Single, MAX_ROWS).is_err());
        assert!(build(&[vec![f64::NAN]], Linkage::Single, MAX_ROWS).is_err());
        assert!(hamming_matrix(&[b"A", b"AA"], MAX_ROWS).is_err());
    }
}
