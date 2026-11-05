//! Exact sequence collapse preserves unweighted single-linkage components.
pub mod threshold;
pub mod v2;
use clonotrace_airr::Record;
use clonotrace_group::{Method, UnionFind};
use std::collections::HashMap;
pub fn assign(
    records: &[Record],
    ids: &[usize],
    method: Method,
    threshold: f64,
) -> Vec<Vec<usize>> {
    if method == Method::Changeo && threshold == 1.0 {
        return vec![ids.to_vec()];
    }
    let mut by_sequence: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for &i in ids {
        by_sequence.entry(&records[i].junction).or_default().push(i);
    }
    let mut unique: Vec<_> = by_sequence.into_values().collect();
    // Canonical representative order makes computation deterministic independent of hash seeds.
    unique.sort_unstable_by_key(|g| g[0]);
    if method == Method::Scoper || threshold == 0.0 {
        return unique;
    }
    let mut uf = UnionFind::new(unique.len());
    for a in 0..unique.len() {
        for b in a + 1..unique.len() {
            if clonotrace_distance::linked(
                &records[unique[a][0]].junction,
                &records[unique[b][0]].junction,
                threshold,
            ) {
                uf.join(a, b);
            }
        }
    }
    uf.groups()
        .into_iter()
        .map(|g| {
            g.into_iter()
                .flat_map(|i| unique[i].iter().copied())
                .collect()
        })
        .collect()
}

/// Independent row-weighted agglomeration. Input IDs must be canonical-ID sorted.
/// Linkage and tie policy are CT hierarchical-v1, not upstream-native behavior.
pub fn hierarchical(
    seqs: &[&[u8]],
    threshold: f64,
    linkage: &str,
) -> Result<Vec<Vec<usize>>, String> {
    if !matches!(linkage, "single" | "average" | "complete") || seqs.len() > 128 {
        return Err("unsupported linkage or hierarchical group bound".into());
    }
    let n = seqs.len();
    if n == 0
        || seqs[0].is_empty()
        || seqs.iter().any(|s| s.len() != seqs[0].len())
        || !threshold.is_finite()
        || !(0.0..=1.0).contains(&threshold)
    {
        return Err("invalid hierarchical input or threshold".into());
    }
    let distances: Vec<Vec<f64>> = seqs
        .iter()
        .map(|a| {
            seqs.iter()
                .map(|b| a.iter().zip(*b).filter(|(x, y)| x != y).count() as f64 / a.len() as f64)
                .collect()
        })
        .collect();
    let mut groups: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    while groups.len() > 1 {
        let mut best: Option<(f64, usize, usize)> = None;
        for i in 0..groups.len() {
            for j in i + 1..groups.len() {
                let mut value = if linkage == "single" {
                    f64::INFINITY
                } else {
                    0.0
                };
                for &a in &groups[i] {
                    for &b in &groups[j] {
                        let d = distances[a][b];
                        value = match linkage {
                            "single" => value.min(d),
                            "complete" => value.max(d),
                            _ => value + d,
                        };
                    }
                }
                if linkage == "average" {
                    value /= (groups[i].len() * groups[j].len()) as f64;
                }
                // Groups remain sorted by their canonical minimum member.
                if best.is_none_or(|(d, _, _)| value < d) {
                    best = Some((value, i, j));
                }
            }
        }
        let (value, i, j) = best.unwrap();
        if value > threshold {
            break;
        }
        let merged = groups.remove(j);
        groups[i].extend(merged);
        groups[i].sort_unstable();
    }
    Ok(groups)
}
