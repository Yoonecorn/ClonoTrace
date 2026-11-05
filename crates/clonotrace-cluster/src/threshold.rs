//! Deterministic nearest-distinct-neighbor gap estimation and exhaustive deletion stability.
//! Stability ranges are not confidence intervals or biological error calibration.
use std::collections::BTreeSet;
#[derive(Clone, Copy, Debug)]
pub struct Config {
    pub max_distinct: usize,
    /// Conservative upper bound counts normalized Hamming evaluations, not nucleotide operations.
    pub max_distance_evaluations: u64,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            max_distinct: 128,
            max_distance_evaluations: 4_000_000,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct Estimate {
    pub threshold: Option<f64>,
    pub reason: Option<&'static str>,
    /// One distance per lexicographically ordered distinct junction.
    pub nearest: Vec<f64>,
    pub distinct_distances: Vec<f64>,
    pub gaps: Vec<f64>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Replicate {
    pub omitted_junction: Vec<u8>,
    pub threshold: Option<f64>,
    pub reason: Option<&'static str>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Sensitivity {
    pub threshold: f64,
    pub clone_count: usize,
    /// Disagreement versus the effective selected threshold on distinct-junction pairs.
    pub pair_disagreements: Option<usize>,
    pub pairs: usize,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Analysis {
    pub input_rows: usize,
    pub unique_junctions: Vec<Vec<u8>>,
    pub base: Estimate,
    pub leave_one_out: Vec<Replicate>,
    pub leave_one_out_accepted: usize,
    /// Range over base and successfully estimated leave-one-distinct-junction-out replicates.
    pub stability_range: Option<[f64; 2]>,
    pub override_threshold: Option<f64>,
    pub selected_threshold: Option<f64>,
    pub sensitivity: Vec<Sensitivity>,
    pub distance_evaluations: u64,
    pub reserved_distance_evaluations: u64,
}
fn hamming(a: &[u8], b: &[u8], evaluations: &mut u64) -> f64 {
    *evaluations += 1;
    a.iter().zip(b).filter(|(x, y)| x != y).count() as f64 / a.len() as f64
}
fn estimate(indices: &[usize], matrix: &[Vec<f64>]) -> Estimate {
    let mut out = Estimate {
        threshold: None,
        reason: None,
        nearest: Vec::new(),
        distinct_distances: Vec::new(),
        gaps: Vec::new(),
    };
    if indices.len() < 2 {
        out.reason = Some("no_distinct_neighbor");
        return out;
    }
    for &i in indices {
        let minimum = indices
            .iter()
            .copied()
            .filter(|&j| j != i)
            .map(|j| matrix[i][j])
            .fold(f64::INFINITY, f64::min);
        out.nearest.push(minimum);
    }
    out.distinct_distances = out.nearest.clone();
    out.distinct_distances.sort_by(f64::total_cmp);
    out.distinct_distances.dedup();
    if out.distinct_distances.len() < 2 {
        out.reason = Some("insufficient_distance_variation");
        return out;
    }
    out.gaps = out
        .distinct_distances
        .windows(2)
        .map(|x| x[1] - x[0])
        .collect();
    let maximum = out.gaps.iter().copied().fold(0.0, f64::max);
    if out.gaps.iter().filter(|&&x| x == maximum).count() != 1 {
        out.reason = Some("tied_largest_gap");
        return out;
    }
    if let Some(i) = out.gaps.iter().position(|&x| x == maximum) {
        out.threshold = Some((out.distinct_distances[i] + out.distinct_distances[i + 1]) / 2.0);
    }
    out
}
fn partition(matrix: &[Vec<f64>], threshold: f64) -> Vec<usize> {
    let mut roots: Vec<usize> = (0..matrix.len()).collect();
    fn root(roots: &[usize], mut i: usize) -> usize {
        while roots[i] != i {
            i = roots[i];
        }
        i
    }
    for (i, row) in matrix.iter().enumerate() {
        for (j, &distance) in row.iter().enumerate().skip(i + 1) {
            if distance <= threshold {
                let a = root(&roots, i);
                let b = root(&roots, j);
                if a != b {
                    roots[a.max(b)] = a.min(b);
                }
            }
        }
    }
    (0..matrix.len()).map(|i| root(&roots, i)).collect()
}
/// Analyze one exact donor/sample/locus/V/J/length work group.
/// Input errors and resource bounds return Err; scientific refusal returns an Analysis
/// with selected_threshold=None unless an explicit finite [0,1] override was supplied.
/// No leave-one-out estimate is substituted for a refused base estimate.
pub fn analyze(
    seqs: &[&[u8]],
    override_threshold: Option<f64>,
    config: Config,
) -> Result<Analysis, &'static str> {
    if override_threshold.is_some_and(|x| !x.is_finite() || !(0.0..=1.0).contains(&x)) {
        return Err("invalid_threshold_override");
    }
    if config.max_distinct == 0 {
        return Err("max_distinct_must_be_positive");
    }
    if config.max_distinct > 128 {
        return Err("max_distinct_exceeds_128");
    }
    let input_rows = seqs.len();
    let length = seqs.first().map_or(0, |s| s.len());
    if seqs
        .iter()
        .any(|s| s.is_empty() || s.len() != length || !s.iter().all(|b| b"ACGT".contains(b)))
    {
        return Err("invalid_equal_length_nucleotide_input");
    }
    let unique_junctions: Vec<Vec<u8>> = seqs
        .iter()
        .copied()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .map(<[u8]>::to_vec)
        .collect();
    let n = unique_junctions.len();
    if n > config.max_distinct {
        return Err("max_distinct_exceeded");
    }
    let n64 = u64::try_from(n).map_err(|_| "distance_budget_overflow")?;
    let pairs = n64
        .checked_mul(n64.saturating_sub(1))
        .ok_or("distance_budget_overflow")?
        / 2;
    let base_cost = n64
        .checked_mul(n64.saturating_sub(1))
        .ok_or("distance_budget_overflow")?;
    let loo_cost = base_cost
        .checked_mul(n64.saturating_sub(2))
        .ok_or("distance_budget_overflow")?;
    let sensitivity_cost = n64
        .checked_add(2)
        .and_then(|m| m.checked_mul(pairs))
        .ok_or("distance_budget_overflow")?;
    let reserved = base_cost
        .checked_add(loo_cost)
        .and_then(|c| c.checked_add(sensitivity_cost))
        .ok_or("distance_budget_overflow")?;
    if reserved > config.max_distance_evaluations {
        return Err("max_distance_evaluations_exceeded");
    }
    // Exact symmetric distance cache; the conservative scalar reservation and failures remain unchanged.
    let mut evaluations = 0;
    let mut matrix = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let value = hamming(&unique_junctions[i], &unique_junctions[j], &mut evaluations);
            matrix[i][j] = value;
            matrix[j][i] = value;
        }
    }
    let indices: Vec<usize> = (0..n).collect();
    let base = estimate(&indices, &matrix);
    let mut leave_one_out = Vec::with_capacity(n);
    for (i, omitted) in unique_junctions.iter().enumerate() {
        let reduced: Vec<usize> = indices.iter().copied().filter(|&j| i != j).collect();
        let e = estimate(&reduced, &matrix);
        leave_one_out.push(Replicate {
            omitted_junction: omitted.clone(),
            threshold: e.threshold,
            reason: e.reason,
        });
    }
    let leave_one_out_accepted = leave_one_out
        .iter()
        .filter(|r| r.threshold.is_some())
        .count();
    let mut estimates: Vec<f64> = base
        .threshold
        .into_iter()
        .chain(leave_one_out.iter().filter_map(|r| r.threshold))
        .collect();
    estimates.sort_by(f64::total_cmp);
    estimates.dedup();
    let stability_range = estimates
        .first()
        .zip(estimates.last())
        .map(|(&a, &b)| [a, b]);
    let selected_threshold = override_threshold.or(base.threshold);
    let mut thresholds = estimates;
    if let Some(t) = selected_threshold {
        thresholds.push(t);
    }
    thresholds.sort_by(f64::total_cmp);
    thresholds.dedup();
    let partitions: Vec<Vec<usize>> = thresholds.iter().map(|&t| partition(&matrix, t)).collect();
    let selected_labels = selected_threshold
        .and_then(|t| thresholds.iter().position(|&x| x == t))
        .map(|i| &partitions[i]);
    let mut sensitivity = Vec::new();
    for (&threshold, labels) in thresholds.iter().zip(&partitions) {
        let pair_disagreements = selected_labels.map(|selected| {
            let mut changes = 0;
            for i in 0..n {
                for j in i + 1..n {
                    if (selected[i] == selected[j]) != (labels[i] == labels[j]) {
                        changes += 1;
                    }
                }
            }
            changes
        });
        sensitivity.push(Sensitivity {
            threshold,
            clone_count: labels.iter().copied().collect::<BTreeSet<_>>().len(),
            pair_disagreements,
            pairs: usize::try_from(pairs).map_err(|_| "distance_budget_overflow")?,
        });
    }
    if evaluations > reserved {
        return Err("internal_distance_budget_violation");
    }
    Ok(Analysis {
        input_rows,
        unique_junctions,
        base,
        leave_one_out,
        leave_one_out_accepted,
        stability_range,
        override_threshold,
        selected_threshold,
        sensitivity,
        distance_evaluations: evaluations,
        reserved_distance_evaluations: reserved,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    fn known() -> Vec<&'static [u8]> {
        vec![b"AAAAAAAA", b"CAAAAAAA", b"AAAACCAA", b"CCCCCCCC"]
    }
    #[test]
    fn absolute_library_bound_precedes_allocation() {
        assert_eq!(
            analyze(
                &known(),
                None,
                Config {
                    max_distinct: 129,
                    ..Config::default()
                }
            )
            .unwrap_err(),
            "max_distinct_exceeds_128"
        );
    }
    #[test]
    fn exhaustive_deletion_stability_and_pair_sensitivity() {
        let x = analyze(&known(), None, Config::default()).unwrap();
        assert_eq!(x.base.threshold, Some(0.5));
        assert_eq!(x.stability_range, Some([0.1875, 0.5625]));
        assert_eq!(x.leave_one_out_accepted, 4);
        assert_eq!(
            x.leave_one_out
                .iter()
                .map(|r| r.threshold)
                .collect::<Vec<_>>(),
            vec![Some(0.5625), Some(0.5), Some(0.5), Some(0.1875)]
        );
        assert_eq!(
            x.sensitivity
                .iter()
                .map(|s| (s.threshold, s.clone_count, s.pair_disagreements, s.pairs))
                .collect::<Vec<_>>(),
            vec![
                (0.1875, 3, Some(2), 6),
                (0.5, 2, Some(0), 6),
                (0.5625, 2, Some(0), 6)
            ]
        );
    }
    #[test]
    fn no_fallback_and_override_keeps_refusal_diagnostics() {
        let x = analyze(&[b"AAAAAAAA", b"CAAAAAAA"], None, Config::default()).unwrap();
        assert_eq!(x.selected_threshold, None);
        assert_eq!(x.base.reason, Some("insufficient_distance_variation"));
        let x = analyze(&[b"AAAAAAAA", b"CAAAAAAA"], Some(0.15), Config::default()).unwrap();
        assert_eq!(x.selected_threshold, Some(0.15));
        assert_eq!(x.base.threshold, None);
        assert_eq!(x.leave_one_out_accepted, 0);
        assert_eq!(x.sensitivity[0].clone_count, 1);
    }
    #[test]
    fn duplicates_order_and_bounds_are_explicit() {
        let x = analyze(&known(), None, Config::default()).unwrap();
        let mut seqs = known();
        seqs.reverse();
        seqs.extend([b"AAAAAAAA".as_slice(); 20]);
        let y = analyze(&seqs, None, Config::default()).unwrap();
        assert_eq!(x.base, y.base);
        assert_eq!(x.leave_one_out, y.leave_one_out);
        assert_eq!(x.sensitivity, y.sensitivity);
        assert_eq!(
            analyze(
                &known(),
                None,
                Config {
                    max_distinct: 3,
                    ..Config::default()
                }
            )
            .unwrap_err(),
            "max_distinct_exceeded"
        );
        assert_eq!(
            analyze(
                &known(),
                None,
                Config {
                    max_distance_evaluations: 0,
                    ..Config::default()
                }
            )
            .unwrap_err(),
            "max_distance_evaluations_exceeded"
        );
        assert!(analyze(&[b"AA", b"A"], None, Config::default()).is_err());
        assert!(analyze(&[b"AN"], None, Config::default()).is_err());
        assert!(analyze(&known(), Some(f64::NAN), Config::default()).is_err());
    }
}
