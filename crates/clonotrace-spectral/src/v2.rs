//! CT spectral-v2: bounded deterministic normalized spectral clustering.
//! Matrix rows/sequence rows MUST be preordered by canonical byte-order record ID.
//! Numerical eigenspace ties and k-means ties are explicit method policy.
pub const MAX_ROWS: usize = 512;
pub const EIGEN_BLOCK_TOLERANCE: f64 = 1e-10;
pub const DISTANCE_TIE_TOLERANCE: f64 = 1e-12;
#[derive(Clone, Debug)]
pub struct Options {
    pub sigma: f64,
    pub edge_cutoff: Option<f64>,
    pub clusters_per_component: usize,
    pub max_group_rows: usize,
    pub max_iterations: usize,
}
#[derive(Clone, Debug)]
pub struct Component {
    pub members: Vec<usize>,
    pub effective_clusters: usize,
    pub eigenvalues: Vec<f64>,
    pub eigen_blocks: Vec<(usize, usize)>,
    pub iterations: usize,
    pub inertia: f64,
    pub residual: f64,
    pub projector: Vec<Vec<f64>>,
}
#[derive(Clone, Debug)]
pub struct Clustering {
    pub groups: Vec<Vec<usize>>,
    pub components: Vec<Component>,
    pub zero_edges: usize,
    pub underflow_edges: usize,
}
fn norm(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |n, &x| n.hypot(x))
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn squared(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum()
}
fn validate(w: &[Vec<f64>], limit: usize) -> Result<(), String> {
    let n = w.len();
    if limit == 0 || limit > MAX_ROWS || n == 0 || n > limit {
        return Err("spectral dense row bound exceeded".into());
    }
    for (i, row) in w.iter().enumerate() {
        if row.len() != n {
            return Err("affinity matrix must be square".into());
        }
        for (j, &x) in row.iter().enumerate() {
            if !x.is_finite() || x < 0.0 || (i == j && x != 0.0) {
                return Err("affinity requires finite nonnegative values and zero diagonal".into());
            }
        }
    }
    for (i, row) in w.iter().enumerate() {
        for (j, other) in w.iter().enumerate().skip(i + 1) {
            if row[j] != other[i] {
                return Err("affinity matrix must be exactly symmetric".into());
            }
        }
    }
    Ok(())
}
fn connected(w: &[Vec<f64>]) -> Vec<Vec<usize>> {
    let n = w.len();
    let mut seen = vec![false; n];
    let mut result = Vec::new();
    for first in 0..n {
        if seen[first] {
            continue;
        }
        let mut group = vec![first];
        seen[first] = true;
        let mut cursor = 0;
        while cursor < group.len() {
            let i = group[cursor];
            for (j, visited) in seen.iter_mut().enumerate() {
                if !*visited && w[i][j] > 0.0 {
                    *visited = true;
                    group.push(j);
                }
            }
            cursor += 1;
        }
        group.sort_unstable();
        result.push(group);
    }
    result
}
type EigenDecomposition = (Vec<f64>, Vec<Vec<f64>>, f64);
type CanonicalBasis = (Vec<Vec<f64>>, Vec<(usize, usize)>);
// Update the authoritative upper triangle and transposed eigenvector columns.
// Each logical cell receives the identical two scalar operations as the accepted
// full symmetric update, with unchanged cyclic pivots. Independent cells can be
// traversed contiguously; no reductions, tolerance or arithmetic order changes.
fn rotate_dense(a: &mut [Vec<f64>], u: &mut [Vec<f64>], p: usize, q: usize, c: f64, s: f64) {
    let n = a.len();
    let (before_q, from_q) = a.split_at_mut(q);
    let row_q = &mut from_q[0];
    let (before_p, from_p) = before_q.split_at_mut(p);
    let (row_p, between) = from_p.split_first_mut().unwrap();
    for row in before_p {
        let x = row[p];
        let y = row[q];
        row[p] = c * x - s * y;
        row[q] = s * x + c * y;
    }
    for (offset, row) in between.iter_mut().enumerate() {
        let k = p + 1 + offset;
        let x = row_p[k];
        let y = row[q];
        row_p[k] = c * x - s * y;
        row[q] = s * x + c * y;
    }
    for (x, y) in row_p[q + 1..n].iter_mut().zip(&mut row_q[q + 1..n]) {
        let old_x = *x;
        let old_y = *y;
        *x = c * old_x - s * old_y;
        *y = s * old_x + c * old_y;
    }
    let (before_q, from_q) = u.split_at_mut(q);
    for (x, y) in before_q[p][..n].iter_mut().zip(&mut from_q[0][..n]) {
        let old_x = *x;
        let old_y = *y;
        *x = c * old_x - s * old_y;
        *y = s * old_x + c * old_y;
    }
}
fn eigen(l: &[Vec<f64>]) -> Result<EigenDecomposition, String> {
    let n = l.len();
    // Extra trailing storage changes cache-set placement, never logical matrix cells.
    // Keep the exact accepted cyclic pivot schedule and every binary64 operation.
    // Small components retain their original allocation path.
    let (mut a, mut u) = if n >= 64 {
        let width = n + 8;
        let rows = l
            .iter()
            .map(|row| {
                let mut copy = Vec::with_capacity(width);
                copy.extend_from_slice(row);
                copy.resize(width, 0.0);
                copy
            })
            .collect::<Vec<_>>();
        (rows, vec![vec![0.0; width]; n])
    } else {
        (l.to_vec(), vec![vec![0.0; n]; n])
    };
    for (i, row) in u.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let mut converged = false;
    for _ in 0..100 {
        let mut largest = 0.0f64;
        for p in 0..n {
            for q in p + 1..n {
                let off = a[p][q];
                largest = largest.max(off.abs());
                if off.abs() <= 1e-13 {
                    continue;
                }
                let tau = (a[q][q] - a[p][p]) / (2.0 * off);
                let t = if tau >= 0.0 { 1.0 } else { -1.0 } / (tau.abs() + tau.hypot(1.0));
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                a[p][p] -= t * off;
                a[q][q] += t * off;
                a[p][q] = 0.0;
                a[q][p] = 0.0;
                if n >= 64 {
                    rotate_dense(&mut a, &mut u, p, q, c, s);
                } else {
                    for k in 0..n {
                        if k != p && k != q {
                            let x = a[k][p];
                            let y = a[k][q];
                            a[k][p] = c * x - s * y;
                            a[p][k] = a[k][p];
                            a[k][q] = s * x + c * y;
                            a[q][k] = a[k][q];
                        }
                        let x = u[k][p];
                        let y = u[k][q];
                        u[k][p] = c * x - s * y;
                        u[k][q] = s * x + c * y;
                    }
                }
            }
        }
        if largest <= 1e-13 {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("spectral eigensolver did not converge".into());
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]).then(i.cmp(&j)));
    let values: Vec<f64> = order.iter().map(|&i| a[i][i]).collect();
    let vectors: Vec<Vec<f64>> = (0..n)
        .map(|r| {
            order
                .iter()
                .map(|&c| if n >= 64 { u[c][r] } else { u[r][c] })
                .collect()
        })
        .collect();
    let mut maximum = 0.0f64;
    for j in 0..n {
        let v: Vec<f64> = vectors.iter().map(|r| r[j]).collect();
        let error: Vec<f64> = l
            .iter()
            .zip(&v)
            .map(|(r, x)| dot(r, &v) - values[j] * x)
            .collect();
        maximum = maximum.max(norm(&error));
    }
    if !maximum.is_finite() || maximum > 1e-9 || values.iter().any(|x| !x.is_finite()) {
        return Err("spectral eigen residual exceeds bound".into());
    }
    Ok((values, vectors, maximum))
}
fn canonical_basis(
    values: &[f64],
    u: &[Vec<f64>],
    stationary: &[f64],
) -> Result<CanonicalBasis, String> {
    let n = values.len();
    let mut columns: Vec<Vec<f64>> = Vec::new();
    let mut blocks = Vec::new();
    let mut start = 0;
    while start < n {
        let mut end = start + 1;
        while end < n && values[end] - values[start] <= EIGEN_BLOCK_TOLERANCE {
            end += 1;
        }
        let mut basis: Vec<Vec<f64>> = if start == 0 {
            vec![stationary.to_vec()]
        } else {
            Vec::new()
        };
        for axis in 0..n {
            if basis.len() == end - start {
                break;
            }
            let mut v: Vec<f64> = (0..n)
                .map(|r| (start..end).map(|c| u[r][c] * u[axis][c]).sum())
                .collect();
            for _ in 0..2 {
                for q in columns.iter().chain(&basis) {
                    let projection = dot(q, &v);
                    for (x, y) in v.iter_mut().zip(q) {
                        *x -= projection * y;
                    }
                }
            }
            let magnitude = norm(&v);
            if magnitude > 1e-8 {
                for x in &mut v {
                    *x /= magnitude;
                }
                basis.push(v);
            }
        }
        if basis.len() != end - start {
            return Err("spectral canonical eigenspace rank failure".into());
        }
        columns.extend(basis);
        blocks.push((start, end));
        start = end;
    }
    let rows = (0..n)
        .map(|r| columns.iter().map(|c| c[r]).collect())
        .collect();
    Ok((rows, blocks))
}
fn lloyd(
    y: &[Vec<f64>],
    k: usize,
    max_iterations: usize,
) -> Result<(Vec<Vec<usize>>, usize, f64), String> {
    let n = y.len();
    let mut selected = vec![0];
    while selected.len() < k {
        let mut best = None;
        let mut value = -1.0;
        for i in 0..n {
            if selected.contains(&i) {
                continue;
            }
            let d = selected
                .iter()
                .map(|&c| squared(&y[i], &y[c]))
                .fold(f64::INFINITY, f64::min);
            if d > value + DISTANCE_TIE_TOLERANCE {
                best = Some(i);
                value = d;
            }
        }
        selected.push(best.ok_or("spectral center initialization failure")?);
    }
    let mut centers: Vec<Vec<f64>> = selected.iter().map(|&i| y[i].clone()).collect();
    let mut previous = None;
    for iteration in 1..=max_iterations {
        let mut labels = Vec::with_capacity(n);
        let mut errors = Vec::with_capacity(n);
        let mut counts = vec![0usize; k];
        for row in y {
            let ds: Vec<f64> = centers.iter().map(|c| squared(row, c)).collect();
            let mut best = 0;
            for j in 1..k {
                if ds[j] < ds[best] - DISTANCE_TIE_TOLERANCE {
                    best = j;
                }
            }
            labels.push(best);
            errors.push(ds[best]);
            counts[best] += 1;
        }
        for j in 0..k {
            if counts[j] > 0 {
                continue;
            }
            let mut pick = None;
            let mut error = -1.0;
            for i in 0..n {
                if counts[labels[i]] > 1 && errors[i] > error + DISTANCE_TIE_TOLERANCE {
                    pick = Some(i);
                    error = errors[i];
                }
            }
            let i = pick.ok_or("spectral empty center repair failure")?;
            counts[labels[i]] -= 1;
            labels[i] = j;
            counts[j] += 1;
        }
        centers = vec![vec![0.0; k]; k];
        for (i, row) in y.iter().enumerate() {
            for (q, &v) in row.iter().enumerate() {
                centers[labels[i]][q] += v;
            }
        }
        for (j, c) in centers.iter_mut().enumerate() {
            for x in c {
                *x /= counts[j] as f64;
            }
        }
        if previous.as_ref() == Some(&labels) {
            let inertia = y
                .iter()
                .enumerate()
                .map(|(i, row)| squared(row, &centers[labels[i]]))
                .sum();
            let mut groups: Vec<Vec<usize>> = (0..k)
                .map(|j| (0..n).filter(|&i| labels[i] == j).collect())
                .collect();
            groups.sort_by_key(|g| g[0]);
            return Ok((groups, iteration, inertia));
        }
        previous = Some(labels);
    }
    Err("spectral kmeans iteration bound exceeded".into())
}
/// Generic affinity API: accepts an explicit finite symmetric nonnegative zero-diagonal graph.
/// Weights may be zero. Components never merge. k saturates at component size.
pub fn cluster_affinity(
    w: &[Vec<f64>],
    clusters_per_component: usize,
    max_group_rows: usize,
    max_iterations: usize,
) -> Result<Clustering, String> {
    validate(w, max_group_rows)?;
    if clusters_per_component == 0 || max_iterations == 0 || max_iterations > 10000 {
        return Err("spectral cluster count/iteration bounds invalid".into());
    }
    let zero_edges = (0..w.len())
        .map(|i| (i + 1..w.len()).filter(|&j| w[i][j] == 0.0).count())
        .sum();
    let mut groups = Vec::new();
    let mut diagnostics = Vec::new();
    for members in connected(w) {
        let n = members.len();
        let k = clusters_per_component.min(n);
        if k == 1 || k == n {
            if k == 1 {
                groups.push(members.clone());
            } else {
                groups.extend(members.iter().map(|&i| vec![i]));
            }
            diagnostics.push(Component {
                members,
                effective_clusters: k,
                eigenvalues: vec![],
                eigen_blocks: vec![],
                iterations: 0,
                inertia: 0.0,
                residual: 0.0,
                projector: vec![],
            });
            continue;
        }
        let degree: Vec<f64> = members
            .iter()
            .map(|&i| members.iter().map(|&j| w[i][j]).sum())
            .collect();
        if degree.iter().any(|d| !d.is_finite() || *d <= 0.0) {
            return Err("spectral degree arithmetic overflow".into());
        }
        let roots: Vec<f64> = degree.iter().map(|x| x.sqrt()).collect();
        let total: f64 = degree.iter().sum();
        if !total.is_finite() {
            return Err("spectral total degree overflow".into());
        }
        let stationary: Vec<f64> = roots.iter().map(|x| x / total.sqrt()).collect();
        let l: Vec<Vec<f64>> = members
            .iter()
            .enumerate()
            .map(|(a, &i)| {
                members
                    .iter()
                    .enumerate()
                    .map(|(b, &j)| {
                        if i == j {
                            1.0
                        } else {
                            -w[i][j] / (roots[a] * roots[b])
                        }
                    })
                    .collect()
            })
            .collect();
        if l.iter().flatten().any(|x| !x.is_finite()) {
            return Err("spectral normalization arithmetic failure".into());
        }
        let (values, u, residual) = eigen(&l)?;
        let (basis, blocks) = canonical_basis(&values, &u, &stationary)?;
        let projector: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| (0..k).map(|q| basis[i][q] * basis[j][q]).sum())
                    .collect()
            })
            .collect();
        let mut embedding: Vec<Vec<f64>> = basis.into_iter().map(|r| r[..k].to_vec()).collect();
        for row in &mut embedding {
            let magnitude = norm(row);
            if !magnitude.is_finite() || magnitude <= 1e-8 {
                return Err("spectral embedding is numerically ill-conditioned".into());
            }
            for x in row {
                *x /= magnitude;
            }
        }
        let (partition, iterations, inertia) = lloyd(&embedding, k, max_iterations)?;
        groups.extend(
            partition
                .into_iter()
                .map(|g| g.into_iter().map(|i| members[i]).collect::<Vec<_>>()),
        );
        diagnostics.push(Component {
            members,
            effective_clusters: k,
            eigenvalues: values,
            eigen_blocks: blocks,
            iterations,
            inertia,
            residual,
            projector,
        });
    }
    groups.sort_by_key(|g| g[0]);
    Ok(Clustering {
        groups,
        components: diagnostics,
        zero_edges,
        underflow_edges: 0,
    })
}
pub fn cluster(seqs: &[&[u8]], options: &Options) -> Result<Clustering, String> {
    if !options.sigma.is_finite()
        || options.sigma <= 0.0
        || options
            .edge_cutoff
            .is_some_and(|t| !t.is_finite() || !(0.0..=1.0).contains(&t))
        || options.max_group_rows == 0
        || options.max_group_rows > MAX_ROWS
        || seqs.is_empty()
        || seqs.len() > options.max_group_rows
        || seqs[0].is_empty()
        || seqs.iter().any(|s| s.len() != seqs[0].len())
    {
        return Err("invalid spectral sequence/options profile".into());
    }
    let n = seqs.len();
    let mut w = vec![vec![0.0; n]; n];
    let mut underflow_edges = 0;
    for i in 0..n {
        for j in i + 1..n {
            let d = seqs[i].iter().zip(seqs[j]).filter(|(a, b)| a != b).count() as f64
                / seqs[i].len() as f64;
            if options.edge_cutoff.is_some_and(|t| d > t) {
                continue;
            }
            let weight = (-(d / options.sigma).powi(2) / 2.0).exp();
            if weight == 0.0 {
                underflow_edges += 1;
            }
            w[i][j] = weight;
            w[j][i] = weight;
        }
    }
    let mut result = cluster_affinity(
        &w,
        options.clusters_per_component,
        options.max_group_rows,
        options.max_iterations,
    )?;
    result.underflow_edges = underflow_edges;
    Ok(result)
}
#[cfg(test)]
mod tests {
    use super::*;
    fn opts(k: usize) -> Options {
        Options {
            sigma: 0.25,
            edge_cutoff: None,
            clusters_per_component: k,
            max_group_rows: MAX_ROWS,
            max_iterations: 256,
        }
    }
    #[test]
    fn repeated_eigenvalues_and_zero_cut() {
        for ss in [
            vec![b"AAAA".as_slice(), b"AAAC", b"AAAG"],
            vec![b"AAAA".as_slice(), b"AAAC", b"AACC"],
        ] {
            let x = cluster(&ss, &opts(2)).unwrap();
            assert_eq!(x.groups.len(), 2);
            assert_eq!(x.groups.iter().map(Vec::len).sum::<usize>(), 3);
            assert!(x.components[0].residual < 1e-9);
        }
    }
    #[test]
    fn cutoff_components_and_underflow() {
        let ss: Vec<&[u8]> = vec![b"AAAA", b"AAAC", b"CCCC"];
        let mut o = opts(1);
        o.edge_cutoff = Some(0.25);
        assert_eq!(cluster(&ss, &o).unwrap().groups, vec![vec![0, 1], vec![2]]);
        o.sigma = 1e-12;
        o.edge_cutoff = None;
        let x = cluster(&ss, &o).unwrap();
        assert_eq!(x.groups.len(), 3);
        assert_eq!(x.underflow_edges, 3);
    }
    #[test]
    fn validate_graph() {
        assert!(cluster_affinity(&[vec![0.0, 1.0], vec![0.0, 0.0]], 2, MAX_ROWS, 256).is_err());
        assert!(cluster_affinity(&[vec![f64::NAN]], 1, MAX_ROWS, 256).is_err());
        assert!(cluster_affinity(&[vec![0.0]], 0, MAX_ROWS, 256).is_err());
    }
}
