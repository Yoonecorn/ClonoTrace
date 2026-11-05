//! Independent experimental normalized-Laplacian binary separation, CT spectral-v1.
pub mod v2;
#[derive(Debug)]
pub struct Split {
    pub groups: Vec<Vec<usize>>,
    pub eigenvalues: Vec<f64>,
    pub fiedler: Vec<f64>,
    pub residual: f64,
}
pub fn split(seqs: &[&[u8]], sigma: f64) -> Result<Split, String> {
    let n = seqs.len();
    if !sigma.is_finite() || sigma <= 0.0 || !(1..=128).contains(&n) {
        return Err("invalid sigma or spectral group bound".into());
    }
    if seqs[0].is_empty() || seqs.iter().any(|s| s.len() != seqs[0].len()) {
        return Err("spectral requires nonempty equal-length sequences".into());
    }
    if seqs.iter().all(|s| *s == seqs[0]) {
        return Ok(Split {
            groups: vec![(0..n).collect()],
            eigenvalues: vec![],
            fiedler: vec![],
            residual: 0.0,
        });
    }
    let mut w = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let d = seqs[i].iter().zip(seqs[j]).filter(|(a, b)| a != b).count() as f64
                / seqs[i].len() as f64;
            let affinity = (-(d / sigma).powi(2) / 2.0).exp();
            if !affinity.is_finite() || affinity <= 0.0 {
                return Err("underflow or disconnected graph".into());
            }
            w[i][j] = affinity;
            w[j][i] = affinity;
        }
    }
    let degree: Vec<f64> = w.iter().map(|r| r.iter().sum()).collect();
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            l[i][j] = if i == j {
                1.0
            } else {
                -w[i][j] / (degree[i].sqrt() * degree[j].sqrt())
            };
        }
    }
    let mut a = l.clone();
    let mut u = vec![vec![0.0; n]; n];
    for (i, row) in u.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let mut converged = false;
    for _ in 0..100 {
        let mut maximum = 0.0_f64;
        for p in 0..n {
            for q in p + 1..n {
                let off = a[p][q];
                maximum = maximum.max(off.abs());
                if off.abs() <= 1e-13 {
                    continue;
                }
                let tau = (a[q][q] - a[p][p]) / (2.0 * off);
                let t =
                    if tau >= 0.0 { 1.0 } else { -1.0 } / (tau.abs() + (1.0 + tau * tau).sqrt());
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                a[p][p] -= t * off;
                a[q][q] += t * off;
                a[p][q] = 0.0;
                a[q][p] = 0.0;
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
        if maximum <= 1e-13 {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("eigensolver did not converge".into());
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]).then(i.cmp(&j)));
    let values: Vec<f64> = order.iter().map(|&i| a[i][i]).collect();
    if values[1] - values[0] <= 1e-8 || (n > 2 && values[2] - values[1] <= 1e-8) {
        return Err("nonunique Fiedler eigenvalue".into());
    }
    let mut v: Vec<f64> = u.iter().map(|row| row[order[1]]).collect();
    let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    for x in &mut v {
        *x /= norm;
    }
    let residual = l
        .iter()
        .zip(&v)
        .map(|(r, x)| (r.iter().zip(&v).map(|(a, b)| a * b).sum::<f64>() - values[1] * x).powi(2))
        .sum::<f64>()
        .sqrt();
    if !residual.is_finite() || residual > 1e-9 {
        return Err("eigen residual".into());
    }
    if v.iter().any(|x| !x.is_finite() || x.abs() <= 1e-10) {
        return Err("unstable zero cut".into());
    }
    let mut pivot = 0;
    for i in 1..n {
        if v[i].abs() > v[pivot].abs() {
            pivot = i;
        }
    }
    if v[pivot] < 0.0 {
        for x in &mut v {
            *x = -*x;
        }
    }
    let groups = vec![
        (0..n).filter(|&i| v[i] < 0.0).collect::<Vec<_>>(),
        (0..n).filter(|&i| v[i] >= 0.0).collect::<Vec<_>>(),
    ];
    if groups.iter().any(Vec::is_empty) {
        return Err("empty sign partition".into());
    }
    Ok(Split {
        groups,
        eigenvalues: values,
        fiedler: v,
        residual,
    })
}
