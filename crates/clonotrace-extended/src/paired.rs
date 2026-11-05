//! Deterministic heavy-first threshold inference with exact light refinement.
use std::collections::{BTreeMap, BTreeSet};
#[derive(Clone, Debug)]
pub struct ChainInput {
    pub sequence_id: String,
    pub donor_id: String,
    pub sample_id: String,
    pub cell_id: String,
    pub locus: String,
    pub v_call: String,
    pub j_call: String,
    pub junction: Vec<u8>,
    /// A prior parser exclusion; never salvage other rows of this cell.
    pub invalid: bool,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissingLightPolicy {
    Exclude,
    HeavyOnlySeparate,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MultiChainPolicy {
    Exclude,
    CollapseIdentical,
}
#[derive(Clone, Copy, Debug)]
pub struct Config {
    pub threshold: f64,
    pub missing_light: MissingLightPolicy,
    pub multiple: MultiChainPolicy,
    pub max_pairs: u64,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CellDiagnostic {
    pub donor_id: String,
    pub sample_id: String,
    pub cell_id: String,
    pub rows: Vec<usize>,
    pub reason: Option<&'static str>,
    pub heavy_only: bool,
    pub collapsed_identical_rows: usize,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Inference {
    pub groups: Vec<Vec<usize>>,
    pub work_units: Vec<Vec<usize>>,
    pub exclusions: Vec<Option<&'static str>>,
    pub cells: Vec<CellDiagnostic>,
    pub retained_cells: usize,
    pub paired_cells: usize,
    pub heavy_only_cells: usize,
    pub heavy_components: usize,
    pub collapsed_identical_rows: usize,
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Signature {
    locus: String,
    v: String,
    j: String,
    seq: Vec<u8>,
}
#[derive(Clone, Debug)]
struct Cell {
    rows: Vec<usize>,
    heavy: Signature,
    light: Option<Signature>,
}
#[derive(Debug)]
struct Unit {
    cells: Vec<Cell>,
}
fn gene(call: &str, locus: &str, segment: &str) -> Option<String> {
    let (g, a) = call
        .split_once('*')
        .map_or((call, None), |(g, a)| (g, Some(a)));
    let prefix = format!("{locus}{segment}");
    let suffix = g.strip_prefix(&prefix)?;
    if suffix.is_empty()
        || !suffix
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || b"/-".contains(&c))
        || a.is_some_and(|a| a.is_empty() || !a.bytes().all(|c| c.is_ascii_digit()))
    {
        return None;
    }
    Some(g.to_string())
}
fn signature(r: &ChainInput) -> Option<Signature> {
    if r.invalid
        || !matches!(r.locus.as_str(), "IGH" | "IGK" | "IGL")
        || r.junction.is_empty()
        || !r.junction.iter().all(|b| b"ACGT".contains(b))
    {
        return None;
    }
    Some(Signature {
        locus: r.locus.clone(),
        v: gene(&r.v_call, &r.locus, "V")?,
        j: gene(&r.j_call, &r.locus, "J")?,
        seq: r.junction.clone(),
    })
}
fn check(config: Config) -> Result<(), String> {
    if !config.threshold.is_finite() || !(0.0..=1.0).contains(&config.threshold) {
        Err("threshold must be finite decimal in [0,1]".into())
    } else {
        Ok(())
    }
}
fn prepare(rows: &[ChainInput], config: Config) -> Result<(Vec<Unit>, Inference), String> {
    check(config)?;
    let mut ids = BTreeSet::new();
    let mut cells = BTreeMap::new();
    let mut counts = BTreeMap::new();
    for (i, r) in rows.iter().enumerate() {
        if r.sequence_id.trim().is_empty() || !ids.insert(r.sequence_id.as_str()) {
            return Err("blank or duplicate sequence_id".into());
        }
        if r.donor_id.trim().is_empty() || r.sample_id.trim().is_empty() {
            return Err("blank scope boundary".into());
        }
        if r.cell_id.trim().is_empty() {
            return Err("blank cell_id".into());
        }
        cells
            .entry((&r.donor_id, &r.sample_id, &r.cell_id))
            .or_insert_with(Vec::new)
            .push(i);
        *counts.entry((&r.donor_id, &r.sample_id)).or_insert(0u64) += 1;
    }
    let mut pairs = 0u64;
    for &n in counts.values() {
        pairs = pairs
            .checked_add(n.checked_mul(n.saturating_sub(1)).ok_or("pair overflow")? / 2)
            .ok_or("pair overflow")?;
    }
    if pairs > config.max_pairs {
        return Err("max-pairs exceeded".into());
    }
    let mut result = Inference {
        groups: vec![],
        work_units: vec![],
        exclusions: vec![None; rows.len()],
        cells: vec![],
        retained_cells: 0,
        paired_cells: 0,
        heavy_only_cells: 0,
        heavy_components: 0,
        collapsed_identical_rows: 0,
    };
    let mut buckets = BTreeMap::new();
    for ((donor, sample, cell), mut indices) in cells {
        indices.sort_by(|&i, &j| rows[i].sequence_id.cmp(&rows[j].sequence_id));
        let mut h = BTreeSet::new();
        let mut l = BTreeSet::new();
        let mut nh = 0;
        let mut nl = 0;
        let mut invalid = false;
        for &i in &indices {
            match signature(&rows[i]) {
                None => invalid = true,
                Some(s) => {
                    if s.locus == "IGH" {
                        nh += 1;
                        h.insert(s);
                    } else {
                        nl += 1;
                        l.insert(s);
                    }
                }
            }
        }
        let reason = if invalid {
            Some("invalid_chain_in_cell")
        } else if nh == 0 {
            Some("missing_heavy_chain")
        } else if config.multiple == MultiChainPolicy::Exclude && nh > 1 {
            Some("multiple_heavy_chains")
        } else if config.multiple == MultiChainPolicy::Exclude && nl > 1 {
            Some("multiple_light_chains")
        } else if h.len() > 1 {
            Some("conflicting_heavy_chains")
        } else if l.len() > 1 {
            Some("conflicting_light_chains")
        } else if nl == 0 && config.missing_light == MissingLightPolicy::Exclude {
            Some("missing_light_chain")
        } else {
            None
        };
        let heavy_only = reason.is_none() && nl == 0;
        let collapsed = if reason.is_none() {
            indices.len() - 1 - usize::from(nl > 0)
        } else {
            0
        };
        result.cells.push(CellDiagnostic {
            donor_id: donor.clone(),
            sample_id: sample.clone(),
            cell_id: cell.clone(),
            rows: indices.clone(),
            reason,
            heavy_only,
            collapsed_identical_rows: collapsed,
        });
        if let Some(reason) = reason {
            for &i in &indices {
                result.exclusions[i] = Some(reason);
            }
            continue;
        }
        let heavy = h.into_iter().next().unwrap();
        let light = l.into_iter().next();
        let key = (
            donor.clone(),
            sample.clone(),
            heavy.v.clone(),
            heavy.j.clone(),
            heavy.seq.len(),
            heavy_only,
        );
        buckets.entry(key).or_insert_with(Vec::new).push(Cell {
            rows: indices,
            heavy,
            light,
        });
        result.retained_cells += 1;
        result.heavy_only_cells += usize::from(heavy_only);
        result.paired_cells += usize::from(!heavy_only);
        result.collapsed_identical_rows += collapsed;
    }
    let units: Vec<Unit> = buckets.into_values().map(|cells| Unit { cells }).collect();
    result.work_units = units
        .iter()
        .map(|u| {
            u.cells
                .iter()
                .flat_map(|c| c.rows.iter().copied())
                .collect()
        })
        .collect();
    Ok((units, result))
}
fn unit(unit: &Unit, threshold: f64) -> (Vec<Vec<usize>>, usize) {
    let mut parent: Vec<_> = (0..unit.cells.len()).collect();
    fn root(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            p[i] = p[p[i]];
            i = p[i];
        }
        i
    }
    for i in 0..unit.cells.len() {
        for j in i + 1..unit.cells.len() {
            let a = &unit.cells[i].heavy.seq;
            let b = &unit.cells[j].heavy.seq;
            let d = a.iter().zip(b).filter(|(x, y)| x != y).count() as f64 / a.len() as f64;
            if d <= threshold {
                let a = root(&mut parent, i);
                let b = root(&mut parent, j);
                parent[b] = a;
            }
        }
    }
    let mut components = BTreeSet::new();
    let mut groups = BTreeMap::new();
    for (i, c) in unit.cells.iter().enumerate() {
        let r = root(&mut parent, i);
        components.insert(r);
        groups
            .entry((r, &c.light))
            .or_insert_with(Vec::new)
            .extend(c.rows.iter().copied());
    }
    (groups.into_values().collect(), components.len())
}
/// Scalar scientific authority. Row-index outputs map directly to input rows.
pub fn infer(rows: &[ChainInput], config: Config) -> Result<Inference, String> {
    infer_threads(rows, config, 1)
}
/// Independent unit parallelism; no shared numerical reductions or hidden fallback.
pub fn infer_threads(
    rows: &[ChainInput],
    config: Config,
    threads: usize,
) -> Result<Inference, String> {
    if !(1..=1024).contains(&threads) {
        return Err("threads must be in 1..=1024".into());
    }
    let (units, mut result) = prepare(rows, config)?;
    let outputs = if threads == 1 || units.len() < 2 {
        units
            .iter()
            .map(|u| unit(u, config.threshold))
            .collect::<Vec<_>>()
    } else {
        let workers = threads.min(units.len());
        let chunk = units.len().div_ceil(workers);
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for chunk in units.chunks(chunk) {
                handles.push(
                    std::thread::Builder::new()
                        .spawn_scoped(scope, move || {
                            chunk
                                .iter()
                                .map(|u| unit(u, config.threshold))
                                .collect::<Vec<_>>()
                        })
                        .map_err(|e| format!("paired worker spawn failed: {e}"))?,
                );
            }
            let mut outputs = Vec::new();
            for h in handles {
                outputs.extend(h.join().map_err(|_| "paired worker panicked".to_string())?);
            }
            Ok::<_, String>(outputs)
        })?
    };
    for (groups, components) in outputs {
        result.groups.extend(groups);
        result.heavy_components += components;
    }
    for g in &mut result.groups {
        g.sort_by(|&i, &j| rows[i].sequence_id.cmp(&rows[j].sequence_id));
    }
    result
        .groups
        .sort_by(|a, b| rows[a[0]].sequence_id.cmp(&rows[b[0]].sequence_id));
    Ok(result)
}
#[cfg(test)]
mod tests {
    use super::*;
    fn row(id: &str, cell: &str, locus: &str, seq: &str) -> ChainInput {
        ChainInput {
            sequence_id: id.into(),
            donor_id: "d".into(),
            sample_id: "s".into(),
            cell_id: cell.into(),
            locus: locus.into(),
            v_call: format!("{locus}V1*01"),
            j_call: format!("{locus}J1*01"),
            junction: seq.as_bytes().to_vec(),
            invalid: false,
        }
    }
    fn cfg() -> Config {
        Config {
            threshold: 0.25,
            missing_light: MissingLightPolicy::HeavyOnlySeparate,
            multiple: MultiChainPolicy::CollapseIdentical,
            max_pairs: 10000,
        }
    }
    #[test]
    fn known_light_bridge_then_refine() {
        let r = vec![
            row("a_h", "a", "IGH", "AAAA"),
            row("a_l", "a", "IGK", "CCCC"),
            row("b_h", "b", "IGH", "AAAC"),
            row("b_l", "b", "IGK", "GGGG"),
            row("c_h", "c", "IGH", "AACC"),
            row("c_l", "c", "IGK", "CCCC"),
        ];
        assert_eq!(
            infer(&r, cfg()).unwrap().groups,
            vec![vec![0, 1, 4, 5], vec![2, 3]]
        );
    }
    #[test]
    fn missing_light_cannot_bridge() {
        let r = vec![
            row("a_h", "a", "IGH", "AAAA"),
            row("a_l", "a", "IGK", "CCCC"),
            row("b_h", "b", "IGH", "AAAC"),
            row("c_h", "c", "IGH", "AACC"),
            row("c_l", "c", "IGK", "CCCC"),
        ];
        assert_eq!(
            infer(&r, cfg()).unwrap().groups,
            vec![vec![0, 1], vec![2], vec![3, 4]]
        );
    }
    #[test]
    fn collapse_retains_original_rows() {
        let r = vec![
            row("h1", "a", "IGH", "AAAA"),
            row("h2", "a", "IGH", "AAAA"),
            row("l", "a", "IGK", "CCCC"),
        ];
        let x = infer(&r, cfg()).unwrap();
        assert_eq!(x.groups, vec![vec![0, 1, 2]]);
        assert_eq!(x.collapsed_identical_rows, 1);
        let mut c = cfg();
        c.multiple = MultiChainPolicy::Exclude;
        assert_eq!(
            infer(&r, c).unwrap().exclusions,
            vec![Some("multiple_heavy_chains"); 3]
        );
    }
    #[test]
    fn scopes_threads_and_structural() {
        let mut r = vec![];
        for d in 0..16 {
            let mut a = row(&format!("h{d}"), "a", "IGH", "AAAA");
            a.donor_id = d.to_string();
            let mut b = row(&format!("l{d}"), "a", "IGK", "CCCC");
            b.donor_id = d.to_string();
            r.extend([a, b]);
        }
        assert_eq!(
            infer_threads(&r, cfg(), 8).unwrap(),
            infer(&r, cfg()).unwrap()
        );
        r[1].sequence_id = r[0].sequence_id.clone();
        assert!(infer(&r, cfg()).is_err());
    }
}
