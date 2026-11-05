//! Separate observable grouping models. SCOPer work units are not internal-group claims.
use clonotrace_airr::{Record, Result};
use std::collections::{BTreeMap, BTreeSet, HashMap};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    Changeo,
    Scoper,
}
impl Method {
    pub fn name(self) -> &'static str {
        match self {
            Self::Changeo => "changeo",
            Self::Scoper => "scoper",
        }
    }
}
#[derive(Clone, Debug)]
pub struct UnionFind {
    parent: Vec<usize>,
}
impl UnionFind {
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }
    pub fn root(&mut self, mut a: usize) -> usize {
        while self.parent[a] != a {
            self.parent[a] = self.parent[self.parent[a]];
            a = self.parent[a];
        }
        a
    }
    pub fn join(&mut self, a: usize, b: usize) {
        let a = self.root(a);
        let b = self.root(b);
        if a != b {
            self.parent[b] = a;
        }
    }
    pub fn groups(mut self) -> Vec<Vec<usize>> {
        let mut g = BTreeMap::new();
        for i in 0..self.parent.len() {
            g.entry(self.root(i)).or_insert_with(Vec::new).push(i);
        }
        g.into_values().collect()
    }
}
fn set(x: &[String], first: bool) -> BTreeSet<String> {
    x.iter()
        .take(if first { 1 } else { x.len() })
        .cloned()
        .collect()
}
struct Group {
    v: BTreeSet<String>,
    j: BTreeSet<String>,
    members: Vec<usize>,
}
pub fn build(
    records: &[Record],
    method: Method,
    first: bool,
    max_pairs: u64,
) -> Result<Vec<Vec<usize>>> {
    let mut buckets = BTreeMap::new();
    for (i, r) in records.iter().enumerate() {
        let key = (
            &r.scope,
            if method == Method::Changeo {
                r.junction.len()
            } else {
                0
            },
        );
        buckets.entry(key).or_insert_with(Vec::new).push(i);
    }
    let mut bound = 0u64;
    for b in buckets.values() {
        let n = b.len() as u64;
        bound = bound
            .checked_add(n.saturating_mul(n.saturating_sub(1)) / 2)
            .ok_or("pair bound overflow")?;
    }
    if bound > max_pairs {
        return Err(format!("max-pairs exceeded: {bound} > {max_pairs}"));
    }
    let mut output = Vec::new();
    for ids in buckets.into_values() {
        if method == Method::Scoper {
            let mut uf = UnionFind::new(ids.len());
            // Keep the index linear in compact annotation size. Very large call
            // products use the exact scalar overlap path without changing acceptance.
            let linear_budget = ids
                .iter()
                .map(|&i| records[i].v.len().saturating_add(records[i].j.len()))
                .fold(0usize, usize::saturating_add)
                .saturating_mul(4);
            let products = ids.iter().try_fold(0usize, |sum, &i| {
                let r = &records[i];
                let v = if first { 1 } else { r.v.len() };
                let j = if first { 1 } else { r.j.len() };
                sum.checked_add(v.checked_mul(j)?)
            });
            if products.is_some_and(|n| n <= linear_budget) {
                let mut pairs: HashMap<(&str, &str), usize> = HashMap::new();
                for (local, &i) in ids.iter().enumerate() {
                    let r = &records[i];
                    for v in r.v.iter().take(if first { 1 } else { r.v.len() }) {
                        for j in r.j.iter().take(if first { 1 } else { r.j.len() }) {
                            if let Some(&previous) = pairs.get(&(v.as_str(), j.as_str())) {
                                uf.join(local, previous);
                            } else {
                                pairs.insert((v.as_str(), j.as_str()), local);
                            }
                        }
                    }
                }
            } else {
                let v: Vec<_> = ids.iter().map(|&i| set(&records[i].v, first)).collect();
                let j: Vec<_> = ids.iter().map(|&i| set(&records[i].j, first)).collect();
                for a in 0..ids.len() {
                    for b in a + 1..ids.len() {
                        if !v[a].is_disjoint(&v[b]) && !j[a].is_disjoint(&j[b]) {
                            uf.join(a, b);
                        }
                    }
                }
            }
            let mut groups = uf.groups();
            groups.sort_unstable_by_key(|g| g[0]);
            output.extend(
                groups
                    .into_iter()
                    .map(|g| g.into_iter().map(|i| ids[i]).collect::<Vec<_>>()),
            );
        } else if first
            || ids.iter().all(|&i| {
                let r = &records[i];
                r.v.iter().all(|v| v == &r.v[0]) && r.j.iter().all(|j| j == &r.j[0])
            })
        {
            // In this exact profile, sets never expand. Preserve native last-use order.
            let mut by_key: HashMap<(&str, &str), Vec<usize>> = HashMap::new();
            for i in ids {
                let r = &records[i];
                by_key.entry((&r.v[0], &r.j[0])).or_default().push(i);
            }
            let mut groups: Vec<_> = by_key.into_values().collect();
            groups.sort_unstable_by_key(|g| *g.last().unwrap());
            output.extend(groups);
        } else {
            let mut groups: Vec<Group> = Vec::new();
            for i in ids {
                let mut v = set(&records[i].v, first);
                let mut j = set(&records[i].j, first);
                let mut members = vec![i];
                let mut matched = Vec::new();
                for (k, g) in groups.iter().enumerate() {
                    if !j.is_disjoint(&g.j) {
                        if !first {
                            j.extend(g.j.iter().cloned());
                        }
                        if !v.is_disjoint(&g.v) {
                            matched.push(k);
                            members.extend_from_slice(&g.members);
                            if !first {
                                v.extend(g.v.iter().cloned());
                            }
                        }
                    }
                }
                for k in matched.into_iter().rev() {
                    groups.remove(k);
                }
                groups.push(Group { v, j, members });
            }
            output.extend(groups.into_iter().map(|g| g.members));
        }
    }
    Ok(output)
}
