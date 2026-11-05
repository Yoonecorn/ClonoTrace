//! Strict declared-matrix schema and transactional CLI/library workflow.
use clonotrace_airr::Result;
use serde::Deserialize;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::HashSet,
    fs::File,
    io::{Read, Write},
    path::PathBuf,
};
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Scope {
    pub donor_id: String,
    pub sample_id: String,
    pub locus: String,
}
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MatrixInput {
    pub schema: String,
    pub kind: String,
    pub ids: Vec<String>,
    pub scope: Scope,
    pub provenance: String,
    pub matrix: Vec<Vec<f64>>,
}
#[derive(Clone, Debug)]
pub struct Config {
    pub input: PathBuf,
    pub out: PathBuf,
    pub method: super::Method,
    pub linkage: Option<String>,
    pub threshold: Option<f64>,
    pub clusters: Option<usize>,
    pub clusters_per_component: Option<usize>,
    pub max_iterations: Option<usize>,
    pub max_rows: usize,
    pub max_input_bytes: u64,
}
fn field(s: &str) -> bool {
    !s.trim().is_empty() && !s.bytes().any(|b| b"\t\r\n\0".contains(&b))
}
impl MatrixInput {
    /// Reorders both matrix axes with stable record IDs; malformed dimensions fail first.
    pub fn canonicalize(&mut self, limit: usize) -> Result<()> {
        if self.schema != "CT-matrix-1"
            || !matches!(self.kind.as_str(), "distance" | "affinity")
            || self.provenance.trim().is_empty()
        {
            return Err("invalid matrix schema/kind/provenance".into());
        }
        let n = self.ids.len();
        if !(1..=512).contains(&limit)
            || n == 0
            || n > limit
            || self.matrix.len() != n
            || self.matrix.iter().any(|r| r.len() != n)
        {
            return Err("invalid matrix row bound or dimensions".into());
        }
        if self.ids.iter().any(|s| !field(s))
            || self.ids.iter().collect::<HashSet<_>>().len() != n
            || ![
                &self.scope.donor_id,
                &self.scope.sample_id,
                &self.scope.locus,
            ]
            .iter()
            .all(|s| field(s))
            || !matches!(self.scope.locus.as_str(), "IGH" | "IGK" | "IGL")
        {
            return Err("invalid matrix IDs or explicit B-chain scope".into());
        }
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&i, &j| self.ids[i].cmp(&self.ids[j]));
        self.matrix = order
            .iter()
            .map(|&i| order.iter().map(|&j| self.matrix[i][j]).collect())
            .collect();
        self.ids = order.into_iter().map(|i| self.ids[i].clone()).collect();
        Ok(())
    }
}
pub fn run(a: Config) -> Result<Value> {
    if std::fs::symlink_metadata(&a.out).is_ok() {
        return Err("output already exists; replacement is prohibited".into());
    }
    if a.max_input_bytes == 0 || a.max_input_bytes > 256_000_000 {
        return Err("matrix byte limit must be in 1..=256000000".into());
    }
    let mut raw = Vec::new();
    File::open(&a.input)
        .map_err(|e| e.to_string())?
        .take(a.max_input_bytes + 1)
        .read_to_end(&mut raw)
        .map_err(|e| e.to_string())?;
    if raw.len() as u64 > a.max_input_bytes {
        return Err("matrix max-input-bytes exceeded".into());
    }
    let digest = format!("{:x}", Sha256::digest(&raw));
    let mut input: MatrixInput = serde_json::from_slice(&raw).map_err(|e| e.to_string())?;
    input.canonicalize(a.max_rows)?;
    let (groups, diagnostics) = match a.method {
        super::Method::HierarchicalV2 => {
            if input.kind != "distance"
                || a.threshold.is_some() == a.clusters.is_some()
                || a.clusters == Some(0)
                || a.clusters_per_component.is_some()
                || a.max_iterations.is_some()
            {
                return Err(
                    "matrix hierarchy requires distance kind, linkage and threshold XOR clusters"
                        .into(),
                );
            }
            let link = clonotrace_cluster::v2::Linkage::parse(
                a.linkage
                    .as_deref()
                    .ok_or("matrix hierarchy requires linkage")?,
            )?;
            let tree = clonotrace_cluster::v2::build(&input.matrix, link, a.max_rows)?;
            let groups = if let Some(t) = a.threshold {
                tree.cut_threshold(t)?
            } else {
                tree.cut_count(a.clusters.unwrap())?
            };
            let merges: Vec<_> = tree
                .merges
                .iter()
                .map(|m| json!({"left":m.left,"right":m.right,"height":m.height,"size":m.size}))
                .collect();
            (groups, json!({"merges":merges}))
        }
        super::Method::SpectralV2 => {
            if input.kind != "affinity"
                || a.clusters_per_component.is_none_or(|n| n == 0)
                || a.threshold.is_some()
                || a.clusters.is_some()
                || a.linkage.is_some()
            {
                return Err(
                    "matrix spectral requires affinity kind and positive clusters-per-component"
                        .into(),
                );
            }
            let result = clonotrace_spectral::v2::cluster_affinity(
                &input.matrix,
                a.clusters_per_component.unwrap(),
                a.max_rows,
                a.max_iterations.unwrap_or(1000),
            )?;
            let components:Vec<_>=result.components.iter().map(|c|json!({"members":c.members,"effective_clusters":c.effective_clusters,"eigenvalues":c.eigenvalues,"eigen_blocks":c.eigen_blocks,"iterations":c.iterations,"inertia":c.inertia,"residual":c.residual,"projector":c.projector})).collect();
            (
                result.groups,
                json!({"components":components,"zero_edges":result.zero_edges}),
            )
        }
        _ => return Err("matrix workflow supports hierarchical-v2 or spectral-v2".into()),
    };
    let mut labels = vec![None; input.ids.len()];
    for (c, g) in groups.iter().enumerate() {
        for &i in g {
            if i >= labels.len() || labels[i].replace(c + 1).is_some() {
                return Err("matrix output membership violation".into());
            }
        }
    }
    if labels.iter().any(Option::is_none) {
        return Err("matrix row conservation failed".into());
    }
    let mut table = String::from(
        "sequence_id\tdonor_id\tsample_id\tlocus\tclone_id\tclonotrace_status\tclonotrace_reason\n",
    );
    for (i, id) in input.ids.iter().enumerate() {
        table.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\tretained\t\n",
            id,
            input.scope.donor_id,
            input.scope.sample_id,
            input.scope.locus,
            labels[i].unwrap()
        ));
    }
    let manifest = json!({"program":"clonotrace","contract":"CT-complete-2","manifest_schema":"CT-run-2","output_schema":"CT-matrix-assignment-1","source_schema":input.schema,"source_kind":input.kind,"scope":{"donor_id":input.scope.donor_id,"sample_id":input.scope.sample_id,"locus":input.scope.locus},"provenance":input.provenance,"input":a.input,"input_sha256":digest,"method":a.method.name(),"execution":"scalar matrix work unit","input_rows":input.ids.len(),"retained":input.ids.len(),"excluded":0,"failed":0,"clones":groups.len(),"configuration":{"linkage":a.linkage,"threshold":a.threshold,"clusters":a.clusters,"clusters_per_component":a.clusters_per_component},"limits":{"max_rows":a.max_rows,"max_input_bytes":a.max_input_bytes,"max_iterations":a.max_iterations.unwrap_or(1000)}});
    let mut recheck = Vec::new();
    File::open(&a.input)
        .map_err(|e| e.to_string())?
        .take(a.max_input_bytes + 1)
        .read_to_end(&mut recheck)
        .map_err(|e| e.to_string())?;
    if recheck.len() as u64 > a.max_input_bytes
        || format!("{:x}", Sha256::digest(&recheck)) != digest
    {
        return Err("matrix input changed before publication".into());
    }
    clonotrace_report::transaction(&a.out, |folder| {
        for (name, bytes) in [
            ("assignments.tsv", table.into_bytes()),
            (
                "diagnostics.json",
                serde_json::to_vec_pretty(&diagnostics).map_err(|e| e.to_string())?,
            ),
            (
                "manifest.json",
                serde_json::to_vec_pretty(&manifest).map_err(|e| e.to_string())?,
            ),
        ] {
            let mut file = File::create(folder.join(name)).map_err(|e| e.to_string())?;
            file.write_all(&bytes).map_err(|e| e.to_string())?;
            file.sync_all().map_err(|e| e.to_string())?;
        }
        Ok(())
    })?;
    Ok(
        json!({"status":"PASS","input_rows":input.ids.len(),"retained":input.ids.len(),"excluded":0,"failed":0,"clones":groups.len()}),
    )
}
