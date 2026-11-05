//! Reusable versioned inference workflows; native CT-0.1-native-1 remains separate.
pub mod amino;
pub mod matrix;
pub mod paired;
mod workflow_v2;
use clonotrace_airr::{line, open, Input, Record, Result, Row};
use clonotrace_group::UnionFind;
use rayon::prelude::*;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    path::PathBuf,
    time::Instant,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    HierarchicalV1,
    SpectralV1,
    PairedExactV1,
    AaIdentityV1,
    AaHammingV1,
    AutoGapV1,
    GermlineMaskV1,
    HierarchicalV2,
    SpectralV2,
    PairedThresholdV2,
    AaIdentityV2,
    AaHammingV2,
    GermlineV2,
    AutoGapStabilityV2,
}
impl Method {
    pub fn name(self) -> &'static str {
        match self {
            Self::HierarchicalV1 => "hierarchical-v1",
            Self::SpectralV1 => "spectral-v1",
            Self::PairedExactV1 => "paired-exact-v1",
            Self::AaIdentityV1 => "aa-identity-v1",
            Self::AaHammingV1 => "aa-hamming-v1",
            Self::AutoGapV1 => "auto-gap-v1",
            Self::GermlineMaskV1 => "germline-mask-v1",
            Self::HierarchicalV2 => "hierarchical-v2",
            Self::SpectralV2 => "spectral-v2",
            Self::PairedThresholdV2 => "paired-threshold-v2",
            Self::AaIdentityV2 => "aa-identity-v2",
            Self::AaHammingV2 => "aa-hamming-v2",
            Self::GermlineV2 => "germline-v2",
            Self::AutoGapStabilityV2 => "auto-gap-stability-v2",
        }
    }
    fn aa(self) -> bool {
        matches!(
            self,
            Self::AaIdentityV1 | Self::AaHammingV1 | Self::AaIdentityV2 | Self::AaHammingV2
        )
    }
}
#[derive(Clone, Debug)]
pub struct Config {
    pub input: PathBuf,
    pub out: PathBuf,
    pub method: Method,
    pub threshold: Option<String>,
    pub linkage: Option<String>,
    pub sigma: Option<f64>,
    /// AA source: direct junction_aa or translate junction with explicit frame.
    pub aa_input: Option<String>,
    pub frame: Option<usize>,
    pub reference_fasta: Option<PathBuf>,
    pub reference_id: Option<String>,
    pub reference_version: Option<String>,
    pub reference_source: Option<String>,
    pub threads: usize,
    pub max_rows: usize,
    pub max_junction_length: usize,
    pub max_pairs: u64,
    pub alphabet: Option<String>,
    pub ambiguity: Option<String>,
    pub clusters: Option<usize>,
    pub clusters_per_component: Option<usize>,
    pub affinity_cutoff: Option<f64>,
    pub missing_light: Option<String>,
    pub multi_chain: Option<String>,
    pub germline_mode: Option<String>,
    pub max_group_rows: Option<usize>,
    pub max_iterations: Option<usize>,
    pub max_alignment_length: Option<usize>,
    pub max_distinct: Option<usize>,
    pub max_distance_evaluations: Option<u64>,
    pub max_reference_bytes: Option<u64>,
}
impl Method {
    pub fn is_v2(self) -> bool {
        matches!(
            self,
            Self::HierarchicalV2
                | Self::SpectralV2
                | Self::PairedThresholdV2
                | Self::AaIdentityV2
                | Self::AaHammingV2
                | Self::GermlineV2
                | Self::AutoGapStabilityV2
        )
    }
    pub fn is_paired(self) -> bool {
        matches!(self, Self::PairedExactV1 | Self::PairedThresholdV2)
    }
    pub fn is_germline(self) -> bool {
        matches!(self, Self::GermlineMaskV1 | Self::GermlineV2)
    }
}
impl std::str::FromStr for Method {
    type Err = String;
    fn from_str(value: &str) -> Result<Self> {
        [
            Self::HierarchicalV1,
            Self::SpectralV1,
            Self::PairedExactV1,
            Self::AaIdentityV1,
            Self::AaHammingV1,
            Self::AutoGapV1,
            Self::GermlineMaskV1,
            Self::HierarchicalV2,
            Self::SpectralV2,
            Self::PairedThresholdV2,
            Self::AaIdentityV2,
            Self::AaHammingV2,
            Self::GermlineV2,
            Self::AutoGapStabilityV2,
        ]
        .into_iter()
        .find(|m| m.name() == value)
        .ok_or_else(|| "unsupported method".into())
    }
}
impl Config {
    pub fn new(input: impl Into<PathBuf>, out: impl Into<PathBuf>, method: Method) -> Self {
        Self {
            input: input.into(),
            out: out.into(),
            method,
            threshold: None,
            linkage: None,
            sigma: None,
            aa_input: None,
            frame: None,
            reference_fasta: None,
            reference_id: None,
            reference_version: None,
            reference_source: None,
            threads: 1,
            max_rows: 1_000_000,
            max_junction_length: 10_000,
            max_pairs: 100_000_000,
            alphabet: None,
            ambiguity: None,
            clusters: None,
            clusters_per_component: None,
            affinity_cutoff: None,
            missing_light: None,
            multi_chain: None,
            germline_mode: None,
            max_group_rows: None,
            max_iterations: None,
            max_alignment_length: None,
            max_distinct: None,
            max_distance_evaluations: None,
            max_reference_bytes: None,
        }
    }
    pub fn analyzes_aa(&self) -> bool {
        self.method.aa()
            || matches!(self.method, Method::HierarchicalV2 | Method::SpectralV2)
                && self.alphabet.as_deref() == Some("aa")
    }
    pub fn validate(&self) -> Result<()> {
        configuration(self).map(|_| ())
    }
}
fn decimal(s: &str) -> Result<f64> {
    if s.is_empty()
        || !s
            .bytes()
            .all(|b| b.is_ascii_digit() || b".-+eE".contains(&b))
    {
        return Err("threshold must be finite decimal in [0,1]".into());
    }
    let value: f64 = s.parse().map_err(|_| "invalid threshold")?;
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err("threshold must be finite decimal in [0,1]".into());
    }
    Ok(value)
}
fn configuration(a: &Config) -> Result<Option<f64>> {
    if a.method.is_v2() {
        workflow_v2::validate(a)
    } else {
        if a.alphabet.is_some()
            || a.ambiguity.is_some()
            || a.clusters.is_some()
            || a.clusters_per_component.is_some()
            || a.affinity_cutoff.is_some()
            || a.missing_light.is_some()
            || a.multi_chain.is_some()
            || a.germline_mode.is_some()
            || a.max_group_rows.is_some()
            || a.max_iterations.is_some()
            || a.max_alignment_length.is_some()
            || a.max_distinct.is_some()
            || a.max_distance_evaluations.is_some()
            || a.max_reference_bytes.is_some()
        {
            return Err("v2 options require an explicit v2 method".into());
        }
        legacy_configuration(a)
    }
}
fn legacy_configuration(a: &Config) -> Result<Option<f64>> {
    if !(1..=1024).contains(&a.threads) {
        return Err("threads must be in 1..=1024".into());
    }
    let t = a.threshold.as_deref().map(decimal).transpose()?;
    if matches!(a.method, Method::HierarchicalV1 | Method::AaHammingV1) && t.is_none() {
        return Err("method requires explicit --threshold".into());
    }
    if !matches!(
        a.method,
        Method::HierarchicalV1 | Method::AaHammingV1 | Method::AutoGapV1
    ) && t.is_some()
    {
        return Err("threshold is unsupported for selected method".into());
    }
    if a.method == Method::HierarchicalV1 {
        if !a
            .linkage
            .as_deref()
            .is_some_and(|x| matches!(x, "single" | "average" | "complete"))
        {
            return Err("hierarchical-v1 requires --linkage single/average/complete".into());
        }
    } else if a.linkage.is_some() {
        return Err("linkage only supports hierarchical-v1".into());
    }
    if a.method == Method::SpectralV1 {
        if !a.sigma.is_some_and(|s| s.is_finite() && s > 0.0) {
            return Err("spectral-v1 requires finite --sigma > 0".into());
        }
    } else if a.sigma.is_some() {
        return Err("sigma only supports spectral-v1".into());
    }
    if a.analyzes_aa() {
        match a.aa_input.as_deref() {
            Some("direct") if a.frame.is_none() => {}
            Some("translate") if a.frame.is_some_and(|f| f <= 2) => {}
            _ => {
                return Err(
                    "AA requires --aa-input direct or --aa-input translate --frame 0/1/2".into(),
                )
            }
        }
    } else if a.aa_input.is_some() || a.frame.is_some() {
        return Err("AA settings only support AA methods".into());
    }
    let refs = [
        a.reference_id.as_deref(),
        a.reference_version.as_deref(),
        a.reference_source.as_deref(),
    ];
    if a.method == Method::GermlineMaskV1 {
        if a.reference_fasta.is_none() || refs.iter().any(|x| x.is_none_or(|s| s.trim().is_empty()))
        {
            return Err("germline requires reference-fasta/id/version/source".into());
        }
    } else if a.reference_fasta.is_some() || refs.iter().any(Option::is_some) {
        return Err("reference settings only support germline-mask-v1".into());
    }
    Ok(t)
}
fn gene(call: &str, locus: &str, kind: &str) -> Option<String> {
    let genes = clonotrace_airr::calls(call, &format!("{locus}{kind}"))?;
    if genes.len() != 1 {
        return None;
    }
    genes.into_iter().next()
}
fn translate(s: &[u8], frame: usize) -> std::result::Result<Vec<u8>, &'static str> {
    if s.is_empty() || !s.iter().all(|b| b"ACGT".contains(b)) {
        return Err("invalid_nucleotide_sequence");
    }
    if frame > 2 || s.len() <= frame || !(s.len() - frame).is_multiple_of(3) {
        return Err("empty_or_partial_codon");
    }
    // NCBI standard genetic code in T,C,A,G lexical codon order; independently frozen in fixtures.
    let code = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    let base = |b: u8| b"TCAG".iter().position(|&x| x == b).unwrap();
    let mut aa = Vec::new();
    for c in s[frame..].as_chunks::<3>().0 {
        let x = code[base(c[0]) * 16 + base(c[1]) * 4 + base(c[2])];
        if x == b'*' {
            return Err("stop_codon");
        }
        aa.push(x);
    }
    Ok(aa)
}
#[derive(Clone, Debug)]
struct GermlineInput {
    v_call: String,
    j_call: String,
    d_call: String,
    v_cigar: String,
    j_cigar: String,
    d_cigar: String,
}
struct Loaded {
    input: Input,
    cells: Vec<(String, String, String)>,
    projections: Vec<Option<Vec<u8>>>,
    chains: Vec<paired::ChainInput>,
    germline_inputs: Vec<Option<GermlineInput>>,
    reconstructions: Vec<Option<clonotrace_germline::v2::Reconstruction>>,
    cell_report: Option<String>,
}
fn load(a: &Config, refs: Option<&BTreeMap<String, Vec<u8>>>) -> Result<Loaded> {
    let mut stream = open(&a.input)?;
    let mut buf = Vec::new();
    if !line(&mut *stream.reader, &mut buf)? {
        return Err("missing header".into());
    }
    let header = String::from_utf8(buf.clone()).unwrap();
    let fields: Vec<&str> = header.split('\t').collect();
    if fields.iter().any(|s| s.trim().is_empty())
        || fields.iter().copied().collect::<HashSet<_>>().len() != fields.len()
    {
        return Err("blank or duplicate header".into());
    }
    if fields.iter().any(|s| {
        [
            "clone_id",
            "CLONE",
            "clonotrace_status",
            "clonotrace_reason",
            "clonotrace_germline",
        ]
        .contains(s)
    }) {
        return Err("reserved output field".into());
    }
    let paired = a.method.is_paired();
    let germline = a.method.is_germline();
    if !paired && !germline && fields.iter().any(|s| ["cell_id", "CELL"].contains(s)) {
        return Err("bulk extension requires explicit source_cell_id projection".into());
    }
    let mut required = vec![
        "sequence_id",
        "v_call",
        "j_call",
        "donor_id",
        "sample_id",
        "locus",
    ];
    if germline {
        required.extend(["sequence", "v_cigar", "j_cigar"]);
    } else if a.analyzes_aa() && a.aa_input.as_deref() == Some("direct") {
        required.push("junction_aa");
    } else {
        required.push("junction");
    }
    if paired {
        required.push("cell_id");
    }
    if a.method.is_v2()
        && fields.iter().any(|s| {
            [
                "clonotrace_query_alignment",
                "clonotrace_germline_alignment",
            ]
            .contains(s)
        })
    {
        return Err("reserved v2 output field".into());
    }
    let indices: BTreeMap<&str, usize> = fields.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    for key in required {
        if !indices.contains_key(key) {
            return Err(format!("missing required column {key}"));
        }
    }
    let mut records = Vec::new();
    let mut rows = Vec::new();
    let mut cells = Vec::new();
    let mut projections = Vec::new();
    let mut chains = Vec::new();
    let mut germline_inputs = Vec::new();
    let mut reconstructions = Vec::new();
    let mut ids = HashSet::new();
    while line(&mut *stream.reader, &mut buf)? {
        if rows.len() >= a.max_rows {
            return Err("max-rows exceeded".into());
        }
        let text = std::str::from_utf8(&buf).unwrap();
        let values: Vec<&str> = text.split('\t').collect();
        if values.len() != fields.len() {
            return Err("ragged row".into());
        }
        let get = |key: &str| indices.get(key).map_or("", |&i| values[i]);
        let id = get("sequence_id");
        let donor = get("donor_id");
        let sample = get("sample_id");
        let locus = get("locus");
        if id.trim().is_empty() || !ids.insert(id.to_owned()) {
            return Err("blank or duplicate sequence_id".into());
        }
        if donor.trim().is_empty() || sample.trim().is_empty() {
            return Err("blank scope boundary".into());
        }
        if paired && get("cell_id").trim().is_empty() {
            return Err("blank cell_id".into());
        }
        cells.push((
            donor.to_owned(),
            sample.to_owned(),
            get("cell_id").to_owned(),
        ));
        let v = gene(get("v_call"), locus, "V");
        let j = gene(get("j_call"), locus, "J");
        let source = if germline {
            get("sequence")
        } else if a.analyzes_aa() && a.aa_input.as_deref() == Some("direct") {
            get("junction_aa")
        } else {
            get("junction")
        };
        if source.len() > a.max_junction_length {
            return Err("max-junction-length exceeded".into());
        }
        let mut sequence = source.as_bytes().to_vec();
        let mut projection = None;
        let reason = if (!paired && !germline && locus != "IGH")
            || ((paired || germline) && !matches!(locus, "IGH" | "IGK" | "IGL"))
        {
            Some("outside_locus_scope")
        } else if v.is_none() || j.is_none() {
            Some("invalid_or_ambiguous_VJ_call")
        } else if a.method == Method::GermlineV2 {
            if sequence.is_empty() || !sequence.iter().all(|b| b"ACGT".contains(b)) {
                Some("invalid_sequence")
            } else {
                None
            }
        } else if germline {
            match clonotrace_germline::project(
                source.as_bytes(),
                get("v_call"),
                get("j_call"),
                get("v_cigar"),
                get("j_cigar"),
                refs.unwrap(),
            ) {
                Ok(p) => {
                    projection = Some(p);
                    None
                }
                Err(e) => Some(e),
            }
        } else if a.analyzes_aa() && a.aa_input.as_deref() == Some("direct") {
            if sequence.is_empty() || !sequence.iter().all(|b| b"ACDEFGHIKLMNPQRSTVWY".contains(b))
            {
                Some("invalid_amino_acid_sequence")
            } else {
                None
            }
        } else if a.method.is_v2()
            && a.analyzes_aa()
            && a.ambiguity.as_deref() == Some("resolve-unanimous")
        {
            if indices.contains_key("junction_length")
                && get("junction_length") != source.len().to_string()
            {
                Some("inconsistent_declared_junction_length")
            } else {
                match amino::translate(
                    &sequence,
                    amino::TranslationConfig {
                        frame: a.frame.unwrap(),
                        ambiguity: amino::AmbiguityPolicy::ResolveUnanimous,
                    },
                ) {
                    Ok(aa) => {
                        sequence = aa.into_bytes();
                        None
                    }
                    Err(e) => Some(e),
                }
            }
        } else if source.is_empty() || !source.bytes().all(|b| b"ACGT".contains(&b)) {
            Some("invalid_nucleotide_sequence")
        } else if indices.contains_key("junction_length")
            && get("junction_length") != source.len().to_string()
        {
            Some("inconsistent_declared_junction_length")
        } else if a.analyzes_aa() {
            match translate(&sequence, a.frame.unwrap()) {
                Ok(aa) => {
                    sequence = aa;
                    None
                }
                Err(e) => Some(e),
            }
        } else {
            None
        };
        if a.method == Method::PairedThresholdV2 {
            chains.push(paired::ChainInput {
                sequence_id: id.into(),
                donor_id: donor.into(),
                sample_id: sample.into(),
                cell_id: get("cell_id").into(),
                locus: locus.into(),
                v_call: get("v_call").into(),
                j_call: get("j_call").into(),
                junction: source.as_bytes().to_vec(),
                invalid: reason.is_some(),
            });
        }
        if a.method == Method::GermlineV2 {
            germline_inputs.push(if reason.is_none() {
                Some(GermlineInput {
                    v_call: get("v_call").into(),
                    j_call: get("j_call").into(),
                    d_call: get("d_call").into(),
                    v_cigar: get("v_cigar").into(),
                    j_cigar: get("j_cigar").into(),
                    d_cigar: get("d_cigar").into(),
                })
            } else {
                None
            });
            reconstructions.push(None);
        }
        projections.push(projection);
        if reason.is_none() {
            records.push(Record {
                row: rows.len(),
                id: id.into(),
                scope: (donor.into(), sample.into(), locus.into()),
                v: vec![v.unwrap()],
                j: vec![j.unwrap()],
                junction: sequence,
            });
        }
        rows.push(Row {
            id: id.into(),
            reason,
        });
    }
    Ok(Loaded {
        input: Input {
            header,
            records,
            rows,
            sha256: stream.digest(),
        },
        cells,
        projections,
        chains,
        germline_inputs,
        reconstructions,
        cell_report: None,
    })
}
fn hamming(a: &[u8], b: &[u8]) -> f64 {
    a.iter().zip(b).filter(|(x, y)| x != y).count() as f64 / a.len() as f64
}
fn singles(seqs: &[&[u8]], t: f64) -> Vec<Vec<usize>> {
    let mut uf = UnionFind::new(seqs.len());
    for i in 0..seqs.len() {
        for j in i + 1..seqs.len() {
            if hamming(seqs[i], seqs[j]) <= t {
                uf.join(i, j);
            }
        }
    }
    uf.groups()
}
fn auto(seqs: &[&[u8]]) -> (Option<f64>, Value) {
    let unique: Vec<&[u8]> = seqs
        .iter()
        .copied()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    if unique.len() < 2 {
        return (
            None,
            json!({"reason":"no_distinct_neighbor","distinct_sequences":unique.len()}),
        );
    }
    let nearest: Vec<f64> = unique
        .iter()
        .enumerate()
        .map(|(i, a)| {
            unique
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, b)| hamming(a, b))
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    let mut distances = nearest.clone();
    distances.sort_by(f64::total_cmp);
    distances.dedup();
    if distances.len() < 2 {
        return (
            None,
            json!({"reason":"insufficient_distance_variation","nearest":nearest,"distinct_distances":distances}),
        );
    }
    let gaps: Vec<f64> = distances.windows(2).map(|p| p[1] - p[0]).collect();
    let maximum = gaps.iter().copied().fold(0.0, f64::max);
    if gaps.iter().filter(|&&g| g == maximum).count() != 1 {
        return (
            None,
            json!({"reason":"tied_largest_gap","nearest":nearest,"distinct_distances":distances,"gaps":gaps}),
        );
    }
    let i = gaps.iter().position(|&g| g == maximum).unwrap();
    let t = (distances[i] + distances[i + 1]) / 2.0;
    (
        Some(t),
        json!({"reason":null,"nearest":nearest,"distinct_distances":distances,"gaps":gaps,"estimated_threshold":t}),
    )
}
fn canonical(records: &[Record], groups: &mut [Vec<usize>]) {
    for g in groups.iter_mut() {
        g.sort_by(|&i, &j| records[i].id.cmp(&records[j].id));
    }
    groups.sort_by(|a, b| records[a[0]].id.cmp(&records[b[0]].id));
}
fn work(a: &Config, input: &Input) -> Result<Vec<Vec<usize>>> {
    let mut buckets = BTreeMap::new();
    let mut scopes: BTreeMap<(&_, usize), u64> = BTreeMap::new();
    for (i, r) in input.records.iter().enumerate() {
        buckets
            .entry((&r.scope, &r.v[0], &r.j[0], r.junction.len()))
            .or_insert_with(Vec::new)
            .push(i);
        *scopes.entry((&r.scope, r.junction.len())).or_default() += 1;
    }
    let mut pairs = 0u64;
    for &n in scopes.values() {
        pairs = pairs
            .checked_add(n.checked_mul(n.saturating_sub(1)).ok_or("pair overflow")? / 2)
            .ok_or("pair overflow")?;
    }
    if pairs > a.max_pairs {
        return Err("max-pairs exceeded".into());
    }
    let mut output: Vec<Vec<usize>> = buckets.into_values().collect();
    for g in &mut output {
        g.sort_by(|&i, &j| input.records[i].id.cmp(&input.records[j].id));
    }
    if a.method.is_v2() {
        if matches!(a.method, Method::HierarchicalV2 | Method::SpectralV2)
            && output
                .iter()
                .any(|g| g.len() > a.max_group_rows.unwrap_or(512))
        {
            return Err("v2 dense group bound exceeded".into());
        }
    } else if !a.analyzes_aa() && output.iter().any(|g| g.len() > 128) {
        return Err("group bound exceeded: maximum 128 rows".into());
    }
    Ok(output)
}
fn cluster(
    a: &Config,
    input: &Input,
    ids: &[usize],
    threshold: Option<f64>,
) -> Result<(Vec<Vec<usize>>, Value)> {
    if a.method.is_v2() {
        return workflow_v2::cluster(a, input, ids, threshold);
    }
    let seqs: Vec<&[u8]> = ids
        .iter()
        .map(|&i| input.records[i].junction.as_slice())
        .collect();
    let (groups, diagnostic) = match a.method {
        Method::HierarchicalV1 => (
            clonotrace_cluster::hierarchical(
                &seqs,
                threshold.unwrap(),
                a.linkage.as_deref().unwrap(),
            )?,
            json!({"linkage":a.linkage,"threshold":threshold}),
        ),
        Method::SpectralV1 => {
            let x = clonotrace_spectral::split(&seqs, a.sigma.unwrap())?;
            (
                x.groups,
                json!({"lower_gap":if x.eigenvalues.len()>1 {Some(x.eigenvalues[1]-x.eigenvalues[0])} else {None},"upper_gap":if x.eigenvalues.len()>2 {Some(x.eigenvalues[2]-x.eigenvalues[1])} else {None},"eigenvalues":x.eigenvalues,"fiedler":x.fiedler,"residual":x.residual}),
            )
        }
        Method::AutoGapV1 => {
            let (estimated, mut diag) = auto(&seqs);
            let t = threshold.or(estimated).ok_or_else(|| {
                format!(
                    "auto-gap-v1 refused group {}: {}",
                    input.records[ids[0]].id, diag
                )
            })?;
            diag["override_threshold"] = json!(threshold);
            diag["selected_threshold"] = json!(t);
            (singles(&seqs, t), diag)
        }
        Method::AaIdentityV1 => {
            // Profiled exact identity collapse: avoids all-pair comparisons without changing
            // any grouping, resource acceptance, row weights, labels or diagnostics.
            let mut by_sequence: BTreeMap<&[u8], Vec<usize>> = BTreeMap::new();
            for (i, &sequence) in seqs.iter().enumerate() {
                by_sequence.entry(sequence).or_default().push(i);
            }
            (
                by_sequence.into_values().collect(),
                json!({"threshold":0.0}),
            )
        }
        Method::AaHammingV1 => (
            singles(&seqs, threshold.unwrap()),
            json!({"threshold":threshold.unwrap()}),
        ),
        _ => unreachable!(),
    };
    Ok((
        groups
            .into_iter()
            .map(|g| g.into_iter().map(|i| ids[i]).collect())
            .collect(),
        diagnostic,
    ))
}
fn paired(loaded: &mut Loaded) -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Value) {
    let mut by_cell: BTreeMap<&(String, String, String), Vec<usize>> = BTreeMap::new();
    for (row, cell) in loaded.cells.iter().enumerate() {
        by_cell.entry(cell).or_default().push(row);
    }
    let before: BTreeMap<usize, &Record> =
        loaded.input.records.iter().map(|r| (r.row, r)).collect();
    let mut valid_cells = 0;
    for rows in by_cell.values() {
        let rs: Vec<&Record> = rows
            .iter()
            .filter_map(|row| before.get(row).copied())
            .collect();
        let valid = rows.len() == 2
            && rs.len() == 2
            && rs.iter().filter(|r| r.scope.2 == "IGH").count() == 1
            && rs
                .iter()
                .filter(|r| matches!(r.scope.2.as_str(), "IGK" | "IGL"))
                .count()
                == 1;
        if valid {
            valid_cells += 1;
        } else {
            for &row in rows {
                loaded.input.rows[row].reason = Some("unsupported_or_invalid_cell");
            }
        }
    }
    loaded
        .input
        .records
        .retain(|r| loaded.input.rows[r.row].reason.is_none());
    let mut cells = BTreeMap::new();
    for (i, r) in loaded.input.records.iter().enumerate() {
        cells
            .entry(&loaded.cells[r.row])
            .or_insert_with(Vec::new)
            .push(i);
    }
    let mut signatures = BTreeMap::new();
    let mut units = Vec::new();
    for (key, mut ids) in cells {
        ids.sort_by(|&i, &j| {
            loaded.input.records[i]
                .scope
                .2
                .cmp(&loaded.input.records[j].scope.2)
        });
        let chain: Vec<_> = ids
            .iter()
            .map(|&i| {
                let r = &loaded.input.records[i];
                (
                    r.scope.2.clone(),
                    r.v[0].clone(),
                    r.j[0].clone(),
                    r.junction.clone(),
                )
            })
            .collect();
        signatures
            .entry((&key.0, &key.1, chain))
            .or_insert_with(Vec::new)
            .extend(ids.iter().copied());
        units.push(ids);
    }
    (
        signatures.into_values().collect(),
        units,
        json!({"input_cells":by_cell.len(),"retained_cells":valid_cells,"excluded_cells":by_cell.len()-valid_cells,"unit":"sequence rows; exactly two rows per retained cell"}),
    )
}
fn reconstruction(a: &Config, loaded: &Loaded, mut manifest: Value) -> Result<()> {
    use std::{
        fs::File,
        io::{BufWriter, Write},
    };
    clonotrace_report::transaction(&a.out, |folder| {
        let err = |e: std::io::Error| e.to_string();
        let started = Instant::now();
        let mut source = open(&a.input)?;
        let mut buf = Vec::new();
        if !line(&mut *source.reader, &mut buf)? || buf != loaded.input.header.as_bytes() {
            return Err("input changed between passes".into());
        }
        let mut w = BufWriter::new(File::create(folder.join("reconstruction.tsv")).map_err(err)?);
        if a.method == Method::GermlineV2 {
            writeln!(w,"{}\tclonotrace_germline\tclonotrace_query_alignment\tclonotrace_germline_alignment\tclonotrace_status\tclonotrace_reason",loaded.input.header).map_err(err)?;
        } else {
            writeln!(
                w,
                "{}\tclonotrace_germline\tclonotrace_status\tclonotrace_reason",
                loaded.input.header
            )
            .map_err(err)?;
        }
        let mut i = 0;
        while line(&mut *source.reader, &mut buf)? {
            if i >= loaded.input.rows.len() {
                return Err("input changed between passes".into());
            }
            w.write_all(&buf).map_err(err)?;
            if a.method == Method::GermlineV2 {
                if let Some(r) = &loaded.reconstructions[i] {
                    for value in [
                        &r.query_projection,
                        &r.query_alignment,
                        &r.germline_alignment,
                    ] {
                        w.write_all(b"\t").map_err(err)?;
                        w.write_all(value).map_err(err)?;
                    }
                    w.write_all(b"\tretained\t\n").map_err(err)?;
                } else {
                    writeln!(
                        w,
                        "\t\t\t\texcluded\t{}",
                        loaded.input.rows[i]
                            .reason
                            .ok_or("unaccounted reconstruction row")?
                    )
                    .map_err(err)?;
                }
            } else if let Some(p) = &loaded.projections[i] {
                w.write_all(b"\t").map_err(err)?;
                w.write_all(p).map_err(err)?;
                w.write_all(b"\tretained\t\n").map_err(err)?;
            } else {
                writeln!(
                    w,
                    "\t\texcluded\t{}",
                    loaded.input.rows[i].reason.ok_or("unaccounted row")?
                )
                .map_err(err)?;
            }
            i += 1;
        }
        if i != loaded.input.rows.len() || source.digest() != loaded.input.sha256 {
            return Err("input changed between passes".into());
        }
        w.flush().map_err(err)?;
        w.get_ref().sync_all().map_err(err)?;
        if a.method == Method::GermlineV2 {
            let mut segments =
                BufWriter::new(File::create(folder.join("segments.tsv")).map_err(err)?);
            writeln!(segments,"sequence_id\tsegment\tcall\tcigar\tquery_start_0\tquery_end_0_exclusive\treference_start_0\treference_end_0_exclusive\taligned_bases\tinserted_bases\tdeleted_bases").map_err(err)?;
            for (i, r) in loaded.reconstructions.iter().enumerate() {
                if let Some(r) = r {
                    for d in &r.segments {
                        writeln!(
                            segments,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            loaded.input.rows[i].id,
                            d.segment,
                            d.call,
                            d.cigar,
                            d.query_start,
                            d.query_end,
                            d.reference_start,
                            d.reference_end,
                            d.aligned_bases,
                            d.inserted_bases,
                            d.deleted_bases
                        )
                        .map_err(err)?;
                    }
                }
            }
            segments.flush().map_err(err)?;
            segments.get_ref().sync_all().map_err(err)?;
        }
        let mut diagnostics = BTreeMap::new();
        for row in &loaded.input.rows {
            *diagnostics
                .entry(row.reason.unwrap_or("retained"))
                .or_insert(0usize) += 1;
        }
        let mut diag = String::from("status_or_reason\trows\n");
        for (reason, count) in diagnostics {
            diag.push_str(&format!("{reason}\t{count}\n"));
        }
        manifest["profile_seconds"]["emit"] = json!(started.elapsed().as_secs_f64());
        for (name, contents) in [
            ("diagnostics.tsv", diag),
            (
                "manifest.json",
                serde_json::to_string_pretty(&manifest).map_err(|e| e.to_string())? + "\n",
            ),
        ] {
            let mut file = File::create(folder.join(name)).map_err(err)?;
            file.write_all(contents.as_bytes()).map_err(err)?;
            file.sync_all().map_err(err)?;
        }
        Ok(())
    })
}
pub fn run(a: Config) -> Result<Value> {
    let threshold = configuration(&a)?;
    if std::fs::symlink_metadata(&a.out).is_ok() {
        return Err("output already exists; replacement is prohibited".into());
    }
    let reference_bytes = if a.method == Method::GermlineV2 {
        use std::io::Read;
        let limit = a.max_reference_bytes.unwrap_or(100_000_000);
        let file =
            std::fs::File::open(a.reference_fasta.as_ref().unwrap()).map_err(|e| e.to_string())?;
        let mut bytes = Vec::new();
        file.take(
            limit
                .checked_add(1)
                .ok_or("reference byte bound overflow")?,
        )
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
        if bytes.len() as u64 > limit {
            return Err("max-reference-bytes exceeded".into());
        }
        Some(bytes)
    } else {
        a.reference_fasta
            .as_ref()
            .map(std::fs::read)
            .transpose()
            .map_err(|e| e.to_string())?
    };
    let refs = reference_bytes
        .as_ref()
        .map(|bytes| {
            let text =
                std::str::from_utf8(bytes).map_err(|_| "reference is not UTF-8".to_string())?;
            clonotrace_germline::fasta(text)
        })
        .transpose()?;
    let started = Instant::now();
    let mut loaded = load(&a, refs.as_ref())?;
    let parse = started.elapsed().as_secs_f64();
    let mut manifest = json!({"program":"clonotrace","version":env!("CARGO_PKG_VERSION"),"contract":"CT-extended-1","method":a.method.name(),"disposition":"experimental","reference":"independent mathematical contract; no upstream-native equivalence claim","input":a.input,"input_sha256":loaded.input.sha256,"threads":a.threads,"limits":{"max_rows":a.max_rows,"max_junction_length":a.max_junction_length,"max_pairs":a.max_pairs,"dense_group_rows":if matches!(a.method,Method::HierarchicalV1|Method::SpectralV1|Method::AutoGapV1) {Some(128usize)} else {None}},"configuration":{"threshold":a.threshold,"threshold_binary64":threshold,"linkage":a.linkage,"sigma":a.sigma,"aa_input":a.aa_input,"frame":a.frame},"platform":{"os":std::env::consts::OS,"architecture":std::env::consts::ARCH},"profile_seconds":{"parse":parse},"warnings":["Experimental independent method; computational agreement is not biological validation."]});
    if a.method.is_v2() {
        workflow_v2::metadata(&a, &mut manifest);
    }
    if a.method == Method::PairedExactV1 {
        let mut counts = BTreeMap::new();
        for (donor, sample, _) in &loaded.cells {
            *counts.entry((donor, sample)).or_insert(0u64) += 1;
        }
        let mut total = 0u64;
        for &n in counts.values() {
            total = total
                .checked_add(n.checked_mul(n.saturating_sub(1)).ok_or("pair overflow")? / 2)
                .ok_or("pair overflow")?;
        }
        if total > a.max_pairs {
            return Err("max-pairs exceeded".into());
        }
    }
    manifest["execution"] = json!(if matches!(
        a.method,
        Method::PairedExactV1 | Method::GermlineMaskV1
    ) || a.threads == 1
    {
        "scalar"
    } else {
        "parallel independent work units"
    });
    let started = Instant::now();
    if a.method.is_germline() {
        if a.method == Method::GermlineV2 {
            let pool = workflow_v2::reconstruct_rows(&a, &mut loaded, refs.as_ref().unwrap())?;
            manifest["execution"] = json!(if pool > 1 {
                "parallel independent row reconstructions"
            } else {
                "scalar"
            });
            manifest["effective_pool_threads"] = json!(pool);
        }
        manifest["reference_database"] = json!({"id":a.reference_id,"version":a.reference_version,"source":a.reference_source,"fasta":a.reference_fasta,"sha256":format!("{:x}",Sha256::digest(reference_bytes.as_ref().unwrap())),"sequences":refs.as_ref().unwrap().len(),"limitations":"Query-coordinate V/J projection only; D/N/uncovered masked; exact reference IDs; no indel reconstruction or inferred ancestry."});
        if a.method == Method::GermlineV2 {
            manifest["reference_database"]["d_policy"] = json!(a.germline_mode);
            manifest["reference_database"]["limitations"]=json!("Conditional annotated V(D)J reference reconstruction with indels; unknown junction masked; exact supplied reference IDs; not de novo alignment, allele inference or ancestral consensus.");
        }
        manifest["input_rows"] = json!(loaded.input.rows.len());
        manifest["retained"] = json!(loaded.input.records.len());
        manifest["excluded"] = json!(loaded.input.rows.len() - loaded.input.records.len());
        manifest["failed"] = json!(0);
        manifest["profile_seconds"]["group_cluster"] = json!(started.elapsed().as_secs_f64());
        reconstruction(&a, &loaded, manifest)?;
        return Ok(
            json!({"status":"PASS","method":a.method.name(),"input_rows":loaded.input.rows.len(),"retained":loaded.input.records.len(),"excluded":loaded.input.rows.len()-loaded.input.records.len(),"failed":0}),
        );
    }
    let (mut clones, units, diagnostics) = if a.method == Method::PairedExactV1 {
        paired(&mut loaded)
    } else if a.method == Method::PairedThresholdV2 {
        workflow_v2::infer_paired(&a, &mut loaded, threshold.unwrap())?
    } else {
        let units = work(&a, &loaded.input)?;
        let compute = |ids: &Vec<usize>| cluster(&a, &loaded.input, ids, threshold);
        // Collect in deterministic work-unit order, including all errors; pick first error in order.
        let results: Vec<_> = if a.threads == 1 {
            units.iter().map(compute).collect()
        } else {
            rayon::ThreadPoolBuilder::new()
                .num_threads(a.threads)
                .build()
                .map_err(|e| e.to_string())?
                .install(|| units.par_iter().map(compute).collect())
        };
        let mut clones = Vec::new();
        let mut diagnostics = Vec::new();
        for (i, result) in results.into_iter().enumerate() {
            let (groups, mut diag) = result?;
            diag["work_unit"] = json!(i + 1);
            diag["minimum_sequence_id"] = json!(loaded.input.records[units[i][0]].id);
            diag["rows"] = json!(units[i].len());
            clones.extend(groups);
            diagnostics.push(diag);
        }
        (clones, units, json!(diagnostics))
    };
    canonical(&loaded.input.records, &mut clones);
    manifest["profile_seconds"]["group_cluster"] = json!(started.elapsed().as_secs_f64());
    if a.method == Method::PairedThresholdV2 {
        manifest["execution"] = diagnostics["execution"].clone();
    }
    manifest["diagnostics"] = diagnostics;
    manifest["input_rows"] = json!(loaded.input.rows.len());
    manifest["retained"] = json!(loaded.input.records.len());
    manifest["excluded"] = json!(loaded.input.rows.len() - loaded.input.records.len());
    manifest["failed"] = json!(0);
    manifest["clones"] = json!(clones.len());
    manifest["work_units"] = json!(units.len());
    if a.method == Method::AutoGapV1 {
        manifest["warnings"]=json!(["Uncalibrated nearest-distinct-neighbor gap heuristic; inspect sample counts, ordered distances, gaps and threshold stability. Explicit override may be needed."]);
    }
    let mut extras = BTreeMap::new();
    if a.analyzes_aa() {
        let mut text = String::from("sequence_id\tclonotrace_junction_aa\n");
        for r in &loaded.input.records {
            text.push_str(&r.id);
            text.push('\t');
            text.push_str(std::str::from_utf8(&r.junction).unwrap());
            text.push('\n');
        }
        extras.insert("derived.tsv".into(), text);
    }
    if let Some(cells) = loaded.cell_report {
        extras.insert("cells.tsv".into(), cells);
    }
    clonotrace_report::write_with_files(
        &a.input,
        &a.out,
        &loaded.input,
        &clones,
        &units,
        manifest,
        &extras,
    )?;
    Ok(
        json!({"status":"PASS","method":a.method.name(),"input_rows":loaded.input.rows.len(),"retained":loaded.input.records.len(),"excluded":loaded.input.rows.len()-loaded.input.records.len(),"failed":0,"clones":clones.len()}),
    )
}
