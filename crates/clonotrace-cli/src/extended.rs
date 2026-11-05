//! CLI adapter; library owns validation, inference and atomic publication.
use clap::{Args as ClapArgs, ValueEnum};
use clonotrace_airr::Result;
use std::path::PathBuf;
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
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
#[derive(ClapArgs)]
pub struct Args {
    #[arg(long)]
    input: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, value_enum)]
    method: Method,
    #[arg(long, allow_hyphen_values = true)]
    threshold: Option<String>,
    #[arg(long)]
    linkage: Option<String>,
    #[arg(long, allow_hyphen_values = true)]
    sigma: Option<f64>,
    /// AA source: direct junction_aa or translate junction with explicit frame.
    #[arg(long)]
    aa_input: Option<String>,
    #[arg(long)]
    frame: Option<usize>,
    #[arg(long)]
    reference_fasta: Option<PathBuf>,
    #[arg(long)]
    reference_id: Option<String>,
    #[arg(long)]
    reference_version: Option<String>,
    #[arg(long)]
    reference_source: Option<String>,
    #[arg(long, default_value_t = 1)]
    threads: usize,
    #[arg(long, default_value_t = 1_000_000)]
    max_rows: usize,
    #[arg(long, default_value_t = 10_000)]
    max_junction_length: usize,
    #[arg(long, default_value_t = 100_000_000)]
    max_pairs: u64,
    #[arg(long)]
    alphabet: Option<String>,
    #[arg(long)]
    ambiguity: Option<String>,
    #[arg(long)]
    clusters: Option<usize>,
    #[arg(long)]
    clusters_per_component: Option<usize>,
    #[arg(long)]
    affinity_cutoff: Option<f64>,
    #[arg(long)]
    missing_light: Option<String>,
    #[arg(long)]
    multi_chain: Option<String>,
    #[arg(long)]
    germline_mode: Option<String>,
    #[arg(long)]
    max_group_rows: Option<usize>,
    #[arg(long)]
    max_iterations: Option<usize>,
    #[arg(long)]
    max_alignment_length: Option<usize>,
    #[arg(long)]
    max_distinct: Option<usize>,
    #[arg(long)]
    max_distance_evaluations: Option<u64>,
    #[arg(long)]
    max_reference_bytes: Option<u64>,
}
pub fn run(a: Args) -> Result<()> {
    let method = match a.method {
        Method::HierarchicalV1 => clonotrace_extended::Method::HierarchicalV1,
        Method::SpectralV1 => clonotrace_extended::Method::SpectralV1,
        Method::PairedExactV1 => clonotrace_extended::Method::PairedExactV1,
        Method::AaIdentityV1 => clonotrace_extended::Method::AaIdentityV1,
        Method::AaHammingV1 => clonotrace_extended::Method::AaHammingV1,
        Method::AutoGapV1 => clonotrace_extended::Method::AutoGapV1,
        Method::GermlineMaskV1 => clonotrace_extended::Method::GermlineMaskV1,
        Method::HierarchicalV2 => clonotrace_extended::Method::HierarchicalV2,
        Method::SpectralV2 => clonotrace_extended::Method::SpectralV2,
        Method::PairedThresholdV2 => clonotrace_extended::Method::PairedThresholdV2,
        Method::AaIdentityV2 => clonotrace_extended::Method::AaIdentityV2,
        Method::AaHammingV2 => clonotrace_extended::Method::AaHammingV2,
        Method::GermlineV2 => clonotrace_extended::Method::GermlineV2,
        Method::AutoGapStabilityV2 => clonotrace_extended::Method::AutoGapStabilityV2,
    };
    let result = clonotrace_extended::run(clonotrace_extended::Config {
        input: a.input,
        out: a.out,
        method,
        threshold: a.threshold,
        linkage: a.linkage,
        sigma: a.sigma,
        aa_input: a.aa_input,
        frame: a.frame,
        reference_fasta: a.reference_fasta,
        reference_id: a.reference_id,
        reference_version: a.reference_version,
        reference_source: a.reference_source,
        threads: a.threads,
        max_rows: a.max_rows,
        max_junction_length: a.max_junction_length,
        max_pairs: a.max_pairs,
        alphabet: a.alphabet,
        ambiguity: a.ambiguity,
        clusters: a.clusters,
        clusters_per_component: a.clusters_per_component,
        affinity_cutoff: a.affinity_cutoff,
        missing_light: a.missing_light,
        multi_chain: a.multi_chain,
        germline_mode: a.germline_mode,
        max_group_rows: a.max_group_rows,
        max_iterations: a.max_iterations,
        max_alignment_length: a.max_alignment_length,
        max_distinct: a.max_distinct,
        max_distance_evaluations: a.max_distance_evaluations,
        max_reference_bytes: a.max_reference_bytes,
    })?;
    println!("{result}");
    Ok(())
}
