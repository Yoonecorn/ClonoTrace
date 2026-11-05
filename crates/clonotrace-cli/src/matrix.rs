//! Declared matrix argument adapter.
use clap::Args;
use std::{path::PathBuf, str::FromStr};
#[derive(Args)]
pub struct MatrixArgs {
    #[arg(long)]
    input: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    method: String,
    #[arg(long)]
    linkage: Option<String>,
    #[arg(long, allow_hyphen_values = true)]
    threshold: Option<f64>,
    #[arg(long)]
    clusters: Option<usize>,
    #[arg(long)]
    clusters_per_component: Option<usize>,
    #[arg(long)]
    max_iterations: Option<usize>,
    #[arg(long, default_value_t = 512)]
    max_rows: usize,
    #[arg(long, default_value_t = 32_000_000)]
    max_input_bytes: u64,
}
pub fn run(a: MatrixArgs) -> clonotrace_airr::Result<()> {
    let result = clonotrace_extended::matrix::run(clonotrace_extended::matrix::Config {
        input: a.input,
        out: a.out,
        method: clonotrace_extended::Method::from_str(&a.method)?,
        linkage: a.linkage,
        threshold: a.threshold,
        clusters: a.clusters,
        clusters_per_component: a.clusters_per_component,
        max_iterations: a.max_iterations,
        max_rows: a.max_rows,
        max_input_bytes: a.max_input_bytes,
    })?;
    println!("{result}");
    Ok(())
}
