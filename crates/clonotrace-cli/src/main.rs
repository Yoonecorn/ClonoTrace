//! Bounded native-compatible bulk IGH assignment.
mod extended;
mod matrix;
use clap::{Parser, Subcommand, ValueEnum};
use clonotrace_airr::{Policy, Result};
use clonotrace_group::Method;
use serde_json::json;
use std::{collections::BTreeMap, path::PathBuf, time::Instant};
#[derive(Parser)]
#[command(
    version,
    about = "ClonoTrace: qualified bulk IGH nucleotide clone assignment"
)]
struct Cli {
    #[arg(long)]
    reference_version: bool,
    #[command(subcommand)]
    command: Option<Command>,
}
#[derive(Subcommand)]
enum Command {
    /// Validate the bounded input profile without clone inference.
    Validate(InputArgs),
    /// Assign clones and atomically publish a new output directory.
    Assign(AssignArgs),
    /// Explicit versioned inference workflows (native compatibility is separate).
    Extend(Box<extended::Args>),
    /// Cluster an explicit provenance-bearing distance or affinity matrix.
    ClusterMatrix(matrix::MatrixArgs),
}
#[derive(clap::Args)]
struct InputArgs {
    #[arg(long)]
    input: PathBuf,
    #[arg(long, value_enum, default_value = "airr")]
    schema: Schema,
    /// Explicit canonical=source column mapping; repeat for multiple columns.
    #[arg(long = "map")]
    mappings: Vec<String>,
    #[arg(long, default_value_t = 10_000_000)]
    max_rows: usize,
    #[arg(long, default_value_t = 10_000)]
    max_junction_length: usize,
}
#[derive(Clone, Copy, ValueEnum)]
enum Schema {
    Airr,
    Changeo,
}
#[derive(Clone, Copy, ValueEnum)]
enum Algorithm {
    Changeo,
    Scoper,
}
#[derive(Clone, Copy, ValueEnum)]
enum Action {
    Set,
    First,
}
#[derive(clap::Args)]
struct AssignArgs {
    #[command(flatten)]
    input: InputArgs,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, value_enum)]
    method: Algorithm,
    #[arg(long, value_enum, default_value = "set")]
    action: Action,
    #[arg(long, allow_hyphen_values = true)]
    threshold: Option<String>,
    #[arg(long, default_value_t = 1)]
    threads: usize,
    /// Upper bound on all row pairs in each method's scope/length buckets.
    #[arg(long, default_value_t = 100_000_000)]
    max_pairs: u64,
}
fn policy(a: &InputArgs) -> Result<Policy> {
    let mut mappings = BTreeMap::new();
    for entry in &a.mappings {
        let (k, v) = entry
            .split_once('=')
            .ok_or("mapping requires canonical=column")?;
        if v.is_empty() || mappings.insert(k.to_string(), v.to_string()).is_some() {
            return Err("empty or duplicate mapping".into());
        }
    }
    Ok(Policy {
        changeo: matches!(a.schema, Schema::Changeo),
        mappings,
        max_rows: a.max_rows,
        max_junction_length: a.max_junction_length,
    })
}
fn run(cli: Cli) -> Result<()> {
    if cli.reference_version {
        println!("Change-O 1.3.5; SCOPer 1.6.0; CT-0.1-native-1");
        return Ok(());
    }
    match cli.command.ok_or("a subcommand is required; use --help")? {
        Command::Extend(a) => extended::run(*a)?,
        Command::ClusterMatrix(a) => matrix::run(a)?,
        Command::Validate(a) => {
            let d = clonotrace_airr::load(&a.input, &policy(&a)?)?;
            println!(
                "{}",
                json!({"status":"PASS","input_rows":d.rows.len(),"retained":d.records.len(),"excluded":d.rows.len()-d.records.len(),"failed":0,"input_sha256":d.sha256})
            );
        }
        Command::Assign(a) => {
            if a.threads == 0 || a.threads > 1024 {
                return Err("threads must be in 1..=1024".into());
            }
            if std::fs::symlink_metadata(&a.out).is_ok() {
                return Err("output already exists; replacement is prohibited".into());
            }
            let method = match a.method {
                Algorithm::Changeo => Method::Changeo,
                Algorithm::Scoper => Method::Scoper,
            };
            let threshold = match (&a.threshold, method) {
                (Some(t), Method::Changeo) => {
                    if !t
                        .bytes()
                        .all(|b| b.is_ascii_digit() || b".-+eE".contains(&b))
                    {
                        return Err("threshold must be finite decimal in [0,1]".into());
                    }
                    let n: f64 = t.parse().map_err(|_| "invalid threshold")?;
                    if !n.is_finite() || !(0.0..=1.0).contains(&n) {
                        return Err("threshold must be finite decimal in [0,1]".into());
                    }
                    n
                }
                (None, Method::Scoper) => 0.,
                (None, _) => return Err("changeo requires an explicit --threshold".into()),
                (Some(_), _) => return Err("scoper identity does not accept --threshold".into()),
            };
            let p = policy(&a.input)?;
            let start = Instant::now();
            let input = clonotrace_airr::load(&a.input.input, &p)?;
            let parse = start.elapsed().as_secs_f64();
            let start = Instant::now();
            let work = clonotrace_group::build(
                &input.records,
                method,
                matches!(a.action, Action::First),
                a.max_pairs,
            )?;
            let group = start.elapsed().as_secs_f64();
            let start = Instant::now();
            let clones =
                clonotrace_exec::assign(&input.records, &work, method, threshold, a.threads)?;
            let cluster = start.elapsed().as_secs_f64();
            let manifest = json!({"program":"clonotrace","version":env!("CARGO_PKG_VERSION"),"contract":"CT-0.1-native-1","method":method.name(),"reference":if method==Method::Changeo{"Change-O 1.3.5"}else{"SCOPer 1.6.0"},"action":if matches!(a.action,Action::Set){"set"}else{"first"},"threshold":a.threshold,"threshold_binary64":threshold,"threads":a.threads,"schema":if p.changeo{"changeo"}else{"airr"},"mappings":p.mappings,"input":a.input.input,"input_sha256":input.sha256,"input_rows":input.rows.len(),"retained":input.records.len(),"excluded":input.rows.len()-input.records.len(),"failed":0,"clones":clones.len(),"work_units":work.len(),"limits":{"max_rows":a.input.max_rows,"max_junction_length":a.input.max_junction_length,"max_pairs":a.max_pairs},"platform":{"os":std::env::consts::OS,"architecture":std::env::consts::ARCH},"profile_seconds":{"parse":parse,"group":group,"cluster":cluster}});
            clonotrace_report::write(&a.input.input, &a.out, &input, &clones, &work, manifest)?;
            println!(
                "{}",
                json!({"status":"PASS","output":a.out,"input_rows":input.rows.len(),"retained":input.records.len(),"excluded":input.rows.len()-input.records.len(),"failed":0,"clones":clones.len()})
            );
        }
    }
    Ok(())
}
fn main() {
    if let Err(e) = run(Cli::parse()) {
        eprintln!("clonotrace: {e}");
        std::process::exit(2)
    }
}
