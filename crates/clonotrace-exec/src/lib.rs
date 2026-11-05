//! Explicit thread pools, independent whole work units, canonical reduction.
use clonotrace_airr::{Record, Result};
use clonotrace_group::Method;
use rayon::prelude::*;
pub fn assign(
    records: &[Record],
    work: &[Vec<usize>],
    method: Method,
    threshold: f64,
    threads: usize,
) -> Result<Vec<Vec<usize>>> {
    if threads == 0 || threads > 1024 {
        return Err("threads must be in 1..=1024".into());
    }
    let run = || {
        work.par_iter()
            .flat_map_iter(|g| clonotrace_cluster::assign(records, g, method, threshold))
            .collect::<Vec<_>>()
    };
    let mut groups = if threads == 1 {
        work.iter()
            .flat_map(|g| clonotrace_cluster::assign(records, g, method, threshold))
            .collect()
    } else {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| e.to_string())?
            .install(run)
    };
    for g in &mut groups {
        g.sort_unstable_by(|&a, &b| records[a].id.cmp(&records[b].id));
    }
    groups.sort_unstable_by(|a, b| records[a[0]].id.cmp(&records[b[0]].id));
    Ok(groups)
}
