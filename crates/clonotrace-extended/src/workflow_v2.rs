//! Checked v2 workflow adapters. Scientific algorithms remain typed reusable modules.
use super::*;

pub(super) fn validate(a: &Config) -> Result<Option<f64>> {
    if !(1..=1024).contains(&a.threads) {
        return Err("threads must be in 1..=1024".into());
    }
    let threshold = a.threshold.as_deref().map(decimal).transpose()?;
    let hier = a.method == Method::HierarchicalV2;
    let spectral = a.method == Method::SpectralV2;
    let paired = a.method == Method::PairedThresholdV2;
    let germline = a.method == Method::GermlineV2;
    let automatic = a.method == Method::AutoGapStabilityV2;
    if hier {
        if threshold.is_some() == a.clusters.is_some() || a.clusters == Some(0) {
            return Err("hierarchical-v2 requires threshold XOR positive clusters".into());
        }
        clonotrace_cluster::v2::Linkage::parse(
            a.linkage
                .as_deref()
                .ok_or("hierarchical-v2 requires linkage")?,
        )?;
    } else if a.clusters.is_some() || a.linkage.is_some() {
        return Err("clusters/linkage require hierarchical-v2".into());
    }
    if matches!(a.method, Method::AaHammingV2 | Method::PairedThresholdV2) && threshold.is_none() {
        return Err("method requires explicit threshold".into());
    }
    if !matches!(
        a.method,
        Method::HierarchicalV2
            | Method::AaHammingV2
            | Method::PairedThresholdV2
            | Method::AutoGapStabilityV2
    ) && threshold.is_some()
    {
        return Err("threshold unsupported for method".into());
    }
    if spectral {
        if !a.sigma.is_some_and(|s| s.is_finite() && s > 0.0)
            || a.clusters_per_component.is_none_or(|k| k == 0)
        {
            return Err("spectral-v2 requires positive sigma and clusters-per-component".into());
        }
        if a.affinity_cutoff
            .is_some_and(|t| !t.is_finite() || !(0.0..=1.0).contains(&t))
        {
            return Err("affinity-cutoff must be finite in [0,1]".into());
        }
        if a.max_iterations.is_some_and(|n| n == 0 || n > 10000) {
            return Err("max-iterations must be in 1..=10000".into());
        }
    } else if a.sigma.is_some()
        || a.clusters_per_component.is_some()
        || a.affinity_cutoff.is_some()
        || a.max_iterations.is_some()
    {
        return Err("spectral options require spectral-v2".into());
    }
    if hier || spectral {
        if a.alphabet
            .as_deref()
            .is_some_and(|s| !matches!(s, "nt" | "aa"))
        {
            return Err("alphabet requires nt or aa".into());
        }
        if a.max_group_rows.is_some_and(|n| n == 0 || n > 512) {
            return Err("max-group-rows must be in 1..=512".into());
        }
    } else if a.alphabet.is_some() || a.max_group_rows.is_some() {
        return Err("alphabet/group rows require hierarchical-v2 or spectral-v2".into());
    }
    if a.analyzes_aa() {
        match a.aa_input.as_deref() {
            Some("direct") if a.frame.is_none() && a.ambiguity.is_none() => {}
            Some("translate") if a.frame.is_some_and(|f| f <= 2) => {
                if a.ambiguity
                    .as_deref()
                    .is_some_and(|p| !matches!(p, "reject" | "resolve-unanimous"))
                {
                    return Err("ambiguity requires reject or resolve-unanimous".into());
                }
            }
            _ => {
                return Err(
                    "AA requires direct without frame/ambiguity or translate with frame 0/1/2"
                        .into(),
                )
            }
        }
    } else if a.aa_input.is_some() || a.frame.is_some() || a.ambiguity.is_some() {
        return Err("AA options require AA analysis".into());
    }
    if paired {
        if !a
            .missing_light
            .as_deref()
            .is_some_and(|s| matches!(s, "exclude" | "heavy-only-separate"))
            || !a
                .multi_chain
                .as_deref()
                .is_some_and(|s| matches!(s, "exclude" | "collapse-identical"))
        {
            return Err(
                "paired-threshold-v2 requires explicit missing-light and multi-chain policies"
                    .into(),
            );
        }
    } else if a.missing_light.is_some() || a.multi_chain.is_some() {
        return Err("chain policies require paired-threshold-v2".into());
    }
    let references = [
        a.reference_id.as_deref(),
        a.reference_version.as_deref(),
        a.reference_source.as_deref(),
    ];
    if germline {
        if !a
            .germline_mode
            .as_deref()
            .is_some_and(|s| matches!(s, "mask-d" | "include-d"))
            || a.reference_fasta.is_none()
            || references
                .iter()
                .any(|s| s.is_none_or(|v| v.trim().is_empty()))
        {
            return Err(
                "germline-v2 requires germline-mode and reference-fasta/id/version/source".into(),
            );
        }
        if a.max_alignment_length == Some(0) {
            return Err("max-alignment-length must be positive".into());
        }
        if a.max_reference_bytes
            .is_some_and(|n| n == 0 || n > 1_000_000_000)
        {
            return Err("max-reference-bytes must be in 1..=1000000000".into());
        }
    } else if a.reference_fasta.is_some()
        || references.iter().any(Option::is_some)
        || a.germline_mode.is_some()
        || a.max_alignment_length.is_some()
        || a.max_reference_bytes.is_some()
    {
        return Err("reference options require germline-v2".into());
    }
    if automatic {
        if a.max_distinct.is_some_and(|n| n == 0 || n > 128) {
            return Err("max-distinct must be in 1..=128".into());
        }
    } else if a.max_distinct.is_some() || a.max_distance_evaluations.is_some() {
        return Err("stability work bounds require auto-gap-stability-v2".into());
    }
    Ok(threshold)
}
fn estimate_json(x: &clonotrace_cluster::threshold::Estimate) -> Value {
    json!({"threshold":x.threshold,"reason":x.reason,"nearest":x.nearest,"distinct_distances":x.distinct_distances,"gaps":x.gaps})
}
pub(super) fn cluster(
    a: &Config,
    input: &Input,
    ids: &[usize],
    threshold: Option<f64>,
) -> Result<(Vec<Vec<usize>>, Value)> {
    let seqs: Vec<&[u8]> = ids
        .iter()
        .map(|&i| input.records[i].junction.as_slice())
        .collect();
    let (groups, diag) = match a.method {
        Method::HierarchicalV2 => {
            let limit = a.max_group_rows.unwrap_or(512);
            let matrix = clonotrace_cluster::v2::hamming_matrix(&seqs, limit)?;
            let linkage = clonotrace_cluster::v2::Linkage::parse(a.linkage.as_deref().unwrap())?;
            let tree = clonotrace_cluster::v2::build(&matrix, linkage, limit)?;
            let groups = if let Some(t) = threshold {
                tree.cut_threshold(t)?
            } else {
                tree.cut_count(a.clusters.unwrap())?
            };
            let merges: Vec<_> = tree
                .merges
                .iter()
                .map(|m| json!({"left":m.left,"right":m.right,"height":m.height,"size":m.size}))
                .collect();
            let effective = groups.len();
            (
                groups,
                json!({"linkage":a.linkage,"threshold":threshold,"requested_clusters":a.clusters,"effective_clusters":effective,"merges":merges,"tie_policy":"minimum canonical member IDs; original cross-row sum order"}),
            )
        }
        Method::SpectralV2 => {
            let result = clonotrace_spectral::v2::cluster(
                &seqs,
                &clonotrace_spectral::v2::Options {
                    sigma: a.sigma.unwrap(),
                    edge_cutoff: a.affinity_cutoff,
                    clusters_per_component: a.clusters_per_component.unwrap(),
                    max_group_rows: a.max_group_rows.unwrap_or(512),
                    max_iterations: a.max_iterations.unwrap_or(1000),
                },
            )?;
            let components:Vec<_>=result.components.iter().map(|c|json!({"members":c.members.iter().map(|&i|&input.records[ids[i]].id).collect::<Vec<_>>(),"effective_clusters":c.effective_clusters,"eigenvalues":c.eigenvalues,"eigen_blocks":c.eigen_blocks,"iterations":c.iterations,"inertia":c.inertia,"residual":c.residual,"projector":c.projector})).collect();
            (
                result.groups,
                json!({"components":components,"zero_edges":result.zero_edges,"underflow_edges":result.underflow_edges,"eigen_block_tolerance":1e-10,"distance_tie_tolerance":1e-12,"minimum_embedding_norm":1e-8,"requested_clusters_per_component":a.clusters_per_component}),
            )
        }
        Method::AaIdentityV2 | Method::AaHammingV2 => {
            let aa: Vec<amino::AaSequence> = seqs
                .iter()
                .map(|s| amino::AaSequence::new(s))
                .collect::<std::result::Result<_, _>>()
                .map_err(str::to_owned)?;
            let groups = if a.method == Method::AaIdentityV2 {
                amino::identity(&aa)
            } else {
                amino::single_linkage(&aa, threshold.unwrap()).map_err(str::to_owned)?
            };
            (
                groups,
                json!({"threshold":threshold.unwrap_or(0.0),"alphabet":"canonical20","translation_ambiguity":a.ambiguity.as_deref().unwrap_or("reject")}),
            )
        }
        Method::AutoGapStabilityV2 => {
            let analysis = clonotrace_cluster::threshold::analyze(
                &seqs,
                threshold,
                clonotrace_cluster::threshold::Config {
                    max_distinct: a.max_distinct.unwrap_or(128),
                    max_distance_evaluations: a.max_distance_evaluations.unwrap_or(4_000_000),
                },
            )
            .map_err(str::to_owned)?;
            let loo:Vec<_>=analysis.leave_one_out.iter().map(|r|json!({"omitted_junction":String::from_utf8_lossy(&r.omitted_junction),"threshold":r.threshold,"reason":r.reason})).collect();
            let sensitivity:Vec<_>=analysis.sensitivity.iter().map(|s|json!({"threshold":s.threshold,"clone_count":s.clone_count,"pair_disagreements":s.pair_disagreements,"distinct_pairs":s.pairs})).collect();
            let status = if analysis.leave_one_out_accepted == 0 {
                "unavailable"
            } else if analysis.leave_one_out_accepted == analysis.leave_one_out.len() {
                "complete"
            } else {
                "partial"
            };
            let diag = json!({"input_rows":analysis.input_rows,"distinct_junctions":analysis.unique_junctions.len(),"base":estimate_json(&analysis.base),"leave_one_out":loo,"accepted_replicates":analysis.leave_one_out_accepted,"total_replicates":analysis.leave_one_out.len(),"stability_status":status,"stability_range":analysis.stability_range,"override_threshold":analysis.override_threshold,"selected_threshold":analysis.selected_threshold,"sensitivity":sensitivity,"distance_evaluations":analysis.distance_evaluations,"reserved_distance_evaluations":analysis.reserved_distance_evaluations,"interpretation":"Range includes base and successful distinct-junction deletion estimates; failed replicates retained; not a biological confidence interval."});
            let t = analysis.selected_threshold.ok_or_else(|| {
                format!(
                    "auto-gap-stability-v2 refused group {}: {}",
                    input.records[ids[0]].id, diag
                )
            })?;
            (singles(&seqs, t), diag)
        }
        _ => return Err("method is not a bulk v2 cluster method".into()),
    };
    Ok((
        groups
            .into_iter()
            .map(|g| g.into_iter().map(|i| ids[i]).collect())
            .collect(),
        diag,
    ))
}
type PairedWorkflow = (Vec<Vec<usize>>, Vec<Vec<usize>>, Value);
pub(super) fn infer_paired(
    a: &Config,
    loaded: &mut Loaded,
    threshold: f64,
) -> Result<PairedWorkflow> {
    let result = paired::infer_threads(
        &loaded.chains,
        paired::Config {
            threshold,
            missing_light: if a.missing_light.as_deref() == Some("exclude") {
                paired::MissingLightPolicy::Exclude
            } else {
                paired::MissingLightPolicy::HeavyOnlySeparate
            },
            multiple: if a.multi_chain.as_deref() == Some("exclude") {
                paired::MultiChainPolicy::Exclude
            } else {
                paired::MultiChainPolicy::CollapseIdentical
            },
            max_pairs: a.max_pairs,
        },
        a.threads,
    )?;
    for (row, reason) in loaded.input.rows.iter_mut().zip(&result.exclusions) {
        row.reason = *reason;
    }
    loaded
        .input
        .records
        .retain(|r| loaded.input.rows[r.row].reason.is_none());
    let by_row: BTreeMap<usize, usize> = loaded
        .input
        .records
        .iter()
        .enumerate()
        .map(|(i, r)| (r.row, i))
        .collect();
    let remap = |rows: Vec<Vec<usize>>| -> Result<Vec<Vec<usize>>> {
        rows.into_iter()
            .map(|g| {
                g.into_iter()
                    .map(|row| {
                        by_row
                            .get(&row)
                            .copied()
                            .ok_or_else(|| "paired retained row missing from compact input".into())
                    })
                    .collect()
            })
            .collect()
    };
    let groups = remap(result.groups)?;
    let units = remap(result.work_units)?;
    let mut cells=String::from("donor_id\tsample_id\tcell_id\trows\tsequence_ids_json\tstatus\treason\theavy_only\tcollapsed_identical_rows\n");
    for c in &result.cells {
        let sequence_ids: Vec<_> = c
            .rows
            .iter()
            .map(|&r| &loaded.chains[r].sequence_id)
            .collect();
        let json = serde_json::to_string(&sequence_ids).map_err(|e| e.to_string())?;
        cells.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            c.donor_id,
            c.sample_id,
            c.cell_id,
            c.rows.len(),
            json,
            if c.reason.is_none() {
                "retained"
            } else {
                "excluded"
            },
            c.reason.unwrap_or(""),
            c.heavy_only,
            c.collapsed_identical_rows
        ));
    }
    loaded.cell_report = Some(cells);
    let execution = if a.threads > 1 && units.len() > 1 {
        "parallel independent heavy-chain work units"
    } else {
        "scalar"
    };
    let diag = json!({"input_cells":result.cells.len(),"retained_cells":result.retained_cells,"excluded_cells":result.cells.len()-result.retained_cells,"paired_cells":result.paired_cells,"heavy_only_cells":result.heavy_only_cells,"heavy_components":result.heavy_components,"collapsed_identical_rows":result.collapsed_identical_rows,"execution":execution,"cell_diagnostics":"cells.tsv","policies":{"missing_light":a.missing_light,"multi_chain":a.multi_chain},"model":"heavy graph within paired-state classes, then exact known-light refinement; no missing-light bridge"});
    Ok((groups, units, diag))
}
pub(super) fn reconstruct_rows(
    a: &Config,
    loaded: &mut Loaded,
    refs: &BTreeMap<String, Vec<u8>>,
) -> Result<usize> {
    use clonotrace_germline::v2::{reconstruct_bounded, DPolicy, SegmentInput};
    let policy = if a.germline_mode.as_deref() == Some("include-d") {
        DPolicy::Include
    } else {
        DPolicy::Mask
    };
    let compute = |r: &Record| {
        let args = loaded.germline_inputs[r.row]
            .as_ref()
            .ok_or("missing reconstruction input")?;
        let d = if args.d_call.is_empty() && args.d_cigar.is_empty() {
            None
        } else {
            Some(SegmentInput {
                call: &args.d_call,
                cigar: &args.d_cigar,
            })
        };
        if policy == DPolicy::Include
            && r.scope.2 == "IGH"
            && !args.d_call.is_empty()
            && gene(&args.d_call, "IGH", "D").is_none()
        {
            return Err("invalid_D_call");
        }
        reconstruct_bounded(
            &r.junction,
            &r.scope.2,
            SegmentInput {
                call: &args.v_call,
                cigar: &args.v_cigar,
            },
            SegmentInput {
                call: &args.j_call,
                cigar: &args.j_cigar,
            },
            d,
            policy,
            refs,
            a.max_alignment_length.unwrap_or(20_000),
        )
    };
    let threads = if a.threads > 1 && loaded.input.records.len() > 1 {
        a.threads
    } else {
        1
    };
    let results: Vec<_> = if threads == 1 {
        loaded.input.records.iter().map(compute).collect()
    } else {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| e.to_string())?
            .install(|| loaded.input.records.par_iter().map(compute).collect())
    };
    for (r, result) in loaded.input.records.iter().zip(results) {
        match result {
            Ok(output) => loaded.reconstructions[r.row] = Some(output),
            Err("max_alignment_length_exceeded" | "alignment_size_overflow") => {
                return Err(format!(
                    "reconstruction resource bound exceeded for {}",
                    r.id
                ))
            }
            Err(reason) => loaded.input.rows[r.row].reason = Some(reason),
        }
    }
    loaded
        .input
        .records
        .retain(|r| loaded.input.rows[r.row].reason.is_none());
    Ok(threads)
}
pub(super) fn metadata(a: &Config, m: &mut Value) {
    m["contract"] = json!("CT-complete-2");
    m["disposition"] = json!("qualified-implementation");
    m["implementation_status"] = json!(
        "production implementation of explicit bounded computational contract; qualified platform"
    );
    m["scientific_evidence"] = json!(if a.method == Method::AutoGapStabilityV2 {
        "Exploratory uncalibrated threshold model; deletion sensitivity is not biological uncertainty calibration."
    } else {
        "Computational contract qualification, not biological superiority or general upstream-native equivalence."
    });
    m["manifest_schema"] = json!("CT-run-2");
    m["output_schemas"] = if a.method.is_germline() {
        json!({"reconstruction":"CT-germline-2","segments":"CT-germline-segments-2"})
    } else {
        json!({"assignments":"CT-assignment-1","clones":"CT-clone-summary-1","derived_aa":"CT-aa-derived-1","cells":"CT-cells-2 when paired"})
    };
    m["configuration"]["alphabet"] = json!(if a.analyzes_aa() { "aa" } else { "nt" });
    for (key, value) in [
        (
            "ambiguity",
            json!(a.ambiguity.as_deref().unwrap_or("reject")),
        ),
        ("clusters", json!(a.clusters)),
        ("clusters_per_component", json!(a.clusters_per_component)),
        ("affinity_cutoff", json!(a.affinity_cutoff)),
        ("missing_light", json!(a.missing_light)),
        ("multi_chain", json!(a.multi_chain)),
        ("germline_mode", json!(a.germline_mode)),
    ] {
        m["configuration"][key] = value;
    }
    m["limits"]["dense_group_rows"] = json!(if matches!(
        a.method,
        Method::HierarchicalV2 | Method::SpectralV2
    ) {
        Some(a.max_group_rows.unwrap_or(512))
    } else {
        None
    });
    m["limits"]["max_iterations"] = json!(a.max_iterations.unwrap_or(1000));
    m["limits"]["max_alignment_length"] = json!(a.max_alignment_length.unwrap_or(20_000));
    m["limits"]["max_distinct"] = json!(a.max_distinct.unwrap_or(128));
    m["limits"]["max_distance_evaluations"] =
        json!(a.max_distance_evaluations.unwrap_or(4_000_000));
    m["limits"]["max_reference_bytes"] = json!(a.max_reference_bytes.unwrap_or(100_000_000));
    m["warnings"]=json!(["Method assumptions and numerical tie policies are explicit; inspect diagnostics. No automatic production biological cutoff or clone-number guarantee."]);
}
