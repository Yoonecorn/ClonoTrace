# ClonoTrace user manual

**Version 1.0.0**

## What ClonoTrace does

B cells descended from a common rearrangement can carry related, but nonidentical,
receptor sequences. Repertoire analysis therefore needs an explicit rule for
deciding which sequence observations belong together. Gene-call ambiguity,
junction mutations, incomplete chain recovery and sampling boundaries can all
change that decision.

ClonoTrace is a deterministic Rust toolkit for B-cell receptor clonal inference.
It reproduces specified Change-O and SCOPer behavior and provides separate
hierarchical, spectral, paired-chain, amino-acid and annotated germline workflows.
The method names identify different models. Their outputs should be interpreted
with their assumptions, input profiles and diagnostics.

This manual describes the current native and v2 workflows. The older `*-v1`
methods remain available for reproducing previous analyses. Package version
`1.0.0` does not replace the scientific contract or output-schema identifiers.
Qualification is local macOS ARM64 computational qualification. Agreement with a
software reference does not establish biological clonality or superior accuracy.

See the model-selection guidance and scientific contracts included later in this manual before selecting a model.

## Build, install and check

The repository pins Rust 1.98.0 with rustfmt and Clippy. Run these commands from
the extracted source-package or repository root:

```sh
export CARGO_TARGET_DIR="$PWD/target"
cargo build --release --locked --offline
target/release/clonotrace --version
target/release/clonotrace --reference-version
target/release/clonotrace --help
```

The source release uses the exact dependencies in Cargo.lock. Populate a clean
checkout cache once with `cargo fetch --locked` when network access is available,
then use the locked offline commands above. Dependencies and the pinned Rust
toolchain are not bundled in the source archive. ClonoTrace itself does not need
Python, R, Change-O or SCOPer for inference.

The version command prints `clonotrace-cli 1.0.0`. Setting
`CARGO_TARGET_DIR` explicitly keeps build products in this workspace even if the
shell inherited a different target directory.

To install into a chosen local directory:

```sh
cargo install --path crates/clonotrace-cli --root "$PWD/local-install" \
  --locked --offline
local-install/bin/clonotrace --version
```

The executable below is `target/release/clonotrace`. Substitute the installed
path if appropriate. 

Run the released checks from the project root:

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked --offline -- -D warnings
cargo test --workspace --locked --offline
```

The examples are synthetic and do not represent biological truth.

For a first native run using the longstanding small input:

```sh
target/release/clonotrace validate --input examples/repertoire.tsv
target/release/clonotrace assign --input examples/repertoire.tsv \
  --out example-changeo --method changeo --threshold .15
target/release/clonotrace assign --input examples/repertoire.tsv \
  --out example-scoper --method scoper
```

Each output directory must be new. The threshold `.15` is an illustration, not a
recommended biological cutoff. Inspect `manifest.json` and the excluded rows
before interpreting clone counts.

## Prepare a receptor table

Inputs are UTF-8 tab-separated tables, either plain or gzip compressed. Tabs and
newlines delimit fields and records; quoted multiline fields, stdin and zstd are
unsupported. Original fields, values and row order are preserved in assignments.
ClonoTrace does not silently repair case, whitespace, ambiguous calls or IDs.

All receptor-table workflows require `sequence_id`, `donor_id`, `sample_id`,
`locus`, `v_call` and `j_call`. IDs must be globally unique and nonblank. Donor and
sample values must be explicit, nonblank strings supported by the data source.
Donor/sample/locus boundaries prevent cross-sample inference. Do not invent a
donor identity from sequence similarity.

Native assignment and nucleotide inference require `junction`. Direct
amino-acid inference requires `junction_aa`. Paired inference requires active
`cell_id`; germline reconstruction also needs `sequence` and V/J CIGAR annotations.
Columns required only by a selected method are described below.

The native bulk profile accepts IGH with uppercase ACGT junctions. Gene calls may
be gene-only or use numeric allele suffixes such as `IGHV1-2*01`. Gene resolution
removes only the numeric allele suffix; it does not infer an allele. Native
methods can use comma-separated V/J calls. Current bulk v2 methods require one
unambiguous V gene and J gene and group by exact donor/sample/locus/gene/length.
The accepted gene-name grammar is deliberately narrower than all AIRR/IMGT names.

If `junction_length` is present, it must agree with the actual junction length.
Duplicate headers, ragged rows, missing required columns, invalid UTF-8,
blank/duplicate IDs and reserved output-column collisions fail the whole run.
Rows outside the selected biological input profile remain visible as exclusions.
Exclusion precedence and accepted fields are fixed by the relevant contract.

For a bulk projection of a single-cell table, deliberately rename `cell_id` to
`source_cell_id` and prior `clone_id` to `source_clone_id`, preserving the values.
Document that chains are being analyzed as rows. Active `cell_id` or `CELL` in a
bulk input is rejected rather than silently ignored. Keep active cell annotations
when selecting the paired-chain model.

Native `validate` and `assign` accept `--schema changeo` for uppercase Change-O
columns. Repeat `--map canonical=source_column` for explicit field mapping:

```sh
target/release/clonotrace validate --input uppercase.tsv --schema changeo
target/release/clonotrace assign --input mapped.tsv --out mapped-result \
  --method changeo --threshold .15 --map donor_id=donor
```

Mapping does not convert arbitrary annotations or broaden the input profile.
The `extend` command uses the lowercase AIRR-profile fields directly.

## Reproduce native assignment behavior

Use `assign --method changeo` for the qualified Change-O 1.3.5 profile:
gene resolution, normalized nucleotide Hamming distance and inclusive
single-linkage thresholding. An explicit finite threshold in `[0,1]` is required.
The default ambiguity action is `set`; `--action first` uses the first listed
V/J gene call.

```sh
target/release/clonotrace assign --input bulk.tsv.gz --out native-changeo \
  --method changeo --action set --threshold .15 --threads 1
```

Some ambiguous-call patterns produce order-dependent Change-O partitions.
ClonoTrace preserves that qualified native behavior. Preserve the input order
and ambiguity action to reproduce a native run. Canonical numeric labels do not
remove this partition-level effect.

Use `assign --method scoper` for the qualified SCOPer 1.6.0 `identicalClones`
nucleotide profile. A threshold is not accepted:

```sh
target/release/clonotrace assign --input bulk.tsv.gz --out native-scoper \
  --method scoper --action set --threads 1
```

Its observed final membership combines compatible ambiguous V/J annotations with
exact junction identity, including qualified cross-length ambiguity bridges.
These rules differ from Change-O grouping. The target is this specified
SCOPer identity profile, not every SCOPer clustering mode. Compare partitions by
member IDs, not by arbitrary reference clone numbers.

## Choose a hierarchical model

`hierarchical-v2` provides row-weighted single, average and complete linkage on
normalized Hamming distances. Select an inclusive threshold or a positive
cluster count per work unit, but not both:

```sh
target/release/clonotrace extend --input bulk.tsv --out hierarchy \
  --method hierarchical-v2 --linkage average --threshold .15
target/release/clonotrace extend --input bulk.tsv --out hierarchy-count \
  --method hierarchical-v2 --linkage complete --clusters 2
```

All rows contribute to average-linkage weighting. Stable ties use canonical
minimum member IDs. The manifest includes the full merge structure and cut
diagnostics. A requested count larger than the work-unit size saturates at that
size and is diagnosed. Hierarchy and spectral sequence work units have an
absolute limit of 512 rows; `--max-group-rows` can lower it.

## Use an explicit spectral model

`spectral-v2` clusters a declared Gaussian affinity graph. Choose a positive
finite bandwidth and positive cluster count per positive-edge component:

```sh
target/release/clonotrace extend --input bulk.tsv --out spectral \
  --method spectral-v2 --sigma .25 --clusters-per-component 2
```

The off-diagonal affinity is `exp(-0.5 * (d / sigma)^2)` for normalized Hamming
distance `d`; the diagonal is zero. Optional `--affinity-cutoff .4` removes edges
with distance greater than `.4`. Underflow can also produce zero edges and is
counted. Positive-edge components never merge; isolates remain singletons. The
requested count saturates within each component.

The method uses a symmetric normalized Laplacian, deterministic cyclic Jacobi
eigendecomposition, canonical repeated-eigenspace bases and deterministic Lloyd
clustering. Numerical tolerances, tie rules and convergence bounds are fixed.
Ill-conditioned embeddings and nonconvergence fail the run. There is no fallback
to a different solver or automatic biological estimate of the number of clones.

The dense implementation improves workspace traversal while
preserving the qualified arithmetic and results. On the qualified dense-512
workflow, the one-worker median fell from 5.03 to 1.15 seconds, a 4.39-fold
improvement. Dense execution did not show meaningful worker scaling; the
improvement must not be attributed to multithreading. These local measurements
do not predict every repertoire or host. Performance details are in the
[manuscript](../publication/ClonoTrace-Manuscript-v1.0-2026-11-04.pdf).

## Analyze amino-acid junctions

Select direct canonical amino acids or explicit forward translation:

```sh
target/release/clonotrace extend --input amino.tsv --out amino-identity \
  --method aa-identity-v2 --aa-input direct
target/release/clonotrace extend --input bulk.tsv --out amino-hamming \
  --method aa-hamming-v2 --aa-input translate --frame 0 --threshold .15
```

Direct input accepts the 20 uppercase canonical residues in `junction_aa`.
It rejects `--frame` and `--ambiguity`. Translation uses standard genetic code 1
in a caller-selected frame `0`, `1` or `2`. The whole nucleotide sequence is
validated, including skipped prefix bases. Stops, invalid symbols and incomplete
final codons exclude the row. Default ambiguity is `reject`.

`--ambiguity resolve-unanimous` accepts an IUPAC codon only when every concrete
expansion encodes the same canonical residue. It never substitutes `X`.
AA identity groups equal derived sequences. AA Hamming uses inclusive normalized
Hamming single linkage within same-length work units.

Hierarchy and spectral inference also accept AA distances:

```sh
target/release/clonotrace extend --input amino.tsv --out amino-hierarchy \
  --method hierarchical-v2 --alphabet aa --aa-input direct \
  --linkage average --clusters 2
```

For AA spectral inference, supply the same explicit AA source options together
with the required sigma and per-component count. Inspect `derived.tsv` to see
the analyzed AA strings; source annotations are preserved. Synonymous nucleotide
changes can disappear in AA space, so an AA partition addresses a different
biological question.

## Retain heavy-light pairing

`paired-threshold-v2` uses active `cell_id` and IGH plus IGK/IGL chain rows.
Both cell policies are mandatory:

```sh
target/release/clonotrace extend --input paired.tsv --out paired \
  --method paired-threshold-v2 --threshold .15 \
  --missing-light exclude --multi-chain collapse-identical
```

Heavy chains form thresholded single-linkage components within donor/sample,
heavy V/J/length and paired-state classes. Exact known light signatures refine
those components. This heavy-first rule does not infer an ancestor or perform
joint heavy-light evolutionary modeling.

`--missing-light exclude` excludes cells lacking light chains.
`--missing-light heavy-only-separate` retains them in a separate class that
cannot bridge paired cells. A missing heavy chain always excludes the cell.
`--multi-chain exclude` rejects multiple-chain cells;
`--multi-chain collapse-identical` permits only identical
locus/gene/junction signatures. Conflicting chains still exclude the whole cell.
Identical observations remain in the output and do not prove molecular identity.
Inspect cell-level diagnostics and keep row counts separate from cell counts.

## Reconstruct supplied germline annotations

`germline-v2` reconstructs reference sequence from supplied segment annotations.
It does not assign clones. In addition to the shared fields, supply `sequence`,
`v_cigar`, `j_cigar` and a canonical reference FASTA with exact call IDs:

```sh
target/release/clonotrace extend --input annotated.tsv --out germline \
  --method germline-v2 --germline-mode mask-d \
  --reference-fasta supplied-vdj.fa --reference-id supplied-reference \
  --reference-version exact-version --reference-source exact-source
```

The reference ID, version, source and raw SHA-256 identify the reference. Headers
must be unique bare IDs and reference sequences must be nonempty uppercase ACGT.
No database or allele is inferred. Biological reference assets retain their own
terms; the source-package license does not grant rights to a chosen database.

`mask-d` ignores D annotations deliberately. `include-d` requires `d_call` and
`d_cigar` for IGH and forbids D annotations on light chains. Unknown rearrangement
positions remain `N`. V/J are required in both modes.

AIRR CIGAR allows terminal `S`/`N` clips around an aligned `M`, `=`, `X`, `I`, `D`
body. `M` cannot mix with `=`/`X`, which must agree with the supplied sequences.
Insertions consume query and yield a gap in aligned germline and `N` in the
query-coordinate projection. Deletions consume reference and retain the deleted
reference bases opposite query gaps. Unsupported operators, overlapping segments,
invalid bounds and missing/ambiguous references are diagnosed and refused or
excluded according to the fixed contract. Invalid reference files and expanded
alignment resource overflow fail the whole run.

Read the query-coordinate projection, paired gapped alignment and segment
diagnostics together. Reconstruction does not align reads de novo, discover
alleles, reconstruct clone ancestors or establish general CreateGermlines parity.

## Inspect an exploratory automatic threshold

`auto-gap-stability-v2` estimates a threshold from nearest-distinct-junction
distances and a unique largest gap. It reruns estimation after deleting each
distinct junction:

```sh
target/release/clonotrace extend --input bulk.tsv --out threshold-analysis \
  --method auto-gap-stability-v2
target/release/clonotrace extend --input bulk.tsv --out threshold-override \
  --method auto-gap-stability-v2 --threshold .15
```

Insufficient, constant or tied-gap evidence causes a diagnostic refusal. Without
an explicit override, a refusal in any work unit fails the whole run. An override
does not suppress the failed estimate or stability diagnostics. Inspect every
deletion result, acceptance fraction, supported range and partition sensitivity.
The deletion range is not a confidence interval or a calibrated biological error
rate. This method remains biologically exploratory.

The absolute limit is 128 distinct junctions per work unit. The conservative
distance-evaluation reservation is checked before computation; exact caching can
reduce actual evaluations without changing that acceptance boundary.

## Supply a declared matrix

`cluster-matrix` accepts one declared distance or affinity matrix in strict JSON
schema `CT-matrix-1`. Provide `kind`, unique nonblank `ids`, a donor/sample/locus
`scope`, nonblank `provenance` and a finite nonnegative exactly symmetric
zero-diagonal square `matrix`. Unknown and duplicate fields fail. Both axes are
reordered with canonical IDs.

```json
{"schema":"CT-matrix-1","kind":"distance","ids":["b","a"],"scope":{"donor_id":"d","sample_id":"s","locus":"IGH"},"provenance":"synthetic distance example, version 1","matrix":[[0,2],[2,0]]}
```

```sh
target/release/clonotrace cluster-matrix --input distance.json \
  --out matrix-hierarchy --method hierarchical-v2 \
  --linkage average --threshold 2
target/release/clonotrace cluster-matrix --input affinity.json \
  --out matrix-spectral --method spectral-v2 --clusters-per-component 2
```

Distance thresholds use the caller's distance units and may exceed one; count
cuts are also supported. Supplied affinities are used directly, without Gaussian
conversion or sigma. Matrix inference is a single scalar work unit and has no
`--threads` option. The 512-row limit remains absolute. A provenance field
identifies a model; it does not validate that model biologically.

## Read the results

Successful CLI workflows print a JSON result and exit 0. Invalid commands and
operational, structural, resource or numerical failures exit 2 and describe the
problem on stderr. Scientific exclusions are reported within a successful run.
For receptor tables, `input_rows = retained + excluded` and `failed = 0`.

| Product | Interpretation |
| --- | --- |
| `assignments.tsv` | Original rows and fields plus run-local clone ID, status and reason; exclusions have blank clone ID. |
| `clones.tsv` | Clone summaries; row counts are sequence observations. |
| `diagnostics.tsv` | Work-unit or method diagnostics; these are not asserted SCOPer internal candidate groups. |
| `manifest.json` | Contract/schema, method, configuration, input hash, limits, accounting, numerical diagnostics, platform and measured profile. |
| `derived.tsv` | Derived AA sequences for AA workflows. |
| `cells.tsv` | Cell-level status, policies and chain accounting under `CT-cells-2`. |
| `reconstruction.tsv` | Germline query projection and paired alignments; no clone assignment. |
| `segments.tsv` | Germline segment-level calls, coordinates and alignment diagnostics. |

The method determines which products are present. Germline emits only
reconstruction, segments, diagnostics and a manifest; it does not emit assignment
or clone-summary tables. Matrix assignment uses
`CT-matrix-assignment-1`. v2 manifests use `CT-complete-2` / `CT-run-2`, while
assignment and clone-summary schemas remain `CT-assignment-1` and
`CT-clone-summary-1`. AA derived output retains `CT-aa-derived-1`; germline
products identify `CT-germline-2` and `CT-germline-segments-2`. These scientific
identifiers remain distinct from the software package version.

Clone numbers are ordered by the smallest member `sequence_id` in bytewise order.
They are run-local labels, not persistent biological identifiers. A new input can
change labels even when some memberships remain identical.

## Output safety and reproducibility

An output path must be new, with an existing parent directory. Existing files,
directories and symlinks are protected. The completed result becomes visible by
exclusive atomic directory rename. Precommit failure publishes no result
directory. A postcommit parent-directory synchronization warning retains the
published result; atomic visibility does not guarantee survival after every crash.

Keep the input unchanged for the entire run. A second raw-input SHA-256 check
detects between-pass mutation. Gzip truncation and checksum failures fail the run.
An abruptly killed process can leave a temporary sibling transaction directory;
inspect it and check that no running process owns it before removing it.

For fixed input bytes and configuration, qualified scientific tables agree
byte-for-byte across repeated runs and tested worker counts. Profile timings,
platform metadata and requested worker counts can change manifest bytes.
v2 methods canonically order IDs before inference. Native Change-O additionally
preserves its qualified input-order dependence for ambiguous calls.

Archive the input hash, executable/version, full command, manifest, assignment
tables, reference hash where relevant and provenance needed to interpret scopes.
Comparing only clone counts cannot establish partition agreement.

## Set resource limits deliberately

| Option | Default | Applies to |
| --- | ---: | --- |
| `--threads` | 1 | Native assign and extend; 1 through 1024 workers. |
| `--max-rows` | 10000000 native; 1000000 extend | Whole receptor table. |
| `--max-junction-length` | 10000 | Analyzed source sequence/query byte length. |
| `--max-pairs` | 100000000 | Conservative pair bound where applicable. |
| `--max-group-rows` | 512 | v2 hierarchy/spectral; can only lower the absolute limit. |
| `--max-iterations` | 1000 | Spectral; maximum configurable value 10000. |
| `--max-alignment-length` | 20000 | Germline expanded output, checked before allocation. |
| `--max-reference-bytes` | 100000000 | Germline reference; maximum 1000000000. |
| `--max-distinct` | 128 | Automatic v2; absolute maximum 128. |
| `--max-distance-evaluations` | 4000000 | Automatic v2 conservative per-unit reservation. |
| Matrix `--max-rows` | 512 | Absolute matrix row limit. |
| Matrix `--max-input-bytes` | 32000000 | Maximum configurable value 256000000. |

Pair reservations can exceed the number of comparisons actually performed by an
exact optimization. Native grouping and v2 grouping use different work-unit
definitions. Paired inference reserves input row-pair work before cell filtering;
germline reconstruction has no pair computation.

Workers distribute complete independent work units where implemented. They do
not make one dense spectral component an internally parallel eigensolver. Small
groups, input/output and startup can dominate. Measure the intended workflow
before choosing more workers. Resource limits reject excessive declared work;
they are not an operating-system RAM cap. Compact records, IDs, sequences and
labels still grow with input size.

## Use the Rust libraries

The workspace exposes checked APIs as well as the CLI. `clonotrace-extended`
provides `Config`, `Method` and `run` for complete table workflows. The scientific
crates expose hierarchy, spectral, threshold and germline modules; AA and paired
APIs are in `clonotrace-extended`. Native APIs remain in their original crates.

```rust
use clonotrace_extended::{run, Config, Method};

fn analyze() -> Result<serde_json::Value, String> {
    let mut config = Config::new("bulk.tsv", "result", Method::HierarchicalV2);
    config.linkage = Some("average".into());
    config.threshold = Some("0.15".into());
    config.validate()?;
    run(config)
}
```

The workflow returns a machine-readable result without printing or exiting.
Algorithm-only callers must establish their own scope units and canonical
bytewise ID order; the file workflow does this for them. Checked APIs return
errors rather than approximate fallbacks. The threshold API can return a
scientific refusal as an `Analysis` with no selected threshold; the no-override
CLI turns that refusal into a failed run. See the build and API reference included below for entry points.

## Reproduce older analyses

The CLI retains `hierarchical-v1`, `spectral-v1`, `paired-exact-v1`,
`aa-identity-v1`, `aa-hamming-v1`, `germline-mask-v1` and `auto-gap-v1`.
They keep their original `CT-extended-1` semantics and bounds. A newer package
does not silently substitute v2 behavior for a named v1 method. Use the historical extension contract included below and the original
manifest when reproducing these methods. Its historical experimental scope is
specific to those versions and should not be read as the current v2 feature list.

## Troubleshooting

| Symptom | Next step |
| --- | --- |
| Offline build cannot find a dependency | Populate the locked cache, then repeat the offline build. |
| Existing output path | Choose a fresh path; replacement is intentionally unavailable. |
| Missing boundary or duplicate ID | Repair the source-supported metadata upstream. |
| Active cell column rejected | Select paired inference or document an explicit bulk projection. |
| Many exclusions | Read status/reason columns and check locus, calls, alphabet and annotations. |
| Pair/group/distinct limit exceeded | Inspect work-unit composition; use supported explicit bounds. The dense maximum remains 512. |
| Numerical refusal | Inspect bandwidth, graph components and diagnostics; preserve the failure rather than changing tolerances. |
| Automatic estimate refuses | Inspect distance evidence; any override must be explicit and scientifically justified. |
| Different Change-O result after sorting | Restore original input order and ambiguity action. |
| Different reference clone numbers | Compare complete member partitions. |
| More workers are slower | Measure startup, I/O and work-unit sizes; one worker remains the default. |
| Changed-input or gzip failure | Restore the intended immutable input and verify its checksum. |
| Output-parent error | Create an accessible existing parent and choose a new child path. |

## Interpretation, licenses and release status

The toolkit does not infer cross-sample clones, genotype alleles, align raw reads,
calibrate biological false-discovery rates or establish a universal clone model.
Unqualified platforms, arbitrary repertoire sizes and biological accuracy require
additional evidence. The model-selection guidance below describes the major scientific
limitations without treating computational acceptance as biological validation.

Original ClonoTrace code is licensed `Apache-2.0 OR MIT`. Dependencies, upstream
tools, reference databases, datasets and retained third-party materials keep their
own terms; see [NOTICE](../NOTICE) and [LICENSE](../LICENSE).



---

# Choosing and interpreting a ClonoTrace method

**ClonoTrace 1.0.0 | Method-selection guide**

A clone is a biological lineage, while an inferred clone is a set of observations
joined by a model. Similarity alone does not prove common descent. The right
method depends on the information available and the question being asked: native
software reproducibility, junction similarity, heavy-light consistency,
amino-acid identity or reconstruction of supplied reference annotations.

Start with explicit donor/sample boundaries, supported chain annotations and an
auditable input projection. Preserve cell pairing when it is part of the
question. Choose parameters before inspecting a favorable clone count, and
report exclusions as part of the result.

## Match the question to a model

| Question | Method | Required choice | Main interpretive limit |
| --- | --- | --- | --- |
| Reproduce qualified Change-O assignment | `assign --method changeo` | Threshold; ambiguity action | Specified 1.3.5 profile only; ambiguous-call grouping can depend on row order. |
| Reproduce qualified SCOPer nucleotide identity | `assign --method scoper` | Ambiguity action | Specified 1.6.0 `identicalClones` profile only. |
| Compare junctions under a linkage rule | `hierarchical-v2` | Linkage plus threshold or count | Fixed work units; linkage and tie policy define the partition. |
| Partition a declared affinity graph | `spectral-v2` | Sigma and per-component count | No automatic biological clone count; numerical refusals remain possible. |
| Retain heavy-light information | `paired-threshold-v2` | Heavy threshold and both cell policies | Heavy-first graph with exact light refinement; incomplete/conflicting cells need explicit handling. |
| Compare protein-level junction identity | `aa-identity-v2` | Direct AA or explicit translation | Synonymous nucleotide distinctions disappear. |
| Allow AA junction differences | `aa-hamming-v2` | AA source and threshold | Hamming substitution costs are uniform. |
| Reconstruct supplied reference annotations | `germline-v2` | Reference identity and D policy | Reconstruction, not de novo alignment or ancestor inference. |
| Explore threshold sensitivity | `auto-gap-stability-v2` | Optional explicit override | Uncalibrated gap estimate and deletion sensitivity. |
| Use a caller-defined distance/affinity model | `cluster-matrix` | Declared matrix, provenance and cut | Matrix validity does not establish biological validity. |

This manual provides complete commands, input profiles,
resource limits and output interpretation. The retained synthetic example inputs
show expected behavior; example parameters are not biological recommendations.

## Native compatibility is method-specific

The native Change-O and SCOPer targets intentionally use different ambiguity and
grouping rules. Do not force both into a common candidate graph to make their
outputs look alike. The compatibility claim concerns member partitions after
ClonoTrace's declared input policy, not all upstream inputs, internal candidate
groups or numeric clone labels.

Change-O uses gene-resolution normalized nucleotide Hamming single linkage with
an explicit inclusive threshold. Some ambiguous V/J patterns are order-sensitive;
reproducing that behavior requires preserving source order. SCOPer nucleotide
identity uses its qualified observed final-membership model, including tested
cross-length ambiguity bridges. Neither native profile establishes a universal
biological definition of clonality.

Current v2 bulk methods instead use exact donor/sample/locus/single V gene/single
J gene/analyzed-length work units. Ambiguous calls are excluded. These methods
have their own computational contracts; a v2 name does not broaden a native
compatibility claim.

## Choose a distance and a cut deliberately

Single linkage can connect observations through intermediate sequences even when
the endpoints exceed the threshold. Complete linkage limits the largest
cross-cluster distance at each merge. Average linkage weights original rows, so
duplicate observations affect its merge criterion. A stable tie policy can differ
from another library's equally valid tied choice.

Normalized Hamming requires equal analyzed lengths and assigns equal cost to
each unequal position. It does not model context-dependent somatic hypermutation,
sequencing quality, indels or evolutionary substitution probabilities. AA analysis
changes what counts as the same sequence. Forward translation requires an explicit
frame and standard code; optional IUPAC resolution accepts only unanimously
translated codons. It does not invent unknown residues.

An explicit cluster count is a requested partition resolution, not evidence that
the repertoire contains that many lineages. A threshold needs independent
scientific justification or a declared sensitivity analysis. No packaged value of
sigma, threshold or cluster count is a validated universal default.

## Interpret spectral inference as a graph model

Spectral v2 uses Gaussian affinities from normalized Hamming distances unless a
matrix is supplied directly. Sigma controls edge strength; an explicit affinity
cutoff and binary64 underflow can disconnect the graph. Positive-edge components
never merge, isolates stay singleton and requested counts saturate per component.

Canonical repeated-eigenspace and tie rules make the implementation deterministic.
They do not reveal a uniquely correct biological split in a symmetric graph.
Ill-conditioned embeddings and convergence failures stop the run. Numerical
tolerances and solver behavior are frozen, and no alternate solver is selected
after failure. Dense work units remain limited to 512 rows, adjustable downward.

The optimization evidence is computational: one-worker dense-512 median wall
time decreased from 5.03 to 1.15 seconds, or 4.39-fold, on the qualified host.
The secondary eight-worker comparison was 5.54-fold versus its own baseline.
Dense execution did not show meaningful worker scaling. More workers did not
cause the optimization, and these timings do not establish biological accuracy,
large-repertoire scaling or performance on another platform.

## Preserve the meaning of cells and references

Paired v2 first separates paired from heavy-only cells, builds heavy-chain
threshold components and refines them with exact known light signatures.
Heavy-only cells cannot bridge paired cells. Missing-heavy cells and conflicting
chains exclude the whole cell. Collapsing identical chain signatures is an
explicit computational policy, not proof of one molecule or one expressed allele.
Sequence-row counts, cell counts, reads and molecules are different quantities.

Germline v2 uses supplied reference IDs and annotations. It can retain reference
bases deleted from a query in a paired gapped alignment, while its query-coordinate
projection has a different purpose. Unknown rearrangement regions remain masked.
Neither representation is an inferred ancestral B-cell sequence. Reference choice,
allelic completeness and annotation correctness remain external scientific inputs.

## Treat automatic thresholding as exploratory

The gap estimator examines nearest-distinct-junction distances and requires a
unique largest gap. It reports every distinct-junction deletion estimate or
refusal, its acceptance fraction, a supported range and partition sensitivity.
This range describes deletion stability. It is not a confidence interval,
posterior interval, false-discovery rate or calibrated uncertainty in clonality.

A refusal without an override fails the CLI run; it does not quietly select a
different estimate. An explicit override retains all estimate/refusal diagnostics.
The absolute 128-distinct-junction bound and conservative work reservation remain
part of the method's computational boundary.

## Report what the evidence supports

Computational qualification covers recorded synthetic and public-data workflows
on macOS ARM64. Frozen Change-O/SCOPer observations, mathematical oracles,
numerical differential checks and repeated workflow results answer different
questions. Their agreement supports implementation behavior within those tests.
It does not supply a biological truth set, population validation, universal
upstream equivalence or superior predictive accuracy.

The retained public B-cell repertoires are finite case studies. Subsets selected
because a bounded method can run are conditional computational analyses; they do
not validate all rows of the original repertoire. Synthetic many-group stress
tests do not establish success on a single very large dense group.

Qualification does not extend automatically to other operating systems,
architectures, compilers or arbitrary input sizes. Memory grows with compact
records and work units; conservative limits are not a hard operating-system RAM
guarantee. The reference R restoration is a same-host content snapshot, not a
cross-platform hermetic source build.

For each analysis, report method and contract, package/executable identity,
input and reference hashes, scope/projection, parameter values, row and cell
accounting, exclusions, numerical diagnostics and relevant sensitivity checks.
Use complete member partitions when comparing methods. Canonical numeric clone
IDs are run-local labels and should not be carried across changed datasets as
biological identifiers.

## Remaining limitations

Cross-sample clone matching, allele inference, de novo alignment, general
ancestral reconstruction, flexible conflicting-chain resolution, calibrated
biological error rates and universal SHM-aware distances are unsupported claims.
Older `*-v1` methods remain separate historical contracts. Their limitations are described in the versioned contracts below.

ClonoTrace code is `Apache-2.0 OR MIT`; third-party software, datasets and references retain their own terms.


---

# ClonoTrace bounded v0.1 scientific contract

The macOS ARM64 qualification preserved native method-specific behavior. Contract identifier: `CT-0.1-native-1`. No biological truth follows from software partition parity.

## Input and accounting

Bulk IGH rows, UTF-8 plain TSV or gzip TSV, exact uppercase ACGT junctions.
Required columns: sequence_id, v_call, j_call, junction, donor_id, sample_id, locus.
Optional junction_length, when present, must equal the decimal actual length.
Calls are comma-separated `IGHV[A-Za-z0-9/-]+` or `IGHJ[A-Za-z0-9/-]+`, with
optional `*[0-9]+` allele suffix. Gene resolution removes only that suffix.
No allele is invented for a gene-only call. No whitespace/case normalization.
The profile is deliberately narrower than all AIRR/IMGT names.

Blank/duplicate headers, missing required columns, ragged rows, invalid UTF-8,
blank/duplicate IDs, blank donor/sample boundaries, active cell_id/CELL columns,
and existing output columns are whole-run errors. No results are published.
Quoted multiline TSV and stdin are unsupported; tabs and newlines are framing,
not escaped content. Literal quotes in annotations are preserved. CRLF and LF
are accepted; emitted framing is LF. A final row need not end in newline.
Original field bytes, order, and values are preserved in the assigned table.

Row exclusions in priority order: locus != IGH; empty/non-ACGT junction;
missing/out-of-profile V/J call; inconsistent supplied junction_length.
Excluded rows remain in output with blank clone_id and a stable reason.
Every successful run accounts for input_rows = retained + excluded, failed=0.
A structural, resource or pre-publication I/O failure fails the whole transaction
and emits no result directory. Exclusive rename is the commit point. A failed
parent durability sync after commit reports a warning and retains successful
publication; crash durability is not guaranteed by atomic visibility. Upstream launch/nonzero failures are computational failures,
never biological exclusions. The upstream tools accept some inputs excluded by
this explicit project preflight; native parity is claimed only after preflight.

Donor/sample/locus tuples are hard scope boundaries. No cross-sample matching.
Scope values are literal strings. Sequence-row counts do not mean cells or reads.
Original single-cell data require an explicit documented bulk projection; retain
cell_id as source_cell_id and clone_id as source_clone_id before analysis.

Change-O schema maps SEQUENCE_ID,V_CALL,J_CALL,JUNCTION,JUNCTION_LENGTH,LOCUS,
DONOR_ID,SAMPLE_ID to their lowercase names. Explicit `--map canonical=column`
can override each field. Input headers and values remain unchanged on output.
This is declared column mapping, not a complete Change-O annotation converter.

## Supported methods and exact options

| CLI method | Reference | Fixed settings | Exposed settings |
| --- | --- | --- | --- |
| changeo | Change-O 1.3.5 DefineClones | mode=gene, model=ham, norm=len, sym=avg, link=single, maxmiss=0 | action=set/first; mandatory threshold |
| scoper | SCOPer 1.6.0 identicalClones | method=nt, cdr3=FALSE, mod3=FALSE, max_n=0, fields=NULL, cell_id=NULL, summarize_clones=FALSE | action=set/first maps first=FALSE/TRUE |

Threshold is a finite decimal parsed as IEEE-754 binary64 in [0,1], with no
default. Compare `(mismatches as f64)/(length as f64) <= threshold`.
This binary64 boundary rule reproduces the observed oracle, including 1/3 and
adjacent representable thresholds; replacing it with real-rational decimal
comparison changes behavior. Hamming counts unequal ACGT positions. Single linkage
is connected components of all qualifying pair distances inside a Change-O group.
Zero threshold preserves exact identity, one connects every same-length group.

## Native Change-O observable grouping

The following independently inferred input/output model passed an
exhaustive, randomized, repeat and process-count corpus. It is not upstream code.
Within each scope and actual junction length, process rows in input order.
With first action, keep the first V and J gene; group matching keys exactly.
With set action, maintain existing groups in insertion order and mutable incoming
V/J gene sets. Traverse groups once: when incoming J intersects a group's J,
union that group's J into incoming J, even without V overlap. If incoming V also
intersects, collect that group and immediately union its V into incoming V.
Do not revisit earlier groups. Remove collected groups; append their records and
the incoming row as one group at the end. Unmatched groups keep their order.

This behavior is order-sensitive for some ambiguous annotations. Preserve it.
Do not sort input before grouping or replace it with static graph components.
Canonical numeric clone labels do not eliminate native partition order effects.

## Native SCOPer observable final partitions

Within each scope, build connections between row pairs with both intersecting
V-gene and J-gene sets, without a length restriction. With first action use only
the first V/J genes. Final partitions are connected components split by exact
junction string. This models observed final membership, including cross-length
ambiguity bridges. It is not a claim about SCOPer's internal candidate groups.
Internal execution work units are implementation details; diagnostics do not
label them as proven native SCOPer candidate groups.

## Determinism, output and resources

Clone IDs are canonical positive integers ordered by the smallest sequence_id
in each clone (UTF-8 byte order). They are run-local, not persistent biology IDs.
Assignment rows retain input order. Output is a new directory published atomically
with no replacement. Files: assignments.tsv, clones.tsv, diagnostics.tsv,
manifest.json. The manifest records input SHA-256, method, contract, configuration,
accounting, stage timings and environment. Scientific tables are byte-identical
across repeated runs and thread counts for a fixed input/configuration. Timings,
thread configuration and resulting manifest bytes may differ.

`--threads` explicitly bounds worker count; default 1. Workers process complete
method work units; grouping stays deterministic before parallel assignment.
`--max-rows`, `--max-junction-length`, and `--max-pairs` reject work exceeding
explicit bounds without partial publication. Memory is O(compact input records)
plus union/find and output labels; full annotation rows are re-streamed in pass 2.
No hard operating-system memory guarantee or universal scalability claim is made.
Input SHA-256 is verified again during emission to detect between-pass mutation.

## Evidence and deferred scope

Live upstream runtimes are not bundled. Cross-platform or hermetic upstream source rebuilds are not claimed.

Deferred: allele-level modes, average/complete/hierarchical/spectral clustering,
SHM/custom distance models, amino acids, paired-chain refinement, germlines,
non-IGH inference, missing/ambiguous nucleotides, stdin, zstd, automatic threshold
selection and cross-sample clone matching. No automatic biological cutoff.


---

# ClonoTrace Extended v1 contracts

The v1 extension methods are distinct from CT-0.1-native-1 and the native Change-O and SCOPer paths. They do not claim native equivalence or biological superiority.

## Shared profile and transaction

New `extend` command, explicit versioned `--method`, new output directory.
AIRR lowercase columns only; no implicit schema mapping. UTF-8 TSV/gzip framing,
unique nonblank sequence_id, donor_id/sample_id, no preexisting output fields.
Original field payload and input order are retained. Structural/resource/numerical
failure publishes nothing. Retained + excluded = all input rows; failed=0 on
success. Whole-run failure is never recoded as a biological exclusion. Bounded
rows, sequence length, conservative scope pair budget and per-group size are
checked before parallel work. No hard OS memory cap is claimed.

Bulk methods require locus=IGH; one unambiguous valid V/J call, numeric allele
suffix removed; uppercase ACGT junctions (AA direct mode uses junction_aa instead).
Work units are exact (donor,sample,locus,V gene,J gene,analyzed sequence length).
These deliberately stricter units are not upstream-native ambiguity groups.
Bulk active cell_id is refused: rename to source_cell_id for an explicit bulk
projection. Paired mode requires active cell_id. Gzip two-pass input SHA verification
and the native exclusive directory publication mechanism protect outputs.
IDs and members sort by UTF-8 sequence_id; labels follow minimum member ID.
Scientific tables repeat exactly across thread counts. Profiles and requested
configuration make manifest bytes variable. New diagnostics are in the manifest.

## hierarchical-v1

Purpose: sensitivity to single, complete or row-weighted average linkage for
same-gene, same-length junctions. Require explicit threshold binary64 in [0,1]
and explicit linkage. Hamming = unequal positions / length, binary64 inclusive
cut. Start one cluster per row; retain duplicate row weights. At every step
choose minimum (distance,min-ID-of-A,min-ID-of-B), orient A/B by their minimum
IDs, merge and sort members. Average sums original cross-row distances in sorted
member order then divides by cross-pair count; single=min, complete=max. Stop
when minimum distance exceeds threshold. Scalar bound: 128 rows/group.

References: SciPy public linkage mathematical documentation and independent
hand-calculated/tie fixtures. SciPy's unspecified alternative tie policy can
produce a different partition. Observed SCOPer average duplicate weighting also
differs: AAAA,AAAA,AAAC,AACC at .375 is one clone in SCOPer versus two here.
No sparse, approximate or duplicate-collapse optimization is qualified.

## spectral-v1

Purpose: deterministic exploratory binary graph separation, one bisection per
work unit, not automatic biological clone number. Explicit finite sigma >0.
Off-diagonal W_ij=exp(-((hamming(i,j)/sigma)^2)/2), W_ii=0,
with exactly that binary64 evaluation order. Normalize using separate square roots
of degrees so multiplying very small positive degrees cannot underflow first. Require every edge
strictly positive in binary64; underflow fails. Symmetric normalized Laplacian
L=I-D^(-1/2)WD^(-1/2). Singleton/all-identical groups remain one clone. Otherwise
compute all eigenpairs using deterministic cyclic Jacobi rotations, binary64,
canonical record order. Absolute off-diagonal stopping tolerance 1e-13; bounded
100 sweeps. Sort eigenvalues ascending. Require adjacent Fiedler eigenvalue gaps
>1e-8, residual Euclidean norm <=1e-9 and all vector components abs>1e-10.
Canonicalize vector sign by first maximum-absolute component. Split signs into
negative/nonnegative; refuse empty cuts. No k-means, randomness or recursion.
Group bound128. Eigensolver diagnostics record eigenvalues, vector, residual and
gap. NumPy LAPACK eigh is the independent numerical oracle; projectors remove
sign ambiguity. Full public data include genuine refused groups; qualification
of a derivative explicitly excludes whole unsupported groups and reports this.

## paired-exact-v1

Purpose: conservative exact paired receptor clonotypes, not lineage inference.
Cell identity is (donor_id,sample_id,cell_id), locus/sequence_id identify chains.
Every cell must have exactly one IGH and exactly one IGK or IGL row; every chain
must pass uppercase ACGT and single-VJ validation. Missing/multiple/invalid chains
exclude every row of that cell, without salvaging a partial pair. Blank cell ID
is a structural run failure. No implicit productive filtering or chain selection.
Group cells by donor/sample plus both (locus,V gene,J gene,junction) signatures;
each clone contains both chain rows of every member cell. Original cell, donor,
sample, locus and chain row identity remain in output. Manifest counts cells
separately from rows; clone row count is never called a biological cell count.

## aa-identity-v1 and aa-hamming-v1

Purpose: AA exact identity or unweighted single-linkage normalized AA Hamming.
Select direct junction_aa or explicitly request translation with frame 0/1/2
(number of leading nucleotide bases skipped). Standard NCBI code table1, forward
strand only, ordinary codon translation without start recoding. Validate whole
junction uppercase ACGT before offset; remaining sequence nonempty and divisible
by3; any stop excludes row. No partial-codon truncation, ambiguous-codon imputation,
reverse complement or inferred frame. Direct AA accepts only uppercase canonical
20 amino acids; *, X, B, Z, J, U, O, gaps and lowercase exclude. Equal analyzed AA
length groups only. AA Hamming needs explicit threshold inclusive binary64;
identity rejects threshold. Bio.Seq.translate public black-box observations of
all64 codons plus provider junction_aa are independent translation evidence.

## automatic threshold and germline reconstruction

Their additional contracts and oracle commit precede integration below. Native
threshold behavior remains frozen; no upstream-native reconstruction is claimed.

## Disposition and release

The historical v1 extensions are distinct from the native methods. Full scientific validity, arbitrary scale, cross-platform behavior and general upstream option coverage are unclaimed.

## auto-gap-v1 (additional frozen contract)

Purpose: diagnostic exploratory cutoff from distinct-junction nearest-neighbor
separation, never a universal or calibrated biological threshold. In each exact
bulk work unit, sort unique junctions lexicographically, compute each one's minimum
normalized Hamming distance to another distinct junction. Sort unique minimum
distances; require at least two. Choose the unique largest adjacent gap and its
binary64 midpoint (lo+hi)/2. Equal largest gaps refuse. No distinct neighbor or
constant distances refuse. All per-group selected thresholds, nearest distances,
gaps and reasons are recorded. `--threshold` explicitly overrides every group,
while still computing/reporting its diagnostic estimate/refusal. The selected
threshold drives ordinary unweighted single-linkage row assignment. This is an
independent deterministic heuristic; no SHazaM density/mixture or native claim.
Bound128/group. Full public input without override is expected to refuse due to
singleton groups. Explicit override remains available without hiding that fact.

## germline-mask-v1 (additional frozen contract)

Purpose: conservative reference-derived V/J sequence in query coordinates,
masking all other positions with N, with no inference of D or ancestral lineage.
Require sequence, v_call,j_call,v_cigar,j_cigar, donor/sample/locus/sequence_id.
Active cell_id is accepted as annotation; no clones are assigned. IGH/IGK/IGL
only, exact FASTA identifiers (no allele guessing). User supplies canonical
FASTA (unique bare headers, uppercase ACGT), explicit reference ID/version and
source/provenance text. Manifest records file path, raw SHA256 and reference
configuration. The reference is external; it is not embedded in the executable.

Interpret positive CIGAR operations: S advances query only, N advances reference
only, M/= /X advance both and project reference bases, including mismatches.
No I/D/H/P or zero/empty operations. Require one contiguous M/= /X block
with optional S then N at each end; internal/out-of-order clipping is excluded. Bounds are checked, and each segment must
align at least one base. Optional unconsumed query tail stays N; overlap between
V/J projections excludes the row even if bases agree. Missing/ambiguous FASTA
calls, unsupported CIGAR, invalid sequence, bounds and overlap are explicit row
exclusions. Corrupt/duplicate/empty reference files are whole-run errors. Output
is reconstruction.tsv plus diagnostics and manifest, not a clone assignment.
No alignment, IMGT gap placement, junction/D reconstruction, consensus ancestor
or Change-O CreateGermlines compatibility is claimed.

Public reference validation uses external 10x GRCh38-alts-ensembl5.0.0 V/J FASTA;
an explicit adapter omits genes with multiple distinct reference sequences.
Original reference SHA, adapter SHA and source headers are recorded. Reference
terms differ from the public dataset CC-BY4 license; reference FASTA/derived
biological sequence assets are not included in source or recommended release.


## Resource and execution details

The 128-row bound applies only to hierarchical-v1, spectral-v1 and auto-gap-v1;
AA graph inference is bounded by conservative retained scope/length pair counts.
Paired applies max-pairs to all source rows within each donor/sample, before
whole-cell filtering. Germline projection has no pair budget because it computes
no pairs; its max-pairs manifest configuration is inapplicable. Paired signature
hashing and germline projection are scalar even if a larger thread count is
requested; manifests say execution=scalar. Independent bulk work units can use
1..1024 explicit workers. The only new retained optimization is exact sequence-map collapse in AA identity,
selected after scalar profiling and qualified by repeated exact ablation; native
measured optimizations remain confined to native paths.


---

# CT-complete-2 technical scientific contract

CT-0.1-native-1 and the older method names and semantics remain distinct. CT-extended-1 v1 outputs remain separately identified; nothing silently falls back between methods.

## Library and workflow integration

clonotrace-extended exposes checked AA/paired scientific APIs and a reusable
Config/run workflow API returning a machine-readable result, without printing or
process exit. Native APIs stay in their existing crates. Hierarchy/spectral/
germline/threshold v2 expose typed public modules in the existing scientific
crates. The CLI is an argument adapter. Structural and resource errors return
errors; all CLI failures exit2 and publish nothing. Scientific row exclusions
remain explicit rows. A numerical failure is a whole-run failure.

Plain/gzip strict UTF-8 AIRR framing, unique sequence_id, literal donor/sample,
per-row locus, immutable-input SHA and exclusive atomic output persist. All
members canonically order by bytewise sequence_id before algorithms; clone IDs
follow minimum member ID. Every successful row is retained or excluded, failed0.
Bulk groups are exact donor/sample/locus/single V gene/single J gene/analyzed
sequence length; no ambiguous-call propagation or upstream-native grouping claim.
Allele suffix removal follows the established numeric-only gene profile. Active
cell_id requires paired or reconstruction mode; bulk needs deliberate source
annotation renaming. Source fields/order remain unchanged.

New manifests use CT-complete-2 with explicit schema identifiers and separate
implementation qualification versus biological evidence. Unchanged assignment
columns remain CT-assignment-1; derived canonical AA remains CT-aa-derived-1.
Germline-v2 appends distinct alignment/projection columns, schema CT-germline-2.
Existing method manifests retain CT-extended-1 and scientific tables byte parity.

## hierarchical-v2

Original row-weighted agglomeration, single/complete/average, Hamming fractions in
binary64, NT or explicitly selected canonical AA alphabet. Complete dendrogram
records left/right nodes, merge height and row size. Stable ties use canonical
minimum member IDs; average sums original cross-row distances in sorted member
order. Choose threshold inclusive [0,1] OR positive requested cluster count per
work unit; requested count saturates at group size, explicitly diagnosed. Public
matrix API accepts any finite nonnegative exactly symmetric zero-diagonal matrix;
its threshold is finite nonnegative (CLI Hamming threshold remains [0,1]). Invalid
or forged dendrograms fail even when an early cut could hide the invalid node.
Default and absolute dense row limit512; user can choose a lower limit. No
approximate nearest-neighbor pruning, duplicate weighting change or sparse merge
substitute. SciPy documentation and independent direct agglomeration define the
oracle; eight historical tied choices differ legitimately from SciPy defaults.

## spectral-v2

Purpose: explicit-k normalized spectral clustering for a declared graph, not
biological auto clone-number inference or SCOPer novj/vj equivalence. NT or checked
AA normalized Hamming. Required positive finite sigma and positive requested
clusters-per-component. Gaussian off-diagonal exp(-((d/sigma)^2)/2), zero diagonal.
Optional explicit affinity distance cutoff [0,1] creates zero edges where d>cutoff;
binary64 exponential underflow also produces zero edges and is counted. This is
an explicit versioned affinity contract; v1 retains its original underflow refusal.

Positive-edge connected components never merge. Isolates are singleton groups.
Effective k=min(requested,component rows), with saturation diagnostics. k1 or k=n
has the exact trivial result and needs no eigensolver. Otherwise symmetric
normalized Laplacian uses separate degree square roots. Deterministic cyclic
Jacobi,100sweeps, off-diagonal tolerance1e-13, eigen residual<=1e-9. Eigenvalues
sort ascending. Blocks within1e-10 of block start use projected canonical
coordinate axes and two-pass Gram-Schmidt; first block pins the analytic
sqrt(degree)/sqrt(total degree) stationary vector. This explicit numerical
symmetry-breaking policy is part of the method, not an unobserved upstream rule.
First k basis columns, pre-row-normalization norm>1e-8, then unit row embeddings.
Ill-conditioned rows fail the run rather than amplify eigensolver noise.

Deterministic farthest-point centers start at canonical first row; squared-distance
ties use tolerance1e-12 and earliest canonical row/center. Lloyd updates sum in
canonical order. Empty centers take the largest-error row from a nonsingleton
cluster, ties earliest canonical ID; no row is lost. Bounded iterations(default1000,
max10000), stable assignments stop; failure to converge is a run failure. Record
components, k, eigenvalues/eigenspace blocks, residual, projector, iterations,
inertia, zero/underflow edges, conditioning and all numerical tolerances. Dense
limit512, adjustable downward. NumPy LAPACK plus an independent projector/Lloyd
oracle establishes computational concordance; no universal biological guarantee.

## paired-threshold-v2

ChainInput preserves row/sequence/cell/donor/sample/locus/calls/junction identity.
Structural IDs/scope/cell errors fail the run. Invalid chain biology excludes the
whole cell. Missing heavy always excludes. Conflicting heavy or light signatures
always exclude the whole cell. No abundance-based or arbitrary chain selection.
Require explicit missing-light policy exclude or heavy-only-separate and explicit
multiple-chain policy exclude or collapse-identical. Under collapse-identical,
only equal (locus,V gene,J gene,NT junction) signatures collapse computationally;
all original chain rows remain in the resulting cell clone. This is not proof
that redundant rows represent one physical molecule or one expressed allele.

Partition paired versus missing-light heavy-only cells BEFORE heavy graph building;
missing-light cells cannot bridge paired clones. Within each class and exact
heavy donor/sample/V/J/length work unit, normalized NT Hamming <=explicit threshold
defines single-linkage heavy components. Refine paired heavy components by exact
light locus/V gene/J gene/NT junction signature. Different known lights can be
heavy-graph intermediates before refinement; that is the explicit heavy-first
model. Heavy-only cells form separately flagged clones. No hidden fallback to a
heavy-only result. Preserve cell and row counts, exclusions and per-cell policy
outcomes. Worker-equivalent independent work units; max-pairs checked before work.

## AA v2 and reusable AA APIs

aa-identity-v2 and aa-hamming-v2 retain explicit direct versus translation/frame.
Direct requires canonical20 uppercase AA; no inferred costs for B/Z/X/J/U/O/*.
Forward standard code1 translation, frame0/1/2, whole-string validation, no start
recoding or partial-codon truncation. Stops exclude. Default ambiguity=reject.
Explicit resolve-unanimous enumerates IUPAC DNA codon resolutions: emit a canonical
AA only if every resolution yields that same residue; all-stop=>stop exclusion,
mixed or nonunanimous=>ambiguous-codon exclusion. No X substitution or silent
fallback. Invalid/lowercase symbols exclude. Unanimous policy is not accepted for
direct AA input. Same checked AA source also serves hierarchical/spectral v2 when
--alphabet aa is explicit. Identity and equal-length Hamming/single linkage have
validated public library APIs; invalid thresholds/lengths return errors.

## germline-v2

Explicit --germline-mode mask-d/include-d. Exact supplied canonical FASTA IDs,
reference ID/version/source/raw SHA, no allele inference. V/J required; include-d
requires D call/CIGAR for IGH and forbids D annotations for light loci. mask-d
ignores D fields explicitly, preserving source annotations. Unknown junction,
unassigned and masked D positions remain N; no invented ancestor or junction.

AIRR CIGAR: optional S then N ends around nonempty aligned body(M/=X/I/D); no
internal clips, no zero operators, no mixing M with =/X. S consumes query, N
consumes reference outside alignment; M/=X consume both and copy reference,
= and X must agree with supplied query/reference meaning. I consumes query,
query projection N and aligned germline '-'; D consumes reference, aligned query
'-' with retained reference base. Every segment must align bases, bounds hold,
segments cannot overlap even through query insertion slots. Fail closed for
unsupported H/P, missing/ambiguous references, conflicts or malformed annotations.

Emit query-coordinate projection plus paired query/germline gapped alignment,
retaining reference bases deleted from query. This is annotated reference
reconstruction, not de novo alignment, genotype discovery, clone consensus,
ancestral reconstruction or general CreateGermlines parity. Segment diagnostics
record exact call/CIGAR/query/ref ranges and match/insertion/deletion counts.
--max-alignment-length(default20000) bounds deletion-expanded output BEFORE
allocation; overflow/limit is whole-run resource failure, not row exclusion.
Invalid reference file is whole-run failure. Direct public map APIs validate the
referenced sequence alphabet. Reference biological assets remain external under
reference-specific rights, not assumed to inherit dataset CC-BY4.

## auto-gap-stability-v2

Exact distinct-junction nearest-neighbor largest-gap estimate remains transparent
and uncalibrated. At least two distinct nearest distances and a unique largest
gap are required; otherwise diagnostic refusal. Compute deterministic deletion of
each DISTINCT junction, re-estimate fully, report every omitted junction,
threshold/refusal, acceptance fraction and finite minimum/maximum stability range.
This is deletion sensitivity, NOT a confidence/credible interval or biological
error rate. Compute resulting clone counts/partition disagreement at supported
range endpoints and base/override threshold. Override remains explicit and does
not suppress any estimate/refusal/stability diagnostic. No alternate estimator.

Bound128distinct junctions and4000000directed distance evaluations by explicit
configuration before computation. An exact symmetric distance cache reuses each
pair across the base estimate, deletions and sensitivity partitions. The original
conservative scalar reservation remains the preflight limit; actual
distance_evaluations is an execution counter and can decrease after optimization.
The checked library and CLI both enforce the absolute128-distinct bound.
Constant/no-neighbor/tied-gap evidence with no
override fails the whole run after diagnostics are constructed. Per-group work
can parallelize; row abundance does not change distinct-junction resampling.

## Qualification boundary

"Production implementation" means complete checked library/CLI implementation of
these explicit computational contracts on the qualified platform. Biological
model validation is a separate evidence dimension. In particular, automatic gap
thresholds remain scientifically exploratory even when engineering is complete.
Unsupported submodes are not exposed as approximate fallbacks. v0.1 and v1 regression behavior remains
subject to qualification. Profile complete workflows before optimization; retain only exact
scientific, row-accounting and failure-semantic parity with measured benefit.

## Declared matrix CLI and schema

The checked generic distance/affinity APIs also have cluster-matrix CLI integration.
Input JSON schema CT-matrix-1 requires kind=distance/affinity, unique nonblank ids,
explicit single donor/sample/locus scope, nonblank provenance and matrix. Unknown
or duplicate JSON fields fail; finite nonnegative exact symmetric zero-diagonal
matrices only. Canonically reorder IDs and both matrix axes before inference.
Distance matrices use hierarchical-v2 and arbitrary finite nonnegative threshold
units (or explicit count), affinity matrices use spectral-v2 with explicit clusters
per component. No Gaussian construction/sigma is applied to supplied affinities.
Input byte(default32000000, max256000000) and row(default/max512) bounds precede
inference. Output CT-matrix-assignment-1 preserves every declared ID and scope;
input JSON checksum and provenance identify the distance source. This is a single
matrix work unit with scalar algorithm execution. It does not infer that externally
supplied costs are a validated biological distance model. The AIRR CLI remains the
qualified direct NT/AA sequence path; this schema makes matrix APIs runnable.


---

# CT-complete-2 build, API and CLI reference

This is the technical interface for the bounded v2 computational contracts in
the complete v2 scientific contract above. CT-0.1-native-1 (`assign`/`validate`)
and explicitly named `*-v1` workflows retain their existing scientific behavior.
Select a versioned method deliberately; there is no fallback between contracts.
Qualification is local macOS ARM64 computational qualification. It does not
establish biological superiority, general upstream equivalence or another
platform's numerical behavior. Automatic gap thresholds remain exploratory.

## Build and validation

The pinned toolchain is Rust 1.98.0 with rustfmt and Clippy. Cargo dependencies
must already be in the local cache for offline commands. From the project root:

```sh
export CARGO_TARGET_DIR="$PWD/target"
cargo build --release --locked --offline
cargo test --workspace --locked --offline
```

The commands above build the release and run the workspace Rust tests.
Live upstream runtime environments are not bundled in this source release. Synthetic reference partitions were tested without upstream software. Historical
runtime qualification is host-specific, not a portable hermetic R rebuild.

`target/release/clonotrace` is the newly built CLI. Pre-optimization executables
are preserved separately and are not required for building or testing this release.
Inspect argument spelling with:

```sh
target/release/clonotrace extend --help
target/release/clonotrace cluster-matrix --help
```

## Input and publication contract

`extend` reads strict UTF-8 tab-separated AIRR records, plain or gzip. Source
fields and row order are preserved. All modes require globally unique nonblank
`sequence_id`, nonblank literal `donor_id`/`sample_id`, and `locus`, `v_call`,
`j_call`. Duplicate headers, ragged records and reserved output columns fail.
Bulk modes require `junction` (NT) or `junction_aa` (direct AA), and support IGH.
Bulk work units are exact donor/sample/locus/V gene/J gene/analyzed length.
Only numeric allele suffixes are removed. Ambiguous V/J calls are excluded.
Rename an intentionally preserved source cell annotation to `source_cell_id`
before a bulk run; active `cell_id`/`CELL` denotes a different input profile.
Paired and germline modes support IGH/IGK/IGL and retain cell annotations.

Every `--out` must be a fresh directory with an existing parent. Publication is
exclusive and atomic, never replacement. Input SHA-256 is checked across passes.
Structural/resource/numerical failures exit 2 and publish no result directory;
scientific row exclusions remain explicit output rows. Every successful run
accounts for all input rows as retained or excluded, with failed=0. Algorithms
use canonical bytewise sequence-ID order; output labels follow minimum member ID.
Worker counts do not change scientific tables. Profile/timing metadata can differ.

## Commands

The examples use prepared input files and fresh output directories. Create the
parent `out/` once. Paths and reference metadata are examples supplied by the
caller, not bundled biological data.

```sh
ct=target/release/clonotrace
mkdir -p out

# All qualified linkage methods; inclusive normalized-Hamming threshold.
for link in single average complete; do
  "$ct" extend --input bulk.tsv --out "out/hier-$link" \
    --method hierarchical-v2 --linkage "$link" --threshold .15 --threads 4
done
# Count cut is per work unit; never combine --clusters with --threshold.
"$ct" extend --input bulk.tsv --out out/hier-count \
  --method hierarchical-v2 --linkage average --clusters 2

# Gaussian affinity, positive-edge components, explicit k per component.
"$ct" extend --input bulk.tsv --out out/spectral \
  --method spectral-v2 --sigma .25 --clusters-per-component 2 \
  --affinity-cutoff .4 --threads 4

# Canonical direct AA, or standard-code translation with explicit frame.
"$ct" extend --input amino.tsv --out out/aa-direct \
  --method aa-identity-v2 --aa-input direct
"$ct" extend --input bulk.tsv --out out/aa-hamming \
  --method aa-hamming-v2 --aa-input translate --frame 0 --threshold .15
"$ct" extend --input iupac.tsv --out out/aa-unanimous \
  --method aa-identity-v2 --aa-input translate --frame 1 \
  --ambiguity resolve-unanimous
# AA can also define hierarchical/spectral Hamming distances.
"$ct" extend --input amino.tsv --out out/hier-aa \
  --method hierarchical-v2 --alphabet aa --aa-input direct \
  --linkage complete --clusters 2

# cell_id is required; both cell policies are mandatory.
"$ct" extend --input paired.tsv --out out/paired \
  --method paired-threshold-v2 --threshold .15 \
  --missing-light exclude --multi-chain collapse-identical --threads 4

# Germline input additionally requires sequence, v_cigar and j_cigar.
for mode in mask-d include-d; do
  "$ct" extend --input annotated.tsv --out "out/germline-$mode" \
    --method germline-v2 --germline-mode "$mode" \
    --reference-fasta supplied-vdj.fa --reference-id supplied-reference \
    --reference-version exact-version --reference-source exact-source --threads 4
done

# Automatic estimate: a refusal in any work unit fails the run.
"$ct" extend --input bulk.tsv --out out/auto \
  --method auto-gap-stability-v2 --threads 4
# An explicit override preserves estimator/refusal/stability diagnostics.
"$ct" extend --input bulk.tsv --out out/auto-override \
  --method auto-gap-stability-v2 --threshold .15 --threads 4
```

Hierarchy supports row-weighted single/average/complete linkage with a full
stable dendrogram; count cuts saturate at work-unit size. NT/AA Hamming thresholds
are finite in [0,1]. Dense groups are exact and bounded, not approximate.

Spectral `--sigma` must be finite and positive. `--affinity-cutoff` is optional
and removes edges whose distance exceeds the supplied value in [0,1]. Positive
components never merge; isolates are singletons and requested k saturates at
component size. Explicit eigenspace/tie/conditioning policies and bounded Lloyd
iterations are recorded. Underflow edges are diagnosed; conditioning or
convergence failures fail the run. There is no biological automatic-k claim.

Direct AA requires uppercase canonical20 residues. Translation uses forward
standard code 1, frame 0/1/2, whole-sequence validation and complete codons. Stops,
partial codons and invalid symbols exclude. Default ambiguity is `reject`.
`resolve-unanimous` accepts an IUPAC codon only when every concrete expansion
encodes the same canonical residue; it never substitutes X. Direct AA accepts
neither `--frame` nor `--ambiguity`. Spectral AA uses the same explicit AA source
options as hierarchy, with its own required sigma/k options.

Paired inference builds heavy single-linkage graphs within donor/sample/heavy
V/J/length and paired-state classes, then refines by exact known light signature.
`--missing-light` is `exclude` or `heavy-only-separate`; missing-heavy cells always
exclude. `--multi-chain` is `exclude` or `collapse-identical`. The latter permits
only identical locus/gene/junction signatures, preserving every original row;
conflicting chains still exclude the whole cell. Heavy-only cells cannot bridge
paired cells. This is an explicit heavy-first model, not an ancestry inference.

Germline reconstruction uses exact supplied FASTA IDs and records the raw
reference SHA with ID/version/source. V/J are required. `include-d` requires
`d_call`/`d_cigar` for IGH and forbids D annotations on light chains; `mask-d`
deliberately ignores D annotations. Unknown junction and masked positions are N.
AIRR CIGAR permits terminal S/N around an M/=X/I/D body; M cannot mix with =/X,
and =/X must agree with query/reference. I consumes query and yields germline
`-` in the alignment, N in the query projection. D consumes reference and yields
query `-` with the deleted reference base retained. Unsupported H/P, bounds,
overlaps and ambiguous/missing references fail or exclude as diagnosed. This
reconstructs supplied annotations; it does not align de novo or infer alleles,
ancestors or clone consensus.

Automatic selection uses nearest-distinct-junction distances and a unique
largest gap. It fully re-estimates after deleting each distinct junction and
reports refusals, acceptance counts, a supported range and partition sensitivity.
This is deletion stability, not a confidence interval or calibrated biological
error rate. Neither a failed base estimate nor a missing override is silently
replaced by a leave-one-out estimate. `distance_evaluations` counts actual Hamming
calls; exact caching can reduce it. The conservative reservation and
`--max-distance-evaluations` rejection boundary remain unchanged.

## Resources and output schemas

| Control | Default | Bound/meaning |
| --- | ---: | --- |
| `--threads` | 1 | 1..1024; independent work units, not a promise of speedup |
| `--max-rows` | 1000000 | Whole input bound |
| `--max-junction-length` | 10000 | Source analyzed sequence/query byte length |
| `--max-pairs` | 100000000 | Conservative pair work bound before inference |
| `--max-group-rows` | 512 | Hierarchy/spectral only; 1..512 |
| `--max-iterations` | 1000 | Spectral only; 1..10000 |
| `--max-alignment-length` | 20000 | Germline expanded alignment, checked before allocation |
| `--max-reference-bytes` | 100000000 | Germline only; 1..1000000000 |
| `--max-distinct` | 128 | Automatic threshold only; absolute 1..128, including library API |
| `--max-distance-evaluations` | 4000000 | Per automatic work-unit conservative reservation; zero allows only zero-work groups |

Successful v2 manifests identify `CT-complete-2` / `CT-run-2`, method/configuration,
input hash, resource limits, row accounting, diagnostics and profile timings.
Assignment/clones schemas remain `CT-assignment-1` / `CT-clone-summary-1`;
AA adds `CT-aa-derived-1`, paired adds `CT-cells-2`. Germline emits
`CT-germline-2` query projection and paired alignment columns plus
`CT-germline-segments-2`. These identifiers separate changed derived products
from the frozen v0.1/v1 outputs. Timing/profile fields are not scientific parity
keys. Read manifest diagnostics rather than inferring status from clone counts.

## Declared matrix input

`cluster-matrix` is a single scalar matrix work unit. Supply JSON schema
`CT-matrix-1` with `kind` (`distance`/`affinity`), unique nonblank `ids`, explicit
`scope` (`donor_id`, `sample_id`, B-chain `locus`), nonblank `provenance`, and a
square finite nonnegative exactly symmetric zero-diagonal `matrix`. Unknown or
duplicate fields fail. Both axes are reordered with canonical IDs before work.
For example, this is a complete distance input:

```json
{"schema":"CT-matrix-1","kind":"distance","ids":["b","a"],"scope":{"donor_id":"d","sample_id":"s","locus":"IGH"},"provenance":"caller-declared costs, version 1","matrix":[[0,2],[2,0]]}
```

```sh
"$ct" cluster-matrix --input distance.json --out out/matrix-hier \
  --method hierarchical-v2 --linkage average --threshold 2
"$ct" cluster-matrix --input affinity.json --out out/matrix-spectral \
  --method spectral-v2 --clusters-per-component 2
```

Distance thresholds are finite nonnegative values in the supplied matrix's
units; count cuts are also supported. Supplied affinity weights are used
directly, without sigma/Gaussian conversion. The CLI has no matrix `--threads`
option. Default/max row bound is 512; default input limit is 32000000 bytes,
adjustable with `--max-input-bytes` through 256000000. Output schema
`CT-matrix-assignment-1` preserves every declared ID/scope and records checksum
and provenance. Caller-supplied costs are not thereby biologically validated.

## Rust library surfaces

Use workspace path dependencies for the relevant crate. `clonotrace-extended`
provides a reusable IO workflow independent of CLI parsing:

```rust
use clonotrace_extended::{run, Config, Method};
fn example() -> Result<serde_json::Value, String> {
    let mut c = Config::new("bulk.tsv", "out/library", Method::HierarchicalV2);
    c.linkage = Some("average".into());
    c.threshold = Some("0.15".into());
    c.threads = 4;
    c.validate()?; // Configuration only; run validates input and output paths.
    run(c)        // Publishes products, returns summary; does not print or exit.
}
```

Checked algorithm APIs return row-index groups/results or errors. Callers of
matrix/sequence hierarchy and spectral APIs must establish canonical bytewise
record-ID order themselves; the CLI/workflow does this. The AA/paired APIs retain
independent duplicate observations. Algorithm-only callers establish their own
scope units; only workflow APIs parse files and publish manifests.

| Crate/module | Public entry points |
| --- | --- |
| `clonotrace_extended::amino` | `AaSequence::new`, `translate(TranslationConfig)`, `normalized_hamming`, `identity`, `single_linkage`; `AmbiguityPolicy::{Reject,ResolveUnanimous}` |
| `clonotrace_extended::paired` | `ChainInput`, `Config`, `infer`, `infer_threads`; `MissingLightPolicy` / `MultiChainPolicy`; checked exclusions, cells and work units in `Inference` |
| `clonotrace_cluster::v2` | `hamming_matrix`, `build(Linkage, limit)`, `Dendrogram::cut_threshold`, `cut_count`; `Linkage::{Single,Average,Complete}` |
| `clonotrace_spectral::v2` | `cluster(sequences, &Options)` or `cluster_affinity(matrix, k, row_limit, iteration_limit)`; component/numerical diagnostics in `Clustering` |
| `clonotrace_cluster::threshold` | `analyze(sequences, override_threshold, Config)`; `Analysis` includes base, deletion replicates, selection, sensitivity and actual/reserved work |
| `clonotrace_germline::v2` | `SegmentInput`, `DPolicy`, `reconstruct` or `reconstruct_bounded`; supplied `BTreeMap<String, Vec<u8>>` reference; projection/alignment/segments in `Reconstruction` |
| `clonotrace_extended::matrix` | `MatrixInput::canonicalize(limit)`, `run(Config)`; strict schema and transactional matrix workflow |

Threshold `analyze` returns scientific refusal as `Analysis` with
`selected_threshold=None`, allowing a library caller to inspect it. Structural
or resource invalidity is an error. The CLI's no-override refusal is a whole-run
failure. The absolute distinct bound remains 128 even when constructing a
library configuration directly. Full mathematical tolerances, scope rules and
unsupported submodes remain in the scientific contract ; these examples do not expand those boundaries.
