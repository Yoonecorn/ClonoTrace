# ClonoTrace

ClonoTrace 1.0.0 is a deterministic Rust toolkit for B-cell receptor clonal inference. It supports specified Change-O and SCOPer compatible behavior and hierarchical, spectral, paired-chain, amino-acid, germline and threshold analyses. The release was qualified on macOS ARM64; computational agreement does not establish biological accuracy.

## Build and quick start

Rust 1.98.0 is pinned. Populate the dependency cache with `cargo fetch --locked` when needed, then run:

```sh
export CARGO_TARGET_DIR="$PWD/target"
cargo build --release --locked --offline
target/release/clonotrace --version
target/release/clonotrace validate --input examples/repertoire.tsv
mkdir -p out
target/release/clonotrace assign --input examples/repertoire.tsv --out out/changeo --method changeo --threshold .25
target/release/clonotrace extend --input examples/v1/bulk.tsv --out out/spectral --method spectral-v2 --sigma .25 --clusters-per-component 2
```

Output paths must be new and have an existing parent. The native example retains four IGH rows and excludes one light-chain row. The spectral example assigns `a,b` and `c,d` to two groups. See the [user manual](docs/user-manual.md) for methods, input requirements, limits, tests and interpretation.

## Citation and license

Yoon Sing Yap, *ClonoTrace: deterministic clonal inference and alternative receptor models for B-cell repertoires*. [Manuscript](publication/ClonoTrace-Manuscript-v1.0-2026-11-04.pdf). Repository: [https://github.com/Yoonecorn/ClonoTrace](https://github.com/Yoonecorn/ClonoTrace). The manuscript refers to supporting data that are preserved separately and are not included in this repository.

Original ClonoTrace code is Apache-2.0 OR MIT; see [LICENSE](LICENSE) and [NOTICE](NOTICE).
