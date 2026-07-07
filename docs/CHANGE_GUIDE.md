# OutliMer Change Guide

This guide explains the changes made during the packaging and extension pass
that moved OutliMer from a prototype-style repository toward a PyPI and
bioconda-ready command line tool.

## Goals

The work had four main goals:

- Make the project installable as a normal Python package.
- Clean the package source so wheels and source distributions do not include
  generated analysis outputs.
- Stabilize the command line behavior and fix known logic errors.
- Add the next practical extensions for bioinformatics users: metadata,
  sourmash signature input, caching, HTML reports, workflow examples, tests,
  release automation, MultiQC output, RO-Crate metadata, profile utilities,
  run comparison, and sparse feature export.

## Repository Layout

The project now uses a cleaner `src/` layout:

```text
.
├── .github/workflows/
│   ├── ci.yml
│   └── release.yml
├── docs/
│   └── CHANGE_GUIDE.md
├── packaging/bioconda/
│   └── meta.yaml
├── src/OutliMer/
│   ├── __init__.py
│   ├── OutliMer.py
│   ├── classification.py
│   └── profiles.py
├── tests/
├── workflows/
├── MANIFEST.in
├── README.md
└── pyproject.toml
```

Generated CSV and PNG outputs were removed from `src/OutliMer/...` so they are
not part of the importable package. The tiny smoke-test dataset under
`examples/toy/` is included in source distributions; bulky generated example
outputs remain ignored and excluded.

## Packaging Changes

`pyproject.toml` was rewritten as the source of truth for packaging metadata.
It now includes:

- PEP 621 project metadata.
- Runtime dependencies:
  `sourmash`, `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`,
  and `seaborn`.
- Development extras:
  `build`, `pytest`, `ruff`, and `twine`.
- Console entry points:
  `outlimer`, `OutliMer`, and `outlimer-classify`.
- Setuptools package discovery configured for a non-namespace package.

`src/OutliMer/__init__.py` was added so setuptools discovers exactly the
`OutliMer` package rather than treating output folders as namespace packages.

`MANIFEST.in` was added to include useful files in sdists while excluding
generated artifacts and local cache files.

The `outlimer-profile` console script was added for profile management.

## Main CLI Changes: `outlimer`

The main CLI in `src/OutliMer/OutliMer.py` now has a more reliable execution
surface.

### Default Outputs

`outlimer` now writes these by default:

- `kmer_report.csv`
- `export_kmers/top_union_summary.csv`
- `run_manifest.json`
- `sketch_cache/` unless disabled

Detailed hash matrices and heavy per-sample exports still require
`--full-exports`.

Optional outputs:

- `--multiqc` writes `outlimer_summary_mqc.json`.
- `--ro-crate` writes `ro-crate-metadata.json`.

### Argument Validation

The CLI now fails early for invalid values such as:

- `--threads 0`
- `--map-threads 0`
- negative `--top-per-sample`
- invalid `--outlier-threshold`
- `--profile-update` without `--profile-name`
- `--write-db-each` without `--incremental`
- `--no-cache` combined with `--cache-dir`

Bad database or profile inputs now stop the command instead of silently falling
back to cohort-vs-cohort comparison.

### Deterministic Sample Discovery

Input files and paired sample prefixes are sorted deterministically. This
matters especially for incremental database mode, where processing order affects
the comparison baseline.

### Fixed or Wired CLI Flags

Previously inert or misleading flags were wired up:

- `--report` now defaults to `kmer_report.csv`.
- `--shared-summary` writes a per-sample shared-hash summary.
- `--top-shared` limits the shared summary.
- `--export-unique-shared` controls per-sample unique/shared hash files.
- `--force-top` now fails if a sample cannot provide the requested number of
  hashes.

The full-export fast path no longer writes and then later overwrites the same
CSV output.

## Sourmash Signature Input

`outlimer` now accepts existing sourmash signatures:

```bash
outlimer \
  --input-dir signatures/ \
  --input-type signature \
  --out-dir outlimer_results
```

Supported extensions:

- `.sig`
- `.sig.gz`

Signature input loads hash abundances from sourmash signatures instead of
reading FASTQ or FASTA. Sequence export options are rejected in signature mode
because `.sig` files do not contain recoverable k-mer strings.

The implementation uses `sourmash.load_file_as_signatures` when available and
falls back to the documented JSON signature loader path.

## Sketch Caching

Sketch caching was added to avoid reprocessing unchanged samples.

By default:

```text
<out-dir>/sketch_cache/
```

Each cache entry stores:

- sample name
- absolute input paths
- file sizes
- file modification times
- input type
- k-mer size
- sourmash scaled value
- hash counts

Use:

```bash
outlimer --no-cache ...
```

to force recomputation, or:

```bash
outlimer --cache-dir shared_cache ...
```

to share cache entries across runs.

## Sample Metadata

Both CLIs now accept sample metadata:

```bash
outlimer \
  --input-dir reads/ \
  --input-type paired-fastq \
  --metadata samples.tsv \
  --metadata-sample-column sample \
  --out-dir outlimer_results
```

The metadata file can be CSV or TSV. The default sample-name column is
`sample`.

In `outlimer`, metadata columns are appended to `kmer_report.csv`.

In `outlimer-classify`, metadata columns are joined into ranked output tables
and the HTML report.

Metadata query support was added for group-aware enrichment:

```bash
outlimer-classify \
  --union-csv top_union_summary.csv \
  --metadata samples.tsv \
  --foreground-query "status == 'case'" \
  --background-query "status == 'control'" \
  --out-dir classify
```

## Classification CLI Changes: `outlimer-classify`

`src/OutliMer/classification.py` was rewritten into a more robust downstream
classifier.

### Modes

Union mode:

```bash
outlimer-classify \
  --union-csv outlimer_results/export_kmers/top_union_summary.csv \
  --out-dir outlimer_results/classify
```

Report mode:

```bash
outlimer-classify \
  --mode report \
  --report-csv outlimer_results/kmer_report.csv \
  --out-dir outlimer_results/classify_report
```

`--mode report` no longer requires `--union-csv`.

### Scoring Fixes

The classifier now handles tiny or tied cohorts more safely. If all scores tie,
samples are not all labeled with an anomaly score of `1.0`.

Heavy dependencies are imported lazily so `--help` and report-mode paths are
lighter and easier to test.

### Enrichment Fixes

Hash enrichment now compares an explicit foreground to a disjoint background.

Example:

```bash
outlimer-classify \
  --union-csv outlimer_results/export_kmers/top_union_summary.csv \
  --foreground-samples sampleA,sampleB \
  --background-samples sampleC,sampleD \
  --out-dir outlimer_results/classify
```

If no foreground is supplied, the top anomalies are used as the foreground and
the remaining samples become background.

The Benjamini-Hochberg multiple-testing correction was fixed.

### HTML Report

`outlimer-classify` now writes:

```text
outlimer_report.html
```

by default. Disable it with:

```bash
outlimer-classify --no-html-report ...
```

The HTML report is static and dependency-free. It links key CSV/JSON artifacts,
embeds diagnostic plot thumbnails, and shows compact summary tables.

### Diagnostic Plots

Union-mode classification now writes a richer plot set by default:

- `pca_samples.png`: PCA of log-count features, colored by anomaly score and
  labeled only for high-priority samples.
- `jaccard_mds.png`: two-dimensional MDS from Jaccard distances.
- `anomaly_qc_scatter.png`: hash depth versus cohort-unique fraction, colored
  by anomaly score.
- `jaccard_cluster_heatmap.png`: clustered sample-by-sample Jaccard similarity
  heatmap.
- `driver_hash_heatmap.png`: log-count heatmap for hashes that distinguish
  top-ranked samples.
- `hash_enrichment_plot.png`: ranked `-log10(p-value)` view of enrichment
  results.

These plots are controlled by `--no-plots`, `--plot-labels`,
`--plot-top-hashes`, and `--plot-top-samples`.

### MultiQC, Sparse Matrix, and Run Comparison

`outlimer-classify --multiqc` writes MultiQC custom-content JSON for anomaly
scores.

`outlimer-classify --sparse-npz` writes a sparse sample-by-hash feature matrix
plus a `.labels.json` sidecar. This is intended for larger cohorts where wide
CSV matrices become cumbersome.

Run comparison mode compares a query union matrix against a baseline union
matrix:

```bash
outlimer-classify \
  --mode compare \
  --baseline-union old/export_kmers/top_union_summary.csv \
  --query-union new/export_kmers/top_union_summary.csv \
  --out-dir new/compare
```

The output ranks query samples by novelty relative to the baseline.

## Provenance

`outlimer` now writes:

```text
run_manifest.json
```

The manifest records:

- input type
- k-mer size
- scaled value
- processed samples
- metadata columns
- cache directory
- primary output paths

This is meant to make pipeline runs easier to audit.

Both `outlimer --ro-crate` and `outlimer-classify --ro-crate` can now emit a
minimal `ro-crate-metadata.json` for workflow-level provenance.

## Profile Utilities

A new command was added:

```bash
outlimer-profile
```

Supported subcommands:

- `describe profile_dir`
- `compare profile_a profile_b`
- `build profile_dir --hash-db db.txt.gz`
- `update profile_dir --hash-db new_hashes.txt.gz`

These operate on the existing profile layout:

- `db.txt.gz`
- `agg_counts.csv`
- `samples.txt`

## Workflow Examples

Minimal workflow templates were added:

- `workflows/snakemake/Snakefile`
- `workflows/snakemake/config.yaml`
- `workflows/nextflow/main.nf`

These run the two-step workflow:

1. `outlimer`
2. `outlimer-classify`

They are intentionally small and should be copied into analysis repositories
before adding cluster resources, containers, conda environments, or project
specific metadata rules.

## Toy Dataset

A tiny example dataset was added under `examples/toy/`. It contains short FASTA
sequences and a small `samples.tsv`. It is intended for smoke tests and docs,
not biological benchmarking.

## Release Automation

Two GitHub Actions workflows were added.

### CI

`.github/workflows/ci.yml` runs on pushes and pull requests:

- installs the package with development dependencies
- runs tests
- runs Ruff
- builds wheel and sdist

### PyPI Release

`.github/workflows/release.yml` builds distributions and publishes to PyPI when
a GitHub release is published.

It uses PyPI Trusted Publishing via OIDC, so the repository does not need a
long-lived PyPI API token. Before the first release, configure the PyPI project
to trust this GitHub repository and workflow.

## Bioconda Preparation

A draft recipe was added:

```text
packaging/bioconda/meta.yaml
```

After publishing a PyPI sdist, fill in:

- the final version
- the PyPI sdist SHA256

Then submit the recipe to `bioconda-recipes`.

The package is pure Python and the draft recipe uses `noarch: python`.

## Tests Added

Tests were added under `tests/`:

- `test_packaging.py`
- `test_outlimer_core.py`
- `test_classification.py`
- `test_extensions.py`
- `test_profiles.py`

Current test coverage includes:

- regular package discovery
- deterministic paired FASTQ discovery
- DB comparison edge cases
- union summary output paths
- argument validation
- tied anomaly scores
- report mode without union CSV
- enrichment foreground/background behavior
- metadata loading
- cache round trips
- sourmash signature input behavior with a fake sourmash loader
- metadata query selection
- run-comparison calculations
- MultiQC custom-content writers
- RO-Crate writers
- profile build/describe/compare helpers

## Verification Performed

The local verification performed so far:

```bash
python -m py_compile ...
python -m unittest discover -s tests -v
python -m pip wheel . --no-deps --no-build-isolation
```

The test suite passed with 19 tests, with one expected skip when SciPy was not
available in the bundled test runtime.

The wheel build passed and the wheel contents were inspected. The wheel contains
only:

- `OutliMer/OutliMer.py`
- `OutliMer/classification.py`
- `OutliMer/profiles.py`
- `OutliMer/__init__.py`
- package metadata
- license metadata

Generated outputs and workflow templates are not included in the wheel payload.

After approval for network access, runtime dependencies were installed into
`/private/tmp/outlimer_runtime_deps` and a real toy workflow was run with
actual `sourmash`, `scikit-learn`, `scipy`, `matplotlib`, `seaborn`, `pandas`,
and `numpy`:

```bash
outlimer \
  --input-dir examples/toy/fasta \
  --input-type fasta \
  --ksize 5 \
  --scaled 1 \
  --metadata examples/toy/samples.tsv \
  --out-dir /private/tmp/outlimer_toy_real \
  --multiqc \
  --ro-crate \
  --plot

outlimer-classify \
  --union-csv /private/tmp/outlimer_toy_real/export_kmers/top_union_summary.csv \
  --metadata examples/toy/samples.tsv \
  --foreground-query "status == 'case'" \
  --background-query "status == 'control'" \
  --out-dir /private/tmp/outlimer_toy_real/classify \
  --multiqc \
  --ro-crate \
  --sparse-npz toy_matrix.npz

outlimer-classify \
  --mode compare \
  --baseline-union /private/tmp/outlimer_toy_real/export_kmers/top_union_summary.csv \
  --query-union /private/tmp/outlimer_toy_real/export_kmers/top_union_summary.csv \
  --metadata examples/toy/samples.tsv \
  --out-dir /private/tmp/outlimer_toy_real/compare \
  --multiqc \
  --ro-crate
```

Those commands produced the expected reports, plots, MultiQC JSON files,
RO-Crate metadata, sparse matrix export, cache entries, and run-comparison CSV.

## Known Caveats

Real end-to-end integration now passes on the toy FASTA dataset locally, but it
still needs to run in CI or a clean release environment with:

```bash
python -m pip install -e ".[dev]"
```

Sourmash signature input is still unit-tested with a fake sourmash loader. A
real `.sig` integration fixture should be added before the first public release.

## Suggested Next Checks

Before publishing:

1. Create a fresh virtual environment with Python 3.10, 3.11, and 3.12.
2. Install with `python -m pip install -e ".[dev]"`.
3. Run `python -m pytest`.
4. Run a tiny real FASTQ end-to-end workflow.
5. Run a tiny real `.sig` end-to-end workflow.
6. Run `python -m build`.
7. Run `twine check dist/*`.
8. Publish a test release or use PyPI Trusted Publishing on a GitHub release.
9. Fill the bioconda SHA256 from the PyPI sdist.

## Common Commands

Run the main workflow:

```bash
outlimer \
  --input-dir reads/ \
  --input-type paired-fastq \
  --metadata samples.tsv \
  --out-dir outlimer_results
```

Classify from the union matrix:

```bash
outlimer-classify \
  --union-csv outlimer_results/export_kmers/top_union_summary.csv \
  --metadata samples.tsv \
  --out-dir outlimer_results/classify
```

Classify from the report:

```bash
outlimer-classify \
  --mode report \
  --report-csv outlimer_results/kmer_report.csv \
  --metadata samples.tsv \
  --out-dir outlimer_results/classify_report
```

Build a wheel locally:

```bash
python -m pip wheel . --no-deps --no-build-isolation -w dist
```
