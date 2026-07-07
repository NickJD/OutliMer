# OutliMer

OutliMer detects outlier genomic or metagenomic samples from sourmash MinHash
k-mer sketches. It scans FASTQ or FASTA inputs, builds per-sample hash
abundance profiles, compares samples against either a supplied hash database or
the rest of the cohort, and writes reports for unusual samples and their
contributing hashes.

The project currently provides two command line tools:

- `outlimer`: sketch input sequence files and write hash/count outputs.
- `outlimer-classify`: run downstream anomaly scoring from an OutliMer union
  matrix or per-sample report.

## Installation

From a checkout:

```bash
python -m pip install .
```

For development:

```bash
python -m pip install -e ".[dev]"
```

OutliMer requires Python 3.10 or later.

## Quick Start

Sketch paired FASTQ files from a directory and write the default report and
union summary:

```bash
outlimer \
  --input-dir reads/ \
  --input-type paired-fastq \
  --out-dir outlimer_results
```

Sketch existing sourmash signatures instead of reads:

```bash
outlimer \
  --input-dir signatures/ \
  --input-type signature \
  --out-dir outlimer_results
```

Include sample metadata in reports:

```bash
outlimer \
  --input-dir reads/ \
  --input-type paired-fastq \
  --metadata samples.tsv \
  --metadata-sample-column sample \
  --out-dir outlimer_results
```

Run downstream anomaly classification from the generated union summary:

```bash
outlimer-classify \
  --union-csv outlimer_results/export_kmers/top_union_summary.csv \
  --out-dir outlimer_results/classify
```

Run classification from the per-sample report instead:

```bash
outlimer-classify \
  --mode report \
  --report-csv outlimer_results/kmer_report.csv \
  --out-dir outlimer_results/classify_report
```

Compare a new run against an existing union matrix:

```bash
outlimer-classify \
  --mode compare \
  --baseline-union old/export_kmers/top_union_summary.csv \
  --query-union new/export_kmers/top_union_summary.csv \
  --out-dir new/compare
```

## Main Outputs

`outlimer` writes:

- `kmer_report.csv`: per-sample hash totals, hashes shared with the comparison
  database or cohort, and hashes new to the sample.
- `export_kmers/top_union_summary.csv`: a wide hash-by-sample count matrix.
- `run_manifest.json`: run parameters, sample list, cache path, and primary
  output paths.
- Optional plot PNGs when `--plot` is supplied.
- Optional per-sample unique/shared hash exports when `--full-exports` and
  `--export-unique-shared` are supplied.
- Optional MultiQC custom-content JSON when `--multiqc` is supplied.
- Optional RO-Crate metadata when `--ro-crate` is supplied.

`outlimer-classify` writes:

- `outliers_report.csv`: combined anomaly scores from Jaccard distance,
  Isolation Forest, and Local Outlier Factor.
- `top_anomalies.csv`: the top 5 percent of samples by anomaly score.
- `explanations_by_sample.csv` and `.json`: top contributing hashes and short
  explanations.
- `hash_enrichment.csv`: Fisher exact enrichment of hashes in a foreground set
  against a background set.
- `outlimer_report.html`: a compact static report linking the main artifacts.
- Diagnostic PNGs unless `--no-plots` is supplied:
  `pca_samples.png`, `jaccard_mds.png`, `anomaly_qc_scatter.png`,
  `dendrogram.png`, `jaccard_cluster_heatmap.png`,
  `driver_hash_heatmap.png`, and `hash_enrichment_plot.png`.
- Plot controls: `--plot-labels`, `--plot-top-hashes`, and
  `--plot-top-samples`.
- Optional MultiQC custom-content JSON with `--multiqc`.
- Optional sparse `.npz` feature matrix with `--sparse-npz`.
- Optional RO-Crate metadata with `--ro-crate`.

## Caching

By default, sketches are cached under `<out-dir>/sketch_cache`. The cache key
includes input file path, size, mtime, k-mer size, scaled value, and input type.
Use `--no-cache` to force recomputation or `--cache-dir` to share a cache across
runs.

## Workflow Examples

Minimal Snakemake and Nextflow examples live under `workflows/`. They are
templates intended to be copied into analysis projects and adapted for local
resources, conda/container setup, and metadata conventions.

## Profile Utilities

`outlimer-profile` provides small profile-management helpers:

```bash
outlimer-profile describe profile_dir
outlimer-profile compare profile_a profile_b
outlimer-profile build profile_dir --hash-db db.txt.gz --sample run_1
outlimer-profile update profile_dir --hash-db new_hashes.txt.gz --sample run_2
```


For PyPI releases:

```bash
python -m build
twine check dist/*
twine upload dist/*
```

For bioconda, use the draft recipe under `packaging/bioconda/meta.yaml` as a
starting point after publishing the PyPI sdist. The recipe should be submitted
to the `bioconda-recipes` repository, where the final source URL and SHA256 can
be filled in from the PyPI release artifact.

## Development Checks

```bash
python -m pytest
python -m ruff check .
```

The tests use small synthetic data and do not require real sequencing files.

For a maintainer-oriented summary of the packaging, CLI, testing, workflow, and
release changes made during the cleanup pass, see `docs/CHANGE_GUIDE.md`.

## License

OutliMer is released under the GNU General Public License v3.0. See
`LICENSE` for details.
