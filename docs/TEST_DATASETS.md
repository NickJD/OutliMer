# OutliMer Test Dataset Guide

OutliMer uses layered test data so fast installation checks remain separate
from scientific validation and performance benchmarking.

## Committed Fixtures

Small fixtures belong under `tests/data/` and `examples/toy/`. They must be
redistributable, deterministic, and small enough to run on every pull request.

The classifier truth fixture contains eight composition controls, one
lower-depth version of the same composition, and one sample with an orthogonal
hash composition. Its expected result is that `composition_outlier` is ranked
first and flagged while `depth_control` is not flagged.

The toy FASTA cohort contains two related controls and one deliberately
different case. CI runs it through both `outlimer` and `outlimer-classify` with
`ksize=5`, `scaled=1`, and the real sourmash implementation.

Add future committed fixtures for these cases:

- valid plain and gzip-compressed FASTQ records;
- common paired-read naming conventions;
- orphan and duplicate mates;
- truncated FASTQ and empty FASTA records;
- abundance signatures at supported and incompatible scales;
- multiple-signature, wrong-seed, wrong-ksize, and non-DNA signature files;
- empty, tied, malformed, duplicate-hash, and negative-count union matrices.

## Generated Synthetic Cohort

Create a deterministic paired-read generator under `validation/` using a fixed
random seed. A useful first cohort contains 24 samples:

- 16 controls generated from one base genome or community with modest
  abundance variation;
- four contamination samples with 1%, 5%, 20%, and 50% unrelated reads;
- two depth controls containing 10% and 1% of normal read depth;
- two technical controls with different read lengths or error profiles.

Evaluate each contamination level in a separate otherwise-identical cohort.
Putting all gradient samples in one cohort lets them cluster with one another
and weakens the intended truth test.

Record the simulator name and version, random seed, exact command, source
sequence accession and checksum, read count, read length, error model, and
expected class. Useful acceptance checks are:

- identical seeds produce identical union matrices with one and four threads;
- gzip and plain-text representations produce identical sketches;
- composition scores increase with contamination fraction;
- 20% and 50% contamination rank at or near the top;
- downsampling alone does not create a composition anomaly;
- malformed or incomplete inputs produce a non-zero exit code;
- `--allow-partial` is required before outputs can omit a failed sample.

## Public Validation Data

Use unrestricted CAMI toy metagenomes for biologically realistic validation.
Do not assume that an original CAMI sample is an outlier. Instead, create a
controlled perturbation or sample-swap experiment and preserve the original
gold-standard metadata.

The existing 80-sample `shallow_MG` run is useful as a case study and runtime
benchmark, but it has no confirmed outlier labels and must not be used to claim
sensitivity or specificity.

Large raw reads should not be committed. Store a `validation/datasets.tsv`
manifest with these fields:

```text
dataset sample source_url accession sha256 license expected_class perturbation seed
```

Commit download and generation scripts, checksums, compact expected summaries,
and the exact OutliMer command. Archive a frozen validation release in a DOI-
backed repository when results are used in a manuscript.

## Scale Benchmarks

Generate sparse matrices at 10, 50, 100, and 500 samples with increasing hash
counts. Capture wall time, peak resident memory, output size, and plot time for
both CLIs. The dense union CSV and pairwise Jaccard calculations are expected
to become the limiting steps; report those limits explicitly rather than
presenting the package as unbounded.

## Test Data Review

Every fixture or external dataset must state its origin, redistribution terms,
truth definition, expected result, and whether it tests installation,
correctness, robustness, or performance. A dataset without an explicit truth
definition is a smoke test or case study, not a validation set.
