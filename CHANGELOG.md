# Changelog

All notable changes to OutliMer are documented here. The project follows
Semantic Versioning.

## [0.1.0] - Unreleased

### Added

- Installable `outlimer`, `outlimer-classify`, and `outlimer-profile` commands.
- FASTQ, FASTA, and sourmash signature input with deterministic discovery.
- Sketch caching, metadata attachment, profiles, MultiQC output, RO-Crate
  metadata, sparse matrices, run comparison, HTML reports, and diagnostic
  plots.
- Strict sequence, pairing, metadata, signature, and union-matrix validation.
- Cohort-relative anomaly ranking with depth-normalised count features.
- Deterministic classifier truth data and a real-sourmash toy integration test.

### Changed

- Failed samples stop a run unless `--allow-partial` is supplied.
- `--contamination` now controls the number of flagged anomaly candidates.
- Effectively constant detector components no longer dominate score ranking.
- Run manifests now record software versions, commands, input fingerprints,
  hash seed, and sample failures.

### Fixed

- Duplicate and orphaned paired reads are no longer silently selected or lost.
- Signature scale, molecule type, seed, abundance, and multiplicity are checked.
- Mean Jaccard distance excludes self-distance.
- Invalid matrix values are no longer silently converted to zero.
- Explanations and driver hashes are restricted to features used for scoring.
