# Contributing To OutliMer

Use Python 3.11 through 3.13 and install the development dependencies:

```bash
python -m pip install -e ".[dev]"
```

Before opening a change, run:

```bash
python -m pytest
python -m ruff check .
python -m build
python -m twine check dist/*
```

Changes to sample discovery, sketch compatibility, scoring, output schemas, or
CLI behaviour should include a regression test. Keep pull-request fixtures
small and deterministic. Place large validation reads outside Git and provide a
download or generation script, checksums, provenance, redistribution terms,
and an explicit expected result. See `docs/TEST_DATASETS.md`.

Report bugs with the OutliMer version, Python and sourmash versions, full
command, relevant log output, input naming pattern, and the smallest
redistributable reproducer possible. Do not attach sensitive sequencing data.
