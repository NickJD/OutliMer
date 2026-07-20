# Toy OutliMer Dataset

This tiny dataset is intended for smoke tests and documentation examples. It is
not biologically meaningful.

Run with:

```bash
outlimer \
  --input-dir examples/toy/fasta \
  --input-type fasta \
  --ksize 5 \
  --scaled 1 \
  --metadata examples/toy/samples.tsv \
  --out-dir outlimer_toy \
  --multiqc \
  --ro-crate

outlimer-classify \
  --union-csv outlimer_toy/export_kmers/top_union_summary.csv \
  --metadata examples/toy/samples.tsv \
  --foreground-query "status == 'case'" \
  --background-query "status == 'control'" \
  --out-dir outlimer_toy/classify \
  --multiqc \
  --ro-crate
```

The sequences are deliberately short so this should be run with a small k-mer
size such as `--ksize 5` if real sourmash integration is being checked.
