import csv
import contextlib
import io
import json
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from OutliMer import OutliMer as outlimer


class OutliMerCoreTests(unittest.TestCase):
    def test_gather_samples_from_dir_is_deterministic_and_pairs_fastq(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for name in [
                "sampleB_R2.fastq",
                "sampleA_R1.fastq",
                "sampleB_R1.fastq",
                "sampleA_R2.fastq",
                "ignored.txt",
            ]:
                with open(os.path.join(tmpdir, name), "w") as fh:
                    fh.write("@r\nACGT\n+\n!!!!\n")

            samples = outlimer.gather_samples_from_dir(tmpdir, "paired-fastq")

        self.assertEqual([row[0] for row in samples], ["sampleA", "sampleB"])
        self.assertTrue(samples[0][1].endswith("sampleA_R1.fastq"))
        self.assertTrue(samples[0][2].endswith("sampleA_R2.fastq"))

    def test_compare_sample_to_db_handles_empty_and_overlap(self):
        self.assertEqual(outlimer.compare_sample_to_db(set(), {1, 2}), (0, 0, 0.0))
        self.assertEqual(
            outlimer.compare_sample_to_db({1, 2, 3}, {2, 3, 4}),
            (2, 3, 2 / 3),
        )

    def test_paired_discovery_rejects_orphans_and_duplicate_mates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            orphan = os.path.join(tmpdir, "orphan_R1.fastq")
            with open(orphan, "w") as fh:
                fh.write("@r\nACGT\n+\n!!!!\n")
            with self.assertRaisesRegex(ValueError, "no matching mate"):
                outlimer.gather_samples_from_dir(tmpdir, "paired-fastq")

        with tempfile.TemporaryDirectory() as tmpdir:
            for name in ("sample_R1.fastq", "sample_R1.fq", "sample_R2.fastq"):
                with open(os.path.join(tmpdir, name), "w") as fh:
                    fh.write("@r\nACGT\n+\n!!!!\n")
            with self.assertRaisesRegex(ValueError, "duplicate R1"):
                outlimer.gather_samples_from_dir(tmpdir, "paired-fastq")

    def test_fastq_parser_rejects_truncated_and_length_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            truncated = os.path.join(tmpdir, "truncated.fastq")
            with open(truncated, "w") as fh:
                fh.write("@r\nACGT\n+\n")
            with self.assertRaisesRegex(ValueError, "truncated FASTQ"):
                list(outlimer.fastq_sequences(truncated))

            mismatch = os.path.join(tmpdir, "mismatch.fastq")
            with open(mismatch, "w") as fh:
                fh.write("@r\nACGT\n+\n!!!\n")
            with self.assertRaisesRegex(ValueError, "lengths differ"):
                list(outlimer.fastq_sequences(mismatch))

    def test_fasta_parser_requires_headers_and_sequences(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "bad.fa")
            with open(path, "w") as fh:
                fh.write("ACGT\n")
            with self.assertRaisesRegex(ValueError, "before first FASTA header"):
                list(outlimer.fasta_sequences(path))

    def test_write_union_summary_honors_explicit_output_path(self):
        counts = {
            "s1": {10: 2, 20: 1},
            "s2": {20: 3},
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "custom", "union.csv")
            written = outlimer.write_union_summary(
                counts,
                tmpdir,
                out_path=out_path,
            )
            with open(written, newline="") as fh:
                rows = list(csv.reader(fh))

        self.assertEqual(written, out_path)
        self.assertEqual(rows[0], ["hash", "s1", "s2"])
        self.assertEqual(rows[1], ["10", "2", "0"])
        self.assertEqual(rows[2], ["20", "1", "3"])

    def test_parser_rejects_bad_threads(self):
        parser = outlimer.build_arg_parser()
        args = parser.parse_args([
            "--input-dir", ".",
            "--input-type", "single-fastq",
            "--out-dir", "out",
            "--threads", "0",
        ])
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                outlimer._validate_args(parser, args)

    def test_partial_sample_failure_requires_explicit_opt_in(self):
        class FakeMinHash:
            def __init__(self, **kwargs):
                self.hashes = {}

            def add_sequence(self, sequence, force=True):
                self.hashes[len(sequence)] = self.hashes.get(len(sequence), 0) + 1

        old_sourmash = outlimer.sourmash
        old_minhash = outlimer.MinHash
        outlimer.sourmash = object()
        outlimer.MinHash = FakeMinHash
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                inputs = os.path.join(tmpdir, "inputs")
                os.makedirs(inputs)
                with open(os.path.join(inputs, "good.fa"), "w") as fh:
                    fh.write(">good\nACGTACGT\n")
                with open(os.path.join(inputs, "bad.fa"), "w") as fh:
                    fh.write("ACGT\n")

                base_args = [
                    "--input-dir", inputs,
                    "--input-type", "fasta",
                    "--out-dir", os.path.join(tmpdir, "out"),
                    "--no-cache",
                ]
                self.assertEqual(outlimer.main(base_args), 1)
                self.assertEqual(outlimer.main(base_args + ["--allow-partial"]), 0)
                with open(os.path.join(tmpdir, "out", "run_manifest.json")) as fh:
                    manifest = json.load(fh)
                self.assertEqual(manifest["schema_version"], 2)
                self.assertEqual(manifest["seed"], 42)
                self.assertIn("OutliMer", manifest["software_versions"])
                self.assertIn("bad", manifest["failed_samples"])
        finally:
            outlimer.sourmash = old_sourmash
            outlimer.MinHash = old_minhash


if __name__ == "__main__":
    unittest.main()
