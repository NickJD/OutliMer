import csv
import contextlib
import io
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


if __name__ == "__main__":
    unittest.main()
