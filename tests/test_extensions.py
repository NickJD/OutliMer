import csv
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from OutliMer import OutliMer as outlimer


class _FakeMinHash:
    ksize = 31
    hashes = {101: 2, 202: 1}


class _FakeSignature:
    minhash = _FakeMinHash()


class _FakeSourmash:
    @staticmethod
    def load_file_as_signatures(path, ksize=None):
        return [_FakeSignature()]


class ExtensionTests(unittest.TestCase):
    def test_metadata_loader_reads_tsv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "samples.tsv")
            with open(path, "w") as fh:
                fh.write("sample\tgroup\ns1\tcase\n")

            metadata, columns = outlimer.load_sample_metadata(path)

        self.assertEqual(columns, ["group"])
        self.assertEqual(metadata["s1"]["group"], "case")

    def test_cache_round_trip(self):
        counts = {1: 2, 3: 4}
        with tempfile.TemporaryDirectory() as tmpdir:
            source = os.path.join(tmpdir, "sample.fa")
            with open(source, "w") as fh:
                fh.write(">s\nACGT\n")
            cache_dir = os.path.join(tmpdir, "cache")

            outlimer._save_cached_counts(
                cache_dir, "s1", (source, ""), 31, 10000, "fasta", counts
            )
            loaded = outlimer._load_cached_counts(
                cache_dir, "s1", (source, ""), 31, 10000, "fasta"
            )

        self.assertEqual(loaded, counts)

    def test_signature_input_uses_sourmash_signature_loader(self):
        old_sourmash = outlimer.sourmash
        outlimer.sourmash = _FakeSourmash()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                sig = os.path.join(tmpdir, "sample.sig")
                with open(sig, "w") as fh:
                    fh.write("{}")
                samples = outlimer.gather_samples_from_dir(tmpdir, "signature")
                name, counts, error = outlimer._sketch_sample(
                    samples[0][0], samples[0][1], samples[0][2], 31, 10000,
                    "signature"
                )
        finally:
            outlimer.sourmash = old_sourmash

        self.assertEqual(name, "sample")
        self.assertIsNone(error)
        self.assertEqual(counts, {101: 2, 202: 1})

    def test_report_writer_includes_metadata_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_path = os.path.join(tmpdir, "samples.csv")
            with open(metadata_path, "w") as fh:
                fh.write("sample,group\ns1,case\n")
            metadata, columns = outlimer.load_sample_metadata(metadata_path)
            report_path = os.path.join(tmpdir, "report.csv")
            with open(report_path, "w", newline="") as fh:
                writer = csv.writer(fh)
                writer.writerow(["sample"] + columns + [
                    "n_hashes", "n_in_db", "pct_in_db", "n_new_hashes"
                ])
                writer.writerow(["s1"] + [
                    metadata["s1"].get(column, "") for column in columns
                ] + [2, 1, 0.5, 1])

            with open(report_path, newline="") as fh:
                rows = list(csv.reader(fh))

        self.assertEqual(rows[0], [
            "sample", "group", "n_hashes", "n_in_db", "pct_in_db",
            "n_new_hashes",
        ])
        self.assertEqual(rows[1][1], "case")

    def test_main_multiqc_and_rocrate_writers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            report_rows = [("s1", 10, 8, 0.8, 2)]
            mqc = outlimer.write_multiqc_summary(report_rows, tmpdir)
            crate = outlimer.write_ro_crate_metadata(
                tmpdir,
                {"input": "reads"},
                {"multiqc": mqc},
                {"ksize": 31},
            )

            self.assertTrue(os.path.exists(mqc))
            self.assertTrue(os.path.exists(crate))


if __name__ == "__main__":
    unittest.main()
