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
from OutliMer import classification


class _FakeMinHash:
    ksize = 31
    scaled = 10000
    seed = 42
    moltype = "DNA"
    track_abundance = True
    hashes = {101: 2, 202: 1}

    def downsample(self, *, scaled):
        self.scaled = scaled
        return self


class _FakeSignature:
    minhash = _FakeMinHash()


class _FakeSignatureWith:
    def __init__(self, minhash):
        self.minhash = minhash


class _FakeSourmash:
    @staticmethod
    def load_file_as_signatures(path, ksize=None, select_moltype=None):
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
                cache_dir, "s1", (source, ""), 31, 10000, 42, "fasta", counts
            )
            loaded = outlimer._load_cached_counts(
                cache_dir, "s1", (source, ""), 31, 10000, 42, "fasta"
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
                    42, "signature"
                )
        finally:
            outlimer.sourmash = old_sourmash

        self.assertEqual(name, "sample")
        self.assertIsNone(error)
        self.assertEqual(counts, {101: 2, 202: 1})

    def test_signature_rejects_sparse_or_multiple_sketches(self):
        class SparseMinHash(_FakeMinHash):
            scaled = 20000

        class SparseSourmash:
            @staticmethod
            def load_file_as_signatures(path, ksize=None, select_moltype=None):
                return [_FakeSignatureWith(SparseMinHash())]

        class MultipleSourmash:
            @staticmethod
            def load_file_as_signatures(path, ksize=None, select_moltype=None):
                return [_FakeSignature(), _FakeSignature()]

        old_sourmash = outlimer.sourmash
        try:
            outlimer.sourmash = SparseSourmash()
            with self.assertRaisesRegex(RuntimeError, "cannot be upsampled"):
                outlimer._load_signature_counts("sparse.sig", 31, 10000, 42)
            outlimer.sourmash = MultipleSourmash()
            with self.assertRaisesRegex(RuntimeError, "one DNA signature per file"):
                outlimer._load_signature_counts("multi.sig", 31, 10000, 42)
        finally:
            outlimer.sourmash = old_sourmash

    def test_signature_validates_seed_molecule_abundance_and_downsampling(self):
        class SignatureSource:
            signatures = []

            @classmethod
            def load_file_as_signatures(
                cls, path, ksize=None, select_moltype=None
            ):
                return cls.signatures

        old_sourmash = outlimer.sourmash
        outlimer.sourmash = SignatureSource
        try:
            wrong_seed = _FakeMinHash()
            wrong_seed.seed = 99
            SignatureSource.signatures = [_FakeSignatureWith(wrong_seed)]
            with self.assertRaisesRegex(RuntimeError, "seed=99"):
                outlimer._load_signature_counts("seed.sig", 31, 10000, 42)

            protein = _FakeMinHash()
            protein.moltype = "protein"
            SignatureSource.signatures = [_FakeSignatureWith(protein)]
            with self.assertRaisesRegex(RuntimeError, "No DNA signature"):
                outlimer._load_signature_counts("protein.sig", 31, 10000, 42)

            flat = _FakeMinHash()
            flat.track_abundance = False
            SignatureSource.signatures = [_FakeSignatureWith(flat)]
            with self.assertRaisesRegex(RuntimeError, "abundance"):
                outlimer._load_signature_counts("flat.sig", 31, 10000, 42)

            dense = _FakeMinHash()
            dense.scaled = 1000
            SignatureSource.signatures = [_FakeSignatureWith(dense)]
            counts = outlimer._load_signature_counts(
                "dense.sig", 31, 10000, 42
            )
            self.assertEqual(dense.scaled, 10000)
            self.assertEqual(counts, {101: 2, 202: 1})
        finally:
            outlimer.sourmash = old_sourmash

    def test_metadata_loader_rejects_duplicate_samples(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "samples.csv")
            with open(path, "w") as fh:
                fh.write("sample,group\ns1,case\ns1,control\n")
            with self.assertRaisesRegex(ValueError, "duplicate metadata sample"):
                outlimer.load_sample_metadata(path)

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

    @unittest.skipUnless(
        outlimer.sourmash is not None and outlimer.MinHash is not None,
        "working sourmash installation is required",
    )
    def test_real_sourmash_toy_workflow_end_to_end(self):
        toy_dir = os.path.join(ROOT, "examples", "toy")
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = os.path.join(tmpdir, "run")
            rc = outlimer.main([
                "--input-dir", os.path.join(toy_dir, "fasta"),
                "--input-type", "fasta",
                "--ksize", "5",
                "--scaled", "1",
                "--metadata", os.path.join(toy_dir, "samples.tsv"),
                "--out-dir", run_dir,
                "--no-cache",
            ])
            self.assertEqual(rc, 0)

            classify_dir = os.path.join(run_dir, "classify")
            rc = classification.main([
                "--union-csv",
                os.path.join(run_dir, "export_kmers", "top_union_summary.csv"),
                "--out-dir", classify_dir,
                "--no-plots",
            ])
            self.assertEqual(rc, 0)

            with open(os.path.join(classify_dir, "top_anomalies.csv")) as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["sample"], "toy_case_1")


if __name__ == "__main__":
    unittest.main()
