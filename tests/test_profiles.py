import json
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from OutliMer import profiles
from OutliMer.OutliMer import save_hash_db


class ProfileTests(unittest.TestCase):
    def test_build_describe_and_compare_profiles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_a = os.path.join(tmpdir, "a.txt")
            db_b = os.path.join(tmpdir, "b.txt")
            save_hash_db(db_a, {1, 2, 3})
            save_hash_db(db_b, {2, 3, 4})
            prof_a = os.path.join(tmpdir, "prof_a")
            prof_b = os.path.join(tmpdir, "prof_b")

            self.assertEqual(profiles.main([
                "build", prof_a, "--hash-db", db_a, "--sample", "s1"
            ]), 0)
            self.assertEqual(profiles.main([
                "build", prof_b, "--hash-db", db_b, "--sample", "s2"
            ]), 0)

            described = profiles.describe_profile(prof_a)
            compared = profiles.compare_profiles(prof_a, prof_b)

        self.assertEqual(described["n_hashes"], 3)
        self.assertEqual(described["samples"], ["s1"])
        self.assertEqual(compared["n_shared_hashes"], 2)
        self.assertAlmostEqual(compared["jaccard"], 0.5)

    def test_describe_json_out(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db = os.path.join(tmpdir, "db.txt")
            out = os.path.join(tmpdir, "describe.json")
            profile_dir = os.path.join(tmpdir, "profile")
            save_hash_db(db, {10})
            profiles.main(["build", profile_dir, "--hash-db", db])

            rc = profiles.main(["describe", profile_dir, "--json-out", out])
            with open(out) as fh:
                payload = json.load(fh)

        self.assertEqual(rc, 0)
        self.assertEqual(payload["n_hashes"], 1)


if __name__ == "__main__":
    unittest.main()
