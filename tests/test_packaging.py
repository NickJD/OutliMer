import os
import sys
import contextlib
import io
import tomllib
import unittest

ROOT = os.path.dirname(os.path.dirname(__file__))
SRC = os.path.join(ROOT, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)


class PackagingTests(unittest.TestCase):
    def test_outlimer_is_regular_package(self):
        import OutliMer

        self.assertTrue(hasattr(OutliMer, "__version__"))
        self.assertTrue(os.path.exists(os.path.join(SRC, "OutliMer", "__init__.py")))

    def test_version_is_single_sourced(self):
        with open(os.path.join(ROOT, "pyproject.toml"), "rb") as fh:
            project = tomllib.load(fh)

        self.assertNotIn("version", project["project"])
        self.assertIn("version", project["project"]["dynamic"])
        self.assertEqual(
            project["tool"]["setuptools"]["dynamic"]["version"]["attr"],
            "OutliMer.__version__",
        )

    def test_all_commands_report_version(self):
        import OutliMer
        from OutliMer import OutliMer as core
        from OutliMer import classification, profiles

        for command in (core.main, classification.main, profiles.main):
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                with self.assertRaises(SystemExit) as raised:
                    command(["--version"])
            self.assertEqual(raised.exception.code, 0)
            self.assertIn(OutliMer.__version__, stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
