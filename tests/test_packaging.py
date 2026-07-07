import os
import sys
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


if __name__ == "__main__":
    unittest.main()
