import importlib
from pathlib import Path
import unittest


class test_all_python_plots(unittest.TestCase):
    def test_run(self):
        directory = Path(__file__).parent

        for path in directory.glob("*.py"):
            # Don't try to import the test file itself
            if path == Path(__file__):
                continue

            module_name = path.stem

            try:
                print(f"Importing {path}")
                importlib.import_module(module_name)
            except Exception as e:
                self.fail(
                    f"Exception importing {path}: {e}"
                )


if __name__ == "__main__":
    unittest.main()

    
