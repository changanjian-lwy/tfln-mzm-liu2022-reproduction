import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_evidence", ROOT/"scripts/check_evidence.py")
check_evidence = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check_evidence)


class TestEvidenceManifest(unittest.TestCase):
    def test_registered_evidence_is_intact(self):
        self.assertEqual(check_evidence.static_check(), [])

    def test_frozen_experiments_are_all_registered(self):
        registered = {e["id"] for e in check_evidence.load_manifest()["experiments"]}
        on_disk = {p.name for p in (ROOT/"experiments/track_A_reproduction").iterdir() if p.is_dir()}
        self.assertEqual(on_disk - registered, set())

    def test_changed_or_missing_files_are_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            folder = root/"exp"
            folder.mkdir()
            for name in ("BOUNDARY.md", "RESULTS.md"):
                (folder/name).write_text(name)
            (folder/"out.json").write_text("{}")
            manifest = {"experiments": [{"id": "X", "dir": "exp",
                                         "outputs": {"exp/out.json": check_evidence.sha256(folder/"out.json"),
                                                     "exp/gone.csv": "0"*64}}]}
            (folder/"out.json").write_text('{"tampered": true}')
            failures = check_evidence.static_check(root, manifest)
        self.assertIn("X: missing INPUTS.md", failures)
        self.assertIn("X: outputs changed: exp/out.json", failures)
        self.assertIn("X: outputs file missing: exp/gone.csv", failures)


if __name__ == "__main__":
    unittest.main()
