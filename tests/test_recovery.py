from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from mempro_repro.recovery import recover_legacy_context, verify_legacy_context


class LegacyRecoveryTests(unittest.TestCase):
    def test_recovery_preserves_provenance_and_does_not_claim_v72_equivalence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repo = root / "repo"
            plan = repo / "config" / "sources" / "legacy_context_recovery_v1.tsv"
            plan.parent.mkdir(parents=True)
            plan.write_text(
                "snapshot_id\tsource_id\tsource_database\tsource_relpath\ttarget_relpath\tartifact_type\trelease_equivalence\tnote\n"
                "test\ttest_source\tTest DB\tdata/input.txt\ttest/input.txt\tSOURCE_EXPORT\tNOT_EQUIVALENT_TO_V72\ttest\n",
                encoding="utf-8",
            )
            source = root / "legacy"
            (source / "data").mkdir(parents=True)
            (source / "data" / "input.txt").write_text("legacy", encoding="utf-8")
            destination = root / "recovered"
            result = recover_legacy_context(repo, source, destination, plan)
            self.assertEqual(result["status"], "PASS")
            self.assertEqual(result["recovered"], 1)
            manifest = (destination / "LEGACY_CONTEXT_MANIFEST.tsv").read_text(encoding="utf-8")
            self.assertIn("NOT_EQUIVALENT_TO_V72", manifest)
            self.assertTrue((destination / "test" / "input.txt").exists())
            self.assertEqual(verify_legacy_context(destination)["status"], "PASS")
            (destination / "test" / "input.txt").write_text("changed", encoding="utf-8")
            self.assertEqual(verify_legacy_context(destination)["status"], "FAIL")


if __name__ == "__main__":
    unittest.main()
