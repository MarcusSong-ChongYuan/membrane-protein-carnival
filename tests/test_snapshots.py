from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from mempro_repro.snapshots import audit_snapshots


class SnapshotAuditTests(unittest.TestCase):
    def _write_tsv(self, path: Path, header: str, rows: list[str]) -> None:
        path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")

    def test_missing_historical_inputs_are_limits_not_integrity_failures(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repo = root / "repo"
            config = repo / "config" / "sources"
            config.mkdir(parents=True)
            data_root = root / "data"
            raw = data_root / "06_scripts" / "legacy_handoff" / "raw"
            raw.mkdir(parents=True)
            snapshot = raw / "retained.txt"
            snapshot.write_text("frozen input", encoding="utf-8")
            digest = hashlib.sha256(snapshot.read_bytes()).hexdigest()

            registry = config / "registry.tsv"
            self._write_tsv(
                registry,
                "source_id\tmodule\tsource_database\tsource_version\trelease_or_retrieval_date\trecord_scope\thistorical_snapshot_status\tlegacy_parser_or_workflow",
                [
                    "retained\tprotein\tRetained DB\tv1\t2026-01\tinput\tSNAPSHOT_PRESENT\tparser.py",
                    "missing\tprotein\tMissing DB\tv1\t2026-01\tinput\tHISTORICAL_SNAPSHOT_NOT_IN_CURRENT_PACKAGE\tNOT_RETAINED",
                ],
            )
            manifest = config / "manifest.tsv"
            self._write_tsv(
                manifest,
                "snapshot_id\tsource_id\tsnapshot_relpath\tbytes\tsha256\tacquisition_mode\tstatus\tnote",
                [f"retained_v1\tretained\traw/retained.txt\t{snapshot.stat().st_size}\t{digest}\tfrozen\tPRESENT\ttest"],
            )

            audit = audit_snapshots(repo, data_root, registry_path=registry, manifest_path=manifest)
            self.assertEqual(audit["status"], "PASS_WITH_LIMITATIONS")
            self.assertEqual(audit["snapshot_counts"]["verified"], 1)
            self.assertEqual(audit["source_counts"]["historical_snapshot_unavailable"], 1)

            snapshot.write_text("changed", encoding="utf-8")
            self.assertEqual(audit_snapshots(repo, data_root, registry_path=registry, manifest_path=manifest)["status"], "FAIL")


if __name__ == "__main__":
    unittest.main()
