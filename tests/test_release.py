from __future__ import annotations

import gzip
import hashlib
import tempfile
import unittest
from pathlib import Path

from mempro_repro.release import tsv_row_count, verify_manifest


class ReleaseUtilityTests(unittest.TestCase):
    def test_tsv_row_count_handles_embedded_newlines(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "records.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                handle.write("id\tname\n1\tplain\n2\t\"line one\nline two\"\n")
            self.assertEqual(tsv_row_count(path), 2)

    def test_manifest_verification_detects_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "05_metadata").mkdir()
            payload = root / "payload.txt"
            payload.write_text("MemPro", encoding="utf-8")
            digest = hashlib.sha256(payload.read_bytes()).hexdigest()
            (root / "05_metadata" / "FORMAL_MANIFEST_SHA256.tsv").write_text(
                f"relative_path\tbytes\tsha256\npayload.txt\t{payload.stat().st_size}\t{digest}\n",
                encoding="utf-8",
            )
            result = verify_manifest(root)
            self.assertEqual(result["status"], "PASS")
            payload.write_text("changed", encoding="utf-8")
            self.assertEqual(verify_manifest(root)["status"], "FAIL")


if __name__ == "__main__":
    unittest.main()

