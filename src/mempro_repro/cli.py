from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .paths import load_project_paths, release_preflight
from .release import build_release_summary, verify_manifest


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def command_doctor(args: argparse.Namespace) -> int:
    paths = load_project_paths(args.data_root)
    missing = release_preflight(paths.data_root)
    payload = {
        "repository_root": str(paths.repo_root),
        "data_root": str(paths.data_root),
        "output_root": str(paths.output_root),
        "website_root": str(paths.website_root),
        "release_preflight": "PASS" if not missing else "FAIL",
        "missing_release_paths": missing,
        "python": sys.version,
    }
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0 if not missing else 2


def command_verify(args: argparse.Namespace) -> int:
    paths = load_project_paths(args.data_root)
    missing = release_preflight(paths.data_root)
    if missing:
        print(json.dumps({"status": "FAIL", "missing_release_paths": missing}, indent=2))
        return 2
    payload = verify_manifest(paths.data_root, None if args.full else args.sample_size)
    _write_json(paths.output_root / "reproducibility" / "release_verification.json", payload)
    print(json.dumps(payload, indent=2))
    return 0 if payload["status"] == "PASS" else 2


def command_summary(args: argparse.Namespace) -> int:
    paths = load_project_paths(args.data_root)
    missing = release_preflight(paths.data_root)
    if missing:
        print(json.dumps({"status": "FAIL", "missing_release_paths": missing}, indent=2))
        return 2
    payload = build_release_summary(paths.data_root)
    _write_json(paths.output_root / "reproducibility" / "release_summary.json", payload)
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0 if payload["release_qa_status"] == "PASS" else 2


def main() -> int:
    parser = argparse.ArgumentParser(description="Portable MemPro FORMAL release tools")
    parser.add_argument("--data-root", help="Path to extracted 01_database_FORMAL directory")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("doctor", help="Check paths and release prerequisites").set_defaults(func=command_doctor)
    verify = subparsers.add_parser("verify", help="Verify formal file hashes")
    verify.add_argument("--full", action="store_true", help="Hash every manifest file; slower")
    verify.add_argument("--sample-size", type=int, default=32, help="Deterministic sample count when not using --full")
    verify.set_defaults(func=command_verify)
    subparsers.add_parser("summary", help="Recount formal tables and emit a release summary").set_defaults(func=command_summary)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

