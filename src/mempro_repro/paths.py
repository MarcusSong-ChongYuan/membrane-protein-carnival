from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass
from pathlib import Path


REQUIRED_RELEASE_PATHS = (
    "01_core_v72/01_release_tables/protein_master_v72.tsv.gz",
    "01_core_v72/01_release_tables/protein_compound_pair_v72.tsv.gz",
    "01_core_v72/01_release_tables/positive_interaction_evidence_v72.tsv.gz",
    "02_companion_v721/02_companion_tables",
    "03_publication_repairs",
    "04_archive_negative_evidence",
    "05_metadata/FORMAL_MANIFEST_SHA256.tsv",
    "07_QA/FORMAL_VALIDATION_REPORT.json",
)


@dataclass(frozen=True)
class ProjectPaths:
    repo_root: Path
    data_root: Path
    output_root: Path
    website_root: Path


def repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _read_config(repo_root: Path) -> dict:
    default_path = repo_root / "config" / "release.example.toml"
    local_path = repo_root / "config" / "release.local.toml"
    path = local_path if local_path.exists() else default_path
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _absolute(root: Path, raw_path: str) -> Path:
    candidate = Path(raw_path).expanduser()
    return candidate if candidate.is_absolute() else (root / candidate).resolve()


def load_project_paths(data_root_override: str | None = None) -> ProjectPaths:
    repo_root = repository_root()
    config = _read_config(repo_root)
    configured = data_root_override or os.environ.get("MEMPRO_DATA_ROOT") or config["paths"]["data_root"]
    data_root = _absolute(repo_root, configured)
    output_root = _absolute(repo_root, os.environ.get("MEMPRO_OUTPUT_ROOT", config["paths"]["output_root"]))
    website_root = _absolute(repo_root, os.environ.get("MEMPRO_WEBSITE_ROOT", config["paths"]["website_root"]))
    return ProjectPaths(repo_root=repo_root, data_root=data_root, output_root=output_root, website_root=website_root)


def release_preflight(data_root: Path) -> list[str]:
    return [entry for entry in REQUIRED_RELEASE_PATHS if not (data_root / entry).exists()]

