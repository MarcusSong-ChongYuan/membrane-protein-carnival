"""Shared utilities for the MemPro protein-language-model clustering workflow."""
from __future__ import annotations
import gzip, json, os, platform, sys
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent

def load_config() -> dict[str, Any]:
    with (ROOT / "config.yaml").open(encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    cfg["output_dir"] = ROOT / cfg["output_dir"]
    cfg["output_dir"].mkdir(exist_ok=True)
    (cfg["output_dir"] / "figures").mkdir(exist_ok=True)
    return cfg

def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    return pd.read_csv(path, sep="\t", compression="gzip" if path.suffix == ".gz" else None,
                       dtype=str, keep_default_na=False)

def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)

def save_json(value: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(value, indent=2, default=str), encoding="utf-8")

def get_unique_proteins(cfg: dict[str, Any]) -> pd.DataFrame:
    master = read_table(cfg["input"]["protein_master"])
    required = {"canonical_uniprot_accession", "canonical_sequence", "sequence_length"}
    missing = required - set(master.columns)
    if missing:
        raise ValueError(f"protein master missing columns: {sorted(missing)}")
    master = master.rename(columns={"canonical_uniprot_accession": "uniprot_id",
                                    "canonical_sequence": "sequence"})
    master["sequence_length"] = pd.to_numeric(master["sequence_length"], errors="coerce")
    master = master.loc[(master["uniprot_id"] != "") & (master["sequence"] != "")].copy()
    if master["uniprot_id"].duplicated().any():
        raise ValueError("canonical protein master contains duplicated UniProt accessions")
    role_path = Path(cfg["input"]["formal_role_overlay"])
    if role_path.exists():
        roles = read_table(role_path).rename(columns={"canonical_uniprot_accession": "uniprot_id",
                                                       "formal_primary_membrane_role": "target_class"})
        roles = roles[["uniprot_id", "target_class"]].drop_duplicates("uniprot_id")
        master = master.merge(roles, on="uniprot_id", how="left", validate="one_to_one")
    else:
        master["target_class"] = ""
    master["target_class"] = master["target_class"].replace("", np.nan).fillna("Unknown")
    cols = ["uniprot_id", "sequence", "sequence_length", "target_class"]
    return master[cols].sort_values("uniprot_id").reset_index(drop=True)

def software_versions() -> str:
    lines = [f"python\t{sys.version.replace(os.linesep, ' ')}", f"platform\t{platform.platform()}"]
    for module in ("numpy", "pandas", "sklearn", "umap", "hdbscan", "torch", "esm", "matplotlib", "seaborn", "plotly"):
        try:
            x = __import__(module)
            lines.append(f"{module}\t{getattr(x, '__version__', 'installed')}")
        except Exception as exc:
            lines.append(f"{module}\tNOT_AVAILABLE ({type(exc).__name__})")
    return "\n".join(lines) + "\n"
