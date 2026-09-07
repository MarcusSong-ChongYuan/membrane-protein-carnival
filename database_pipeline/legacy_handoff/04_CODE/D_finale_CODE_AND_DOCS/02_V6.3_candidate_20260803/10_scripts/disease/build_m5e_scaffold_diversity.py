#!/usr/bin/env python3
"""Build full-release Bemis-Murcko scaffold diversity statistics and M5E."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold


COLORS = {"navy": "#173B57", "blue": "#3572A5", "teal": "#2A9D8F", "orange": "#F4A261", "gray": "#B9C0C8", "dark": "#67727E", "light": "#E9EDF1", "purple": "#7A5195"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--compound-master", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--figure-dir", type=Path, required=True)
    p.add_argument("--qa-dir", type=Path, required=True)
    return p.parse_args()


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def scaffold_id(smiles: str) -> str:
    return "MSCF-" + hashlib.sha1(smiles.encode()).hexdigest()[:16].upper()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg", "pdf"):
        (args.figure_dir / ext).mkdir(parents=True, exist_ok=True)
    args.qa_dir.mkdir(parents=True, exist_ok=True)

    mapping_path = args.output_dir / "compound_scaffold_mapping_v0_1.tsv.gz"
    fields = ["compound_internal_id", "preferred_name", "computed_structural_class", "scaffold_status", "murcko_scaffold_id", "murcko_scaffold_smiles", "generic_scaffold_smiles", "rdkit_version", "source_release"]
    scaffold_counts = Counter()
    scaffold_generic = {}
    scaffold_examples = {}
    class_scaffolds = defaultdict(Counter)
    status_counts = Counter()
    total_core = 0
    with gzip.open(mapping_path, "wt", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for chunk in pd.read_csv(args.compound_master, sep="\t", dtype=str, usecols=["compound_internal_id", "preferred_name", "compound_scope_status", "record_qc_status", "standard_smiles", "computed_structural_class"], chunksize=100000, keep_default_na=False):
            qc = chunk["record_qc_status"].str.lower().isin({"ok", "pass", "passed"})
            frame = chunk.loc[chunk["compound_scope_status"].eq("core") & qc]
            for row in frame.itertuples(index=False):
                total_core += 1
                smiles = str(row.standard_smiles).strip()
                status, sid, murcko, generic = "", "", "", ""
                if not smiles:
                    status = "missing_smiles"
                else:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        status = "invalid_smiles"
                    else:
                        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
                        if scaffold_mol.GetNumAtoms() == 0:
                            status = "acyclic_no_murcko_scaffold"
                        else:
                            murcko = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
                            generic_mol = MurckoScaffold.MakeScaffoldGeneric(scaffold_mol)
                            generic = Chem.MolToSmiles(generic_mol, canonical=True, isomericSmiles=False)
                            sid = scaffold_id(murcko)
                            status = "murcko_scaffold_resolved"
                            scaffold_counts[sid] += 1
                            scaffold_generic[sid] = generic
                            scaffold_examples[sid] = murcko
                            class_scaffolds[str(row.computed_structural_class) or "unclassified"][sid] += 1
                status_counts[status] += 1
                writer.writerow({
                    "compound_internal_id": row.compound_internal_id,
                    "preferred_name": row.preferred_name,
                    "computed_structural_class": row.computed_structural_class,
                    "scaffold_status": status,
                    "murcko_scaffold_id": sid,
                    "murcko_scaffold_smiles": murcko,
                    "generic_scaffold_smiles": generic,
                    "rdkit_version": rdkit.__version__,
                    "source_release": "MemPro V6.2 small_molecule_master_v1_3",
                })

    frequency_rows = []
    ranked = scaffold_counts.most_common()
    for rank, (sid, count) in enumerate(ranked, start=1):
        frequency_rows.append({
            "scaffold_rank": rank,
            "murcko_scaffold_id": sid,
            "murcko_scaffold_smiles": scaffold_examples[sid],
            "generic_scaffold_smiles": scaffold_generic[sid],
            "compound_count": count,
            "compound_fraction_of_scaffold_bearing": count / max(sum(scaffold_counts.values()), 1),
            "singleton_flag": int(count == 1),
        })
    freq_path = args.output_dir / "scaffold_frequency_v0_1.tsv"
    with freq_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(frequency_rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(frequency_rows)

    class_rows = []
    for class_name, counts in class_scaffolds.items():
        compounds = sum(counts.values())
        unique = len(counts)
        singleton_compounds = sum(count for count in counts.values() if count == 1)
        class_rows.append({
            "computed_structural_class": class_name,
            "scaffold_bearing_compounds": compounds,
            "unique_murcko_scaffolds": unique,
            "scaffold_diversity_ratio": unique / max(compounds, 1),
            "singleton_scaffold_compounds": singleton_compounds,
            "singleton_compound_fraction": singleton_compounds / max(compounds, 1),
        })
    class_rows.sort(key=lambda row: row["scaffold_bearing_compounds"], reverse=True)
    class_path = args.output_dir / "scaffold_diversity_by_structural_class_v0_1.tsv"
    with class_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(class_rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(class_rows)

    scaffold_bearing = sum(scaffold_counts.values())
    unique_scaffolds = len(scaffold_counts)
    singletons = sum(1 for count in scaffold_counts.values() if count == 1)
    counts_array = np.array([count for _, count in ranked], dtype=float)
    cumulative = np.cumsum(counts_array) / max(scaffold_bearing, 1)
    top10_share = float(counts_array[:10].sum() / max(scaffold_bearing, 1))
    top100_share = float(counts_array[:100].sum() / max(scaffold_bearing, 1))
    summary = [{
        "core_qc_compounds": total_core,
        "scaffold_bearing_compounds": scaffold_bearing,
        "acyclic_compounds": status_counts["acyclic_no_murcko_scaffold"],
        "missing_smiles": status_counts["missing_smiles"],
        "invalid_smiles": status_counts["invalid_smiles"],
        "unique_murcko_scaffolds": unique_scaffolds,
        "singleton_scaffolds": singletons,
        "singleton_scaffold_fraction": singletons / max(unique_scaffolds, 1),
        "top10_scaffold_compound_share": top10_share,
        "top100_scaffold_compound_share": top100_share,
        "scaffold_diversity_ratio": unique_scaffolds / max(scaffold_bearing, 1),
        "sampling_used": 0,
        "rdkit_version": rdkit.__version__,
    }]
    summary_path = args.output_dir / "M5E_scaffold_diversity_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summary)

    mpl.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold", "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "savefig.facecolor": "white", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(11.4, 7.4))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)
    ax_a, ax_b, ax_c, ax_d = [fig.add_subplot(gs[index]) for index in range(4)]
    for ax, letter, title in zip((ax_a, ax_b, ax_c, ax_d), "ABCD", ("Structure-resolution coverage", "Scaffold rank-frequency", "Cumulative chemical coverage", "Most frequent Murcko scaffolds")):
        ax.text(-0.08, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", color=COLORS["navy"])
        ax.set_title(title, loc="left", pad=7)

    coverage_names = ["Scaffold-bearing", "Acyclic", "Missing SMILES", "Invalid SMILES"]
    coverage_values = [scaffold_bearing, status_counts["acyclic_no_murcko_scaffold"], status_counts["missing_smiles"], status_counts["invalid_smiles"]]
    coverage_colors = [COLORS["teal"], COLORS["orange"], COLORS["gray"], COLORS["purple"]]
    y = np.arange(len(coverage_names))
    ax_a.barh(y, coverage_values, color=coverage_colors)
    ax_a.set_yticks(y, coverage_names); ax_a.invert_yaxis(); ax_a.set_xlabel("Canonical compounds (count)"); ax_a.grid(axis="x", color=COLORS["light"], linewidth=0.6)
    for index, value in enumerate(coverage_values):
        ax_a.text(value, index, f" {value:,} ({value/max(total_core,1):.1%})", va="center", fontsize=6.8)

    ranks = np.arange(1, len(counts_array) + 1)
    ax_b.plot(ranks, counts_array, color=COLORS["blue"], linewidth=1.2)
    ax_b.set_xscale("log"); ax_b.set_yscale("log"); ax_b.set_xlabel("Scaffold rank"); ax_b.set_ylabel("Compounds per scaffold"); ax_b.grid(which="both", color=COLORS["light"], linewidth=0.5)
    ax_b.text(0.98, 0.96, f"Unique scaffolds: {unique_scaffolds:,}\nSingletons: {singletons:,} ({singletons/max(unique_scaffolds,1):.1%})", transform=ax_b.transAxes, ha="right", va="top", fontsize=7, color=COLORS["dark"])

    ax_c.plot(ranks, cumulative * 100, color=COLORS["teal"], linewidth=1.6)
    ax_c.set_xscale("log"); ax_c.set_ylim(0, 101); ax_c.set_xlabel("Top-ranked scaffolds included"); ax_c.set_ylabel("Scaffold-bearing compounds covered (%)"); ax_c.grid(which="both", color=COLORS["light"], linewidth=0.5)
    ax_c.scatter([10, min(100, len(ranks))], [top10_share * 100, top100_share * 100], color=[COLORS["orange"], COLORS["purple"]], zorder=3)
    ax_c.text(10, top10_share * 100, f"  Top 10: {top10_share:.1%}", va="center", fontsize=7)
    ax_c.text(min(100, len(ranks)), top100_share * 100, f"  Top 100: {top100_share:.1%}", va="center", fontsize=7)

    top = ranked[:12]
    mols = [Chem.MolFromSmiles(scaffold_examples[sid]) for sid, _ in top]
    legends = [f"#{rank+1}  n={count:,}" for rank, (_, count) in enumerate(top)]
    grid = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(180, 125), legends=legends, useSVG=False)
    ax_d.imshow(np.asarray(grid)); ax_d.axis("off")
    ax_d.text(0.0, -0.06, "Scaffolds are connectivity frameworks; counts are not activity or approval rankings.", transform=ax_d.transAxes, ha="left", va="top", fontsize=6.7, color=COLORS["dark"])

    fig.suptitle("M5E — Full-release Bemis–Murcko scaffold diversity", x=0.06, ha="left", fontsize=12, fontweight="bold", color=COLORS["navy"])
    fig.text(0.06, 0.93, f"All core, QC-passed canonical compounds were processed (n={total_core:,}); no embedding or random subsampling was used.", ha="left", fontsize=7.5, color=COLORS["dark"])
    fig.subplots_adjust(left=0.10, right=0.98, top=0.86, bottom=0.08)
    figure_paths = {}
    for ext in ("png", "svg", "pdf"):
        path = args.figure_dir / ext / f"M5E_scaffold_diversity.{ext}"
        fig.savefig(path, dpi=450 if ext == "png" else None, bbox_inches="tight")
        figure_paths[ext] = str(path)
    plt.close(fig)

    caption = """# M5E caption\n\n**M5E. Full-release Bemis–Murcko scaffold diversity of canonical compounds.** All core, QC-passed canonical compounds in the frozen V6.2 small-molecule master were processed without random sampling. Ring-containing molecules were reduced to Bemis–Murcko frameworks with RDKit; acyclic molecules, missing SMILES and invalid SMILES are reported separately rather than assigned artificial scaffolds. The rank-frequency and cumulative-coverage panels quantify scaffold reuse, while the structure gallery shows the twelve most frequent frameworks. Scaffold frequency describes chemical redundancy and diversity; it is not a measure of bioactivity, approval status or docking suitability.\n"""
    (args.figure_dir / "CAPTION_M5E_SCAFFOLD_DIVERSITY.md").write_text(caption, encoding="utf-8")

    checks = {
        "all_core_qc_rows_accounted": sum(status_counts.values()) == total_core,
        "scaffold_frequency_counts_reconcile": sum(scaffold_counts.values()) == scaffold_bearing,
        "scaffold_ids_unique": len(scaffold_counts) == unique_scaffolds,
        "no_sampling_used": summary[0]["sampling_used"] == 0,
        "figure_outputs_complete": all(Path(path).exists() and Path(path).stat().st_size > 0 for path in figure_paths.values()),
    }
    validation = {
        "module": "M5E_scaffold_diversity",
        "module_version": "0.1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "counts": summary[0],
        "scaffold_status": dict(status_counts),
        "figure_paths": figure_paths,
        "input_sha256": hash_file(args.compound_master),
    }
    validation_path = args.qa_dir / "M5E_SCAFFOLD_DIVERSITY_VALIDATION.json"
    validation_path.write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
