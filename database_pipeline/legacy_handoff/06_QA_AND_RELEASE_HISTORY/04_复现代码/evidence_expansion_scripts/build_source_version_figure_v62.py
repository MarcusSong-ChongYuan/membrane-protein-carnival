from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figure_style_v62 import COLORS, DATA_DIR, PDF_DIR, PNG_DIR, SVG_DIR, configure_style


ROWS = [
    ("Protein", "UniProtKB", "2026_02\n10 Jun 2026",
     "Canonical sequence/accession; gene; names; TM/intramembrane/lipid-anchor features; localization; families; disease/PDB cross-references",
     "Primary/secondary accession → current primary UniProt"),
    ("Protein", "Human Protein Atlas", "HPA 25.1\nMay 2026",
     "Tissue and cell-type RNA; IHC; MS; subcellular localization; predicted membrane-proteome context",
     "Ensembl gene ↔ unique UniProt; ambiguous mappings reviewed"),
    ("Protein", "Human Transmembrane Proteome", "UniTmp HTP d.2.2\nsnapshot 24 Jul 2026",
     "Predicted TM status, segment count, topology orientation and reliability",
     "UniProt entry name → primary accession"),
    ("Protein", "Membranome / BioMembHub", "Current snapshot\n24 Jul 2026",
     "Membrane family, membrane system, membrane segment and OPM/PDB links",
     "Supplied UniProt accession"),
    ("Structure", "PDBe", "Weekly/current snapshot\nJul 2026",
     "PDB entries, chains, biological assemblies, ligands, SIFTS mappings, residue contacts and validation metadata",
     "PDB chain → SIFTS → UniProt; CCD ligand → structure identity"),
    ("Structure", "OPM", "Current snapshot\n24 Jul 2026",
     "Membrane placement/orientation, membrane thickness, topology, family and integral/peripheral structural classification",
     "PDB entry/chain → UniProt"),
    ("Structure", "PDBTM / UniTmp", "Current snapshot\n24 Jul 2026",
     "Transmembrane PDB chains, segment counts and chain topology",
     "PDB entry/chain → UniProt"),
    ("Structure", "PDBbind", "2020.R1\nreprocessed with 2024 workflow",
     "Processed protein–ligand complexes, affinity index, receptor, ligand and pocket coordinate files",
     "PDB chain/SIFTS → UniProt; ligand structure → full InChIKey"),
    ("Interaction", "ChEMBL", "ChEMBL 37\n1 May 2026",
     "Target, assay, molecule, activity type/relation/value, mechanism and compound development status",
     "Single-protein target/accession → UniProt; molecule structure → InChIKey"),
    ("Interaction", "BindingDB", "Monthly snapshot\n2026-07",
     "Measured protein–small-molecule affinities, including Kd, Ki and IC50, with target and publication provenance",
     "UniProt/target sequence → unique protein; ligand structure → InChIKey"),
    ("Interaction", "PubChem BioAssay", "Current snapshot\n2026-07",
     "AID assay metadata; SID/CID substances and compounds; active/inactive outcomes; concentration and assay measurements",
     "NCBI protein/gene target → unique UniProt; SID/CID → standardized structure"),
    ("Interaction", "BRENDA", "2026.1\nMar 2026",
     "Human enzymes, EC numbers, substrates, products, cofactors, inhibitors, activators and kinetic values",
     "Human + EC + unique UniProt; ligand structure or conservative name review"),
    ("Disease", "Open Targets Platform", "26.06\n24 Jun 2026",
     "Target–disease associations, evidence channels, scores, literature/genetic/clinical evidence and MONDO identifiers",
     "Target Ensembl/UniProt → canonical UniProt; disease → MONDO/ontology ID"),
    ("Disease", "UniProtKB disease", "2026_02\n10 Jun 2026",
     "Expert-curated disease comments and disease cross-references, including MIM/Orphanet where available",
     "Current primary UniProt accession; original disease identifier retained"),
]


def wrap(value: str, width: int) -> str:
    return "\n".join(
        textwrap.fill(part, width=width, break_long_words=False)
        for part in str(value).splitlines()
    )


def main() -> None:
    configure_style()
    frame = pd.DataFrame(
        ROWS,
        columns=["Module", "Database", "Version / snapshot", "Information contributed", "Identity mapping anchor"],
    )
    frame.to_csv(DATA_DIR / "M0_source_versions_contributions.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(16.5, 12.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[4.7, 1.35])
    ax = fig.add_subplot(grid[0])
    ax.axis("off")
    fig.suptitle(
        "M0  Data-source versions, information contributed and membrane-scope benchmark",
        x=0.015, y=1.01, ha="left", color=COLORS["navy"], fontsize=14, fontweight="bold",
    )
    ax.set_title(
        "A  Frozen source inventory used by MemPro V6.2",
        loc="left", pad=10, fontsize=10.5, fontweight="bold",
    )
    display = frame.copy()
    display["Module"] = display["Module"].map(lambda x: wrap(x, 11))
    display["Database"] = display["Database"].map(lambda x: wrap(x, 23))
    display["Version / snapshot"] = display["Version / snapshot"].map(lambda x: wrap(x, 21))
    display["Information contributed"] = display["Information contributed"].map(lambda x: wrap(x, 54))
    display["Identity mapping anchor"] = display["Identity mapping anchor"].map(lambda x: wrap(x, 42))
    table = ax.table(
        cellText=display.values,
        colLabels=display.columns,
        cellLoc="left",
        colLoc="left",
        colWidths=[0.075, 0.14, 0.13, 0.385, 0.27],
        bbox=[0, 0, 1, 0.965],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(6.6)
    module_colors = {
        "Protein": "#D9EAF7",
        "Structure": "#E2F0D9",
        "Interaction": "#FCE4D6",
        "Disease": "#E4DFEC",
    }
    for (row, col), cell in table.get_celld().items():
        cell.set_linewidth(0.45)
        cell.set_edgecolor("#C9CED6")
        if row == 0:
            cell.set_facecolor(COLORS["navy"])
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        else:
            module = frame.iloc[row - 1]["Module"]
            cell.set_facecolor(module_colors[module] if col == 0 else ("#F7F9FB" if row % 2 == 0 else "white"))
            if col == 0:
                cell.get_text().set_weight("bold")
        cell.PAD = 0.035

    ax2 = fig.add_subplot(grid[1])
    ax2.set_title(
        "B  Human membrane-scope comparison (definitions are not identical)",
        loc="left", pad=8, fontsize=10.5, fontweight="bold",
    )
    labels = [
        "HPA predicted TM genes",
        "UniProt reviewed human\nwith TM feature",
        "MemPro A all",
        "MemPro A default E1/E2",
        "UniProt non-TM\nintramembrane/lipid feature",
        "MemPro B default E1/E2",
        "UniProt explicit peripheral,\nnon-TM/intramembrane",
        "MemPro C default E1/E2",
        "MemPro C candidate E3",
    ]
    values = [5573, 5230, 5608, 5490, 629, 440, 1028, 1974, 2532]
    colors = [
        "#5B9BD5", "#2F5597", COLORS["A"], COLORS["A"],
        "#70AD47", COLORS["B"], "#B39DDB", COLORS["C"], COLORS["review"],
    ]
    y = np.arange(len(labels))[::-1]
    bars = ax2.barh(y, values, color=colors)
    ax2.set_yticks(y, labels)
    ax2.bar_label(bars, labels=[f"{value:,}" for value in values], padding=3, fontsize=7)
    ax2.set_xlabel("Genes or UniProt accessions (scope-specific; not additive)")
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.grid(axis="x", color="#D9D9D9", linewidth=0.45)
    ax2.text(
        0.995, 0.03,
        "MemPro default = website_default_v52 = 1. E3 candidates and E0 audit records are not included in default counts.",
        transform=ax2.transAxes, ha="right", va="bottom", fontsize=6.8, color="#555555",
    )

    stem = "M0_source_versions_contributions"
    fig.savefig(SVG_DIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(PDF_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(PNG_DIR / f"{stem}.png", bbox_inches="tight", dpi=600)
    plt.close(fig)


if __name__ == "__main__":
    main()
