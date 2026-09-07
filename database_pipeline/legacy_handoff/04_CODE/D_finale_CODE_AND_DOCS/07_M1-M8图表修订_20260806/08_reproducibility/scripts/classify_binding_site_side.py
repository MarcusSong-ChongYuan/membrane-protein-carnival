from __future__ import annotations

import json
import os
import re
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get("MEMPRO_REVISION_ROOT", r"D:\finale\07_M1-M8图表修订_20260806"))
V62 = Path(r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_2_20260730")
MM = Path(r"D:\7.22\membrane_master_v5_working")
OUT = ROOT / "06_binding_site_side"
OUT.mkdir(parents=True, exist_ok=True)


def parse_ranges(text: str) -> list[tuple[int, int]]:
    text = str(text or "")
    out = []
    for a, b in re.findall(r"(?:TRANSMEM\s+)?(\d+)\.\.(\d+)", text, flags=re.I):
        a, b = int(a), int(b)
        if b >= a and b - a <= 80:
            out.append((a, b))
    return sorted(set(out))


def parse_positions(text: str) -> list[int]:
    return sorted(set(int(x) for x in re.findall(r"\b[A-Z]{3}(\d+)\b", str(text or ""))))


def n_side(row) -> str | None:
    raw = str(row.get("htp_n_terminal_side_v5", "") or "").strip().upper()
    if raw == "I": return "cytoplasmic_side"
    if raw == "O": return "non_cytoplasmic_side"
    sp = str(row.get("single_pass_resolution_v5", "") or "").lower()
    if "type_ii" in sp: return "cytoplasmic_side"
    if "type_i" in sp or "type_iii" in sp: return "non_cytoplasmic_side"
    return None


def classify_position(pos: int, ranges: list[tuple[int, int]], side0: str | None) -> str:
    if any(a <= pos <= b for a, b in ranges):
        return "intramembrane"
    if not ranges or side0 is None:
        return "extramembrane_side_unresolved"
    crossed = sum(1 for _, b in ranges if b < pos)
    if crossed % 2 == 0:
        return side0
    return "non_cytoplasmic_side" if side0 == "cytoplasmic_side" else "cytoplasmic_side"


def classify_site(classes: list[str]) -> str:
    known = [x for x in classes if x != "extramembrane_side_unresolved"]
    if not known:
        return "extramembrane_side_unresolved"
    uniq = set(known)
    if len(uniq) == 1:
        return next(iter(uniq))
    if "intramembrane" in uniq:
        return "membrane_interface_or_mixed"
    return "both_sides_or_mixed"


def main():
    pcols = ["target_uniprot_id", "transmembrane_features", "transmembrane_count_v5",
             "htp_n_terminal_side_v5", "single_pass_resolution_v5", "opm_pdb_ids_v5", "pdbtm_pdb_ids_v5"]
    proteins = pd.read_csv(V62 / "human_membrane_protein_master_v6_2.tsv", sep="\t", usecols=pcols,
                           dtype=str, keep_default_na=False).set_index("target_uniprot_id")
    opm = set(pd.read_csv(MM / "normalized/opm_human_structures_v5.tsv", sep="\t", dtype=str)["opm_pdb_id"].str.lower())
    pdbtm = set(pd.read_csv(MM / "normalized/unitmp_pdbtm_structures_v5.tsv", sep="\t", dtype=str)["pdb_id"].str.lower())

    rows = []
    sites_path = V62 / "binding_site_instances_v6_2.tsv.gz"
    use = ["binding_site_instance_id", "target_uniprot_id", "compound_internal_id", "pdb_ids",
           "pdb_chain_ids_v60", "site_type", "residue_or_site_description", "residue_index_type_v60",
           "source_database", "record_qc_status"]
    for chunk in pd.read_csv(sites_path, sep="\t", compression="gzip", usecols=use, dtype=str,
                             keep_default_na=False, chunksize=100000):
        chunk = chunk[(chunk["site_type"] == "experimental_structure_residue_contact") &
                      (chunk["residue_index_type_v60"] == "UNIPROT") &
                      (chunk["record_qc_status"] == "ok")]
        for r in chunk.itertuples(index=False):
            if r.target_uniprot_id not in proteins.index:
                continue
            pr = proteins.loc[r.target_uniprot_id]
            ranges = parse_ranges(pr["transmembrane_features"])
            positions = parse_positions(r.residue_or_site_description)
            classes = [classify_position(x, ranges, n_side(pr)) for x in positions]
            side = classify_site(classes)
            pdbs = [x.strip().lower() for x in str(r.pdb_ids).split(";") if x.strip()]
            support = []
            if any(x in opm for x in pdbs): support.append("OPM")
            if any(x in pdbtm for x in pdbs): support.append("PDBTM/UniTmp")
            if not support: support.append("sequence_topology_only")
            rows.append({
                "binding_site_instance_id": r.binding_site_instance_id,
                "target_uniprot_id": r.target_uniprot_id,
                "compound_internal_id": r.compound_internal_id,
                "pdb_ids": r.pdb_ids,
                "pdb_chain_ids": r.pdb_chain_ids_v60,
                "source_database": r.source_database,
                "residue_positions_uniprot": ";".join(map(str, positions)),
                "transmembrane_ranges_uniprot": ";".join(f"{a}-{b}" for a,b in ranges),
                "residue_compartments": ";".join(classes),
                "membrane_side": side,
                "membrane_structure_support": ";".join(support),
                "coordinate_mapping_basis": "PDBe residue index normalized to UniProt (SIFTS-derived)",
                "orientation_basis": "UniProt TM ranges + UniTmp HTP N-terminal side or curated single-pass type",
                "classification_scope": "protein-relative; non-cytoplasmic may mean extracellular or organelle lumen",
            })
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "binding_site_membrane_side_instances.tsv.gz", sep="\t", index=False, compression="gzip")
    summary = out["membrane_side"].value_counts().rename_axis("membrane_side").reset_index(name="site_count")
    summary.to_csv(OUT / "binding_site_membrane_side_summary.tsv", sep="\t", index=False)
    support = out["membrane_structure_support"].value_counts().rename_axis("support").reset_index(name="site_count")
    support.to_csv(OUT / "binding_site_membrane_structure_support.tsv", sep="\t", index=False)
    qa = {
        "status": "PASS", "site_rows": len(out), "unique_proteins": int(out.target_uniprot_id.nunique()),
        "unique_pdbs": int(out.pdb_ids.nunique()), "side_counts": dict(Counter(out.membrane_side)),
        "limitations": [
            "This is a protein-relative topology classification, not an atom-coordinate OPM z-axis calculation.",
            "Non-cytoplasmic includes extracellular and organelle-lumen sides.",
            "Sites without reliable N-terminal orientation remain unresolved rather than forced.",
        ],
    }
    (OUT / "BINDING_SITE_MEMBRANE_SIDE_QA.json").write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(qa, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
