import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813")
SITE = ROOT / "02_frozen_review" / "binding_site_coordinate_frozen_v7.tsv.gz"

chain_pattern = re.compile(r"@([0-9A-Za-z]{4}):([A-Za-z0-9]+)(?:=|\b)")
residue_pattern = re.compile(r"\b([A-Z]{1,3})(-?\d+)([A-Za-z]?)\b")
pdb_pattern = re.compile(r"(?i)\b[0-9][A-Za-z0-9]{3}\b")

total = 0
by_reason = Counter()
by_source = Counter()
by_type = Counter()
recover = Counter()
cross = Counter()
examples = defaultdict(list)

with gzip.open(SITE, "rt", encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        total += 1
        reason = row["coordinate_mapping_reason_v7"]
        source = row["source_database"] or "<EMPTY>"
        typ = row["site_type"] or "<EMPTY>"
        desc = row["residue_or_site_description"] or ""
        pdbs = [x.strip().lower() for x in re.split(r"[;,|]", row["pdb_ids"] or "") if re.fullmatch(r"(?i)[0-9][a-z0-9]{3}", x.strip())]
        chain_hits = chain_pattern.findall(desc)
        residues = residue_pattern.findall(desc)
        by_reason[reason] += 1; by_source[source] += 1; by_type[typ] += 1
        if reason == "CURRENT_PDB_CHAIN_MISSING":
            if chain_hits and residues:
                cls = "R1_EXPLICIT_PDB_CHAIN_AND_RESIDUES_IN_DESCRIPTION"
            elif chain_hits:
                cls = "R2_EXPLICIT_PDB_CHAIN_ONLY_IN_DESCRIPTION"
            elif len(pdbs) == 1 and residues:
                cls = "R3_SINGLE_PDB_RESIDUES_CHAIN_FROM_SIFTS_ENTITY"
            elif len(pdbs) > 1 and residues:
                cls = "R4_MULTI_PDB_RESIDUES_NEED_STRUCTURE_RANKING"
            elif pdbs:
                cls = "R5_PDB_ONLY_NO_RESIDUES_USE_COCRYSTAL_LIGAND_POCKET"
            else:
                cls = "R6_NO_USABLE_STRUCTURE_CONTEXT"
        else:
            requested = int(row.get("requested_residue_count_v7") or 0)
            mapped = int(row.get("mapped_residue_count_v7") or 0)
            if requested and not mapped:
                cls = "R7_CHAIN_PRESENT_ZERO_SIFTS_MATCH_SEQUENCE_ALIGN"
            else:
                cls = "R8_CHAIN_PRESENT_AMBIGUOUS_OR_PARTIAL"
        recover[cls] += 1
        cross[(cls, source)] += 1
        if len(examples[cls]) < 5:
            examples[cls].append({k: row[k] for k in ["binding_site_instance_id","evidence_id","target_uniprot_id","source_database","site_type","pdb_ids","residue_or_site_description","pdb_chain_ids_v60","requested_uniprot_residue_tokens_v7","coordinate_mapping_reason_v7"]})

report = {
    "total": total,
    "by_reason": dict(by_reason),
    "by_source": dict(by_source.most_common()),
    "by_site_type": dict(by_type.most_common()),
    "recoverability_classes": dict(recover),
    "recoverability_by_source": {f"{k[0]} | {k[1]}": v for k, v in sorted(cross.items(), key=lambda x: (-x[1], x[0]))},
    "examples": dict(examples),
}
out = Path(r"D:\finale\12_V7_final_completion_20260813\07_qa\FROZEN_SITE_RECOVERABILITY.json")
out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({k: report[k] for k in ["total","by_reason","by_source","by_site_type","recoverability_classes"]}, ensure_ascii=False, indent=2))
print(out)
