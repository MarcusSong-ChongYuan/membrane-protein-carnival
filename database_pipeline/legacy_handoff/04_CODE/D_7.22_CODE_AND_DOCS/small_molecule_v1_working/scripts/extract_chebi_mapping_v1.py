from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\7.22")
WORK = ROOT / "small_molecule_v1_working"
CHUNKS = WORK / "intermediate" / "standardized_chunks"
RAW = WORK / "raw" / "chebi"
REPORTS = WORK / "reports"
OUT = RAW / "chebi_202607_selected_mapping.tsv"
REPORT = REPORTS / "CHEBI_202607_MAPPING_REPORT.json"

SDF = RAW / "chebi_202607.sdf.gz"
OBO = RAW / "chebi_lite_202607.obo"
PROPERTY_PATTERN = re.compile(r"^>\s*<([^>]+)>")
SELECTED_PROPERTIES = {
    "ChEBI ID",
    "ChEBI NAME",
    "STAR",
    "INCHI",
    "INCHIKEY",
    "IUPAC_NAME",
    "SMILES",
    "FORMULA",
    "MASS",
    "MONOISOTOPIC_MASS",
    "SYNONYM",
    "SECONDARY_ID",
    "PubChem Compound Database Links",
    "ChEMBL Database Links",
    "DrugCentral Database Links",
    "HMDB Database Links",
    "CAS Registry Numbers",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_wanted_inchikeys() -> set[str]:
    wanted = set()
    for path in sorted(CHUNKS.glob("standardized_chunk_*.tsv")):
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                for field in (
                    "exact_standard_inchikey",
                    "parent_standard_inchikey",
                ):
                    value = row.get(field, "").strip()
                    if value:
                        wanted.add(value)
    return wanted


def read_selected_sdf_records(wanted: set[str]) -> list[dict]:
    matches = []
    props: dict[str, list[str]] = defaultdict(list)
    current_property = ""
    with gzip.open(SDF, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line == "$$$$":
                inchikey = " ".join(props.get("INCHIKEY", [])).strip()
                if inchikey in wanted:
                    record = {}
                    for key in SELECTED_PROPERTIES:
                        values = []
                        seen = set()
                        for value in props.get(key, []):
                            value = value.strip()
                            if value and value.casefold() not in seen:
                                seen.add(value.casefold())
                                values.append(value)
                        record[key] = "|".join(values)
                    matches.append(record)
                props = defaultdict(list)
                current_property = ""
                continue
            header = PROPERTY_PATTERN.match(line)
            if header:
                candidate = header.group(1)
                current_property = candidate if candidate in SELECTED_PROPERTIES else ""
                continue
            if current_property:
                if line == "":
                    current_property = ""
                else:
                    props[current_property].append(line)
    return matches


def parse_obo() -> dict[str, dict]:
    terms: dict[str, dict] = {}
    current: dict[str, object] | None = None

    def commit() -> None:
        nonlocal current
        if current and current.get("id"):
            terms[str(current["id"])] = current
        current = None

    with OBO.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line == "[Term]":
                commit()
                current = {"parents": [], "roles": []}
                continue
            if line.startswith("[") and line != "[Term]":
                commit()
                continue
            if current is None:
                continue
            if line.startswith("id: "):
                current["id"] = line[4:].strip()
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("is_a: "):
                current["parents"].append(line[6:].split(" ! ", 1)[0].strip())
            elif line.startswith("relationship: has_role "):
                current["roles"].append(
                    line[len("relationship: has_role ") :].split(" ! ", 1)[0].strip()
                )
            elif line.startswith("relationship: RO:0000087 "):
                current["roles"].append(
                    line[len("relationship: RO:0000087 ") :]
                    .split(" ! ", 1)[0]
                    .strip()
                )
        commit()
    return terms


def main() -> None:
    wanted = load_wanted_inchikeys()
    matches = read_selected_sdf_records(wanted)
    terms = parse_obo()

    fields = [
        "chebi_id",
        "chebi_name",
        "chebi_star",
        "inchikey",
        "inchi",
        "smiles",
        "formula",
        "mass",
        "monoisotopic_mass",
        "iupac_name",
        "synonyms",
        "secondary_ids",
        "pubchem_cids",
        "chembl_ids",
        "drugcentral_ids",
        "hmdb_ids",
        "cas_numbers",
        "direct_parent_ids",
        "direct_parent_names",
        "direct_role_ids",
        "direct_role_names",
    ]
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for record in sorted(matches, key=lambda item: item.get("ChEBI ID", "")):
            chebi_id = record.get("ChEBI ID", "")
            term = terms.get(chebi_id, {})
            parents = list(term.get("parents", []))
            roles = list(term.get("roles", []))
            writer.writerow(
                {
                    "chebi_id": chebi_id,
                    "chebi_name": record.get("ChEBI NAME", ""),
                    "chebi_star": record.get("STAR", ""),
                    "inchikey": record.get("INCHIKEY", ""),
                    "inchi": record.get("INCHI", ""),
                    "smiles": record.get("SMILES", ""),
                    "formula": record.get("FORMULA", ""),
                    "mass": record.get("MASS", ""),
                    "monoisotopic_mass": record.get("MONOISOTOPIC_MASS", ""),
                    "iupac_name": record.get("IUPAC_NAME", ""),
                    "synonyms": record.get("SYNONYM", ""),
                    "secondary_ids": record.get("SECONDARY_ID", ""),
                    "pubchem_cids": record.get(
                        "PubChem Compound Database Links", ""
                    ),
                    "chembl_ids": record.get("ChEMBL Database Links", ""),
                    "drugcentral_ids": record.get("DrugCentral Database Links", ""),
                    "hmdb_ids": record.get("HMDB Database Links", ""),
                    "cas_numbers": record.get("CAS Registry Numbers", ""),
                    "direct_parent_ids": "|".join(parents),
                    "direct_parent_names": "|".join(
                        terms.get(parent, {}).get("name", parent) for parent in parents
                    ),
                    "direct_role_ids": "|".join(roles),
                    "direct_role_names": "|".join(
                        terms.get(role, {}).get("name", role) for role in roles
                    ),
                }
            )

    matched_keys = {record.get("INCHIKEY", "") for record in matches}
    report = {
        "wanted_structure_inchikeys": len(wanted),
        "matched_chebi_records": len(matches),
        "matched_distinct_inchikeys": len(matched_keys),
        "match_rate": len(matched_keys) / len(wanted) if wanted else 0,
        "ontology_term_count": len(terms),
        "sdf_file": str(SDF),
        "sdf_sha256": sha256(SDF),
        "ontology_file": str(OBO),
        "ontology_sha256": sha256(OBO),
        "source_version": "ChEBI July 2026",
        "source_url": "https://www.ebi.ac.uk/chebi/downloads/",
        "output_file": str(OUT),
    }
    REPORT.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
