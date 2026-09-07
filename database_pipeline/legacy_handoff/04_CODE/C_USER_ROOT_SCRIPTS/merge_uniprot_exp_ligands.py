from __future__ import annotations

import csv
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")
REPORT = BASE / "uniprot_exp_ligand_merge_report.json"
EXCLUDED = BASE / "uniprot_exp_ligand_excluded_terms.tsv"
CACHE = BASE / "uniprot_exp_ligand_pubchem_cache.json"

NOT_FOUND = "not_found_in_current_sources"

GENERIC_EXCLUDE_EXACT = {
    "",
    "substrate",
    "a substrate",
    "a divalent metal cation",
    "a metal cation",
    "a carbohydrate",
    "a lipid",
    "a peptide",
    "a peptide antigen",
    "a protein",
    "an alcohol",
    "a nucleoside 3',5'-cyclic phosphate",
    "a carbohydrate group",
    "a lipid group",
    "a metal ion",
}

GENERIC_EXCLUDE_PREFIXES = (
    "a 1,2-diacyl-",
    "an n-acyl-",
    "a 1-acyl-",
    "a divalent ",
    "a metal ",
    "a peptide ",
    "a carbohydrate",
    "a lipid",
    "a protein",
    "an alcohol",
)

ION_TERMS = {
    "Ca(2+)",
    "Zn(2+)",
    "Mg(2+)",
    "K(+)",
    "Na(+)",
    "Mn(2+)",
    "Cu(2+)",
    "H(+)",
    "chloride",
    "Fe cation",
    "Cu cation",
}

COFACTOR_TERMS = {
    "ATP",
    "ADP",
    "AMP",
    "GTP",
    "GDP",
    "CTP",
    "IMP",
    "NAD",
    "NAD(+)",
    "NADH",
    "NADP(+)",
    "FAD",
    "FMN",
    "CoA",
    "heme",
    "heme b",
    "[4Fe-4S] cluster",
    "[2Fe-2S] cluster",
    "pyridoxal 5'-phosphate",
    "S-adenosyl-L-methionine",
    "(6R)-L-erythro-5,6,7,8-tetrahydrobiopterin",
    "(6R)-5,10-methylene-5,6,7,8-tetrahydrofolate",
    "folate",
}

LIPID_TERMS_EXACT = {
    "cholesterol",
    "sphing-4-enine 1-phosphate",
    "1,2-dioctanoyl-sn-glycero-3-phospho-(1D-myo-inositol-4,5-bisphosphate)",
    "phylloquinone",
    "hexadecanoate",
    "octanoate",
    "3-hydroxyoctanoate",
}

PUBCHEM_ALIASES = {
    "ATP": ("5957", "5'-Atp"),
    "ADP": ("6022", "5'-Adp"),
    "AMP": ("6083", "Adenosine monophosphate"),
    "GTP": ("6830", "Guanosine triphosphate"),
    "GDP": ("8977", "Guanosine diphosphate"),
    "CTP": ("6176", "Cytidine triphosphate"),
    "IMP": ("8582", "Inosine monophosphate"),
    "NAD": ("5892", "NAD"),
    "NAD(+)": ("5893", "Nicotinamide adenine dinucleotide"),
    "NADH": ("439153", "NADH"),
    "NADP(+)": ("5886", "NADP"),
    "FAD": ("643975", "FAD"),
    "FMN": ("643976", "FMN"),
    "CoA": ("87642", "Coenzyme A"),
    "cholesterol": ("5997", "Cholesterol"),
    "serotonin": ("5202", "Serotonin"),
    "L-glutamate": ("33032", "L-glutamic acid"),
    "glycine": ("750", "Glycine"),
    "4-aminobutanoate": ("119", "Gamma-aminobutyric acid"),
    "D-serine": ("71077", "D-serine"),
    "L-arginine": ("6322", "L-arginine"),
    "L-tryptophan": ("6305", "L-tryptophan"),
    "adenosine": ("60961", "Adenosine"),
    "dopamine": ("681", "Dopamine"),
    "histamine": ("774", "Histamine"),
    "melatonin": ("896", "Melatonin"),
    "progesterone": ("5994", "Progesterone"),
    "aldosterone": ("5839", "Aldosterone"),
    "cortisol": ("5754", "Cortisol"),
    "17beta-estradiol": ("5757", "Estradiol"),
    "caffeine": ("2519", "Caffeine"),
    "phosphate": ("1061", "Phosphate"),
    "spermine": ("1103", "Spermine"),
    "phylloquinone": ("5280483", "Vitamin K1"),
    "D-glucose": ("5793", "D-glucose"),
    "D-galactose": ("6036", "D-galactose"),
    "deoxycholate": ("222528", "Deoxycholic acid"),
    "Ca(2+)": ("271", "Calcium ion"),
    "Zn(2+)": ("32051", "Zinc ion"),
    "Mg(2+)": ("888", "Magnesium ion"),
    "K(+)": ("813", "Potassium ion"),
    "Na(+)": ("923", "Sodium ion"),
    "Mn(2+)": ("27854", "Manganese(2+)"),
    "Cu(2+)": ("27099", "Copper(2+)"),
    "H(+)": ("1038", "Hydron"),
    "chloride": ("312", "Chloride"),
}

PRIORITY = [
    "structure_affinity_ligand",
    "approved_or_clinical_target_ligand",
    "cofactor_or_functional_ligand",
    "membrane_lipid_or_sterol",
    "uniprot_experimental_site_ligand",
    "curated_target_ligand",
    "structure_bound_ligand",
    "bioassay_active_ligand",
    "ion_or_metal",
    "buffer_salt_solvent",
    "unknown_structure_ligand",
]
RANK = {value: idx for idx, value in enumerate(PRIORITY)}


def safe_id(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").upper()[:80] or "UNKNOWN"


def is_generic(ligand: str) -> bool:
    s = ligand.strip()
    low = s.lower()
    if low in GENERIC_EXCLUDE_EXACT:
        return True
    if low.startswith(GENERIC_EXCLUDE_PREFIXES):
        return True
    if low in {"dna", "rna"}:
        return True
    return False


def role_for_ligand(ligand: str) -> str:
    if ligand in ION_TERMS:
        return "ion_or_metal"
    if ligand in COFACTOR_TERMS:
        return "endogenous_ligand_or_cofactor"
    if ligand in LIPID_TERMS_EXACT or "cholesterol" in ligand.lower() or "sphing" in ligand.lower():
        return "membrane_lipid_or_sterol"
    return "bioactive_research_ligand"


def interpretation_for_role(role: str) -> str:
    if role == "ion_or_metal":
        return "ion_or_metal"
    if role == "endogenous_ligand_or_cofactor":
        return "cofactor_or_functional_ligand"
    if role == "membrane_lipid_or_sterol":
        return "membrane_lipid_or_sterol"
    return "uniprot_experimental_site_ligand"


def parse_site_entry(entry: str) -> tuple[str, str, str] | None:
    if "@" not in entry:
        return None
    ligand, rest = entry.split("@", 1)
    ligand = ligand.strip()
    rest = rest.strip()
    m = re.match(r"(.+?)\(([^()]*)\)$", rest)
    if m:
        return ligand, m.group(1), m.group(2)
    return ligand, rest, ""


def load_existing_small_molecules() -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    path = BASE / "small_molecule_database.tsv"
    rows: dict[str, dict[str, str]] = {}
    name_to_id: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows[row["drug_id"]] = row
            name_to_id.setdefault(row["drug_name"].lower(), row["drug_id"])
    return rows, name_to_id


def load_cache() -> dict[str, dict[str, str]]:
    if CACHE.exists():
        return json.loads(CACHE.read_text(encoding="utf-8"))
    return {}


def save_cache(cache: dict[str, dict[str, str]]) -> None:
    CACHE.write_text(json.dumps(cache, ensure_ascii=False, indent=2), encoding="utf-8")


def pubchem_lookup(name: str, cache: dict[str, dict[str, str]]) -> dict[str, str] | None:
    if name in cache:
        return cache[name] or None
    encoded = urllib.parse.quote(name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/cids/JSON"
    try:
        with urllib.request.urlopen(url, timeout=12) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        cid = str(data["IdentifierList"]["CID"][0])
        cache[name] = {"cid": cid, "name": name, "method": "pubchem_name_lookup"}
        time.sleep(0.15)
        return cache[name]
    except Exception:
        cache[name] = {}
        time.sleep(0.15)
        return None


def resolve_ligand(ligand: str, small_rows: dict[str, dict[str, str]], name_to_id: dict[str, str], cache: dict[str, dict[str, str]]) -> dict[str, str]:
    if ligand in PUBCHEM_ALIASES:
        cid, canonical = PUBCHEM_ALIASES[ligand]
        return {"drug_id": cid, "compound_id_type": "PubChem CID", "drug_name": canonical, "mapping_status": "manual_pubchem_alias"}
    if ligand.lower() in name_to_id:
        drug_id = name_to_id[ligand.lower()]
        existing = small_rows[drug_id]
        return {
            "drug_id": drug_id,
            "compound_id_type": existing["compound_id_type"],
            "drug_name": existing["drug_name"],
            "mapping_status": "existing_small_molecule_name",
        }
    lookup = pubchem_lookup(ligand, cache)
    if lookup:
        return {"drug_id": lookup["cid"], "compound_id_type": "PubChem CID", "drug_name": ligand, "mapping_status": lookup["method"]}
    return {
        "drug_id": f"UNMAPPED_UNIPROT_LIGAND:{safe_id(ligand)}",
        "compound_id_type": "unmapped_uniprot_ligand",
        "drug_name": ligand,
        "mapping_status": "unmapped_uniprot_ligand",
    }


def build_uniprot_records() -> tuple[list[dict[str, str]], list[dict[str, str]], Counter]:
    small_rows, name_to_id = load_existing_small_molecules()
    cache = load_cache()
    records: dict[tuple[str, str], dict[str, str]] = {}
    excluded: list[dict[str, str]] = []
    map_status = Counter()

    with (BASE / "protein_database.tsv").open("r", encoding="utf-8", newline="") as f:
        for protein in csv.DictReader(f, delimiter="\t"):
            sites_by_lig: dict[str, list[tuple[str, str]]] = defaultdict(list)
            for entry in [x.strip() for x in protein.get("exp_binding_sites", "").split(";") if x.strip()]:
                parsed = parse_site_entry(entry)
                if not parsed:
                    continue
                ligand, site, quality = parsed
                if not ligand:
                    continue
                sites_by_lig[ligand].append((site, quality))

            for ligand in [x.strip() for x in protein.get("exp_ligands", "").split(";") if x.strip()]:
                if is_generic(ligand):
                    excluded.append(
                        {
                            "target_uniprot_id": protein["target_uniprot_id"],
                            "approved_symbol": protein["approved_symbol"],
                            "protein_name": protein["protein_name"],
                            "excluded_ligand": ligand,
                            "reason": "generic_or_non_specific_ligand_term",
                            "exp_binding_sites": ";".join(
                                f"{ligand}@{site}({quality})" if quality else f"{ligand}@{site}"
                                for site, quality in sites_by_lig.get(ligand, [])
                            ),
                            "exp_pubmed": protein.get("exp_pubmed", ""),
                        }
                    )
                    continue
                resolved = resolve_ligand(ligand, small_rows, name_to_id, cache)
                map_status[resolved["mapping_status"]] += 1
                role = role_for_ligand(ligand)
                interp = interpretation_for_role(role)
                site_entries = [
                    f"{ligand}@{site}({quality})" if quality else f"{ligand}@{site}"
                    for site, quality in sites_by_lig.get(ligand, [])
                ]
                evidence_quality = [
                    f"{site}:{quality}" for site, quality in sites_by_lig.get(ligand, []) if quality
                ]
                key = (protein["target_uniprot_id"], resolved["drug_id"])
                records[key] = {
                    "target_uniprot_id": protein["target_uniprot_id"],
                    "approved_symbol": protein["approved_symbol"],
                    "ncbi_gene_id": protein.get("ncbi_gene_id", ""),
                    "protein_name": protein["protein_name"],
                    "drug_id": resolved["drug_id"],
                    "compound_id_type": resolved["compound_id_type"],
                    "drug_name": resolved["drug_name"],
                    "compound_name_category": "common_or_drug_name" if resolved["compound_id_type"] == "PubChem CID" else "uniprot_ligand_name",
                    "compound_biological_role": role,
                    "ligand_interpretation_class": interp,
                    "compound_source": "uniprot_binding_site_annotation",
                    "source_database": "UniProt",
                    "assay_or_mechanism": "UniProt experimental binding site annotation",
                    "activity_type": "",
                    "activity_value_uM": "",
                    "clinical_or_approval_status": "",
                    "has_pdb_biolip_matched": "0",
                    "pdb_biolip_matched_sites": "",
                    "has_pdb_scpdb_matched": "0",
                    "pdb_scpdb_matched_sites": "",
                    "has_pdbbind_matched": "0",
                    "pdbbind_matched_sites": "",
                    "has_stitch_matched": "0",
                    "stitch_matched_compounds": "",
                    "has_exp_binding_site": "1",
                    "exp_binding_sites": ";".join(site_entries),
                    "exp_ligands": ligand,
                    "exp_evidence_quality": ";".join(evidence_quality),
                    "exp_pubmed": protein.get("exp_pubmed", ""),
                    "exp_pdb": protein.get("exp_pdb", ""),
                    "target_group_id": "",
                }
    save_cache(cache)
    return list(records.values()), excluded, map_status


def rewrite_binary(new_records: list[dict[str, str]]) -> None:
    path = BASE / "drug_protein_binary_relationships.tsv"
    temp = path.with_suffix(".tmp")
    exp_cols = {"exp_binding_sites", "exp_ligands", "exp_evidence_quality", "exp_pubmed", "exp_pdb"}
    with path.open("r", encoding="utf-8", newline="") as src, temp.open("w", encoding="utf-8", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        fields = reader.fieldnames or []
        writer = csv.DictWriter(dst, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            for col in exp_cols:
                row[col] = ""
            writer.writerow(row)
        for row in new_records:
            writer.writerow({field: row.get(field, "") for field in fields})
    temp.replace(path)


def aggregate_binary() -> dict[str, dict[str, object]]:
    agg: dict[str, dict[str, object]] = {}
    with (BASE / "drug_protein_binary_relationships.tsv").open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            did = row["drug_id"]
            item = agg.setdefault(
                did,
                {
                    "source_databases": set(),
                    "compound_source": set(),
                    "targets": set(),
                    "classes": set(),
                    "statuses": set(),
                    "activity_types": set(),
                    "activity_values": [],
                    "row_count": 0,
                    "name": row["drug_name"],
                    "id_type": row["compound_id_type"],
                    "name_category": row["compound_name_category"],
                    "role": row["compound_biological_role"],
                },
            )
            item["row_count"] = int(item["row_count"]) + 1
            item["targets"].add(row["target_uniprot_id"])
            item["classes"].add(row["ligand_interpretation_class"])
            item["compound_source"].update(x.strip() for x in row["compound_source"].split(";") if x.strip())
            item["source_databases"].update(x.strip() for x in row["source_database"].split(";") if x.strip())
            item["statuses"].update(x.strip() for x in row["clinical_or_approval_status"].split(";") if x.strip())
            if row["activity_type"]:
                item["activity_types"].add(row["activity_type"])
            if row["activity_value_uM"]:
                try:
                    item["activity_values"].append(float(row["activity_value_uM"]))
                except ValueError:
                    pass
    return agg


def sort_classes(values: set[str]) -> list[str]:
    return sorted(values, key=lambda value: RANK.get(value, len(RANK)))


def median(values: list[float]) -> str:
    if not values:
        return ""
    vals = sorted(values)
    n = len(vals)
    if n % 2:
        return f"{vals[n//2]:g}"
    return f"{(vals[n//2-1] + vals[n//2]) / 2:g}"


def evidence_status(statuses: set[str]) -> str:
    if "approved_or_listed" in statuses:
        return "approved_drug"
    if "clinical_trial" in statuses:
        return "clinical_trial_compound"
    return "binding_evidence_only"


def rewrite_small_molecules() -> tuple[int, int]:
    path = BASE / "small_molecule_database.tsv"
    temp = path.with_suffix(".tmp")
    agg = aggregate_binary()
    existing: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        for row in reader:
            existing[row["drug_id"]] = row

    added = 0
    for did, item in agg.items():
        if did not in existing:
            added += 1
            existing[did] = {
                field: "" for field in fields
            }
            existing[did].update(
                {
                    "drug_id": did,
                    "compound_id_type": str(item["id_type"]),
                    "drug_name": str(item["name"]),
                    "compound_name_category": str(item["name_category"]),
                    "compound_biological_role": str(item["role"]),
                    "ghs_hazard_classes": NOT_FOUND,
                    "ghs_signal_words": NOT_FOUND,
                    "toxicity_summary": NOT_FOUND,
                    "livertox": NOT_FOUND,
                    "drug_classes": NOT_FOUND,
                    "pharmacodynamics": NOT_FOUND,
                    "drug_indication": NOT_FOUND,
                }
            )

    for did, row in existing.items():
        if did not in agg:
            continue
        item = agg[did]
        classes = sort_classes(item["classes"])
        values = item["activity_values"]
        row["source_databases"] = ";".join(sorted(item["source_databases"]))
        row["compound_source"] = ";".join(sorted(item["compound_source"], key=lambda x: ["curated_target_relation", "bioassay_active_relation", "structure_ligand_context", "uniprot_binding_site_annotation"].index(x) if x in ["curated_target_relation", "bioassay_active_relation", "structure_ligand_context", "uniprot_binding_site_annotation"] else 99))
        row["source_row_count"] = str(item["row_count"])
        row["unique_target_count"] = str(len(item["targets"]))
        row["ligand_interpretation_classes"] = ";".join(classes)
        row["best_ligand_interpretation_class"] = classes[0] if classes else ""
        row["compound_evidence_status"] = evidence_status(item["statuses"])
        row["activity_types"] = ";".join(sorted(item["activity_types"]))
        row["activity_value_uM_min"] = f"{min(values):g}" if values else ""
        row["activity_value_uM_median"] = median(values)

    with temp.open("w", encoding="utf-8", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in existing.values():
            writer.writerow({field: row.get(field, "") for field in fields})
    temp.replace(path)
    return len(existing), added


def write_excluded(excluded: list[dict[str, str]]) -> None:
    fields = ["target_uniprot_id", "approved_symbol", "protein_name", "excluded_ligand", "reason", "exp_binding_sites", "exp_pubmed"]
    with EXCLUDED.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(excluded)


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def refresh_manifest() -> None:
    path = BASE / "MANIFEST.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    tracked = [Path(item["path"]).name for item in data.get("files", [])]
    for name in [
        "drug_protein_binary_relationships.tsv",
        "small_molecule_database.tsv",
        "uniprot_exp_ligand_excluded_terms.tsv",
        "uniprot_exp_ligand_merge_report.json",
    ]:
        if name not in tracked:
            tracked.append(name)
    items = []
    for name in tracked:
        p = BASE / name
        if p.exists():
            items.append({"path": f"upload/normalized_tables/{name}", "sha256": hash_file(p), "bytes": p.stat().st_size})
    data["files"] = items
    data["outputs"] = [item for item in items if Path(item["path"]).name in {
        "drug_protein_binary_relationships.tsv",
        "protein_database.tsv",
        "small_molecule_database.tsv",
        "pocket_instances.tsv",
        "normalized_tables_v3.xlsx",
    }]
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    csv.field_size_limit(100_000_000)
    new_records, excluded, map_status = build_uniprot_records()
    rewrite_binary(new_records)
    small_total, small_added = rewrite_small_molecules()
    write_excluded(excluded)
    report = {
        "new_uniprot_binary_rows": len(new_records),
        "excluded_generic_ligand_mentions": len(excluded),
        "small_molecule_rows_after": small_total,
        "new_small_molecule_entities_added": small_added,
        "mapping_status_counts": dict(map_status),
        "excluded_terms": dict(Counter(row["excluded_ligand"] for row in excluded)),
    }
    REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    refresh_manifest()
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
