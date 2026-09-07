import csv
import hashlib
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import openpyxl
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

BASE = Path(r"C:\github-repos\upload\normalized_tables")
PROTEIN = BASE / "protein_database.tsv"
SMALL = BASE / "small_molecule_database.tsv"
OUT_TSV = BASE / "protein_structure_ligand_context_review.tsv"
OUT_XLSX = BASE / "protein_structure_ligand_context_review.xlsx"
CACHE = BASE / "xlsx_build" / "pubchem_ligand_lookup_cache.json"
SUMMARY = BASE / "protein_structure_ligand_context_review_summary.json"
MANIFEST = BASE / "MANIFEST.json"

csv.field_size_limit(100_000_000)

SOURCE_COLUMNS = {
    "pdb_biolip_other": "BioLiP",
    "pdb_scpdb_other": "sc-PDB",
    "pdbbind_other": "PDBbind",
}

WATER = {"HOH", "WAT", "DOD"}
IONS = {
    "NA", "K", "CL", "CA", "MG", "MN", "ZN", "FE", "CU", "CO", "NI", "CD", "HG",
    "BR", "IOD", "F", "LI", "RB", "CS", "SR", "BA", "AL",
}
BUFFER_SALT_SOLVENT = {
    "SO4", "PO4", "NO3", "CO3", "ACT", "ACE", "FMT", "GOL", "EDO", "DMS", "PEG",
    "PG4", "PGE", "MPD", "TRS", "MES", "HEP", "BME", "IPA", "EOH", "MOH",
}
FUNCTIONAL_CONTEXT = {
    "ATP", "ADP", "AMP", "ANP", "ACP", "GDP", "GTP", "GNP", "GMP", "NAD", "NAP",
    "FAD", "FMN", "HEM", "HEC", "PLP", "SAM", "SAH", "COA", "NCT", "ADN", "ACH",
}


def load_tsv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def split_entries(value):
    return [part.strip() for part in value.split(";") if part.strip()]


def ligand_code(label):
    if label in {"()", ""}:
        return ""
    matches = re.findall(r"\(([A-Za-z0-9]{1,8})\)", label)
    return matches[-1].upper() if matches else ""


def ligand_name(label, code, chembl_id):
    if label.strip() == "()":
        return ""
    name = label.strip()
    if chembl_id:
        name = re.sub(r"CHEMBL\d+", "", name, flags=re.I).strip()
    if code:
        name = re.sub(rf"\({re.escape(code)}\)\s*$", "", name, flags=re.I).strip()
    return name.strip()


def parse_entry(target, symbol, source, raw_entry):
    if "@" not in raw_entry:
        return None
    label, rest = raw_entry.split("@", 1)
    label = label.strip()
    if not label or re.match(r"^[A-Za-z]:[A-Z]{3}\d+", label) or label.startswith(":"):
        return None
    code = ligand_code(label)
    chembl_match = re.search(r"(CHEMBL\d+)", label, flags=re.I)
    chembl_id = chembl_match.group(1).upper() if chembl_match else ""
    name = ligand_name(label, code, chembl_id)

    pdb_id = ""
    chain = ""
    residues = ""
    affinity_type = ""
    affinity_value = ""
    affinity_unit = ""

    if source == "BioLiP":
        m = re.match(r"([^:=]+):([^=]+)=(.*)", rest)
        if m:
            pdb_id, chain, residues = m.group(1), m.group(2), m.group(3)
    elif source == "sc-PDB":
        m = re.match(r"([^:]+):(.*)", rest)
        if m:
            pdb_id, residues = m.group(1), m.group(2)
    elif source == "PDBbind":
        m = re.match(r"([^:]+):([^=]+)=([^=]+)=([^:]*):(.*)", rest)
        if m:
            pdb_id, affinity_type, affinity_raw, chain, residues = m.groups()
            am = re.match(r"([-+0-9.eE]+)\s*([A-Za-zµμ]*)", affinity_raw)
            if am:
                affinity_value, affinity_unit = am.group(1), am.group(2).replace("μ", "u").replace("µ", "u")
            else:
                affinity_value = affinity_raw
        else:
            m = re.match(r"([^:]+):(.*)", rest)
            if m:
                pdb_id, residues = m.group(1), m.group(2)

    pdb_id = pdb_id.upper()
    code = code.upper()
    category, keep_small, keep_context, reason = classify_ligand(code, name, chembl_id, affinity_type)

    return {
        "target_uniprot_id": target,
        "approved_symbol": symbol,
        "structure_source": source,
        "raw_ligand_label": label,
        "ligand_name_raw": name,
        "pdb_ligand_code": code,
        "chembl_id": chembl_id,
        "pdb_id": pdb_id,
        "chain": chain,
        "residues": residues,
        "affinity_type": affinity_type,
        "affinity_value": affinity_value,
        "affinity_unit": affinity_unit,
        "ligand_category": category,
        "keep_for_small_molecule_binding_site": keep_small,
        "keep_for_binding_site_context": keep_context,
        "exclusion_reason": reason,
        "raw_entry": raw_entry,
    }


def classify_ligand(code, name, chembl_id, affinity_type):
    if not code and not name and not chembl_id:
        return "empty_or_unparsed", "0", "0", "empty ligand label"
    if code in WATER or name.upper() in WATER:
        return "water", "0", "0", "water molecule"
    if code in BUFFER_SALT_SOLVENT:
        return "buffer_salt_or_solvent", "0", "0", "crystallization solvent/buffer/salt"
    if code in IONS:
        return "ion_or_metal", "0", "1", "ion/metal context only"
    if code in FUNCTIONAL_CONTEXT:
        return "cofactor_or_functional_ligand", "0", "1", ""
    if chembl_id:
        return "chembl_mapped_small_molecule", "1", "1", ""
    if affinity_type:
        return "pdbbind_affinity_ligand", "1", "1", ""
    if name and len(name) > 3:
        return "organic_or_named_small_molecule", "1", "1", ""
    if code:
        return "pdb_ligand_code_unknown", "1", "1", "review PDB ligand code"
    return "unknown", "0", "0", "unmapped unknown"


def normalize_key(value):
    return re.sub(r"\s+", " ", value.strip().lower())


def build_local_indexes(small_rows):
    by_name = {}
    by_iupac = {}
    cids = set()
    for row in small_rows:
        drug_id = row["drug_id"]
        if row.get("compound_id_type") == "PubChem CID" and drug_id.isdigit():
            cids.add(drug_id)
        if row.get("drug_name"):
            by_name.setdefault(normalize_key(row["drug_name"]), drug_id)
        if row.get("iupac_name"):
            by_iupac.setdefault(normalize_key(row["iupac_name"]), drug_id)
    return by_name, by_iupac, cids


def pubchem_lookup(query):
    if not query:
        return ""
    encoded = urllib.parse.quote(query, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/cids/TXT"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as response:
                text = response.read().decode("utf-8", errors="ignore").strip()
            for line in text.splitlines():
                line = line.strip()
                if line.isdigit():
                    return line
            return ""
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return ""
        except Exception:
            if attempt == 2:
                return ""
        time.sleep(0.5 * (attempt + 1))
    return ""


def resolve_cids(records, by_name, by_iupac, existing_cids):
    if CACHE.exists():
        cache = json.loads(CACHE.read_text(encoding="utf-8"))
    else:
        cache = {}

    skip_pubchem = os.environ.get("SKIP_PUBCHEM_LOOKUP", "0") == "1"
    need_queries = {}
    for rec in records:
        candidates = []
        if rec["chembl_id"]:
            candidates.append(("pubchem_name_chembl", rec["chembl_id"]))
        if rec["ligand_name_raw"]:
            key = normalize_key(rec["ligand_name_raw"])
            if key in by_name:
                rec["ligand_pubchem_cid"] = by_name[key]
                rec["cid_mapping_method"] = "local_small_molecule_drug_name"
                rec["cid_mapping_status"] = "mapped"
                continue
            if key in by_iupac:
                rec["ligand_pubchem_cid"] = by_iupac[key]
                rec["cid_mapping_method"] = "local_small_molecule_iupac"
                rec["cid_mapping_status"] = "mapped"
                continue
            candidates.append(("pubchem_name_ligand_name", rec["ligand_name_raw"]))
        if rec["pdb_ligand_code"]:
            candidates.append(("unichem_inchikey_pubchem", rec["pdb_ligand_code"]))
            candidates.append(("rcsb_pubchem_xref", rec["pdb_ligand_code"]))
            candidates.append(("pubchem_name_pdb_ligand_code", rec["pdb_ligand_code"]))

        rec["_candidates"] = candidates
        if not skip_pubchem:
            for method, query in candidates:
                key = f"{method}:{query}"
                if key not in cache:
                    need_queries[key] = query

    if need_queries:
        with ThreadPoolExecutor(max_workers=8) as executor:
            futures = {executor.submit(pubchem_lookup, query): key for key, query in need_queries.items()}
            for future in as_completed(futures):
                key = futures[future]
                cache[key] = future.result()
        CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")

    for rec in records:
        if rec.get("cid_mapping_status") == "mapped":
            pass
        else:
            rec["ligand_pubchem_cid"] = ""
            rec["cid_mapping_method"] = ""
            rec["cid_mapping_status"] = "unmapped"
            for method, query in rec.get("_candidates", []):
                cid = cache.get(f"{method}:{query}", "")
                if cid:
                    rec["ligand_pubchem_cid"] = cid
                    rec["cid_mapping_method"] = method
                    rec["cid_mapping_status"] = "mapped"
                    break
        cid = rec.get("ligand_pubchem_cid", "")
        rec["in_small_molecule_database"] = "1" if cid in existing_cids else "0"
        rec["matched_drug_id"] = cid if cid in existing_cids else ""
        rec.pop("_candidates", None)


def write_outputs(records):
    fields = [
        "structure_ligand_id",
        "target_uniprot_id",
        "approved_symbol",
        "ligand_pubchem_cid",
        "cid_mapping_status",
        "cid_mapping_method",
        "in_small_molecule_database",
        "matched_drug_id",
        "ligand_category",
        "keep_for_small_molecule_binding_site",
        "keep_for_binding_site_context",
        "exclusion_reason",
        "structure_source",
        "pdb_ligand_code",
        "chembl_id",
        "ligand_name_raw",
        "pdb_id",
        "chain",
        "residues",
        "affinity_type",
        "affinity_value",
        "affinity_unit",
        "raw_ligand_label",
        "raw_entry",
    ]
    with OUT_TSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    wb = Workbook()
    ws = wb.active
    ws.title = "Structure_Ligand_Context"
    ws.append(fields)
    for rec in records:
        ws.append([rec.get(field, "") for field in fields])
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(fields))}{len(records) + 1}"
    fill = PatternFill("solid", fgColor="1F4E79")
    font = Font(bold=True, color="FFFFFF")
    for cell in ws[1]:
        cell.fill = fill
        cell.font = font
        cell.alignment = Alignment(wrap_text=True, vertical="top")
    for col_idx, field in enumerate(fields, 1):
        max_len = len(field)
        for row_idx in range(2, min(len(records) + 2, 250)):
            max_len = max(max_len, len(str(ws.cell(row_idx, col_idx).value or "")))
        ws.column_dimensions[get_column_letter(col_idx)].width = max(10, min(max_len + 2, 45))
    wb.save(OUT_XLSX)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def update_manifest():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    files = manifest.get("files", [])
    by_name = {Path(item["path"]).name: item for item in files}
    for path in [OUT_TSV, OUT_XLSX, SUMMARY]:
        item = by_name.get(path.name)
        if item is None:
            item = {"path": f"upload/normalized_tables/{path.name}"}
            files.append(item)
        item["sha256"] = sha256(path)
        item["bytes"] = path.stat().st_size
    manifest["files"] = files
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main():
    protein_rows = load_tsv(PROTEIN)
    small_rows = load_tsv(SMALL)
    by_name, by_iupac, existing_cids = build_local_indexes(small_rows)

    records = []
    seen = set()
    for row in protein_rows:
        for column, source in SOURCE_COLUMNS.items():
            for raw_entry in split_entries(row.get(column, "")):
                rec = parse_entry(row["target_uniprot_id"], row["approved_symbol"], source, raw_entry)
                if not rec:
                    continue
                key = (
                    rec["target_uniprot_id"],
                    rec["structure_source"],
                    rec["raw_ligand_label"],
                    rec["pdb_id"],
                    rec["chain"],
                    rec["residues"],
                    rec["affinity_type"],
                    rec["affinity_value"],
                )
                if key in seen:
                    continue
                seen.add(key)
                records.append(rec)

    resolve_cids(records, by_name, by_iupac, existing_cids)
    for index, rec in enumerate(records, 1):
        rec["structure_ligand_id"] = f"SLCTX_{index:06d}"

    write_outputs(records)

    summary = {
        "records": len(records),
        "mapped_to_pubchem_cid": sum(1 for r in records if r["cid_mapping_status"] == "mapped"),
        "in_small_molecule_database": sum(1 for r in records if r["in_small_molecule_database"] == "1"),
        "keep_for_small_molecule_binding_site": sum(1 for r in records if r["keep_for_small_molecule_binding_site"] == "1"),
        "keep_for_binding_site_context": sum(1 for r in records if r["keep_for_binding_site_context"] == "1"),
        "by_source": Counter(r["structure_source"] for r in records),
        "by_category": Counter(r["ligand_category"] for r in records),
        "by_mapping_method": Counter(r["cid_mapping_method"] or "unmapped" for r in records),
        "top_unmapped_labels": Counter(
            r["raw_ligand_label"] for r in records if r["cid_mapping_status"] != "mapped"
        ).most_common(50),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=dict) + "\n", encoding="utf-8")
    update_manifest()
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=dict))


if __name__ == "__main__":
    main()
