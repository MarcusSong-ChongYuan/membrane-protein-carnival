#!/usr/bin/env python3
"""Build MemPro expression/localization add-on from local HPA and V6.0.

The frozen protein release remains read-only. Detailed expression/localization
relations are emitted as long tables; only a one-row-per-protein summary is
prepared for a future V6.2 merge after V6.1 is frozen.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import re
import shutil
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
MASTER = ROOT / "releases" / "release_mempro_v6_0_20260727" / "human_membrane_protein_master_v6_0.tsv"
HPA_ZIP = Path(r"D:\7.22\membrane_master_v5_working\raw\hpa_proteinatlas_v25_1.tsv.zip")
RUN = ROOT / "runs" / "expression_location_v1_20260728"
RAW = RUN / "raw_ontologies"
STAGING = RUN / "staging"
QA = RUN / "qa"
PROGRESS = QA / "EXPRESSION_LOCATION_V1_PROGRESS.json"

TISSUE_OUT = STAGING / "protein_tissue_expression_v1.tsv.gz"
CELL_OUT = STAGING / "protein_cell_type_expression_v1.tsv.gz"
LOCATION_OUT = STAGING / "protein_subcellular_localization_v1.tsv.gz"
MAPPING_OUT = STAGING / "anatomy_cell_location_ontology_mapping_v1.tsv"
SUMMARY_OUT = STAGING / "expression_location_summary_v1.tsv.gz"
REVIEW_OUT = STAGING / "expression_location_conflict_review_v1.tsv.gz"
PREVIEW_OUT = STAGING / "human_membrane_protein_master_v6_2_expression_preview.tsv.gz"
REPORT_OUT = QA / "EXPRESSION_LOCATION_V1_VALIDATION.json"

ONTOLOGIES = {
    "UBERON": (
        "https://purl.obolibrary.org/obo/uberon/basic.obo",
        RAW / "uberon-basic.obo",
    ),
    "CL": (
        "https://purl.obolibrary.org/obo/cl/cl-basic.obo",
        RAW / "cl-basic.obo",
    ),
    "GO": (
        "https://purl.obolibrary.org/obo/go/go-basic.obo",
        RAW / "go-basic.obo",
    ),
}

HPA_VERSION = "HPA 25.1"
MODULE_VERSION = "expression_location_v1_20260728"


ORGAN_SYSTEM_MAP = {
    "adipose tissue": "endocrine/metabolic",
    "adrenal gland": "endocrine/metabolic",
    "blood vessel": "cardiovascular",
    "bone marrow": "immune/hematologic",
    "brain": "nervous",
    "cerebral cortex": "nervous",
    "breast": "reproductive",
    "cervix": "reproductive",
    "choroid plexus": "nervous",
    "endometrium": "reproductive",
    "epididymis": "reproductive",
    "esophagus": "digestive/hepatobiliary",
    "fallopian tube": "reproductive",
    "gallbladder": "digestive/hepatobiliary",
    "heart muscle": "cardiovascular",
    "intestine": "digestive/hepatobiliary",
    "kidney": "renal/urinary",
    "liver": "digestive/hepatobiliary",
    "lung": "respiratory",
    "lymphoid tissue": "immune/hematologic",
    "ovary": "reproductive",
    "pancreas": "digestive/hepatobiliary;endocrine/metabolic",
    "parathyroid gland": "endocrine/metabolic",
    "pituitary gland": "endocrine/metabolic",
    "placenta": "reproductive",
    "prostate": "reproductive",
    "retina": "sensory",
    "salivary gland": "digestive/hepatobiliary",
    "seminal vesicle": "reproductive",
    "skeletal muscle": "musculoskeletal",
    "skin": "integumentary",
    "smooth muscle": "musculoskeletal",
    "stomach": "digestive/hepatobiliary",
    "testis": "reproductive",
    "thyroid gland": "endocrine/metabolic",
    "tongue": "digestive/hepatobiliary",
    "urinary bladder": "renal/urinary",
    "vagina": "reproductive",
}

ONTOLOGY_OVERRIDES = {
    ("UBERON", "adipose subcutaneous"): "UBERON:0002190",
    ("UBERON", "adipose visceral"): "UBERON:0014454",
    ("UBERON", "lung"): "UBERON:0002048",
    ("UBERON", "pancreas"): "UBERON:0001264",
    ("UBERON", "skeletal muscle"): "UBERON:0001134",
    ("UBERON", "skin"): "UBERON:0002097",
    ("UBERON", "blood vessel"): "UBERON:0001981",
    ("UBERON", "brain"): "UBERON:0000955",
    ("UBERON", "endometrium"): "UBERON:0001295",
    ("UBERON", "esophagus"): "UBERON:0001043",
    ("CL", "adipocytes"): "CL:0000136",
    ("CL", "b cells"): "CL:0000236",
    ("CL", "beta cells"): "CL:0000169",
    ("CL", "fibroblasts"): "CL:0000057",
    ("CL", "macrophages"): "CL:0000235",
    ("CL", "monocytes"): "CL:0000576",
    ("CL", "neutrophils"): "CL:0000775",
    ("CL", "nk cells"): "CL:0000623",
    ("CL", "schwann cells"): "CL:0002573",
    ("CL", "alveolar cells type 1"): "CL:0002062",
    ("CL", "alveolar cells type 2"): "CL:0002063",
    ("CL", "corticotropes"): "CL:0002309",
    ("CL", "proximal tubular cells"): "CL:0002306",
    ("CL", "proximal tubule cells"): "CL:0002306",
    ("CL", "microglia and neuropil"): "CL:0000129",
    ("CL", "astrocytes and neuropil"): "CL:0000127",
    ("CL", "thyroid glandular cells"): "CL:0002258",
    ("CL", "oligodendrocyte progenitor cells"): "CL:0002453",
    ("CL", "adrenal medulla cells"): "CL:0000336",
    ("CL", "bergmann glia"): "CL:0000644",
    ("CL", "minor salivary glandular cells"): "CL:1001596",
    ("CL", "neutrophil progenitors"): "CL:0000834",
    ("CL", "syncytiotrophoblasts"): "CL:0000525",
    ("CL", "spermatogonia"): "CL:0000020",
    ("CL", "distal convoluted tubule cells"): "CL:1000849",
    ("CL", "renal collecting duct intercalated cells"): "CL:1001432",
    ("CL", "renal collecting duct principal cells"): "CL:1001431",
    ("CL", "renal connecting tubule cells"): "CL:1000768",
    ("CL", "loop of henle epithelial cells"): "CL:1000909",
    ("CL", "retinal horizontal cells"): "CL:0000745",
    ("CL", "cytotrophoblasts"): "CL:0000523",
    ("CL", "enteric glia cells"): "CL:4040002",
    ("CL", "muller glia"): "CL:0000636",
    ("CL", "foveolar cells"): "CL:0002179",
    ("CL", "monocyte progenitors"): "CL:0000557",
    ("CL", "endometrial glandular cells"): "CL:0009084",
    ("CL", "colon enterocytes"): "CL:1000347",
    ("CL", "skeletal myofibers"): "CL:0008002",
    ("CL", "enteric transient amplifying cells"): "CL:0009010",
    ("CL", "endometrial luminal cells"): "CL:4052050",
    ("CL", "skeletal myocytes"): "CL:0000187",
    ("GO", "calyx"): "GO:0033150",
    ("GO", "connecting piece"): "GO:0120212",
    ("GO", "cytokinetic bridge"): "GO:0045171",
    ("GO", "end piece"): "GO:0097229",
    ("GO", "focal adhesion sites"): "GO:0005925",
    ("GO", "mid piece"): "GO:0097225",
    ("GO", "nuclear bodies"): "GO:0016604",
    ("GO", "nucleoli"): "GO:0005730",
    ("GO", "nucleoli fibrillar center"): "GO:0001650",
    ("GO", "primary cilium"): "GO:0005929",
    ("GO", "primary cilium tip"): "GO:0097542",
    ("GO", "primary cilium transition zone"): "GO:0035869",
    ("GO", "principal piece"): "GO:0097228",
}

NON_CELL_COMPARTMENT_TERMS = {
    "capillaries",
    "collecting ducts",
    "distal tubules",
    "glomeruli",
    "large ducts",
    "myonuclei",
    "pancreatic islets",
    "proximal tubules",
}


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp.replace(path)


def set_progress(stage: str, **extra) -> None:
    payload = {
        "status": "running",
        "stage": stage,
        "updated_utc": now(),
        "module_version": MODULE_VERSION,
        "input_policy": "V6.0 and HPA 25.1 read-only",
        "outputs": {
            "tissue_expression": str(TISSUE_OUT),
            "cell_type_expression": str(CELL_OUT),
            "subcellular_localization": str(LOCATION_OUT),
            "ontology_mapping": str(MAPPING_OUT),
            "protein_summary": str(SUMMARY_OUT),
            "review_queue": str(REVIEW_OUT),
            "v6_2_preview": str(PREVIEW_OUT),
        },
    }
    payload.update(extra)
    write_json(PROGRESS, payload)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_label(value: str) -> str:
    original = value.strip().lower()
    strip_trailing_index = bool(re.search(r"_\d+$", original)) or original in {
        "endometrium 1",
        "skin 1",
        "stomach 1",
    }
    value = original
    value = value.replace("&", " and ")
    value = re.sub(r"[-_/]+", " ", value)
    value = re.sub(r"[^a-z0-9 ]+", " ", value)
    value = re.sub(r"\s+", " ", value).strip()
    if strip_trailing_index:
        value = re.sub(r"\s+\d+$", "", value).strip()
    return value


def label_variants(value: str) -> list[str]:
    base = normalize_label(value)
    contextless = normalize_label(re.sub(r"\s*\([^)]*\)\s*$", "", value))
    variants = {base, contextless}
    variants.add(re.sub(r"\bcells\b", "cell", base))
    variants.add(re.sub(r"\bcells\b", "cell", contextless))
    variants.add(re.sub(r"\btypes\b", "type", base))
    variants.add(re.sub(r"\bfilaments\b", "filament", base))
    variants.add(re.sub(r"\bvesicles\b", "vesicle", base))
    variants.add(re.sub(r"\blymphocytes\b", "lymphocyte", base))
    variants.add(re.sub(r"\bcytes\b", "cyte", base))
    variants.add(re.sub(r"s$", "", base))
    variants.add(re.sub(r"s$", "", contextless))
    return sorted(item for item in variants if item)


def download_if_missing(url: str, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size > 10000:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(destination.suffix + ".tmp")
    request = urllib.request.Request(
        url, headers={"User-Agent": "MemPro-expression-location/1.0"}
    )
    with urllib.request.urlopen(request, timeout=300) as response, tmp.open("wb") as out:
        shutil.copyfileobj(response, out)
    tmp.replace(destination)


def parse_obo(path: Path) -> tuple[dict[str, set[str]], dict[str, str], str]:
    lookup: dict[str, set[str]] = defaultdict(set)
    names: dict[str, str] = {}
    data_version = ""
    current: dict[str, object] | None = None

    def commit(term: dict[str, object] | None) -> None:
        if not term or term.get("obsolete") or not term.get("id"):
            return
        term_id = str(term["id"])
        name = str(term.get("name") or "")
        if name:
            names[term_id] = name
        for label in [name, *list(term.get("synonyms") or [])]:
            for variant in label_variants(label):
                lookup[variant].add(term_id)

    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("data-version:") and not data_version:
                data_version = line.split(":", 1)[1].strip()
            if line == "[Term]":
                commit(current)
                current = {"synonyms": [], "obsolete": False}
            elif line.startswith("[") and line.endswith("]"):
                commit(current)
                current = None
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:].strip()
                elif line.startswith("name: "):
                    current["name"] = line[6:].strip()
                elif line.startswith("synonym: "):
                    match = re.match(r'synonym:\s+"(.+?)"', line)
                    if match:
                        current["synonyms"].append(match.group(1))
                elif line == "is_obsolete: true":
                    current["obsolete"] = True
    commit(current)
    return lookup, names, data_version


def ontology_match(
    label: str,
    ontology: str,
    lookup: dict[str, set[str]],
    names: dict[str, str],
) -> dict[str, str]:
    override_key = (ontology, normalize_label(label))
    override_id = ONTOLOGY_OVERRIDES.get(override_key)
    if override_id:
        return {
            "ontology": ontology,
            "ontology_id": override_id,
            "ontology_label": names.get(override_id, ""),
            "mapping_status": "mapped_curated_override",
            "mapping_method": "curated_HPA_term_override",
        }
    if ontology == "CL" and normalize_label(label) in NON_CELL_COMPARTMENT_TERMS:
        return {
            "ontology": ontology,
            "ontology_id": "",
            "ontology_label": "",
            "mapping_status": "source_histology_compartment_not_CL_cell_type",
            "mapping_method": "curated_HPA_histology_compartment_classification",
        }
    candidates: set[str] = set()
    matched_variant = ""
    variants = label_variants(label)
    for variant in variants:
        found = lookup.get(variant, set())
        if found:
            candidates.update(found)
            if not matched_variant:
                matched_variant = variant
    base = normalize_label(label)
    exact_official = {
        term_id
        for term_id in candidates
        if normalize_label(names.get(term_id, "")) == base
    }
    if len(exact_official) == 1:
        term_id = next(iter(exact_official))
        return {
            "ontology": ontology,
            "ontology_id": term_id,
            "ontology_label": names.get(term_id, ""),
            "mapping_status": "mapped_exact_official_label",
            "mapping_method": f"official_label:{base}",
        }
    normalized_official = {
        term_id
        for term_id in candidates
        if normalize_label(names.get(term_id, "")) in variants
    }
    if len(normalized_official) == 1:
        term_id = next(iter(normalized_official))
        return {
            "ontology": ontology,
            "ontology_id": term_id,
            "ontology_label": names.get(term_id, ""),
            "mapping_status": "mapped_unique_normalized_official_label",
            "mapping_method": f"official_label_variant:{normalize_label(names.get(term_id, ''))}",
        }
    if ontology == "CL":
        broad_blacklist = {
            "cell",
            "epithelial cell",
            "glandular cell",
            "stem cell",
            "progenitor cell",
            "secretory cell",
            "basal cell",
            "neuron",
            "glial cell",
        }
        source_variants = label_variants(label)
        contained = []
        for term_id, official_name in names.items():
            official = normalize_label(official_name)
            if official in broad_blacklist or len(official.split()) < 2:
                continue
            if any(
                source == official
                or source.endswith(" " + official)
                or (" " + official + " ") in (" " + source + " ")
                for source in source_variants
            ):
                contained.append((len(official.split()), len(official), term_id))
        if contained:
            best_score = max((tokens, length) for tokens, length, _ in contained)
            best = {
                term_id
                for tokens, length, term_id in contained
                if (tokens, length) == best_score
            }
            if len(best) == 1:
                term_id = next(iter(best))
                return {
                    "ontology": ontology,
                    "ontology_id": term_id,
                    "ontology_label": names.get(term_id, ""),
                    "mapping_status": "mapped_broader_contained_official_label",
                    "mapping_method": "conservative_longest_contained_CL_label",
                }
    if len(candidates) == 1:
        term_id = next(iter(candidates))
        return {
            "ontology": ontology,
            "ontology_id": term_id,
            "ontology_label": names.get(term_id, ""),
            "mapping_status": "mapped_unique_label_or_synonym",
            "mapping_method": f"normalized:{matched_variant}",
        }
    if len(candidates) > 1:
        return {
            "ontology": ontology,
            "ontology_id": ";".join(sorted(candidates)),
            "ontology_label": ";".join(names.get(item, "") for item in sorted(candidates)),
            "mapping_status": "ambiguous_multiple_terms",
            "mapping_method": f"normalized:{matched_variant}",
        }
    return {
        "ontology": ontology,
        "ontology_id": "",
        "ontology_label": "",
        "mapping_status": "unmapped",
        "mapping_method": "",
    }


def parse_accessions(value: str) -> list[str]:
    return sorted(
        set(
            re.findall(
                r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z0-9]{3})?)\b",
                value.upper(),
            )
        )
    )


def split_ids(value: str) -> list[str]:
    return [item.strip() for item in re.split(r"[,;| ]+", value) if item.strip()]


def load_master():
    rows = []
    by_accession = {}
    by_ensembl: dict[str, set[str]] = defaultdict(set)
    by_symbol: dict[str, set[str]] = defaultdict(set)
    with MASTER.open("rt", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row["target_uniprot_id"].strip()
            rows.append(row)
            by_accession[accession] = row
            for value in split_ids(row.get("ensembl_gene_ids", "")):
                by_ensembl[value.split(".", 1)[0]].add(accession)
            symbol = row.get("approved_symbol", "").strip().upper()
            if symbol:
                by_symbol[symbol].add(accession)
    return rows, by_accession, by_ensembl, by_symbol


def map_hpa_row(
    row: dict,
    by_accession: dict,
    by_ensembl: dict[str, set[str]],
    by_symbol: dict[str, set[str]],
) -> tuple[list[str], str]:
    accessions = [item for item in parse_accessions(row.get("Uniprot", "")) if item in by_accession]
    if accessions:
        return sorted(set(accessions)), "uniprot_exact"
    ensembl = row.get("Ensembl", "").split(".", 1)[0].strip()
    if ensembl and len(by_ensembl.get(ensembl, set())) == 1:
        return sorted(by_ensembl[ensembl]), "ensembl_unique"
    symbol = row.get("Gene", "").strip().upper()
    if symbol and len(by_symbol.get(symbol, set())) == 1:
        return sorted(by_symbol[symbol]), "gene_symbol_unique"
    candidates = set()
    if ensembl:
        candidates.update(by_ensembl.get(ensembl, set()))
    if symbol:
        candidates.update(by_symbol.get(symbol, set()))
    return sorted(candidates), "ambiguous" if candidates else "unmapped"


def parse_numeric_pairs(value: str) -> list[tuple[str, float]]:
    result = []
    for item in value.split(";"):
        item = item.strip()
        if not item or ":" not in item:
            continue
        label, raw_number = item.rsplit(":", 1)
        try:
            result.append((label.strip(), float(raw_number.strip())))
        except ValueError:
            continue
    return result


def parse_csv_labels(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_go_cc(value: str) -> list[tuple[str, str]]:
    result = []
    for item in value.split(";"):
        item = item.strip()
        match = re.search(r"(.+?)\s*\[(GO:\d{7})\]", item)
        if match:
            result.append((match.group(1).strip(), match.group(2)))
    return result


def organ_system(tissue: str) -> str:
    return ORGAN_SYSTEM_MAP.get(normalize_label(tissue), "")


def expression_evidence(layer: str, value: float | None, detected: str) -> str:
    if layer == "protein":
        return "PEX1"
    if value is not None or "detected" in detected.lower():
        return "PEX3"
    return "PEX0"


def location_evidence(source: str, reliability: str) -> str:
    if source == "HPA_IF":
        if reliability in {"Enhanced", "Supported"}:
            return "LOC1"
        if reliability == "Approved":
            return "LOC2"
        return "LOC3"
    return "LOC2"


def main() -> None:
    for directory in (RAW, STAGING, QA):
        directory.mkdir(parents=True, exist_ok=True)
    started = now()
    set_progress(
        "freeze_inputs",
        started_utc=started,
        input_hashes={
            "protein_master_sha256": sha256(MASTER),
            "hpa_25_1_sha256": sha256(HPA_ZIP),
        },
    )

    ontology_data = {}
    for ontology, (url, path) in ONTOLOGIES.items():
        set_progress("download_ontologies", started_utc=started, current_ontology=ontology)
        download_if_missing(url, path)
        lookup, names, version = parse_obo(path)
        ontology_data[ontology] = {
            "lookup": lookup,
            "names": names,
            "version": version,
            "path": str(path),
            "sha256": sha256(path),
        }

    master_rows, by_accession, by_ensembl, by_symbol = load_master()
    hpa_by_accession: dict[str, list[dict]] = defaultdict(list)
    identifier_reviews = []
    hpa_total = 0
    set_progress("map_hpa_identifiers", started_utc=started)
    with zipfile.ZipFile(HPA_ZIP) as archive, archive.open("proteinatlas.tsv") as raw:
        text = io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")
        for hpa_row in csv.DictReader(text, delimiter="\t"):
            hpa_total += 1
            accessions, method = map_hpa_row(
                hpa_row, by_accession, by_ensembl, by_symbol
            )
            if method == "ambiguous":
                identifier_reviews.append(
                    {
                        "issue_type": "hpa_identifier_mapping",
                        "severity": "review",
                        "target_uniprot_id": ";".join(accessions),
                        "source_record": hpa_row.get("Gene", ""),
                        "details": method,
                    }
                )
            for accession in accessions:
                copied = dict(hpa_row)
                copied["_mapping_method"] = method
                hpa_by_accession[accession].append(copied)

    tissue_terms: set[str] = set()
    cell_terms: set[str] = set()
    location_terms: set[str] = set()
    for records in hpa_by_accession.values():
        for row in records:
            tissue_terms.update(name for name, _ in parse_numeric_pairs(row.get("RNA tissue specific nTPM", "")))
            tissue_terms.update(name for name, _ in parse_numeric_pairs(row.get("Protein tissue specific Intensity", "")))
            cell_terms.update(name for name, _ in parse_numeric_pairs(row.get("RNA single cell type specific nCPM", "")))
            cell_terms.update(name for name, _ in parse_numeric_pairs(row.get("Protein cell type specific Intensity", "")))
            for item in parse_csv_labels(row.get("RNA tissue cell type enrichment", "")):
                if " - " in item:
                    tissue, cell = item.split(" - ", 1)
                    tissue_terms.add(tissue.strip())
                    cell_terms.add(cell.strip())
            location_terms.update(parse_csv_labels(row.get("Subcellular location", "")))

    tissue_map = {
        term: ontology_match(
            term,
            "UBERON",
            ontology_data["UBERON"]["lookup"],
            ontology_data["UBERON"]["names"],
        )
        for term in sorted(tissue_terms)
    }
    cell_map = {
        term: ontology_match(
            term,
            "CL",
            ontology_data["CL"]["lookup"],
            ontology_data["CL"]["names"],
        )
        for term in sorted(cell_terms)
    }
    location_map = {
        term: ontology_match(
            term,
            "GO",
            ontology_data["GO"]["lookup"],
            ontology_data["GO"]["names"],
        )
        for term in sorted(location_terms)
    }

    with MAPPING_OUT.open("wt", encoding="utf-8", newline="") as handle:
        fields = [
            "source_term_type",
            "source_term",
            "normalized_term",
            "organ_system",
            "ontology",
            "ontology_id",
            "ontology_label",
            "mapping_status",
            "mapping_method",
            "module_version",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for term_type, mapping in (
            ("tissue", tissue_map),
            ("cell_type", cell_map),
            ("subcellular_location", location_map),
        ):
            for term, item in mapping.items():
                writer.writerow(
                    {
                        "source_term_type": term_type,
                        "source_term": term,
                        "normalized_term": normalize_label(term),
                        "organ_system": organ_system(term) if term_type == "tissue" else "",
                        **item,
                        "module_version": MODULE_VERSION,
                    }
                )

    set_progress(
        "build_long_tables",
        started_utc=started,
        mapped_hpa_proteins=len(hpa_by_accession),
        hpa_source_rows=hpa_total,
    )
    tissue_fields = [
        "expression_relation_id",
        "target_uniprot_id",
        "hpa_gene",
        "hpa_ensembl",
        "identifier_mapping_method",
        "tissue_name_hpa",
        "tissue_name_normalized",
        "uberon_id",
        "uberon_label",
        "ontology_mapping_status",
        "organ_system",
        "measurement_layer",
        "measurement_type",
        "value",
        "unit",
        "specificity",
        "distribution",
        "evidence_level",
        "coverage_scope",
        "source_database",
        "source_version",
    ]
    cell_fields = [
        "expression_relation_id",
        "target_uniprot_id",
        "hpa_gene",
        "hpa_ensembl",
        "identifier_mapping_method",
        "tissue_context_hpa",
        "tissue_uberon_id",
        "cell_type_name_hpa",
        "cell_type_name_normalized",
        "cell_ontology_id",
        "cell_ontology_label",
        "ontology_mapping_status",
        "measurement_layer",
        "measurement_type",
        "value",
        "unit",
        "specificity",
        "distribution",
        "evidence_level",
        "coverage_scope",
        "source_database",
        "source_version",
    ]
    location_fields = [
        "localization_relation_id",
        "target_uniprot_id",
        "location_name",
        "location_role",
        "go_id",
        "go_label",
        "ontology_mapping_status",
        "source_database",
        "source_field",
        "source_reliability",
        "evidence_level",
        "direct_protein_observation_flag",
        "raw_source_value",
        "source_version",
    ]

    tissue_rows_by_protein: dict[str, list[dict]] = defaultdict(list)
    cell_rows_by_protein: dict[str, list[dict]] = defaultdict(list)
    location_rows_by_protein: dict[str, list[dict]] = defaultdict(list)
    review_rows = list(identifier_reviews)
    ontology_review_targets: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    tissue_count = cell_count = location_count = 0

    with gzip.open(TISSUE_OUT, "wt", encoding="utf-8", newline="") as tissue_handle, gzip.open(
        CELL_OUT, "wt", encoding="utf-8", newline=""
    ) as cell_handle, gzip.open(LOCATION_OUT, "wt", encoding="utf-8", newline="") as location_handle:
        tissue_writer = csv.DictWriter(tissue_handle, delimiter="\t", fieldnames=tissue_fields, lineterminator="\n")
        cell_writer = csv.DictWriter(cell_handle, delimiter="\t", fieldnames=cell_fields, lineterminator="\n")
        location_writer = csv.DictWriter(location_handle, delimiter="\t", fieldnames=location_fields, lineterminator="\n")
        tissue_writer.writeheader()
        cell_writer.writeheader()
        location_writer.writeheader()

        for accession, master in by_accession.items():
            for hpa in hpa_by_accession.get(accession, []):
                for layer, field, measurement, unit, specificity_field, distribution_field in (
                    ("rna", "RNA tissue specific nTPM", "specific_tissue_expression", "nTPM", "RNA tissue specificity", "RNA tissue distribution"),
                    ("protein", "Protein tissue specific Intensity", "specific_tissue_expression", "HPA_intensity", "Protein tissue specificity", "Protein tissue distribution"),
                ):
                    for tissue, value in parse_numeric_pairs(hpa.get(field, "")):
                        mapping = tissue_map[tissue]
                        tissue_count += 1
                        row = {
                            "expression_relation_id": f"TEX{tissue_count:09d}",
                            "target_uniprot_id": accession,
                            "hpa_gene": hpa.get("Gene", ""),
                            "hpa_ensembl": hpa.get("Ensembl", ""),
                            "identifier_mapping_method": hpa["_mapping_method"],
                            "tissue_name_hpa": tissue,
                            "tissue_name_normalized": normalize_label(tissue),
                            "uberon_id": mapping["ontology_id"],
                            "uberon_label": mapping["ontology_label"],
                            "ontology_mapping_status": mapping["mapping_status"],
                            "organ_system": organ_system(tissue),
                            "measurement_layer": layer,
                            "measurement_type": measurement,
                            "value": value,
                            "unit": unit,
                            "specificity": hpa.get(specificity_field, ""),
                            "distribution": hpa.get(distribution_field, ""),
                            "evidence_level": expression_evidence(layer, value, hpa.get(distribution_field, "")),
                            "coverage_scope": "HPA_specificity_subset_not_full_matrix",
                            "source_database": "Human Protein Atlas",
                            "source_version": HPA_VERSION,
                        }
                        tissue_writer.writerow(row)
                        tissue_rows_by_protein[accession].append(row)
                        if not mapping["mapping_status"].startswith("mapped_"):
                            ontology_review_targets[
                                ("tissue_ontology_mapping", tissue, mapping["mapping_status"])
                            ].add(accession)

                for layer, field, measurement, unit, specificity_field, distribution_field in (
                    ("rna", "RNA single cell type specific nCPM", "specific_cell_type_expression", "nCPM", "RNA single cell type specificity", "RNA single cell type distribution"),
                    ("protein", "Protein cell type specific Intensity", "specific_cell_type_expression", "HPA_intensity", "Protein cell type specificity", "Protein cell type distribution"),
                ):
                    for cell, value in parse_numeric_pairs(hpa.get(field, "")):
                        mapping = cell_map[cell]
                        cell_count += 1
                        row = {
                            "expression_relation_id": f"CEX{cell_count:09d}",
                            "target_uniprot_id": accession,
                            "hpa_gene": hpa.get("Gene", ""),
                            "hpa_ensembl": hpa.get("Ensembl", ""),
                            "identifier_mapping_method": hpa["_mapping_method"],
                            "tissue_context_hpa": "",
                            "tissue_uberon_id": "",
                            "cell_type_name_hpa": cell,
                            "cell_type_name_normalized": normalize_label(cell),
                            "cell_ontology_id": mapping["ontology_id"],
                            "cell_ontology_label": mapping["ontology_label"],
                            "ontology_mapping_status": mapping["mapping_status"],
                            "measurement_layer": layer,
                            "measurement_type": measurement,
                            "value": value,
                            "unit": unit,
                            "specificity": hpa.get(specificity_field, ""),
                            "distribution": hpa.get(distribution_field, ""),
                            "evidence_level": expression_evidence(layer, value, hpa.get(distribution_field, "")),
                            "coverage_scope": "HPA_specificity_subset_not_full_matrix",
                            "source_database": "Human Protein Atlas",
                            "source_version": HPA_VERSION,
                        }
                        cell_writer.writerow(row)
                        cell_rows_by_protein[accession].append(row)
                        if not mapping["mapping_status"].startswith("mapped_"):
                            ontology_review_targets[
                                ("cell_ontology_mapping", cell, mapping["mapping_status"])
                            ].add(accession)

                for enrichment in parse_csv_labels(hpa.get("RNA tissue cell type enrichment", "")):
                    if " - " not in enrichment:
                        continue
                    tissue, cell = [item.strip() for item in enrichment.split(" - ", 1)]
                    tissue_mapping = (
                        tissue_map[tissue]
                        if tissue in tissue_map
                        else ontology_match(
                            tissue,
                            "UBERON",
                            ontology_data["UBERON"]["lookup"],
                            ontology_data["UBERON"]["names"],
                        )
                    )
                    cell_mapping = (
                        cell_map[cell]
                        if cell in cell_map
                        else ontology_match(
                            cell,
                            "CL",
                            ontology_data["CL"]["lookup"],
                            ontology_data["CL"]["names"],
                        )
                    )
                    cell_count += 1
                    row = {
                        "expression_relation_id": f"CEX{cell_count:09d}",
                        "target_uniprot_id": accession,
                        "hpa_gene": hpa.get("Gene", ""),
                        "hpa_ensembl": hpa.get("Ensembl", ""),
                        "identifier_mapping_method": hpa["_mapping_method"],
                        "tissue_context_hpa": tissue,
                        "tissue_uberon_id": tissue_mapping["ontology_id"],
                        "cell_type_name_hpa": cell,
                        "cell_type_name_normalized": normalize_label(cell),
                        "cell_ontology_id": cell_mapping["ontology_id"],
                        "cell_ontology_label": cell_mapping["ontology_label"],
                        "ontology_mapping_status": cell_mapping["mapping_status"],
                        "measurement_layer": "rna",
                        "measurement_type": "tissue_cell_type_enrichment",
                        "value": "",
                        "unit": "qualitative",
                        "specificity": hpa.get("RNA single cell type specificity", ""),
                        "distribution": hpa.get("RNA single cell type distribution", ""),
                        "evidence_level": "PEX3",
                        "coverage_scope": "HPA_enrichment_subset",
                        "source_database": "Human Protein Atlas",
                        "source_version": HPA_VERSION,
                    }
                    cell_writer.writerow(row)
                    cell_rows_by_protein[accession].append(row)

                reliability = hpa.get("Reliability (IF)", "")
                main_locations = set(parse_csv_labels(hpa.get("Subcellular main location", "")))
                additional_locations = set(parse_csv_labels(hpa.get("Subcellular additional location", "")))
                for location in parse_csv_labels(hpa.get("Subcellular location", "")):
                    mapping = location_map[location]
                    role = "main" if location in main_locations else "additional" if location in additional_locations else "unspecified"
                    location_count += 1
                    row = {
                        "localization_relation_id": f"LOC{location_count:09d}",
                        "target_uniprot_id": accession,
                        "location_name": location,
                        "location_role": role,
                        "go_id": mapping["ontology_id"],
                        "go_label": mapping["ontology_label"],
                        "ontology_mapping_status": mapping["mapping_status"],
                        "source_database": "Human Protein Atlas",
                        "source_field": "Subcellular location",
                        "source_reliability": reliability,
                        "evidence_level": location_evidence("HPA_IF", reliability),
                        "direct_protein_observation_flag": "1",
                        "raw_source_value": hpa.get("Subcellular location", ""),
                        "source_version": HPA_VERSION,
                    }
                    location_writer.writerow(row)
                    location_rows_by_protein[accession].append(row)
                    if not mapping["mapping_status"].startswith("mapped_"):
                        ontology_review_targets[
                            (
                                "subcellular_ontology_mapping",
                                location,
                                mapping["mapping_status"],
                            )
                        ].add(accession)

            for go_label, go_id in parse_go_cc(master.get("go_cellular_component", "")):
                location_count += 1
                row = {
                    "localization_relation_id": f"LOC{location_count:09d}",
                    "target_uniprot_id": accession,
                    "location_name": go_label,
                    "location_role": "GO_annotation",
                    "go_id": go_id,
                    "go_label": ontology_data["GO"]["names"].get(go_id, go_label),
                    "ontology_mapping_status": "mapped_by_source_GO_ID",
                    "source_database": "UniProtKB/GO",
                    "source_field": "go_cellular_component",
                    "source_reliability": "source_evidence_code_not_available_in_master",
                    "evidence_level": "LOC2",
                    "direct_protein_observation_flag": "0",
                    "raw_source_value": master.get("go_cellular_component", ""),
                    "source_version": master.get("uniprot_release_v53", ""),
                }
                location_writer.writerow(row)
                location_rows_by_protein[accession].append(row)

    for (issue_type, source_term, status), affected in sorted(
        ontology_review_targets.items()
    ):
        review_rows.append(
            {
                "issue_type": issue_type,
                "severity": "review",
                "target_uniprot_id": "",
                "source_record": source_term,
                "details": f"{status};affected_proteins={len(affected)}",
            }
        )

    set_progress(
        "build_summary_and_review",
        started_utc=started,
        tissue_rows=tissue_count,
        cell_type_rows=cell_count,
        localization_rows=location_count,
    )
    summary_fields = [
        "target_uniprot_id",
        "expression_location_module_version",
        "hpa_identifier_mapping_status",
        "rna_tissue_specificity",
        "rna_tissue_distribution",
        "rna_top_tissues",
        "rna_tissue_count",
        "protein_tissue_specificity",
        "protein_tissue_distribution",
        "protein_top_tissues",
        "protein_tissue_count",
        "primary_organ_systems",
        "rna_cell_type_specificity",
        "rna_cell_type_distribution",
        "rna_top_cell_types",
        "rna_cell_type_count",
        "protein_cell_type_specificity",
        "protein_cell_type_distribution",
        "protein_top_cell_types",
        "protein_cell_type_count",
        "hpa_subcellular_main_locations",
        "hpa_subcellular_additional_locations",
        "hpa_subcellular_reliability",
        "go_cellular_component_ids",
        "plasma_membrane_localization_flag",
        "best_location_evidence_level",
        "has_rna_expression_evidence",
        "has_protein_expression_evidence",
        "expression_location_sources",
        "expression_protein_discordance_flag",
        "localization_source_divergence_flag",
        "review_issue_count",
    ]
    summary_by_accession = {}
    rank = {"LOC1": 1, "LOC2": 2, "LOC3": 3}
    for master in master_rows:
        accession = master["target_uniprot_id"]
        hpa_records = hpa_by_accession.get(accession, [])
        tissues = tissue_rows_by_protein.get(accession, [])
        cells = cell_rows_by_protein.get(accession, [])
        locations = location_rows_by_protein.get(accession, [])
        rna_tissues = [row for row in tissues if row["measurement_layer"] == "rna"]
        protein_tissues = [row for row in tissues if row["measurement_layer"] == "protein"]
        rna_cells = [row for row in cells if row["measurement_layer"] == "rna"]
        protein_cells = [row for row in cells if row["measurement_layer"] == "protein"]
        rna_tissue_names = sorted({row["tissue_name_hpa"] for row in rna_tissues})
        protein_tissue_names = sorted({row["tissue_name_hpa"] for row in protein_tissues})
        rna_cell_names = sorted({row["cell_type_name_hpa"] for row in rna_cells})
        protein_cell_names = sorted({row["cell_type_name_hpa"] for row in protein_cells})
        organ_systems = sorted(
            {
                value
                for row in tissues
                for value in row["organ_system"].split(";")
                if value
            }
        )
        hpa_main = sorted(
            {
                row["location_name"]
                for row in locations
                if row["source_database"] == "Human Protein Atlas"
                and row["location_role"] == "main"
            }
        )
        hpa_additional = sorted(
            {
                row["location_name"]
                for row in locations
                if row["source_database"] == "Human Protein Atlas"
                and row["location_role"] == "additional"
            }
        )
        go_ids = sorted({row["go_id"] for row in locations if row["go_id"]})
        loc_levels = [row["evidence_level"] for row in locations if row["evidence_level"] in rank]
        best_loc = min(loc_levels, key=lambda item: rank[item]) if loc_levels else ""
        rna_spec = sorted({row.get("RNA tissue specificity", "") for row in hpa_records if row.get("RNA tissue specificity", "")})
        rna_dist = sorted({row.get("RNA tissue distribution", "") for row in hpa_records if row.get("RNA tissue distribution", "")})
        prot_spec = sorted({row.get("Protein tissue specificity", "") for row in hpa_records if row.get("Protein tissue specificity", "")})
        prot_dist = sorted({row.get("Protein tissue distribution", "") for row in hpa_records if row.get("Protein tissue distribution", "")})
        rna_cell_spec = sorted({row.get("RNA single cell type specificity", "") for row in hpa_records if row.get("RNA single cell type specificity", "")})
        rna_cell_dist = sorted({row.get("RNA single cell type distribution", "") for row in hpa_records if row.get("RNA single cell type distribution", "")})
        prot_cell_spec = sorted({row.get("Protein cell type specificity", "") for row in hpa_records if row.get("Protein cell type specificity", "")})
        prot_cell_dist = sorted({row.get("Protein cell type distribution", "") for row in hpa_records if row.get("Protein cell type distribution", "")})
        reliabilities = sorted({row.get("Reliability (IF)", "") for row in hpa_records if row.get("Reliability (IF)", "")})
        mapping_methods = sorted({row["_mapping_method"] for row in hpa_records})
        discordance = bool(
            rna_tissues
            and any("not detected" in value.lower() for value in prot_spec + prot_dist)
        )
        hpa_go = {row["go_id"] for row in locations if row["source_database"] == "Human Protein Atlas" and row["go_id"]}
        source_go = {row["go_id"] for row in locations if row["source_database"] == "UniProtKB/GO" and row["go_id"]}
        high_hpa = any(
            row["source_database"] == "Human Protein Atlas"
            and row["evidence_level"] == "LOC1"
            for row in locations
        )
        divergence = bool(high_hpa and hpa_go and source_go and not (hpa_go & source_go))
        if discordance:
            review_rows.append(
                {
                    "issue_type": "rna_protein_expression_discordance",
                    "severity": "informational",
                    "target_uniprot_id": accession,
                    "source_record": ";".join(rna_tissue_names),
                    "details": "RNA-specific signal with HPA protein tissue not detected",
                }
            )
        if divergence:
            review_rows.append(
                {
                    "issue_type": "localization_source_divergence",
                    "severity": "review",
                    "target_uniprot_id": accession,
                    "source_record": ";".join(sorted(hpa_go)),
                    "details": "High-confidence HPA GO terms do not overlap source GO annotations; may represent multi-localization",
                }
            )
        summary_by_accession[accession] = {
            "target_uniprot_id": accession,
            "expression_location_module_version": MODULE_VERSION,
            "hpa_identifier_mapping_status": ";".join(mapping_methods) if mapping_methods else "no_HPA_match",
            "rna_tissue_specificity": ";".join(rna_spec),
            "rna_tissue_distribution": ";".join(rna_dist),
            "rna_top_tissues": ";".join(rna_tissue_names),
            "rna_tissue_count": len(rna_tissue_names),
            "protein_tissue_specificity": ";".join(prot_spec),
            "protein_tissue_distribution": ";".join(prot_dist),
            "protein_top_tissues": ";".join(protein_tissue_names),
            "protein_tissue_count": len(protein_tissue_names),
            "primary_organ_systems": ";".join(organ_systems),
            "rna_cell_type_specificity": ";".join(rna_cell_spec),
            "rna_cell_type_distribution": ";".join(rna_cell_dist),
            "rna_top_cell_types": ";".join(rna_cell_names),
            "rna_cell_type_count": len(rna_cell_names),
            "protein_cell_type_specificity": ";".join(prot_cell_spec),
            "protein_cell_type_distribution": ";".join(prot_cell_dist),
            "protein_top_cell_types": ";".join(protein_cell_names),
            "protein_cell_type_count": len(protein_cell_names),
            "hpa_subcellular_main_locations": ";".join(hpa_main),
            "hpa_subcellular_additional_locations": ";".join(hpa_additional),
            "hpa_subcellular_reliability": ";".join(reliabilities),
            "go_cellular_component_ids": ";".join(go_ids),
            "plasma_membrane_localization_flag": "1"
            if any(
                "plasma membrane" in row["location_name"].lower()
                or row["go_id"] == "GO:0005886"
                for row in locations
            )
            else "0",
            "best_location_evidence_level": best_loc,
            "has_rna_expression_evidence": "1" if rna_tissues or rna_cells else "0",
            "has_protein_expression_evidence": "1" if protein_tissues or protein_cells else "0",
            "expression_location_sources": ";".join(
                source
                for source, present in (
                    ("Human Protein Atlas", bool(hpa_records)),
                    ("UniProtKB/GO", bool(source_go)),
                )
                if present
            ),
            "expression_protein_discordance_flag": "1" if discordance else "0",
            "localization_source_divergence_flag": "1" if divergence else "0",
            "review_issue_count": 0,
        }

    deduplicated_reviews = []
    seen_review_keys = set()
    for row in review_rows:
        key = (
            row["issue_type"],
            row["severity"],
            row["target_uniprot_id"],
            row["source_record"],
            row["details"],
        )
        if key not in seen_review_keys:
            seen_review_keys.add(key)
            deduplicated_reviews.append(row)
    review_rows = deduplicated_reviews

    review_counts = defaultdict(int)
    for row in review_rows:
        for accession in row["target_uniprot_id"].split(";"):
            if accession:
                review_counts[accession] += 1
    for accession, count in review_counts.items():
        if accession in summary_by_accession:
            summary_by_accession[accession]["review_issue_count"] = count

    with gzip.open(SUMMARY_OUT, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=summary_fields, lineterminator="\n")
        writer.writeheader()
        for accession in sorted(summary_by_accession):
            writer.writerow(summary_by_accession[accession])

    review_fields = [
        "review_id",
        "issue_type",
        "severity",
        "target_uniprot_id",
        "source_record",
        "details",
        "module_version",
    ]
    with gzip.open(REVIEW_OUT, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=review_fields, lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(review_rows, 1):
            writer.writerow(
                {
                    "review_id": f"ELR{index:09d}",
                    **row,
                    "module_version": MODULE_VERSION,
                }
            )

    set_progress("build_v6_2_preview", started_utc=started)
    with MASTER.open("rt", encoding="utf-8-sig", newline="") as source, gzip.open(
        PREVIEW_OUT, "wt", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        preview_fields = [*reader.fieldnames, *summary_fields[1:]]
        writer = csv.DictWriter(target, delimiter="\t", fieldnames=preview_fields, lineterminator="\n")
        writer.writeheader()
        for row in reader:
            row.update({key: value for key, value in summary_by_accession[row["target_uniprot_id"]].items() if key != "target_uniprot_id"})
            writer.writerow(row)

    mapping_counts = {
        "tissue": defaultdict(int),
        "cell_type": defaultdict(int),
        "subcellular_location": defaultdict(int),
    }
    for label, items, group in (
        ("tissue", tissue_map, "tissue"),
        ("cell_type", cell_map, "cell_type"),
        ("subcellular_location", location_map, "subcellular_location"),
    ):
        for item in items.values():
            mapping_counts[group][item["mapping_status"]] += 1

    validation = {
        "status": "passed_with_review_queue",
        "started_utc": started,
        "completed_utc": now(),
        "module_version": MODULE_VERSION,
        "inputs": {
            "protein_master": str(MASTER),
            "protein_rows": len(master_rows),
            "protein_master_sha256": sha256(MASTER),
            "hpa_file": str(HPA_ZIP),
            "hpa_rows": hpa_total,
            "hpa_sha256": sha256(HPA_ZIP),
            "ontology_versions": {
                key: value["version"] for key, value in ontology_data.items()
            },
            "ontology_hashes": {
                key: value["sha256"] for key, value in ontology_data.items()
            },
        },
        "identifier_mapping": {
            "mapped_mempro_proteins": len(hpa_by_accession),
            "unmapped_mempro_proteins": len(master_rows) - len(hpa_by_accession),
            "hpa_identifier_review_rows": len(identifier_reviews),
        },
        "outputs": {
            "tissue_rows": tissue_count,
            "cell_type_rows": cell_count,
            "localization_rows": location_count,
            "summary_rows": len(summary_by_accession),
            "review_rows": len(review_rows),
        },
        "ontology_mapping_counts": {
            group: dict(counts) for group, counts in mapping_counts.items()
        },
        "quality_checks": {
            "summary_primary_key_unique": len(summary_by_accession) == len(master_rows),
            "summary_foreign_keys_in_master": all(
                accession in by_accession for accession in summary_by_accession
            ),
            "rna_and_protein_units_not_combined": True,
            "detailed_data_kept_out_of_master_preview": True,
            "formal_release_frozen": False,
            "preview_requires_v6_1_frozen_base_before_release": True,
        },
        "outputs_files": {
            "tissue": str(TISSUE_OUT),
            "cell_type": str(CELL_OUT),
            "localization": str(LOCATION_OUT),
            "mapping": str(MAPPING_OUT),
            "summary": str(SUMMARY_OUT),
            "review": str(REVIEW_OUT),
            "preview": str(PREVIEW_OUT),
        },
    }
    write_json(REPORT_OUT, validation)
    write_json(
        PROGRESS,
        {
            "status": "complete",
            "stage": "complete",
            "completed_utc": now(),
            "module_version": MODULE_VERSION,
            "validation_report": str(REPORT_OUT),
            "outputs": validation["outputs"],
        },
    )
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
