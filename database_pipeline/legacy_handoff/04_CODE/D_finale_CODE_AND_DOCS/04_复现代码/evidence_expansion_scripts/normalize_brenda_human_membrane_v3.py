#!/usr/bin/env python3
"""Stream BRENDA 2026.1 and stage human membrane-enzyme ligand evidence.

This is an incremental supplement.  The frozen protein/compound/evidence
releases are read-only.  BRENDA ligand names are never merged into a canonical
compound by name alone; identity mapping is a later, separately audited step.
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import pathlib
import re
import sqlite3
import sys
import tarfile
from collections import Counter, defaultdict
from typing import BinaryIO, Iterator


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "runs" / "incremental_20260727"
OUT_DIR = RUN_ROOT / "staging" / "brenda"
ARCHIVE = (
    ROOT
    / "raw"
    / "manual_downloads"
    / "brenda_2026_1"
    / "brenda_2026_1.json.tar.gz"
)
INDEX_DB = ROOT / "intermediate" / "integration_index_v2.sqlite"
SCHEMA_PATH = ROOT / "config" / "normalized_source_schema.tsv"

TARGET_FIELDS = {
    "natural_substrates_products": {
        "kind": "reaction",
        "tier": "BE3",
        "type": "natural_enzyme_reaction",
    },
    "substrates_products": {
        "kind": "reaction",
        "tier": "BE3",
        "type": "enzyme_reaction",
    },
    "cofactor": {
        "kind": "named",
        "tier": "BE3",
        "role": "cofactor",
        "type": "cofactor_requirement",
    },
    "activating_compound": {
        "kind": "named",
        "tier": "BE3",
        "role": "activator",
        "type": "activating_compound",
    },
    "inhibitor": {
        "kind": "named",
        "tier": "BE3",
        "role": "inhibitor",
        "type": "inhibitor",
    },
    "metals_ions": {
        "kind": "named",
        "tier": "BE3",
        "role": "metal_or_ion",
        "type": "metal_or_ion_requirement",
    },
    "ki_value": {
        "kind": "kinetic",
        "tier": "BE2",
        "role": "inhibitor",
        "type": "Ki",
        "unit": "mM",
    },
    "ic50_value": {
        "kind": "kinetic",
        "tier": "BE2",
        "role": "inhibitor",
        "type": "IC50",
        "unit": "mM",
    },
    "km_value": {
        "kind": "kinetic",
        "tier": "BE3",
        "role": "substrate",
        "type": "Km",
        "unit": "mM",
    },
    "turnover_number": {
        "kind": "kinetic",
        "tier": "BE3",
        "role": "substrate",
        "type": "kcat",
        "unit": "1/s",
    },
    "kcat_km_value": {
        "kind": "kinetic",
        "tier": "BE3",
        "role": "substrate",
        "type": "kcat/Km",
        "unit": "mM/s",
    },
}

EXTRA_FIELDS = [
    "brenda_ec_number",
    "brenda_protein_id",
    "brenda_role",
    "brenda_raw_value",
    "brenda_comment",
    "protein_mapping_method",
    "compound_name_normalized",
]

COMMON_SPECIES = {
    "?",
    "H+",
    "H2O",
    "D2O",
    "O2",
    "CO2",
    "NH3",
    "NH4+",
    "e-",
    "electron",
    "reduced electron acceptor",
    "oxidized electron acceptor",
}
GENERIC_NAME_PATTERNS = [
    re.compile(pattern, re.I)
    for pattern in [
        r"^an? .+",
        r"^several .+",
        r"^various .+",
        r"^other .+",
        r"^.+ protein$",
        r"^.+ peptide$",
        r"^.+ oligosaccharide",
        r"^poly[- ]",
        r"^\?$",
    ]
]


def stable_id(prefix: str, *parts: object) -> str:
    material = "\x1f".join("" if part is None else str(part) for part in parts)
    return prefix + hashlib.sha256(material.encode("utf-8")).hexdigest()[:24].upper()


def normalize_name(name: str) -> str:
    value = name.replace("\u00ad", "").strip()
    value = re.sub(r"\s+", " ", value)
    value = re.sub(r"^\d+\s+(?=[A-Za-z(\[])", "", value)
    return value


def is_generic_or_common(name: str) -> tuple[bool, str]:
    if name in COMMON_SPECIES:
        return True, "common reaction species"
    if not name or len(name) > 500:
        return True, "empty or implausibly long compound name"
    if any(pattern.search(name) for pattern in GENERIC_NAME_PATTERNS):
        return True, "generic or macromolecular ligand description"
    return False, ""


class IncrementalJSONObjectReader:
    """Iterate key/value pairs from the top-level BRENDA `data` object."""

    def __init__(self, binary: BinaryIO, chunk_size: int = 1024 * 1024):
        self.text = io.TextIOWrapper(binary, encoding="utf-8")
        self.chunk_size = chunk_size
        self.buffer = ""
        self.pos = 0
        self.eof = False
        self.decoder = json.JSONDecoder()

    def _fill(self) -> bool:
        if self.eof:
            return False
        chunk = self.text.read(self.chunk_size)
        if chunk:
            self.buffer += chunk
            return True
        self.eof = True
        return False

    def _compact(self) -> None:
        if self.pos > self.chunk_size:
            self.buffer = self.buffer[self.pos :]
            self.pos = 0

    def _ensure_data_start(self) -> None:
        marker = '"data"'
        while True:
            index = self.buffer.find(marker, self.pos)
            if index >= 0:
                colon = self.buffer.find(":", index + len(marker))
                brace = self.buffer.find("{", colon + 1 if colon >= 0 else index)
                if colon >= 0 and brace >= 0:
                    self.pos = brace + 1
                    return
            if not self._fill():
                raise ValueError("BRENDA JSON does not contain a data object")
            if len(self.buffer) > 4 * self.chunk_size:
                self.buffer = self.buffer[-2 * self.chunk_size :]
                self.pos = 0

    def _skip_ws_and_commas(self) -> None:
        while True:
            while self.pos < len(self.buffer) and (
                self.buffer[self.pos].isspace() or self.buffer[self.pos] == ","
            ):
                self.pos += 1
            if self.pos < len(self.buffer) or not self._fill():
                return

    def _decode_at_pos(self):
        while True:
            try:
                return self.decoder.raw_decode(self.buffer, self.pos)
            except json.JSONDecodeError:
                if not self._fill():
                    raise

    def __iter__(self) -> Iterator[tuple[str, dict[str, object]]]:
        self._ensure_data_start()
        while True:
            self._skip_ws_and_commas()
            if self.pos >= len(self.buffer):
                return
            if self.buffer[self.pos] == "}":
                return
            key, end = self._decode_at_pos()
            if not isinstance(key, str):
                raise ValueError(f"Expected string data key at offset {self.pos}")
            self.pos = end
            while True:
                while self.pos < len(self.buffer) and self.buffer[self.pos].isspace():
                    self.pos += 1
                if self.pos < len(self.buffer):
                    break
                if not self._fill():
                    raise EOFError("Unexpected EOF after BRENDA data key")
            if self.buffer[self.pos] != ":":
                raise ValueError(f"Expected colon after BRENDA data key {key}")
            self.pos += 1
            while True:
                while self.pos < len(self.buffer) and self.buffer[self.pos].isspace():
                    self.pos += 1
                if self.pos < len(self.buffer):
                    break
                if not self._fill():
                    raise EOFError("Unexpected EOF before BRENDA data value")
            value, end = self._decode_at_pos()
            self.pos = end
            self._compact()
            if isinstance(value, dict):
                yield key, value


def load_schema() -> list[str]:
    with SCHEMA_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        return [row["column_name"] for row in csv.DictReader(handle, delimiter="\t")]


def load_protein_maps():
    con = sqlite3.connect(INDEX_DB)
    proteins = {
        row[0]: {
            "target_uniprot_id": row[0],
            "approved_symbol": row[1] or "",
        }
        for row in con.execute("SELECT target_uniprot_id, approved_symbol FROM protein")
    }
    aliases: dict[str, set[str]] = defaultdict(set)
    for alias, target in con.execute(
        "SELECT alias, target_uniprot_id FROM protein_alias"
    ):
        aliases[str(alias).upper()].add(str(target))
    con.close()
    return proteins, aliases


def resolve_protein_record(
    record: dict[str, object],
    proteins: dict[str, dict[str, str]],
    aliases: dict[str, set[str]],
) -> tuple[list[str], str, str]:
    if str(record.get("organism") or "").strip() != "Homo sapiens":
        return [], "nonhuman", "organism is not Homo sapiens"
    accessions = [
        str(value).strip()
        for value in (record.get("accessions") or [])
        if str(value).strip()
    ]
    targets: set[str] = set()
    methods: set[str] = set()
    for accession in accessions:
        base = accession.split("-")[0].upper()
        if accession.upper() in proteins:
            targets.add(accession.upper())
            methods.add("exact_uniprot")
        elif base in proteins:
            targets.add(base)
            methods.add("canonicalized_isoform_accession")
        else:
            resolved = aliases.get(accession.upper(), set()) | aliases.get(base, set())
            targets.update(resolved)
            if resolved:
                methods.add("existing_protein_alias")
    if len(targets) == 1:
        return sorted(targets), "+".join(sorted(methods)), ""
    if not accessions:
        return [], "unmapped", "human BRENDA protein has no accession"
    if not targets:
        return [], "unmapped", "accession not found in frozen ABC protein universe"
    return sorted(targets), "ambiguous", "accessions resolve to multiple membrane proteins"


def protein_ids_for_item(item: dict[str, object]) -> list[str]:
    return [str(value) for value in (item.get("proteins") or [])]


def clean_reaction_markup(value: str) -> str:
    return re.sub(r"\s*\{(?:r|ir)\}\s*$", "", value, flags=re.I).strip()


def split_reaction(value: str) -> list[tuple[str, str]]:
    cleaned = clean_reaction_markup(value)
    if " = " not in cleaned:
        return []
    left, right = cleaned.split(" = ", 1)
    result = []
    for name in left.split(" + "):
        result.append(("substrate", normalize_name(name)))
    for name in right.split(" + "):
        result.append(("product", normalize_name(name)))
    return [(role, name) for role, name in result if name]


KINETIC_RE = re.compile(
    r"^\s*(?P<relation><=|>=|<|>|~|ca\.?)?\s*"
    r"(?P<value>[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?\d+)?)"
    r"\s*\{(?P<compound>.+)\}\s*$"
)


def parse_kinetic(value: str):
    match = KINETIC_RE.match(value)
    if not match:
        return None
    relation_map = {
        None: "=",
        "<": "<",
        ">": ">",
        "<=": "<=",
        ">=": ">=",
        "~": "~",
        "ca.": "~",
        "ca": "~",
    }
    return (
        relation_map.get(match.group("relation"), match.group("relation") or "="),
        match.group("value"),
        normalize_name(match.group("compound")),
    )


def item_compounds(
    field: str, item: dict[str, object], policy: dict[str, str]
) -> list[dict[str, str]]:
    raw = str(item.get("value") or "").strip()
    if policy["kind"] == "reaction":
        return [
            {
                "role": role,
                "compound": compound,
                "activity_relation": "",
                "activity_value": "",
                "activity_unit": "",
                "standard_value_nM": "",
            }
            for role, compound in split_reaction(raw)
        ]
    if policy["kind"] == "named":
        return [
            {
                "role": policy["role"],
                "compound": normalize_name(raw),
                "activity_relation": "",
                "activity_value": "",
                "activity_unit": "",
                "standard_value_nM": "",
            }
        ]
    parsed = parse_kinetic(raw)
    if not parsed:
        return []
    relation, value, compound = parsed
    standard_value = ""
    if policy["unit"] == "mM":
        try:
            standard_value = str(float(value) * 1_000_000.0)
        except ValueError:
            standard_value = ""
    return [
        {
            "role": policy["role"],
            "compound": compound,
            "activity_relation": relation,
            "activity_value": value,
            "activity_unit": policy["unit"],
            "standard_value_nM": standard_value,
        }
    ]


def atomic_gzip_writer(path: pathlib.Path, fields: list[str]):
    temporary = path.with_suffix(path.suffix + ".tmp")
    handle = gzip.open(temporary, "wt", encoding="utf-8", newline="")
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=fields,
        extrasaction="ignore",
        lineterminator="\n",
    )
    writer.writeheader()
    return temporary, handle, writer


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = load_schema() + EXTRA_FIELDS
    proteins, aliases = load_protein_maps()

    evidence_path = OUT_DIR / "brenda_normalized_evidence_v3.tsv.gz"
    unresolved_protein_path = OUT_DIR / "brenda_unresolved_human_proteins_v3.tsv.gz"
    compound_inventory_path = OUT_DIR / "brenda_compound_name_inventory_v3.tsv.gz"
    ev_tmp, ev_handle, ev_writer = atomic_gzip_writer(evidence_path, fields)
    up_fields = [
        "brenda_ec_number",
        "brenda_protein_id",
        "organism",
        "accessions",
        "references",
        "mapping_status",
        "mapping_reason",
    ]
    up_tmp, up_handle, up_writer = atomic_gzip_writer(
        unresolved_protein_path, up_fields
    )

    counts: Counter[str] = Counter()
    targets: set[str] = set()
    ecs: set[str] = set()
    target_compound_names: set[tuple[str, str]] = set()
    compound_stats: dict[str, dict[str, object]] = {}
    parser_error = ""

    try:
        with tarfile.open(ARCHIVE, "r:gz") as archive:
            members = [member for member in archive.getmembers() if member.isfile()]
            if len(members) != 1:
                raise ValueError("Expected one BRENDA JSON member in archive")
            binary = archive.extractfile(members[0])
            if binary is None:
                raise FileNotFoundError("Could not open BRENDA JSON archive member")
            reader = IncrementalJSONObjectReader(binary)
            for ec_number, entry in reader:
                counts["ec_entries_seen"] += 1
                if ec_number == "spontaneous":
                    counts["spontaneous_entries_skipped"] += 1
                    continue
                protein_records = entry.get("protein") or {}
                local_targets: dict[str, tuple[list[str], str, str]] = {}
                for protein_id, record in protein_records.items():
                    if not isinstance(record, dict):
                        continue
                    organism = str(record.get("organism") or "").strip()
                    if organism != "Homo sapiens":
                        continue
                    counts["human_protein_records"] += 1
                    resolution = resolve_protein_record(record, proteins, aliases)
                    local_targets[str(protein_id)] = resolution
                    mapped, method, reason = resolution
                    if mapped:
                        counts["human_protein_records_mapped"] += 1
                    else:
                        counts["human_protein_records_unresolved"] += 1
                        up_writer.writerow(
                            {
                                "brenda_ec_number": ec_number,
                                "brenda_protein_id": protein_id,
                                "organism": organism,
                                "accessions": ";".join(
                                    str(v) for v in (record.get("accessions") or [])
                                ),
                                "references": ";".join(
                                    str(v) for v in (record.get("references") or [])
                                ),
                                "mapping_status": method,
                                "mapping_reason": reason,
                            }
                        )
                mapped_local_targets = {
                    key: value for key, value in local_targets.items() if value[0]
                }
                if not mapped_local_targets:
                    continue
                counts["ec_entries_with_mapped_human_membrane_protein"] += 1
                ecs.add(ec_number)
                references = entry.get("reference") or {}

                for field, policy in TARGET_FIELDS.items():
                    items = entry.get(field) or []
                    for item_index, item in enumerate(items):
                        if not isinstance(item, dict):
                            continue
                        item_protein_ids = protein_ids_for_item(item)
                        selected = [
                            protein_id
                            for protein_id in item_protein_ids
                            if protein_id in mapped_local_targets
                        ]
                        if not selected:
                            counts[f"{field}_items_without_mapped_human_protein"] += 1
                            continue
                        compounds = item_compounds(field, item, policy)
                        if not compounds:
                            counts[f"{field}_items_unparsed"] += 1
                            continue
                        reference_ids = [
                            str(value) for value in (item.get("references") or [])
                        ]
                        pmids = sorted(
                            {
                                str(references[ref].get("pmid"))
                                for ref in reference_ids
                                if ref in references
                                and isinstance(references[ref], dict)
                                and references[ref].get("pmid")
                            }
                        )
                        comment = str(item.get("comment") or "").strip()
                        raw_value = str(item.get("value") or "").strip()

                        for protein_id in selected:
                            mapped_targets, mapping_method, _ = mapped_local_targets[
                                protein_id
                            ]
                            for target in mapped_targets:
                                target_row = proteins[target]
                                for compound_item in compounds:
                                    compound = compound_item["compound"]
                                    generic, generic_reason = is_generic_or_common(
                                        compound
                                    )
                                    compound_normalized = normalize_name(compound).casefold()
                                    source_record_id = (
                                        f"{ec_number}:{field}:{item_index}:"
                                        f"{protein_id}:{target}:{compound_normalized}"
                                    )
                                    evidence_id = stable_id(
                                        "BRENDA-", source_record_id, reference_ids
                                    )
                                    exclusion_reasons = [
                                        "compound identity pending; name alone is not merged"
                                    ]
                                    if generic:
                                        exclusion_reasons.append(generic_reason)
                                    row = {
                                        "source_evidence_id": evidence_id,
                                        "source_database": "BRENDA",
                                        "source_version": "2026.1",
                                        "source_record_id": source_record_id,
                                        "source_url": (
                                            "https://www.brenda-enzymes.org/"
                                            f"enzyme.php?ecno={ec_number}"
                                        ),
                                        "retrieval_date": "2026-07-27",
                                        "target_uniprot_id": target,
                                        "target_mapping_status": "single_human_uniprot",
                                        "approved_symbol": target_row[
                                            "approved_symbol"
                                        ],
                                        "compound_source_id": stable_id(
                                            "BRENDA-NAME-", compound_normalized
                                        ),
                                        "compound_name": compound,
                                        "compound_internal_id": "",
                                        "compound_form_id": "",
                                        "compound_mapping_status": "unresolved",
                                        "evidence_tier": policy["tier"],
                                        "evidence_type": policy["type"],
                                        "activity_type": (
                                            policy.get("type")
                                            if policy["kind"] == "kinetic"
                                            else compound_item["role"]
                                        ),
                                        "activity_relation": compound_item[
                                            "activity_relation"
                                        ],
                                        "activity_value": compound_item[
                                            "activity_value"
                                        ],
                                        "activity_unit": compound_item["activity_unit"],
                                        "standard_value_nM": compound_item[
                                            "standard_value_nM"
                                        ],
                                        "activity_outcome": "observed",
                                        "assay_or_mechanism": comment,
                                        "pdb_ids": "",
                                        "binding_site_residues": "",
                                        "pubmed_ids": ";".join(pmids),
                                        "doi": "",
                                        "record_qc_status": (
                                            "review" if generic else "ok"
                                        ),
                                        "default_release_inclusion": "0",
                                        "exclusion_reason": "; ".join(
                                            exclusion_reasons
                                        ),
                                        "brenda_ec_number": ec_number,
                                        "brenda_protein_id": protein_id,
                                        "brenda_role": compound_item["role"],
                                        "brenda_raw_value": raw_value,
                                        "brenda_comment": comment,
                                        "protein_mapping_method": mapping_method,
                                        "compound_name_normalized": compound_normalized,
                                    }
                                    ev_writer.writerow(row)
                                    counts["evidence_rows"] += 1
                                    counts[f"evidence_tier_{policy['tier']}"] += 1
                                    counts[f"evidence_field_{field}"] += 1
                                    counts[
                                        "generic_or_common_compound_rows"
                                        if generic
                                        else "specific_compound_name_rows"
                                    ] += 1
                                    targets.add(target)
                                    target_compound_names.add(
                                        (target, compound_normalized)
                                    )
                                    stat = compound_stats.setdefault(
                                        compound_normalized,
                                        {
                                            "preferred_source_name": compound,
                                            "roles": set(),
                                            "fields": set(),
                                            "targets": set(),
                                            "evidence_rows": 0,
                                            "generic_or_common": generic,
                                            "generic_reason": generic_reason,
                                        },
                                    )
                                    stat["roles"].add(compound_item["role"])
                                    stat["fields"].add(field)
                                    stat["targets"].add(target)
                                    stat["evidence_rows"] += 1
    except Exception as exc:
        parser_error = repr(exc)
        raise
    finally:
        ev_handle.close()
        up_handle.close()

    ev_tmp.replace(evidence_path)
    up_tmp.replace(unresolved_protein_path)

    ci_fields = [
        "compound_name_normalized",
        "preferred_source_name",
        "roles",
        "brenda_fields",
        "target_count",
        "evidence_row_count",
        "generic_or_common",
        "generic_reason",
        "mapping_status",
    ]
    ci_tmp, ci_handle, ci_writer = atomic_gzip_writer(
        compound_inventory_path, ci_fields
    )
    try:
        for key, stat in sorted(compound_stats.items()):
            ci_writer.writerow(
                {
                    "compound_name_normalized": key,
                    "preferred_source_name": stat["preferred_source_name"],
                    "roles": ";".join(sorted(stat["roles"])),
                    "brenda_fields": ";".join(sorted(stat["fields"])),
                    "target_count": len(stat["targets"]),
                    "evidence_row_count": stat["evidence_rows"],
                    "generic_or_common": int(bool(stat["generic_or_common"])),
                    "generic_reason": stat["generic_reason"],
                    "mapping_status": "pending_structure_identity_mapping",
                }
            )
    finally:
        ci_handle.close()
    ci_tmp.replace(compound_inventory_path)

    report = {
        "status": "passed" if not parser_error else "failed",
        "run_type": "incremental_supplement_not_rebuild",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_archive": str(ARCHIVE),
        "source_version": "2026.1",
        "license": "CC BY 4.0",
        "baseline_policy": "read_only",
        "compound_policy": (
            "names alone are never merged; all BRENDA names remain pending "
            "until structure/identifier identity mapping"
        ),
        "counts": dict(sorted(counts.items())),
        "unique_mapped_targets": len(targets),
        "unique_ec_numbers": len(ecs),
        "unique_compound_names": len(compound_stats),
        "unique_target_compound_name_pairs": len(target_compound_names),
        "parser_error": parser_error,
        "outputs": {
            "normalized_evidence": str(evidence_path),
            "unresolved_human_proteins": str(unresolved_protein_path),
            "compound_name_inventory": str(compound_inventory_path),
        },
        "unit_policy": {
            "Ki": "mM per BRENDA datafield documentation",
            "IC50": "mM per BRENDA datafield documentation",
            "Km": "mM per BRENDA datafield documentation",
            "kcat": "1/s per BRENDA datafield documentation",
            "kcat/Km": "mM/s as labelled by BRENDA datafield documentation",
        },
    }
    qa_path = RUN_ROOT / "qa" / "BRENDA_NORMALIZATION_V3_QA.json"
    temporary = qa_path.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    temporary.replace(qa_path)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
