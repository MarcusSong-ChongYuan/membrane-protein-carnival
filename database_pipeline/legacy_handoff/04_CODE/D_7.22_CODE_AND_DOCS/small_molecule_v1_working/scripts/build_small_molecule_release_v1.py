from __future__ import annotations

import csv
import hashlib
import itertools
import json
import re
import sqlite3
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path


ROOT = Path(r"D:\7.22")
INPUT_RELEASE = (
    ROOT / "membrane_master_v5_working" / "releases" / "release_v5_4_pgd_binding"
)
WORK = ROOT / "small_molecule_v1_working"
CHUNKS = WORK / "intermediate" / "standardized_chunks"
RAW = WORK / "raw"
REPORTS = WORK / "reports"
RELEASE = WORK / "releases" / "release_small_molecule_v1_0"
RELEASE.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

DB_PATH = WORK / "intermediate" / "small_molecule_v1.sqlite"
CHEMBL_METADATA = RAW / "chembl37" / "chembl37_selected_metadata.jsonl"
CHEMBL_STRUCTURES = RAW / "chembl37" / "chembl37_selected_chemreps.tsv"
CHEBI_MAPPING = RAW / "chebi" / "chebi_202607_selected_mapping.tsv"

INDEX_INPUT = INPUT_RELEASE / "small_molecule_index_v4_1.tsv"
EVIDENCE_INPUT = INPUT_RELEASE / "binding_evidence_master_v4_1.tsv"
SITE_INPUT = INPUT_RELEASE / "binding_site_instances_v4_1.tsv"
PROTEIN_MASTER_INPUT = INPUT_RELEASE / "human_membrane_audit_master_v5_4.tsv"

PARENT_MASTER_OUT = RELEASE / "small_molecule_master_v1_0.tsv"
FORM_OUT = RELEASE / "compound_form_hierarchy_v1_0.tsv"
XREF_OUT = RELEASE / "compound_source_xref_v1_0.tsv"
SYNONYM_OUT = RELEASE / "compound_synonyms_v1_0.tsv"
UNRESOLVED_OUT = RELEASE / "compound_unresolved_candidates_v1_0.tsv"
EXCLUSION_OUT = RELEASE / "compound_exclusion_audit_v1_0.tsv"
MERGE_AUDIT_OUT = RELEASE / "compound_merge_audit_v1_0.tsv"
CONNECTIVITY_OUT = RELEASE / "compound_connectivity_families_v1_0.tsv"
EVIDENCE_OUT = RELEASE / "binding_evidence_master_v4_2.tsv"
SITE_OUT = RELEASE / "binding_site_instances_v4_2.tsv"
PAIR_OUT = RELEASE / "protein_compound_summary_v4_2.tsv"
COMPOUND_PROTEIN_OUT = RELEASE / "compound_protein_summary_v1_0.tsv"
PROTEIN_BINDING_OUT = RELEASE / "protein_binding_summary_v4_2.tsv"
PROTEIN_MASTER_V55_OUT = RELEASE / "human_membrane_audit_master_v5_5.tsv"
MASTER_REVIEW_OUT = RELEASE / "small_molecule_excel_review_v1_0.tsv"
FORM_REVIEW_OUT = RELEASE / "compound_form_excel_review_v1_0.tsv"
PAIR_REVIEW_OUT = RELEASE / "compound_protein_excel_review_v1_0.tsv"
AUDIT_REVIEW_OUT = RELEASE / "compound_audit_excel_review_v1_0.tsv"
VALIDATION_OUT = RELEASE / "SMALL_MOLECULE_V1_VALIDATION_REPORT.json"

csv.field_size_limit(100_000_000)
CHEMBL_RE = re.compile(r"CHEMBL\d+", re.IGNORECASE)
ID_LIKE_RE = re.compile(
    r"^(?:CHEMBL[:_]?\d+|CID[:_]?\d+|PUBCHEM[:_]?\d+|\d+|UNMAPPED[_:]|UNKNOWN$)",
    re.IGNORECASE,
)
BIOLOGIC_TYPES = {
    "protein",
    "antibody",
    "antibody drug conjugate",
    "oligonucleotide",
    "gene",
    "cell",
    "enzyme",
}
EXTENDED_TYPES = {"oligosaccharide"}
EXTENDED_STRUCTURAL_CLASSES = {
    "peptidic_or_peptidomimetic",
    "inorganic_chemical_entity",
    "organometallic_or_metal_containing",
}

SOURCE_FIELDS = [
    "source_record_id",
    "input_row_number",
    "source_compound_id",
    "source_compound_id_type",
    "source_compound_name",
    "source_scope_status",
    "source_scope_class",
    "source_is_core_small_molecule",
    "source_databases",
    "source_versions",
    "source_record_qc_status",
    "pubchem_cids",
    "chembl_ids",
    "structure_input_source",
    "input_smiles",
    "input_inchikey",
    "structure_parse_status",
    "exact_standard_smiles",
    "exact_standard_inchi",
    "exact_standard_inchikey",
    "exact_identity_key",
    "parent_standard_smiles",
    "parent_standard_inchi",
    "parent_standard_inchikey",
    "parent_identity_key",
    "exact_connectivity_key",
    "parent_connectivity_key",
    "fragment_count",
    "form_type",
    "exact_molecular_formula",
    "exact_molecular_weight",
    "exact_formal_charge",
    "parent_molecular_formula",
    "parent_molecular_weight",
    "parent_exact_mass",
    "parent_xlogp",
    "parent_tpsa",
    "parent_hbond_donor_count",
    "parent_hbond_acceptor_count",
    "parent_rotatable_bond_count",
    "parent_formal_charge",
    "parent_heavy_atom_count",
    "parent_ring_count",
    "parent_max_ring_size",
    "parent_amide_bond_count",
    "computed_structural_class",
    "contains_metal",
    "standardization_actions",
    "structure_qc_note",
    "primary_chembl_id",
    "chembl_pref_name",
    "chembl_molecule_type",
    "chembl_max_phase",
    "chembl_first_approval",
    "chembl_parent_id",
    "chembl_parent_structure_inchikey",
    "chembl_parent_alignment_status",
    "chembl_natural_product",
    "chembl_chemical_probe",
    "final_scope_status",
    "final_scope_reason",
]

PARENT_FIELDS = [
    "compound_internal_id",
    "preferred_name",
    "preferred_name_source",
    "compound_scope_status",
    "compound_scope_class",
    "identity_confidence",
    "standard_smiles",
    "standard_inchi",
    "standard_inchikey",
    "identity_key",
    "connectivity_key",
    "molecular_formula",
    "molecular_weight",
    "exact_mass",
    "xlogp",
    "tpsa",
    "hbond_donor_count",
    "hbond_acceptor_count",
    "rotatable_bond_count",
    "formal_charge",
    "heavy_atom_count",
    "ring_count",
    "max_ring_size",
    "amide_bond_count",
    "computed_structural_class",
    "molecule_types",
    "form_count",
    "source_record_count",
    "source_database_count",
    "source_databases",
    "pubchem_cids",
    "chembl_ids",
    "chebi_ids",
    "drugcentral_ids",
    "gtopdb_ligand_ids",
    "hmdb_ids",
    "cas_numbers",
    "synonym_count",
    "development_status",
    "max_chembl_phase",
    "first_approval_year",
    "is_approved_drug",
    "is_clinical_candidate",
    "is_endogenous_ligand",
    "is_natural_product",
    "is_chemical_probe",
    "compound_category_tags",
    "chebi_direct_classes",
    "chebi_roles",
    "record_qc_status",
    "qc_notes",
]


def split_values(value: str) -> list[str]:
    return [
        item.strip()
        for item in re.split(r"[;|]", value or "")
        if item.strip()
    ]


def numeric_or_blank(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def stable_id(prefix: str, number: int) -> str:
    return f"{prefix}-{number:07d}"


def hash_id(prefix: str, *parts: str) -> str:
    digest = hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:20].upper()
    return f"{prefix}-{digest}"


def chembl_ids_from(*values: str) -> list[str]:
    result = set()
    for value in values:
        result.update(match.upper() for match in CHEMBL_RE.findall(value or ""))
    return sorted(result, key=lambda item: int(item[6:]))


def load_jsonl_by_id(path: Path, id_field: str) -> dict[str, dict]:
    result = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                row = json.loads(line)
                result[str(row[id_field])] = row
    return result


def load_chembl_structures() -> dict[str, dict]:
    result = {}
    with CHEMBL_STRUCTURES.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["chembl_id"]] = row
    return result


def load_chebi() -> dict[str, list[dict]]:
    result: dict[str, list[dict]] = defaultdict(list)
    with CHEBI_MAPPING.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["inchikey"]].append(row)
    return result


def final_scope(row: dict, meta: dict) -> tuple[str, str]:
    source_scope = row.get("source_scope_status", "")
    parse_status = row.get("structure_parse_status", "")
    molecule_type = str(meta.get("molecule_type") or "").strip()
    molecule_type_lower = molecule_type.casefold()
    structural_class = row.get("computed_structural_class", "")

    if source_scope == "excluded":
        return "excluded", "excluded by the audited v4.1 scope decision"
    if molecule_type_lower in BIOLOGIC_TYPES:
        return "excluded", f"ChEMBL molecule type is {molecule_type}"
    if parse_status != "parsed" or not row.get("parent_identity_key"):
        return "unresolved", "no validated standard chemical structure"
    if (
        source_scope in {"separate_not_core", "review"}
        or molecule_type_lower in EXTENDED_TYPES
        or structural_class in EXTENDED_STRUCTURAL_CLASSES
    ):
        reasons = []
        if source_scope in {"separate_not_core", "review"}:
            reasons.append(f"v4.1 scope={source_scope}")
        if molecule_type_lower in EXTENDED_TYPES:
            reasons.append(f"ChEMBL molecule type={molecule_type}")
        if structural_class in EXTENDED_STRUCTURAL_CLASSES:
            reasons.append(f"computed class={structural_class}")
        return "extended", "; ".join(reasons)
    return "core", "validated defined non-polymeric chemical structure"


def refine_structural_class(row: dict) -> None:
    if row.get("computed_structural_class") == "polycyclic_or_steroid_like":
        row["computed_structural_class"] = "polycyclic_organic_compound"
        return
    if row.get("computed_structural_class") != "lipid_like":
        return
    ring_count = numeric_or_blank(row.get("parent_ring_count", "")) or 0
    rotatable = numeric_or_blank(row.get("parent_rotatable_bond_count", "")) or 0
    if ring_count >= 2:
        row["computed_structural_class"] = "polycyclic_organic_compound"
    elif rotatable < 8:
        row["computed_structural_class"] = "other_defined_organic_compound"


def create_database() -> sqlite3.Connection:
    if DB_PATH.exists():
        DB_PATH.unlink()
    connection = sqlite3.connect(DB_PATH)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA synchronous=NORMAL")
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.execute("PRAGMA cache_size=-500000")
    columns = ", ".join(f'"{field}" TEXT' for field in SOURCE_FIELDS)
    connection.execute(f"CREATE TABLE source_curated ({columns})")
    connection.execute(
        "CREATE UNIQUE INDEX idx_source_record_id ON source_curated(source_record_id)"
    )
    connection.execute(
        "CREATE UNIQUE INDEX idx_source_compound_id ON source_curated(source_compound_id)"
    )
    connection.execute(
        "CREATE INDEX idx_source_parent_key ON source_curated(parent_identity_key)"
    )
    connection.execute(
        "CREATE INDEX idx_source_exact_key ON source_curated(exact_identity_key)"
    )
    return connection


def ingest_source_records(
    connection: sqlite3.Connection,
    metadata: dict[str, dict],
    structures: dict[str, dict],
) -> Counter:
    placeholders = ",".join("?" for _ in SOURCE_FIELDS)
    insert_sql = (
        f"INSERT INTO source_curated ({','.join(SOURCE_FIELDS)}) VALUES ({placeholders})"
    )
    counts = Counter()
    batch = []
    for path in sorted(CHUNKS.glob("standardized_chunk_*.tsv")):
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                refine_structural_class(row)
                ids = chembl_ids_from(
                    row.get("source_compound_id", ""),
                    row.get("chembl_ids", ""),
                )
                primary = ids[0] if ids else ""
                meta = metadata.get(primary, {})
                parent_id = str(meta.get("parent_chembl_id") or "")
                parent_key = (
                    structures.get(parent_id, {}).get("standard_inchi_key", "")
                    if parent_id
                    else ""
                )
                computed_parent_key = row.get("parent_standard_inchikey", "")
                if not parent_id:
                    alignment = "not_applicable"
                elif not parent_key:
                    alignment = "chembl_parent_structure_unavailable"
                elif parent_key == computed_parent_key:
                    alignment = "aligned"
                else:
                    alignment = "conflict"
                scope, reason = final_scope(row, meta)
                row.update(
                    {
                        "primary_chembl_id": primary,
                        "chembl_pref_name": str(meta.get("pref_name") or ""),
                        "chembl_molecule_type": str(meta.get("molecule_type") or ""),
                        "chembl_max_phase": str(meta.get("max_phase") or ""),
                        "chembl_first_approval": str(meta.get("first_approval") or ""),
                        "chembl_parent_id": parent_id,
                        "chembl_parent_structure_inchikey": parent_key,
                        "chembl_parent_alignment_status": alignment,
                        "chembl_natural_product": str(meta.get("natural_product") or ""),
                        "chembl_chemical_probe": str(meta.get("chemical_probe") or ""),
                        "final_scope_status": scope,
                        "final_scope_reason": reason,
                    }
                )
                counts[f"scope_{scope}"] += 1
                counts[f"parse_{row.get('structure_parse_status', '')}"] += 1
                if alignment == "conflict":
                    counts["chembl_parent_alignment_conflict"] += 1
                batch.append(tuple(row.get(field, "") for field in SOURCE_FIELDS))
                if len(batch) >= 5_000:
                    connection.executemany(insert_sql, batch)
                    connection.commit()
                    batch.clear()
    if batch:
        connection.executemany(insert_sql, batch)
        connection.commit()
    return counts


def build_group_maps(
    connection: sqlite3.Connection,
) -> tuple[dict[str, str], dict[str, str], dict[str, str], Counter]:
    parent_keys = [
        row[0]
        for row in connection.execute(
            """
            SELECT DISTINCT parent_identity_key
            FROM source_curated
            WHERE final_scope_status IN ('core','extended')
              AND parent_identity_key <> ''
            ORDER BY parent_identity_key
            """
        )
    ]
    parent_id_map = {
        key: stable_id("HMPD-CMPD", index)
        for index, key in enumerate(parent_keys, start=1)
    }
    connection.execute(
        "CREATE TABLE parent_groups (parent_identity_key TEXT PRIMARY KEY, compound_internal_id TEXT UNIQUE)"
    )
    connection.executemany(
        "INSERT INTO parent_groups VALUES (?,?)",
        parent_id_map.items(),
    )

    exact_parent_counts = connection.execute(
        """
        SELECT exact_identity_key, parent_identity_key, COUNT(*) AS n
        FROM source_curated
        WHERE final_scope_status IN ('core','extended')
          AND exact_identity_key <> ''
          AND parent_identity_key <> ''
        GROUP BY exact_identity_key, parent_identity_key
        ORDER BY exact_identity_key, n DESC, parent_identity_key
        """
    )
    selected_parent: dict[str, str] = {}
    conflict_counts = Counter()
    for exact_key, parent_key, count in exact_parent_counts:
        if exact_key not in selected_parent:
            selected_parent[exact_key] = parent_key
        elif selected_parent[exact_key] != parent_key:
            conflict_counts["exact_form_multiple_parent_groups"] += 1

    form_id_map = {
        key: stable_id("HMPD-FORM", index)
        for index, key in enumerate(sorted(selected_parent), start=1)
    }
    connection.execute(
        """
        CREATE TABLE form_groups (
            exact_identity_key TEXT PRIMARY KEY,
            compound_form_id TEXT UNIQUE,
            parent_identity_key TEXT,
            compound_internal_id TEXT
        )
        """
    )
    connection.executemany(
        "INSERT INTO form_groups VALUES (?,?,?,?)",
        (
            (
                exact_key,
                form_id_map[exact_key],
                selected_parent[exact_key],
                parent_id_map[selected_parent[exact_key]],
            )
            for exact_key in sorted(selected_parent)
        ),
    )
    connection.commit()
    return parent_id_map, form_id_map, selected_parent, conflict_counts


def source_priority(source_type: str, source_databases: set[str]) -> int:
    if "DrugCentral" in source_databases:
        return 2
    if "IUPHAR/BPS Guide to PHARMACOLOGY" in source_databases:
        return 3
    if source_type == "PubChem CID":
        return 4
    if source_type in {"ChEMBL", "ChEMBL ID"}:
        return 5
    return 6


def name_candidate(
    candidates: list[tuple[float, str, str]],
    priority: float,
    name: str,
    source: str,
) -> None:
    name = (name or "").strip()
    if not name:
        return
    candidates.append((priority, name, source))


def choose_name(candidates: list[tuple[float, str, str]], fallback: str) -> tuple[str, str]:
    usable = [
        item
        for item in candidates
        if not ID_LIKE_RE.match(item[1]) and len(item[1]) <= 300
    ]
    pool = usable or candidates
    if not pool:
        return fallback, "internal identifier"
    priority, name, source = min(pool, key=lambda item: (item[0], len(item[1]), item[1]))
    return name, source


def serialize(values: set[str]) -> str:
    return "|".join(sorted(value for value in values if value))


def build_parent_master(
    connection: sqlite3.Connection,
    parent_id_map: dict[str, str],
    metadata: dict[str, dict],
    chebi_by_key: dict[str, list[dict]],
) -> tuple[Counter, dict[str, int]]:
    connection.execute(
        f"CREATE TABLE parent_master ({', '.join(f'{field} TEXT' for field in PARENT_FIELDS)})"
    )
    insert_sql = (
        f"INSERT INTO parent_master ({','.join(PARENT_FIELDS)}) "
        f"VALUES ({','.join('?' for _ in PARENT_FIELDS)})"
    )
    form_count_by_parent = {
        row[0]: row[1]
        for row in connection.execute(
            "SELECT parent_identity_key, COUNT(*) FROM form_groups GROUP BY parent_identity_key"
        )
    }
    counts = Counter()
    synonym_writer_handle = SYNONYM_OUT.open("w", encoding="utf-8", newline="")
    synonym_writer = csv.DictWriter(
        synonym_writer_handle,
        fieldnames=[
            "synonym_id",
            "compound_internal_id",
            "synonym",
            "synonym_type",
            "source",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    synonym_writer.writeheader()

    query = connection.execute(
        """
        SELECT *
        FROM source_curated
        WHERE final_scope_status IN ('core','extended')
          AND parent_identity_key <> ''
        ORDER BY parent_identity_key
        """
    )
    batch = []
    for parent_key, iterator in itertools.groupby(
        query, key=lambda row: row[SOURCE_FIELDS.index("parent_identity_key")]
    ):
        rows = [dict(zip(SOURCE_FIELDS, row)) for row in iterator]
        compound_id = parent_id_map[parent_key]
        rep = min(
            rows,
            key=lambda row: (
                row.get("form_type") != "parent_form",
                source_priority(
                    row.get("source_compound_id_type", ""),
                    set(split_values(row.get("source_databases", ""))),
                ),
                row.get("source_record_id", ""),
            ),
        )
        source_databases = set()
        pubchem_ids = set()
        chembl_ids = set()
        drugcentral_ids = set()
        gtopdb_ids = set()
        molecule_types = set()
        candidates: list[tuple[float, str, str]] = []
        synonyms: dict[str, tuple[str, set[str]]] = {}
        max_phase = -1.0
        first_approvals = []
        natural_product = False
        chemical_probe = False
        core_present = False
        alignment_conflicts = 0

        def add_synonym(value: str, source: str) -> None:
            value = (value or "").strip()
            if not value or len(value) > 500:
                return
            key = value.casefold()
            if key not in synonyms:
                synonyms[key] = (value, {source})
            else:
                synonyms[key][1].add(source)

        for row in rows:
            core_present = core_present or row.get("final_scope_status") == "core"
            sources = set(split_values(row.get("source_databases", "")))
            source_databases.update(sources)
            source_id = row.get("source_compound_id", "")
            source_type = row.get("source_compound_id_type", "")
            if source_type == "PubChem CID" and source_id.isdigit():
                pubchem_ids.add(source_id)
            if source_type in {"ChEMBL", "ChEMBL ID"}:
                chembl_ids.update(chembl_ids_from(source_id))
            if source_type == "DrugCentral":
                drugcentral_ids.add(source_id)
            if source_type == "GtoPdb ligand ID":
                gtopdb_ids.add(source_id)
            pubchem_ids.update(split_values(row.get("pubchem_cids", "")))
            chembl_ids.update(chembl_ids_from(row.get("chembl_ids", "")))
            if row.get("primary_chembl_id"):
                chembl_ids.add(row["primary_chembl_id"])
            parent_chembl = row.get("chembl_parent_id", "")
            if parent_chembl:
                chembl_ids.add(parent_chembl)
            molecule_type = row.get("chembl_molecule_type", "")
            if molecule_type:
                molecule_types.add(molecule_type)
            phase = numeric_or_blank(row.get("chembl_max_phase", ""))
            if phase is not None:
                max_phase = max(max_phase, phase)
            approval = numeric_or_blank(row.get("chembl_first_approval", ""))
            if approval and approval > 0:
                first_approvals.append(int(approval))
            natural_product = natural_product or row.get("chembl_natural_product") in {
                "1",
                "True",
                "true",
            }
            chemical_probe = chemical_probe or row.get("chembl_chemical_probe") in {
                "1",
                "True",
                "true",
            }
            if row.get("chembl_parent_alignment_status") == "conflict":
                alignment_conflicts += 1

            source_name = row.get("source_compound_name", "")
            add_synonym(source_name, row.get("source_databases", "v4.1 source"))
            name_candidate(
                candidates,
                float(source_priority(source_type, sources)),
                source_name,
                row.get("source_databases", "v4.1 source"),
            )
            pref_name = row.get("chembl_pref_name", "")
            if pref_name:
                priority = 1.0 if (phase or -1) >= 1 else 3.0
                name_candidate(candidates, priority, pref_name, "ChEMBL 37")
                add_synonym(pref_name, "ChEMBL 37")
            if parent_chembl:
                parent_meta = metadata.get(parent_chembl, {})
                parent_name = str(parent_meta.get("pref_name") or "")
                if parent_name:
                    parent_phase = numeric_or_blank(parent_meta.get("max_phase"))
                    priority = 0.5 if (parent_phase or -1) >= 1 else 2.5
                    name_candidate(
                        candidates,
                        priority,
                        parent_name,
                        "ChEMBL 37 parent",
                    )
                    add_synonym(parent_name, "ChEMBL 37 parent")

        parent_inchikey = rep.get("parent_standard_inchikey", "")
        chebi_records = chebi_by_key.get(parent_inchikey, [])
        chebi_ids = set()
        hmdb_ids = set()
        cas_numbers = set()
        chebi_classes = set()
        chebi_roles = set()
        for chebi in chebi_records:
            chebi_ids.add(chebi.get("chebi_id", ""))
            hmdb_ids.update(split_values(chebi.get("hmdb_ids", "")))
            cas_numbers.update(split_values(chebi.get("cas_numbers", "")))
            pubchem_ids.update(split_values(chebi.get("pubchem_cids", "")))
            chembl_ids.update(chembl_ids_from(chebi.get("chembl_ids", "")))
            chebi_classes.update(split_values(chebi.get("direct_parent_names", "")))
            chebi_roles.update(split_values(chebi.get("direct_role_names", "")))
            star = numeric_or_blank(chebi.get("chebi_star", "")) or 0
            chebi_name = chebi.get("chebi_name", "")
            name_candidate(
                candidates,
                1.5 if star >= 3 else 3.5,
                chebi_name,
                f"ChEBI {chebi.get('chebi_id', '')}",
            )
            add_synonym(chebi_name, f"ChEBI {chebi.get('chebi_id', '')}")
            for synonym in split_values(chebi.get("synonyms", "")):
                add_synonym(synonym, f"ChEBI {chebi.get('chebi_id', '')}")

        preferred_name, preferred_source = choose_name(candidates, compound_id)
        source_count = len(source_databases)
        identity_confidence = (
            "high"
            if source_count >= 2
            or chebi_ids
            or (pubchem_ids and chembl_ids)
            else "medium"
        )
        roles_lower = " ".join(chebi_roles).casefold()
        endogenous = any(
            token in roles_lower
            for token in (
                "metabolite",
                "neurotransmitter",
                "hormone",
                "vitamin",
                "cofactor",
                "signalling molecule",
                "signaling molecule",
                "endogenous",
            )
        )
        approved = max_phase >= 4
        clinical = 0 < max_phase < 4
        if approved:
            development_status = "approved_drug"
        elif clinical:
            development_status = f"clinical_phase_{max_phase:g}"
        else:
            development_status = "research_or_preclinical"
        tags = set()
        if approved:
            tags.add("approved_drug")
        if clinical:
            tags.add("clinical_candidate")
        if endogenous:
            tags.add("endogenous_ligand")
        if natural_product:
            tags.add("natural_product")
        if chemical_probe:
            tags.add("chemical_probe")
        if not tags:
            tags.add("research_ligand")
        qc_notes = []
        distinct_parent_smiles = {
            row.get("parent_standard_smiles", "") for row in rows
        }
        if len(distinct_parent_smiles) > 1:
            qc_notes.append("same identity key has multiple standardized SMILES")
            counts["parent_inchikey_smiles_conflicts"] += 1
        if alignment_conflicts:
            qc_notes.append(
                f"{alignment_conflicts} ChEMBL parent assignment(s) conflict with structure-derived parent"
            )
        record_qc = "review" if qc_notes else "ok"
        scope_status = "core" if core_present else "extended"
        scope_class = (
            "core_small_molecule"
            if scope_status == "core"
            else "extended_chemical_entity"
        )
        row_out = {
            "compound_internal_id": compound_id,
            "preferred_name": preferred_name,
            "preferred_name_source": preferred_source,
            "compound_scope_status": scope_status,
            "compound_scope_class": scope_class,
            "identity_confidence": identity_confidence,
            "standard_smiles": rep.get("parent_standard_smiles", ""),
            "standard_inchi": rep.get("parent_standard_inchi", ""),
            "standard_inchikey": parent_inchikey,
            "identity_key": parent_key,
            "connectivity_key": rep.get("parent_connectivity_key", ""),
            "molecular_formula": rep.get("parent_molecular_formula", ""),
            "molecular_weight": rep.get("parent_molecular_weight", ""),
            "exact_mass": rep.get("parent_exact_mass", ""),
            "xlogp": rep.get("parent_xlogp", ""),
            "tpsa": rep.get("parent_tpsa", ""),
            "hbond_donor_count": rep.get("parent_hbond_donor_count", ""),
            "hbond_acceptor_count": rep.get("parent_hbond_acceptor_count", ""),
            "rotatable_bond_count": rep.get("parent_rotatable_bond_count", ""),
            "formal_charge": rep.get("parent_formal_charge", ""),
            "heavy_atom_count": rep.get("parent_heavy_atom_count", ""),
            "ring_count": rep.get("parent_ring_count", ""),
            "max_ring_size": rep.get("parent_max_ring_size", ""),
            "amide_bond_count": rep.get("parent_amide_bond_count", ""),
            "computed_structural_class": rep.get("computed_structural_class", ""),
            "molecule_types": serialize(molecule_types),
            "form_count": form_count_by_parent.get(parent_key, 0),
            "source_record_count": len(rows),
            "source_database_count": source_count,
            "source_databases": serialize(source_databases),
            "pubchem_cids": serialize(pubchem_ids),
            "chembl_ids": serialize(chembl_ids),
            "chebi_ids": serialize(chebi_ids),
            "drugcentral_ids": serialize(drugcentral_ids),
            "gtopdb_ligand_ids": serialize(gtopdb_ids),
            "hmdb_ids": serialize(hmdb_ids),
            "cas_numbers": serialize(cas_numbers),
            "synonym_count": len(synonyms),
            "development_status": development_status,
            "max_chembl_phase": "" if max_phase < 0 else max_phase,
            "first_approval_year": min(first_approvals) if first_approvals else "",
            "is_approved_drug": int(approved),
            "is_clinical_candidate": int(clinical),
            "is_endogenous_ligand": int(endogenous),
            "is_natural_product": int(natural_product),
            "is_chemical_probe": int(chemical_probe),
            "compound_category_tags": serialize(tags),
            "chebi_direct_classes": serialize(chebi_classes),
            "chebi_roles": serialize(chebi_roles),
            "record_qc_status": record_qc,
            "qc_notes": "; ".join(qc_notes),
        }
        batch.append(tuple(str(row_out.get(field, "")) for field in PARENT_FIELDS))
        counts[f"parent_scope_{scope_status}"] += 1
        counts[f"identity_confidence_{identity_confidence}"] += 1
        if approved:
            counts["approved_parent_compounds"] += 1
        if endogenous:
            counts["endogenous_parent_compounds"] += 1
        for synonym, sources in synonyms.values():
            synonym_writer.writerow(
                {
                    "synonym_id": hash_id(
                        "SYN",
                        compound_id,
                        synonym.casefold(),
                    ),
                    "compound_internal_id": compound_id,
                    "synonym": synonym,
                    "synonym_type": (
                        "preferred_name"
                        if synonym.casefold() == preferred_name.casefold()
                        else "synonym"
                    ),
                    "source": serialize(sources),
                }
            )
        if len(batch) >= 5_000:
            connection.executemany(insert_sql, batch)
            connection.commit()
            batch.clear()
    if batch:
        connection.executemany(insert_sql, batch)
        connection.commit()
    synonym_writer_handle.close()
    connection.execute(
        "CREATE UNIQUE INDEX idx_parent_master_id ON parent_master(compound_internal_id)"
    )
    connection.execute(
        "CREATE UNIQUE INDEX idx_parent_master_key ON parent_master(identity_key)"
    )
    connection.commit()
    return counts, form_count_by_parent


def build_forms_and_xrefs(
    connection: sqlite3.Connection,
    parent_id_map: dict[str, str],
    form_id_map: dict[str, str],
    selected_parent: dict[str, str],
) -> Counter:
    counts = Counter()
    form_fields = [
        "compound_form_id",
        "compound_internal_id",
        "form_type",
        "exact_smiles",
        "exact_inchi",
        "exact_inchikey",
        "exact_identity_key",
        "molecular_formula",
        "molecular_weight",
        "formal_charge",
        "source_record_count",
        "source_databases",
        "source_compound_ids",
        "chembl_parent_ids",
        "chembl_parent_alignment_status",
        "record_qc_status",
        "qc_notes",
    ]
    with FORM_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=form_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        query = connection.execute(
            """
            SELECT *
            FROM source_curated
            WHERE final_scope_status IN ('core','extended')
              AND exact_identity_key <> ''
            ORDER BY exact_identity_key
            """
        )
        for exact_key, iterator in itertools.groupby(
            query, key=lambda row: row[SOURCE_FIELDS.index("exact_identity_key")]
        ):
            rows = [dict(zip(SOURCE_FIELDS, row)) for row in iterator]
            rep = rows[0]
            parent_keys = {row.get("parent_identity_key", "") for row in rows}
            selected = selected_parent[exact_key]
            qc_notes = []
            distinct_exact_smiles = {
                row.get("exact_standard_smiles", "") for row in rows
            }
            if len(distinct_exact_smiles) > 1:
                qc_notes.append(
                    "standard InChI identity contains multiple standardized SMILES representations"
                )
            if len(parent_keys) > 1:
                qc_notes.append("exact form maps to multiple structure-derived parents")
            alignment = {
                row.get("chembl_parent_alignment_status", "") for row in rows
            }
            if "conflict" in alignment:
                qc_notes.append("ChEMBL parent assignment conflicts with structure parent")
            form_types = {row.get("form_type", "") for row in rows}
            form_type = (
                "salt_or_multicomponent_form"
                if "salt_or_multicomponent_form" in form_types
                else sorted(form_types)[0]
            )
            writer.writerow(
                {
                    "compound_form_id": form_id_map[exact_key],
                    "compound_internal_id": parent_id_map[selected],
                    "form_type": form_type,
                    "exact_smiles": rep.get("exact_standard_smiles", ""),
                    "exact_inchi": rep.get("exact_standard_inchi", ""),
                    "exact_inchikey": rep.get("exact_standard_inchikey", ""),
                    "exact_identity_key": exact_key,
                    "molecular_formula": rep.get("exact_molecular_formula", ""),
                    "molecular_weight": rep.get("exact_molecular_weight", ""),
                    "formal_charge": rep.get("exact_formal_charge", ""),
                    "source_record_count": len(rows),
                    "source_databases": serialize(
                        {
                            source
                            for row in rows
                            for source in split_values(row.get("source_databases", ""))
                        }
                    ),
                    "source_compound_ids": serialize(
                        {row.get("source_compound_id", "") for row in rows}
                    ),
                    "chembl_parent_ids": serialize(
                        {row.get("chembl_parent_id", "") for row in rows}
                    ),
                    "chembl_parent_alignment_status": serialize(alignment),
                    "record_qc_status": "review" if qc_notes else "ok",
                    "qc_notes": "; ".join(qc_notes),
                }
            )
            counts[f"form_type_{form_type}"] += 1

    xref_fields = [
        "source_record_id",
        "compound_internal_id",
        "compound_form_id",
        "mapping_status",
        "final_scope_status",
        "final_scope_reason",
        "source_compound_id",
        "source_compound_id_type",
        "source_compound_name",
        "source_databases",
        "source_versions",
        "pubchem_cids",
        "chembl_ids",
        "structure_parse_status",
        "original_smiles",
        "standard_exact_smiles",
        "standard_exact_inchikey",
        "standard_parent_inchikey",
        "structure_input_source",
        "source_record_qc_status",
        "structure_qc_note",
    ]
    unresolved_fields = xref_fields
    with (
        XREF_OUT.open("w", encoding="utf-8", newline="") as xref_handle,
        UNRESOLVED_OUT.open("w", encoding="utf-8", newline="") as unresolved_handle,
        EXCLUSION_OUT.open("w", encoding="utf-8", newline="") as exclusion_handle,
    ):
        xref_writer = csv.DictWriter(
            xref_handle, fieldnames=xref_fields, delimiter="\t", lineterminator="\n"
        )
        unresolved_writer = csv.DictWriter(
            unresolved_handle,
            fieldnames=unresolved_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        exclusion_writer = csv.DictWriter(
            exclusion_handle,
            fieldnames=unresolved_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        xref_writer.writeheader()
        unresolved_writer.writeheader()
        exclusion_writer.writeheader()
        for raw in connection.execute("SELECT * FROM source_curated ORDER BY input_row_number"):
            row = dict(zip(SOURCE_FIELDS, raw))
            exact_key = row.get("exact_identity_key", "")
            parent_key = row.get("parent_identity_key", "")
            scope = row.get("final_scope_status", "")
            if (
                scope in {"core", "extended"}
                and exact_key in form_id_map
                and parent_key in parent_id_map
            ):
                mapping_status = f"mapped_{scope}"
                compound_id = parent_id_map[parent_key]
                form_id = form_id_map[exact_key]
            else:
                mapping_status = scope or "unresolved"
                compound_id = ""
                form_id = ""
            output = {
                "source_record_id": row.get("source_record_id", ""),
                "compound_internal_id": compound_id,
                "compound_form_id": form_id,
                "mapping_status": mapping_status,
                "final_scope_status": scope,
                "final_scope_reason": row.get("final_scope_reason", ""),
                "source_compound_id": row.get("source_compound_id", ""),
                "source_compound_id_type": row.get("source_compound_id_type", ""),
                "source_compound_name": row.get("source_compound_name", ""),
                "source_databases": row.get("source_databases", ""),
                "source_versions": row.get("source_versions", ""),
                "pubchem_cids": row.get("pubchem_cids", ""),
                "chembl_ids": row.get("chembl_ids", ""),
                "structure_parse_status": row.get("structure_parse_status", ""),
                "original_smiles": row.get("input_smiles", ""),
                "standard_exact_smiles": row.get("exact_standard_smiles", ""),
                "standard_exact_inchikey": row.get("exact_standard_inchikey", ""),
                "standard_parent_inchikey": row.get("parent_standard_inchikey", ""),
                "structure_input_source": row.get("structure_input_source", ""),
                "source_record_qc_status": row.get("source_record_qc_status", ""),
                "structure_qc_note": row.get("structure_qc_note", ""),
            }
            xref_writer.writerow(output)
            counts[f"xref_mapping_{mapping_status}"] += 1
            if scope == "unresolved":
                unresolved_writer.writerow(output)
            elif scope == "excluded":
                exclusion_writer.writerow(output)
    return counts


def build_merge_audits(
    connection: sqlite3.Connection,
    parent_id_map: dict[str, str],
) -> Counter:
    counts = Counter()
    fields = [
        "audit_id",
        "issue_type",
        "issue_key",
        "record_count",
        "compound_internal_ids",
        "status",
        "note",
    ]
    with MERGE_AUDIT_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for name, parent_keys, count in connection.execute(
            """
            SELECT LOWER(TRIM(source_compound_name)),
                   GROUP_CONCAT(DISTINCT parent_identity_key),
                   COUNT(DISTINCT parent_identity_key)
            FROM source_curated
            WHERE final_scope_status IN ('core','extended')
              AND source_compound_name <> ''
              AND parent_identity_key <> ''
            GROUP BY LOWER(TRIM(source_compound_name))
            HAVING COUNT(DISTINCT parent_identity_key) > 1
            ORDER BY COUNT(DISTINCT parent_identity_key) DESC
            """
        ):
            if ID_LIKE_RE.match(name or ""):
                continue
            keys = (parent_keys or "").split(",")
            ids = [parent_id_map[key] for key in keys if key in parent_id_map]
            writer.writerow(
                {
                    "audit_id": hash_id("AUD", "same_name_multiple_structures", name),
                    "issue_type": "same_name_multiple_structures",
                    "issue_key": name,
                    "record_count": count,
                    "compound_internal_ids": "|".join(ids[:100]),
                    "status": "kept_separate",
                    "note": "No automatic name-only merge was performed.",
                }
            )
            counts["same_name_multiple_structures"] += 1
        for exact_key, parent_keys, count in connection.execute(
            """
            SELECT exact_identity_key,
                   GROUP_CONCAT(DISTINCT parent_identity_key),
                   COUNT(DISTINCT exact_standard_smiles)
            FROM source_curated
            WHERE final_scope_status IN ('core','extended')
              AND exact_identity_key <> ''
            GROUP BY exact_identity_key
            HAVING COUNT(DISTINCT exact_standard_smiles) > 1
            ORDER BY COUNT(DISTINCT exact_standard_smiles) DESC
            """
        ):
            keys = (parent_keys or "").split(",")
            writer.writerow(
                {
                    "audit_id": hash_id(
                        "AUD", "standard_inchi_multiple_exact_smiles", exact_key
                    ),
                    "issue_type": "standard_inchi_multiple_exact_smiles",
                    "issue_key": exact_key,
                    "record_count": count,
                    "compound_internal_ids": "|".join(
                        parent_id_map[key] for key in keys if key in parent_id_map
                    ),
                    "status": "grouped_by_standard_inchi_flagged_for_review",
                    "note": (
                        "Standard InChI may normalize tautomeric or protonation representations; "
                        "all source SMILES remain preserved in the cross-reference table."
                    ),
                }
            )
            counts["standard_inchi_multiple_exact_smiles"] += 1
        for parent_key, count in connection.execute(
            """
            SELECT parent_identity_key, COUNT(DISTINCT parent_standard_smiles)
            FROM source_curated
            WHERE final_scope_status IN ('core','extended')
              AND parent_identity_key <> ''
            GROUP BY parent_identity_key
            HAVING COUNT(DISTINCT parent_standard_smiles) > 1
            ORDER BY COUNT(DISTINCT parent_standard_smiles) DESC
            """
        ):
            writer.writerow(
                {
                    "audit_id": hash_id(
                        "AUD", "standard_inchi_multiple_parent_smiles", parent_key
                    ),
                    "issue_type": "standard_inchi_multiple_parent_smiles",
                    "issue_key": parent_key,
                    "record_count": count,
                    "compound_internal_ids": parent_id_map.get(parent_key, ""),
                    "status": "grouped_by_standard_inchi_flagged_for_review",
                    "note": (
                        "The parent Standard InChI is shared by multiple normalized SMILES; "
                        "the group is retained with an explicit QC flag."
                    ),
                }
            )
            counts["standard_inchi_multiple_parent_smiles"] += 1
        for source_id, compound_name, parent_key, chembl_parent, computed_key in connection.execute(
            """
            SELECT source_compound_id, source_compound_name, parent_identity_key,
                   chembl_parent_id, chembl_parent_structure_inchikey
            FROM source_curated
            WHERE chembl_parent_alignment_status='conflict'
            """
        ):
            writer.writerow(
                {
                    "audit_id": hash_id("AUD", "chembl_parent_conflict", source_id),
                    "issue_type": "chembl_parent_conflict",
                    "issue_key": source_id,
                    "record_count": 1,
                    "compound_internal_ids": parent_id_map.get(parent_key, ""),
                    "status": "structure_parent_retained",
                    "note": (
                        f"ChEMBL parent {chembl_parent} has structure key {computed_key}; "
                        "the reproducible structure-derived parent is retained pending review."
                    ),
                }
            )
            counts["chembl_parent_conflict"] += 1

    connectivity_fields = [
        "connectivity_key",
        "distinct_parent_compound_count",
        "compound_internal_ids",
        "interpretation",
    ]
    with CONNECTIVITY_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=connectivity_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for connectivity, parent_keys, count in connection.execute(
            """
            SELECT parent_connectivity_key,
                   GROUP_CONCAT(DISTINCT parent_identity_key),
                   COUNT(DISTINCT parent_identity_key)
            FROM source_curated
            WHERE final_scope_status IN ('core','extended')
              AND parent_connectivity_key <> ''
            GROUP BY parent_connectivity_key
            HAVING COUNT(DISTINCT parent_identity_key) > 1
            ORDER BY COUNT(DISTINCT parent_identity_key) DESC
            """
        ):
            keys = (parent_keys or "").split(",")
            writer.writerow(
                {
                    "connectivity_key": connectivity,
                    "distinct_parent_compound_count": count,
                    "compound_internal_ids": "|".join(
                        parent_id_map[key] for key in keys if key in parent_id_map
                    ),
                    "interpretation": (
                        "Same atom connectivity but distinct full InChIKeys; "
                        "stereoisomers/isotopologues remain separate."
                    ),
                }
            )
            counts["connectivity_families"] += 1
    return counts


def build_binding_outputs(
    connection: sqlite3.Connection,
    parent_id_map: dict[str, str],
    form_id_map: dict[str, str],
) -> tuple[Counter, dict[str, dict]]:
    source_mapping = {}
    for (
        compound_id,
        scope,
        exact_key,
        parent_key,
    ) in connection.execute(
        """
        SELECT source_compound_id, final_scope_status,
               exact_identity_key, parent_identity_key
        FROM source_curated
        """
    ):
        source_mapping[compound_id] = {
            "scope": scope,
            "form_id": form_id_map.get(exact_key, ""),
            "compound_internal_id": parent_id_map.get(parent_key, ""),
        }

    with EVIDENCE_INPUT.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        original_fields = list(reader.fieldnames or [])
        new_fields = [
            "compound_internal_id",
            "compound_form_id",
            "compound_mapping_status",
            "compound_scope_status_v1",
            "default_release_inclusion_v4_2",
            "compound_master_version",
        ]
        with EVIDENCE_OUT.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.DictWriter(
                destination,
                fieldnames=original_fields + new_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            connection.execute(
                """
                CREATE TABLE evidence_mapped (
                    evidence_id TEXT PRIMARY KEY,
                    target_uniprot_id TEXT,
                    approved_symbol TEXT,
                    compound_internal_id TEXT,
                    compound_form_id TEXT,
                    evidence_tier TEXT,
                    tier_rank INTEGER,
                    standard_value_nM REAL,
                    source_database TEXT,
                    pdb_ids TEXT,
                    default_v42 INTEGER
                )
                """
            )
            counts = Counter()
            batch = []
            for row in reader:
                compound_id = row.get("compound_id", "")
                mapping = source_mapping.get(compound_id)
                if mapping is None:
                    status = "not_in_v4_1_compound_index"
                    scope = "unresolved"
                    internal_id = ""
                    form_id = ""
                else:
                    scope = mapping["scope"]
                    internal_id = mapping["compound_internal_id"]
                    form_id = mapping["form_id"]
                    status = (
                        f"mapped_{scope}"
                        if internal_id and scope in {"core", "extended"}
                        else scope
                    )
                default_v42 = int(
                    row.get("default_release_inclusion", "") == "1"
                    and status == "mapped_core"
                )
                row.update(
                    {
                        "compound_internal_id": internal_id,
                        "compound_form_id": form_id,
                        "compound_mapping_status": status,
                        "compound_scope_status_v1": scope,
                        "default_release_inclusion_v4_2": default_v42,
                        "compound_master_version": "v1.0",
                    }
                )
                writer.writerow(row)
                tier = row.get("evidence_tier", "")
                rank = {"BE1": 1, "BE2": 2, "BE3": 3}.get(tier, 9)
                batch.append(
                    (
                        row.get("evidence_id", ""),
                        row.get("target_uniprot_id", ""),
                        row.get("approved_symbol", ""),
                        internal_id,
                        form_id,
                        tier,
                        rank,
                        numeric_or_blank(row.get("standard_value_nM", "")),
                        row.get("source_database", ""),
                        row.get("pdb_ids", ""),
                        default_v42,
                    )
                )
                counts[f"evidence_mapping_{status}"] += 1
                counts[f"evidence_tier_{tier}"] += 1
                counts["binding_evidence_rows"] += 1
                counts["default_v42_rows"] += default_v42
                if len(batch) >= 5_000:
                    connection.executemany(
                        "INSERT INTO evidence_mapped VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                        batch,
                    )
                    connection.commit()
                    batch.clear()
            if batch:
                connection.executemany(
                    "INSERT INTO evidence_mapped VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    batch,
                )
                connection.commit()

    connection.execute(
        "CREATE INDEX idx_evidence_parent ON evidence_mapped(compound_internal_id)"
    )
    connection.execute(
        "CREATE INDEX idx_evidence_target ON evidence_mapped(target_uniprot_id)"
    )
    connection.execute(
        "CREATE INDEX idx_evidence_default ON evidence_mapped(default_v42)"
    )
    connection.commit()

    with SITE_INPUT.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or []) + [
            "compound_internal_id",
            "compound_form_id",
            "compound_mapping_status",
            "compound_master_version",
        ]
        with SITE_OUT.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.DictWriter(
                destination, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for row in reader:
                mapping = source_mapping.get(row.get("compound_id", ""), {})
                scope = mapping.get("scope", "unresolved")
                internal_id = mapping.get("compound_internal_id", "")
                form_id = mapping.get("form_id", "")
                status = (
                    f"mapped_{scope}"
                    if internal_id and scope in {"core", "extended"}
                    else scope
                )
                row.update(
                    {
                        "compound_internal_id": internal_id,
                        "compound_form_id": form_id,
                        "compound_mapping_status": status,
                        "compound_master_version": "v1.0",
                    }
                )
                writer.writerow(row)
                counts["binding_site_rows"] += 1

    parent_names = {
        row[0]: row[1]
        for row in connection.execute(
            "SELECT compound_internal_id, preferred_name FROM parent_master"
        )
    }
    pair_fields = [
        "target_uniprot_id",
        "approved_symbol",
        "compound_internal_id",
        "preferred_name",
        "best_binding_evidence_level",
        "BE1_evidence_count",
        "BE2_evidence_count",
        "BE3_evidence_count",
        "binding_evidence_count",
        "independent_source_count",
        "independent_sources",
        "best_standard_value_nM",
        "compound_form_count",
        "compound_master_version",
    ]
    with PAIR_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=pair_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        query = connection.execute(
            """
            SELECT target_uniprot_id, MAX(approved_symbol), compound_internal_id,
                   MIN(tier_rank),
                   SUM(CASE WHEN evidence_tier='BE1' THEN 1 ELSE 0 END),
                   SUM(CASE WHEN evidence_tier='BE2' THEN 1 ELSE 0 END),
                   SUM(CASE WHEN evidence_tier='BE3' THEN 1 ELSE 0 END),
                   COUNT(*),
                   COUNT(DISTINCT source_database),
                   GROUP_CONCAT(DISTINCT source_database),
                   MIN(CASE WHEN standard_value_nM > 0 THEN standard_value_nM END),
                   COUNT(DISTINCT compound_form_id)
            FROM evidence_mapped
            WHERE default_v42=1 AND compound_internal_id <> ''
            GROUP BY target_uniprot_id, compound_internal_id
            ORDER BY target_uniprot_id, compound_internal_id
            """
        )
        for row in query:
            rank = row[3]
            writer.writerow(
                {
                    "target_uniprot_id": row[0],
                    "approved_symbol": row[1],
                    "compound_internal_id": row[2],
                    "preferred_name": parent_names.get(row[2], ""),
                    "best_binding_evidence_level": {1: "BE1", 2: "BE2", 3: "BE3"}.get(
                        rank, ""
                    ),
                    "BE1_evidence_count": row[4],
                    "BE2_evidence_count": row[5],
                    "BE3_evidence_count": row[6],
                    "binding_evidence_count": row[7],
                    "independent_source_count": row[8],
                    "independent_sources": (row[9] or "").replace(",", "|"),
                    "best_standard_value_nM": row[10] if row[10] is not None else "",
                    "compound_form_count": row[11],
                    "compound_master_version": "v1.0",
                }
            )
            counts["protein_compound_pair_rows"] += 1

    compound_summary = {}
    query = connection.execute(
        """
        SELECT compound_internal_id,
               COUNT(DISTINCT target_uniprot_id),
               MIN(tier_rank),
               SUM(CASE WHEN evidence_tier='BE1' THEN 1 ELSE 0 END),
               SUM(CASE WHEN evidence_tier='BE2' THEN 1 ELSE 0 END),
               SUM(CASE WHEN evidence_tier='BE3' THEN 1 ELSE 0 END),
               COUNT(*),
               COUNT(DISTINCT source_database),
               GROUP_CONCAT(DISTINCT source_database),
               MIN(CASE WHEN standard_value_nM > 0 THEN standard_value_nM END)
        FROM evidence_mapped
        WHERE default_v42=1 AND compound_internal_id <> ''
        GROUP BY compound_internal_id
        """
    )
    for row in query:
        compound_summary[row[0]] = {
            "protein_target_count": row[1],
            "best_binding_evidence_level": {1: "BE1", 2: "BE2", 3: "BE3"}.get(
                row[2], ""
            ),
            "BE1_evidence_count": row[3],
            "BE2_evidence_count": row[4],
            "BE3_evidence_count": row[5],
            "binding_evidence_count": row[6],
            "binding_source_count": row[7],
            "binding_sources": (row[8] or "").replace(",", "|"),
            "best_standard_value_nM": row[9] if row[9] is not None else "",
        }

    compound_protein_fields = [
        "compound_internal_id",
        "preferred_name",
        "protein_target_count",
        "best_binding_evidence_level",
        "BE1_evidence_count",
        "BE2_evidence_count",
        "BE3_evidence_count",
        "binding_evidence_count",
        "binding_source_count",
        "binding_sources",
        "best_standard_value_nM",
        "compound_master_version",
    ]
    with COMPOUND_PROTEIN_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=compound_protein_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for compound_id in sorted(compound_summary):
            writer.writerow(
                {
                    "compound_internal_id": compound_id,
                    "preferred_name": parent_names.get(compound_id, ""),
                    **compound_summary[compound_id],
                    "compound_master_version": "v1.0",
                }
            )

    return counts, compound_summary


def export_final_parent_master(
    connection: sqlite3.Connection,
    compound_summary: dict[str, dict],
) -> Counter:
    evidence_fields = [
        "protein_target_count",
        "best_binding_evidence_level",
        "BE1_evidence_count",
        "BE2_evidence_count",
        "BE3_evidence_count",
        "binding_evidence_count",
        "binding_source_count",
        "binding_sources",
        "best_standard_value_nM",
    ]
    release_fields = [
        "release_version",
        "release_date",
        "standardization_policy",
    ]
    fields = PARENT_FIELDS + evidence_fields + release_fields
    counts = Counter()
    with PARENT_MASTER_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for raw in connection.execute(
            "SELECT * FROM parent_master ORDER BY compound_internal_id"
        ):
            row = dict(zip(PARENT_FIELDS, raw))
            summary = compound_summary.get(row["compound_internal_id"], {})
            for field in evidence_fields:
                row[field] = summary.get(field, 0 if "count" in field else "")
            row.update(
                {
                    "release_version": "small-molecule-v1.0",
                    "release_date": "2026-07-27",
                    "standardization_policy": "MemProDB-compound-standardization-v1",
                }
            )
            writer.writerow(row)
            counts["parent_master_rows"] += 1
            counts[f"parent_master_scope_{row['compound_scope_status']}"] += 1
            if summary:
                counts["parent_compounds_with_default_binding"] += 1
    return counts


def build_protein_outputs(
    connection: sqlite3.Connection,
) -> Counter:
    aggregate = {}
    for row in connection.execute(
        """
        SELECT target_uniprot_id,
               COUNT(DISTINCT compound_internal_id),
               COUNT(DISTINCT CASE WHEN evidence_tier='BE1' THEN compound_internal_id END),
               COUNT(DISTINCT CASE WHEN evidence_tier='BE2' THEN compound_internal_id END),
               COUNT(DISTINCT CASE WHEN evidence_tier='BE3' THEN compound_internal_id END),
               MIN(tier_rank),
               COUNT(*),
               COUNT(DISTINCT source_database),
               GROUP_CONCAT(DISTINCT source_database),
               MIN(CASE WHEN standard_value_nM > 0 THEN standard_value_nM END)
        FROM evidence_mapped
        WHERE default_v42=1 AND compound_internal_id <> ''
        GROUP BY target_uniprot_id
        """
    ):
        aggregate[row[0]] = {
            "canonical_core_small_molecule_count": row[1],
            "BE1_compound_count": row[2],
            "BE2_compound_count": row[3],
            "BE3_compound_count": row[4],
            "best_binding_evidence_level_v4_2": {1: "BE1", 2: "BE2", 3: "BE3"}.get(
                row[5], ""
            ),
            "binding_evidence_count_v4_2": row[6],
            "binding_source_count_v4_2": row[7],
            "binding_sources_v4_2": (row[8] or "").replace(",", "|"),
            "best_standard_value_nM": row[9] if row[9] is not None else "",
        }
    fields = [
        "target_uniprot_id",
        "approved_symbol",
        "canonical_core_small_molecule_count",
        "BE1_compound_count",
        "BE2_compound_count",
        "BE3_compound_count",
        "has_BE1_structural_binding",
        "has_BE2_direct_binding",
        "has_BE3_pharmacology",
        "best_binding_evidence_level_v4_2",
        "binding_evidence_count_v4_2",
        "binding_source_count_v4_2",
        "binding_sources_v4_2",
        "best_standard_value_nM",
        "compound_master_version",
    ]
    counts = Counter()
    with (
        PROTEIN_MASTER_INPUT.open("r", encoding="utf-8", newline="") as source,
        PROTEIN_BINDING_OUT.open("w", encoding="utf-8", newline="") as destination,
        PROTEIN_MASTER_V55_OUT.open("w", encoding="utf-8", newline="") as master_out,
    ):
        reader = csv.DictReader(source, delimiter="\t")
        protein_writer = csv.DictWriter(
            destination, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        master_new_fields = list(reader.fieldnames or []) + [
            "canonical_core_small_molecule_count_v55",
            "canonical_BE1_compound_count_v55",
            "canonical_BE2_compound_count_v55",
            "canonical_BE3_compound_count_v55",
            "canonical_best_binding_evidence_level_v55",
            "canonical_binding_evidence_count_v55",
            "canonical_binding_source_count_v55",
            "canonical_binding_sources_v55",
            "canonical_best_standard_value_nM_v55",
            "compound_master_version_v55",
            "v55_release_date",
        ]
        master_writer = csv.DictWriter(
            master_out,
            fieldnames=master_new_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        protein_writer.writeheader()
        master_writer.writeheader()
        for row in reader:
            target = row.get("target_uniprot_id", "")
            summary = aggregate.get(target, {})
            protein_row = {
                "target_uniprot_id": target,
                "approved_symbol": row.get("approved_symbol", ""),
                "canonical_core_small_molecule_count": summary.get(
                    "canonical_core_small_molecule_count", 0
                ),
                "BE1_compound_count": summary.get("BE1_compound_count", 0),
                "BE2_compound_count": summary.get("BE2_compound_count", 0),
                "BE3_compound_count": summary.get("BE3_compound_count", 0),
                "has_BE1_structural_binding": int(
                    summary.get("BE1_compound_count", 0) > 0
                ),
                "has_BE2_direct_binding": int(
                    summary.get("BE2_compound_count", 0) > 0
                ),
                "has_BE3_pharmacology": int(
                    summary.get("BE3_compound_count", 0) > 0
                ),
                "best_binding_evidence_level_v4_2": summary.get(
                    "best_binding_evidence_level_v4_2", ""
                ),
                "binding_evidence_count_v4_2": summary.get(
                    "binding_evidence_count_v4_2", 0
                ),
                "binding_source_count_v4_2": summary.get(
                    "binding_source_count_v4_2", 0
                ),
                "binding_sources_v4_2": summary.get("binding_sources_v4_2", ""),
                "best_standard_value_nM": summary.get("best_standard_value_nM", ""),
                "compound_master_version": "v1.0",
            }
            protein_writer.writerow(protein_row)
            row.update(
                {
                    "canonical_core_small_molecule_count_v55": protein_row[
                        "canonical_core_small_molecule_count"
                    ],
                    "canonical_BE1_compound_count_v55": protein_row[
                        "BE1_compound_count"
                    ],
                    "canonical_BE2_compound_count_v55": protein_row[
                        "BE2_compound_count"
                    ],
                    "canonical_BE3_compound_count_v55": protein_row[
                        "BE3_compound_count"
                    ],
                    "canonical_best_binding_evidence_level_v55": protein_row[
                        "best_binding_evidence_level_v4_2"
                    ],
                    "canonical_binding_evidence_count_v55": protein_row[
                        "binding_evidence_count_v4_2"
                    ],
                    "canonical_binding_source_count_v55": protein_row[
                        "binding_source_count_v4_2"
                    ],
                    "canonical_binding_sources_v55": protein_row[
                        "binding_sources_v4_2"
                    ],
                    "canonical_best_standard_value_nM_v55": protein_row[
                        "best_standard_value_nM"
                    ],
                    "compound_master_version_v55": "v1.0",
                    "v55_release_date": "2026-07-27",
                }
            )
            master_writer.writerow(row)
            counts["protein_summary_rows"] += 1
            if target in aggregate:
                counts["proteins_with_canonical_core_compounds"] += 1
    return counts


def make_review_views() -> None:
    def copy_selected(
        source_path: Path,
        destination_path: Path,
        score,
        limit: int,
    ) -> None:
        with source_path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
            fields = list(rows[0].keys()) if rows else []
        selected = sorted(rows, key=score)[:limit]
        with destination_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(selected)

    make_num = lambda value: numeric_or_blank(value) or 0
    copy_selected(
        PARENT_MASTER_OUT,
        MASTER_REVIEW_OUT,
        lambda row: (
            -int(row.get("is_approved_drug", "0") or 0),
            -int(row.get("is_endogenous_ligand", "0") or 0),
            -make_num(row.get("BE1_evidence_count", "")),
            -make_num(row.get("protein_target_count", "")),
            -make_num(row.get("binding_evidence_count", "")),
            row.get("preferred_name", ""),
        ),
        10_000,
    )
    copy_selected(
        FORM_OUT,
        FORM_REVIEW_OUT,
        lambda row: (
            row.get("record_qc_status") != "review",
            row.get("form_type") != "salt_or_multicomponent_form",
            -make_num(row.get("source_record_count", "")),
            row.get("compound_form_id", ""),
        ),
        8_000,
    )
    copy_selected(
        COMPOUND_PROTEIN_OUT,
        PAIR_REVIEW_OUT,
        lambda row: (
            row.get("best_binding_evidence_level") != "BE1",
            -make_num(row.get("protein_target_count", "")),
            -make_num(row.get("binding_evidence_count", "")),
            row.get("preferred_name", ""),
        ),
        12_000,
    )
    audit_rows = []
    for path, audit_type in (
        (UNRESOLVED_OUT, "unresolved"),
        (EXCLUSION_OUT, "excluded"),
    ):
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in itertools.islice(csv.DictReader(handle, delimiter="\t"), 2_500):
                row["audit_type"] = audit_type
                audit_rows.append(row)
    fields = ["audit_type"] + [
        field for field in audit_rows[0].keys() if field != "audit_type"
    ] if audit_rows else ["audit_type"]
    with AUDIT_REVIEW_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(audit_rows)


def line_count(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def validate(
    connection: sqlite3.Connection,
    counters: Counter,
) -> dict:
    source_rows = connection.execute("SELECT COUNT(*) FROM source_curated").fetchone()[0]
    source_unique = connection.execute(
        "SELECT COUNT(DISTINCT source_record_id) FROM source_curated"
    ).fetchone()[0]
    parent_rows = connection.execute("SELECT COUNT(*) FROM parent_master").fetchone()[0]
    parent_unique = connection.execute(
        "SELECT COUNT(DISTINCT compound_internal_id) FROM parent_master"
    ).fetchone()[0]
    form_rows = connection.execute("SELECT COUNT(*) FROM form_groups").fetchone()[0]
    form_unique = connection.execute(
        "SELECT COUNT(DISTINCT compound_form_id) FROM form_groups"
    ).fetchone()[0]
    evidence_rows = connection.execute("SELECT COUNT(*) FROM evidence_mapped").fetchone()[0]
    mapped_fk_failures = connection.execute(
        """
        SELECT COUNT(*)
        FROM evidence_mapped e
        LEFT JOIN parent_master p
          ON e.compound_internal_id=p.compound_internal_id
        WHERE e.compound_internal_id <> ''
          AND p.compound_internal_id IS NULL
        """
    ).fetchone()[0]
    validation = {
        "source_index_rows": source_rows,
        "source_record_key_unique": source_rows == source_unique,
        "original_v4_1_index_rows": 302_879,
        "repaired_missing_source_records": source_rows - 302_879,
        "parent_master_rows": parent_rows,
        "parent_master_key_unique": parent_rows == parent_unique,
        "form_rows": form_rows,
        "form_key_unique": form_rows == form_unique,
        "binding_evidence_rows_v4_2": evidence_rows,
        "binding_evidence_row_count_preserved": evidence_rows == 665_214,
        "binding_parent_foreign_key_failures": mapped_fk_failures,
        "binding_parent_foreign_key_closed": mapped_fk_failures == 0,
        "binding_site_rows_v4_2": line_count(SITE_OUT),
        "protein_compound_pair_rows_v4_2": line_count(PAIR_OUT),
        "compound_protein_summary_rows": line_count(COMPOUND_PROTEIN_OUT),
        "protein_binding_summary_rows": line_count(PROTEIN_BINDING_OUT),
        "protein_master_v5_5_rows": line_count(PROTEIN_MASTER_V55_OUT),
        "unresolved_source_records": line_count(UNRESOLVED_OUT),
        "excluded_source_records": line_count(EXCLUSION_OUT),
        "merge_audit_rows": line_count(MERGE_AUDIT_OUT),
        "connectivity_family_rows": line_count(CONNECTIVITY_OUT),
        "counters": dict(counters),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    validation["status"] = (
        "passed"
        if (
            validation["source_record_key_unique"]
            and validation["parent_master_key_unique"]
            and validation["form_key_unique"]
            and validation["binding_evidence_row_count_preserved"]
            and validation["binding_parent_foreign_key_closed"]
            and validation["protein_binding_summary_rows"] == 10_997
            and validation["protein_master_v5_5_rows"] == 10_997
        )
        else "failed"
    )
    VALIDATION_OUT.write_text(
        json.dumps(validation, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return validation


def main() -> None:
    metadata = load_jsonl_by_id(CHEMBL_METADATA, "molecule_chembl_id")
    structures = load_chembl_structures()
    chebi = load_chebi()
    connection = create_database()
    counters = Counter()
    counters.update(ingest_source_records(connection, metadata, structures))
    parent_id_map, form_id_map, selected_parent, group_counts = build_group_maps(
        connection
    )
    counters.update(group_counts)
    parent_counts, _ = build_parent_master(
        connection, parent_id_map, metadata, chebi
    )
    counters.update(parent_counts)
    counters.update(
        build_forms_and_xrefs(
            connection,
            parent_id_map,
            form_id_map,
            selected_parent,
        )
    )
    counters.update(build_merge_audits(connection, parent_id_map))
    binding_counts, compound_summary = build_binding_outputs(
        connection,
        parent_id_map,
        form_id_map,
    )
    counters.update(binding_counts)
    counters.update(export_final_parent_master(connection, compound_summary))
    counters.update(build_protein_outputs(connection))
    make_review_views()
    validation = validate(connection, counters)
    connection.close()
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
