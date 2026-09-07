#!/usr/bin/env python3
import csv
import json
from collections import Counter
from pathlib import Path


WORK = Path(__file__).resolve().parents[1]
RELEASE = WORK / "release"
REPORTS = WORK / "reports"
REVIEW = WORK / "review"
CROSSWALKS = WORK / "crosswalks"
STATS = json.loads((REPORTS / "V5_BUILD_STATS.json").read_text(encoding="utf-8"))


def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


master = read_tsv(RELEASE / "human_membrane_protein_master_v5.tsv")
unreviewed = read_tsv(CROSSWALKS / "unreviewed_canonical_resolution_v5.tsv")
special82 = read_tsv(REVIEW / "unreviewed_82_no_reviewed_gene_symbol_match_v5.tsv")
legacy = read_tsv(CROSSWALKS / "legacy_377_resolution_v5.tsv")
topology = read_tsv(CROSSWALKS / "single_pass_resolution_v5.tsv")

new_by_tier = Counter(
    row["release_tier_v5"]
    for row in master
    if row["record_status_v5"] == "new_reviewed_external_union"
)
core_unclassified = sum(
    row["release_tier_v5"] in {"A", "B"}
    and row["functional_primary_class_v5"] == "unclassified"
    for row in master
)
extended_unclassified = sum(
    row["release_tier_v5"] in {"A", "B", "C"}
    and row["functional_primary_class_v5"] == "unclassified"
    for row in master
)
resolved_exact = sum(
    row["resolved_topology_v5"]
    in {"single_pass_type_i", "single_pass_type_ii", "single_pass_type_iii"}
    for row in topology
)
resolved_partial = sum(
    row["resolved_topology_v5"] == "single_pass_type_i_or_iii"
    for row in topology
)

sources = [
    {
        "source_id": "UniProtKB_Swiss-Prot",
        "version_or_date": "retrieved 2026-07-24",
        "role": "canonical reviewed human reference proteome and primary annotation",
        "url": "https://rest.uniprot.org/",
        "license_status": "UniProt terms/citation required; verify redistribution terms at publication freeze",
        "promotion_rule": "curated transmembrane/intramembrane/lipid-anchor annotations seed v4 tiers",
    },
    {
        "source_id": "Human_Protein_Atlas",
        "version_or_date": "25.1",
        "role": "independent protein-class and subcellular-location evidence",
        "url": "https://www.proteinatlas.org/about/download",
        "license_status": "CC BY 4.0 for copyrightable HPA database content; third-party constraints may apply",
        "promotion_rule": "never promotes location-only evidence to integral; HPA-only records remain Tier R",
    },
    {
        "source_id": "UniTmp_HTP",
        "version_or_date": "d.2.2",
        "role": "human transmembrane topology and evidence class",
        "url": "https://htp.unitmp.org/downloads",
        "license_status": "public bulk download and citation available; explicit redistribution license not confirmed",
        "promotion_rule": "3D or Experiment evidence may promote to Tier A; Exists/Prediction/TOPDOM-only additions remain Tier R",
    },
    {
        "source_id": "UniTmp_PDBTM",
        "version_or_date": "current download 2026-07-24",
        "role": "PDB-level experimentally determined membrane structure orientation",
        "url": "https://pdbtm.unitmp.org/downloads",
        "license_status": "public bulk download; XML carries copyright notice; obtain redistribution confirmation before public mirroring",
        "promotion_rule": "cross-validation only when mapped through a protein PDB cross-reference; not sufficient alone for new Tier A",
    },
    {
        "source_id": "OPM",
        "version_or_date": "current API 2026-07-24",
        "role": "orientation of experimental structures in membranes",
        "url": "https://opm.phar.umich.edu/download",
        "license_status": "public CSV/API and citation available; explicit database redistribution license not confirmed",
        "promotion_rule": "Tier A only when an accession-mapped structure has an integral topology subunit/segment",
    },
    {
        "source_id": "Membranome",
        "version_or_date": "current API 2026-07-24",
        "role": "curated single-pass membrane-protein set and hierarchy",
        "url": "https://biomembhub.org/membranome/download",
        "license_status": "public CSV/API and citation available; inspect linked site license before public redistribution",
        "promotion_rule": "reviewed human accession in Membranome supports Tier A single-pass membership",
    },
    {
        "source_id": "GPCRdb",
        "version_or_date": "current API 2026-07-24",
        "role": "formal GPCR class/family/ligand-type/subfamily hierarchy",
        "url": "https://gpcrdb.org/services/",
        "license_status": "public API and citation available; verify redistribution terms at publication freeze",
        "promotion_rule": "classification only; does not independently determine membrane membership",
    },
    {
        "source_id": "IUPHAR_GtoPdb",
        "version_or_date": "2026.2",
        "role": "expert target class and family hierarchy",
        "url": "https://www.guidetopharmacology.org/download.jsp",
        "license_status": "database ODbL; contents CC BY-SA 4.0",
        "promotion_rule": "classification only; target membership does not independently determine membrane membership",
    },
    {
        "source_id": "TCDB",
        "version_or_date": "human.csv retrieved 2026-07-24",
        "role": "IUBMB-style transport classification IDs for human transport systems",
        "url": "https://www.tcdb.org/public/",
        "license_status": "public download; explicit redistribution license not confirmed",
        "promotion_rule": "classification only because TCDB can include non-membrane transport-system subunits",
    },
]
write_tsv(
    RELEASE / "SOURCE_REGISTRY_v5.tsv",
    sources,
    [
        "source_id",
        "version_or_date",
        "role",
        "url",
        "license_status",
        "promotion_rule",
    ],
)

readme = f"""# Human membrane-protein master database v5

This release completes the first-table rebuild as a **reviewed human membrane-protein union with explicit evidence tiers**, while keeping unreviewed TrEMBL candidates and uncertain external-only entries in separate review layers.

## Main result

- Reviewed human reference proteome examined: **{STATS['reviewed_reference_proteome_rows']:,}**
- v4 membrane/membrane-associated baseline: **{STATS['v4_master_rows']:,}**
- v5 reviewed union: **{STATS['v5_union_rows']:,}**
- New reviewed proteins surfaced by independent sources: **{STATS['v5_new_reviewed_external_union_rows']:,}**
- Tier A integral core: **{STATS['release_tiers']['A']:,}**
- Tier B monotopic/lipid-anchored extension: **{STATS['release_tiers']['B']:,}**
- Tier C peripheral membrane-associated extension: **{STATS['release_tiers']['C']:,}**
- Tier R evidence/review queue: **{STATS['release_tiers']['R']:,}**

The three publication views are:

1. `human_integral_membrane_view_v5.tsv` — Tier A.
2. `human_core_membrane_view_v5.tsv` — Tiers A+B.
3. `human_extended_membrane_view_v5.tsv` — Tiers A+B+C.

Tier R is not mixed into the default website result set. It is preserved in `human_membrane_review_queue_v5.tsv` so that coverage is broad without presenting predictions or location-only evidence as settled integral-membrane biology.

## Previously open issues

- **9,823 TrEMBL candidates:** all received a deterministic canonical-resolution status. {STATS['unreviewed_resolution_status']['absorbed_to_reviewed_canonical']:,} map to one reviewed canonical accession; the remainder are retained in auditable ambiguous or review states.
- **82 special TrEMBL rows:** the exact original 82 `no_reviewed_gene_symbol_match` records are isolated; {sum(r['canonical_resolution_status']=='absorbed_to_reviewed_canonical' for r in special82)} were rescued by HGNC/GeneID/Ensembl evidence, {sum(r['canonical_resolution_status']=='retain_review_unreviewed_external_support' for r in special82)} have external membrane support without a reviewed canonical match, and {sum(r['canonical_resolution_status']=='retain_review_unreviewed_unconfirmed' for r in special82)} remain unconfirmed.
- **377 legacy v3.1 proteins:** all identifiers were resolved: 376 remain current reviewed accessions and one was mapped by approved symbol. Membrane-scope decisions are {STATS['legacy_membrane_scope_status']}.
- **680 single-pass unresolved:** {resolved_exact} received exact I/II/III assignments, {resolved_partial} received orientation but remain I-versus-III ambiguous, and {sum(r['resolved_topology_v5']=='single_pass_unresolved' for r in topology)} remain unresolved.
- **4,423 v4 unclassified:** strict residual unclassified among retained v4 rows is now **{STATS['v5_unclassified_retained_v4']:,}**. Family-level classifications are explicitly labelled and not misrepresented as mature functional classes.
- **GPCR/ion-channel/transporter hierarchy:** GPCRdb, GtoPdb 2026.2 and TCDB identifiers are stored in dedicated columns.
- **Peripheral proteins:** Tier C is a separate extension and is excluded from integral/core defaults.

## Important interpretation

“Complete” here means the union of the named, versioned sources under the documented policy. It does not mean that biology has a permanently closed list. New isoforms, new reviewed accessions, revised topology predictions and database updates will change future releases.
"""
(RELEASE / "README_v5.md").write_text(readme, encoding="utf-8")

methods = """# Methods

## Record unit and species

The public record unit is a canonical UniProtKB accession from the Homo sapiens reference proteome UP000005640. Alternative/unreviewed accessions are not silently merged: every TrEMBL candidate has a crosswalk row that records the identifiers used and whether the match is unique, ambiguous or absent.

## Evidence hierarchy

Evidence is stored source by source. The decision priority is:

1. Experimentally determined membrane structure/orientation (PDBTM or accession-mapped OPM integral segment).
2. Curated UniProt transmembrane/intramembrane annotation.
3. Specialist curated database classification (notably Membranome).
4. Integrated or multiple-source topology evidence (UniTmp HTP).
5. Prediction-only membrane evidence.
6. Subcellular localization evidence.
7. Name/family inference.

Location evidence is never treated as proof of membrane insertion. TCDB and pharmacology target membership are used for classification, not as sole membrane-membership evidence.

## Release tiers

- Tier A: integral membrane proteins supported by the v4 UniProt layer or strong specialist/structural evidence.
- Tier B: integral monotopic or lipid-anchored proteins.
- Tier C: peripheral membrane-associated proteins.
- Tier R: records with prediction-only, location-only, non-integral OPM association, or unresolved conflicts.
- Tier E: excluded/non-membrane decisions are retained in audit crosswalks rather than the published membrane union.

## Single-pass topology

For the 680 v4 `single_pass_unresolved_type` rows, UniTmp HTP terminal orientation was applied only when HTP reported one transmembrane segment. N-in/C-out supports type II. N-out/C-in plus a UniProt signal peptide supports type I; N-out/C-in with an N-terminal signal anchor supports type III. N-out/C-in without enough evidence to distinguish type I from III is labelled `single_pass_type_i_or_iii`, not forced.

## Functional classification

GPCRdb has highest precedence for GPCR hierarchy. GtoPdb supplies expert receptor, ion-channel, transporter and other target families. TCDB supplies TC identifiers for transporters. Existing v4 enzyme/receptor/channel assignments are retained when no specialist source overrides them. Conservative UniProt name/family rules provide a broad browsing category; `other_family_defined` means a family exists but a mature functional primary class is still unresolved.

## Reproducibility

Raw source URLs, retrieval dates, byte sizes and SHA-256 hashes are recorded in `source_download_manifest.json` and `SOURCE_REGISTRY_v5.tsv`. The build and validation scripts are included in the release.
"""
(RELEASE / "METHODS_v5.md").write_text(methods, encoding="utf-8")

policy = """# Inclusion policy

The default website should expose Tier A+B as the **core** membrane database. Tier A alone is the strict integral-membrane view. Tier C should be an opt-in peripheral extension. Tier R should be visible only in curator/review interfaces or explicitly labelled exploratory downloads.

Conflicting sources are not settled by majority vote. Structure and curated topology outrank predictions; location-only evidence cannot promote a record to Tier A. A record may therefore appear in the union and still remain Tier R.

TrEMBL candidate rows are not public canonical protein records unless a reviewed canonical accession is uniquely identified. Ambiguous mappings stay in the review crosswalk.
"""
(RELEASE / "INCLUSION_POLICY_v5.md").write_text(policy, encoding="utf-8")

dictionary = """# Data dictionary (v5 additions)

| Column | Meaning |
|---|---|
| `release_tier_v5` | A integral; B monotopic/lipid-anchored; C peripheral; R review |
| `record_status_v5` | Retained v4 record or new reviewed external-union record |
| `membrane_decision_v5` | Machine-readable reason for tier assignment |
| `membrane_confidence_v5` | High, moderate or review |
| `membrane_scope_v5` | Final v5 scope; preserves the original v4 scope in `membrane_scope` |
| `membrane_topology_v5` | Final v5 topology after conservative single-pass resolution |
| `transmembrane_count_v5` | Final v5 TM count used for browsing/QC |
| `independent_membrane_sources_v5` | Semicolon-separated source evidence present for the record |
| `hpa_*_v5` | HPA v25.1 prediction/localization fields; localization is not insertion proof |
| `htp_*_v5` | UniTmp HTP evidence class, TM count, reliability and terminal sides |
| `membranome_*_v5` | Membranome single-pass family and helix segment |
| `opm_*_v5` | OPM structure and integral-segment status |
| `pdbtm_*_v5` | PDB IDs that overlap PDBTM membrane structures |
| `single_pass_resolution_v5` | Exact, partial or unresolved single-pass topology result |
| `functional_primary_class_v5` | Main website browsing class |
| `functional_subclass_v5` | Specialist hierarchy or family-level fallback |
| `classification_source_v5` | Source/rule that produced the v5 class |
| `classification_status_v5` | Specialist aligned, rule classified, family-level only or unclassified |
| `gpcrdb_*_v5` | GPCRdb class/family/ligand type/subfamily |
| `gtopdb_*_v5` | GtoPdb target type and family |
| `tcdb_tcids_v5` | Semicolon-separated TC classification identifiers |
| `cross_source_conflict_v5` | 1 when the record needs conflict/review attention |
| `review_reason_v5` | Explicit conflict or review reason |
"""
(RELEASE / "DATA_DICTIONARY_v5.md").write_text(dictionary, encoding="utf-8")

completeness = f"""# Completeness and residual review report

## What is complete in this release

- All {STATS['reviewed_reference_proteome_rows']:,} reviewed human reference-proteome entries were eligible for the external-source union.
- All {STATS['unreviewed_candidates']:,} TrEMBL transmembrane candidates have canonical-resolution rows.
- The exact 82 no-reviewed-gene-symbol special rows are isolated and adjudicated by additional identifiers/external evidence.
- All {STATS['legacy_rows']} v3.1 legacy rows have identifier and membrane-scope outcomes.
- All {STATS['single_pass_original_unresolved']} unresolved single-pass rows have a v5 topology outcome.
- All v5 master rows have a release tier, membrane decision, evidence matrix row and classification status.

## Residual review, intentionally not hidden

- Tier R contains **{STATS['release_tiers']['R']:,}** reviewed proteins.
- TrEMBL rows without one unique reviewed canonical accession: **{STATS['unreviewed_without_single_canonical']:,}**.
- Single-pass exact type still unresolved: **{sum(r['resolved_topology_v5']=='single_pass_unresolved' for r in topology):,}**; another **{resolved_partial:,}** have orientation but remain type I/III ambiguous.
- Strict unclassified: **{STATS['v5_unclassified_all_union']:,}** in the full union, **{core_unclassified:,}** in Tier A+B, and **{extended_unclassified:,}** in Tier A+B+C.
- Cross-source conflict/review rows: **{STATS['conflict_review_rows']:,}**.
- OPM, Membranome, UniTmp/PDBTM, GPCRdb and TCDB provide public downloads/APIs, but explicit redistribution terms were not uniformly stated. Confirm those terms before mirroring third-party annotations in a public production website.

## Coverage is versioned

The union is complete only relative to the sources and dates in `SOURCE_REGISTRY_v5.tsv`. A release refresh must rerun downloads, mapping, tiering and QA; counts should not be expected to remain constant.
"""
(RELEASE / "COMPLETENESS_REPORT_v5.md").write_text(
    completeness, encoding="utf-8"
)

print("Wrote v5 documentation and source registry.")
