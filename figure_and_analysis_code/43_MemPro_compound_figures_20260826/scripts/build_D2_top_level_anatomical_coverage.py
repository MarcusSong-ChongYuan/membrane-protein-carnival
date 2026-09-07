"""Rebuild D2 as an ontology-projected top-level anatomical-system coverage chart.

The frozen V7.2 disease-to-anatomy assertions are read only.  A frozen
Uberon visualization snapshot is used solely to create a reproducible,
same-level display projection for this figure and its associated source data.
"""

from __future__ import annotations

import csv
import hashlib
from collections import deque
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


FORMAL = Path(r"D:\finale\FORMAL\01_core_v72\01_release_tables")
ROOT = Path(r"D:\finale\43_MemPro_compound_figures_20260826")
OUT = ROOT / "results" / "context_figure_pool"
SNAPSHOT = OUT / "source_data" / "uberon-basic_visualization_snapshot.obo"
SYSTEM_ROOT = "UBERON:0000467"  # anatomical system

mpl.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


def clean_obo_target(value: str) -> str:
    """Remove OBO comments and provenance qualifiers from a target identifier."""
    return value.split(" ! ", 1)[0].split(" {", 1)[0].strip()


def parse_uberon(path: Path) -> tuple[dict[str, dict], str]:
    terms: dict[str, dict] = {}
    current: dict | None = None
    version = "not found"
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("data-version: "):
                version = line.split(": ", 1)[1]
            if line == "[Term]":
                current = {"id": None, "name": "", "is_a": set(), "parents": set(), "obsolete": False}
                continue
            if current is None:
                continue
            if not line.strip():
                if current["id"]:
                    terms[current["id"]] = current
                current = None
                continue
            if line.startswith("id: "):
                current["id"] = line[4:]
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line.startswith("is_a: "):
                parent = clean_obo_target(line[6:])
                current["is_a"].add(parent)
                current["parents"].add(parent)
            elif line.startswith("relationship: part_of "):
                current["parents"].add(clean_obo_target(line[len("relationship: part_of "):]))
            elif line == "is_obsolete: true":
                current["obsolete"] = True
    if current is not None and current["id"]:
        terms[current["id"]] = current
    return terms, version


def ancestor_distances(term_id: str, terms: dict[str, dict]) -> dict[str, int]:
    result: dict[str, int] = {}
    queue: deque[tuple[str, int]] = deque([(term_id, 0)])
    while queue:
        node, distance = queue.popleft()
        for parent in terms.get(node, {}).get("parents", set()):
            if parent not in result or distance + 1 < result[parent]:
                result[parent] = distance + 1
                queue.append((parent, distance + 1))
    return result


def build_terminal_anchor_map(terms: dict[str, dict]) -> tuple[set[str], dict[str, str]]:
    """Return direct system candidates and their non-nested terminal anchors.

    Candidates are direct *is_a* children of ``anatomical system``.  A candidate
    is collapsed only when Uberon itself places it under another candidate via
    ``is_a`` or ``part_of``.  This avoids any label-based manual merge.
    """
    candidates = {
        term_id
        for term_id, term in terms.items()
        if not term["obsolete"] and SYSTEM_ROOT in term["is_a"]
    }
    terminal: dict[str, str] = {}
    for candidate in candidates:
        ancestors = ancestor_distances(candidate, terms)
        ancestor_candidates = [(distance, ancestor) for ancestor, distance in ancestors.items() if ancestor in candidates]
        if not ancestor_candidates:
            terminal[candidate] = candidate
            continue
        # Project to the nearest candidate ancestor, then follow the chain until
        # no candidate ancestor remains.  This prevents vascular -> cardiovascular
        # from stopping before circulatory when both are nested candidates.
        nearest = min(ancestor_candidates, key=lambda item: (item[0], item[1]))[1]
        terminal[candidate] = nearest

    def resolve(candidate: str) -> str:
        seen: set[str] = set()
        node = candidate
        while terminal[node] != node and node not in seen:
            seen.add(node)
            node = terminal[node]
        return node

    terminal = {candidate: resolve(candidate) for candidate in candidates}
    return candidates, terminal


def display_name(ontology_name: str) -> str:
    # Presentation-only alias; ontology ID and original term remain in source data.
    return {"entire sense organ system": "sense organ system"}.get(ontology_name, ontology_name)


def fmt(value: int) -> str:
    return f"{int(value):,}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not SNAPSHOT.exists():
        raise FileNotFoundError(f"Missing frozen visualization snapshot: {SNAPSHOT}")

    terms, version = parse_uberon(SNAPSHOT)
    candidates, terminal = build_terminal_anchor_map(terms)
    snapshot_sha256 = hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest()

    relation = pd.read_csv(
        FORMAL / "protein_disease_relation_v72.tsv.gz", sep="\t", dtype=str, keep_default_na=False,
        usecols=["target_uniprot_id", "canonical_disease_id"],
    ).drop_duplicates()
    anatomy = pd.read_csv(
        FORMAL / "disease_anatomy_v72.tsv.gz", sep="\t", dtype=str, keep_default_na=False,
        usecols=["canonical_disease_id", "anatomy_id", "anatomy_name", "assertion_status", "source_disease_id",
                 "hierarchy_distance", "mapping_predicate", "axiom_type", "source_ontology", "uberon_source", "review_flag"],
    )
    anatomy = anatomy.drop_duplicates()

    projection_cache: dict[str, tuple[list[str], str]] = {}

    def project(anatomy_id: str) -> tuple[list[str], str]:
        if anatomy_id in projection_cache:
            return projection_cache[anatomy_id]
        if anatomy_id not in terms:
            projection_cache[anatomy_id] = ([], "not_in_uberon_snapshot")
            return projection_cache[anatomy_id]
        ancestors = ancestor_distances(anatomy_id, terms)
        hits: list[tuple[int, str]] = []
        if anatomy_id in candidates:
            hits.append((0, anatomy_id))
        hits.extend((distance, candidate) for candidate, distance in ancestors.items() if candidate in candidates)
        anchors = sorted({terminal[candidate] for _, candidate in hits})
        status = "mapped_to_terminal_top_level_system" if anchors else "no_top_level_system_ancestor"
        projection_cache[anatomy_id] = (anchors, status)
        return projection_cache[anatomy_id]

    mapping_rows: list[dict] = []
    for (anatomy_id, anatomy_name), _ in anatomy.groupby(["anatomy_id", "anatomy_name"], dropna=False):
        projected, status = project(anatomy_id)
        original_term = terms.get(anatomy_id, {}).get("name", anatomy_name)
        original_is_candidate = anatomy_id in candidates
        if projected:
            for target in projected:
                relation_type = "identity_top_level" if anatomy_id == target else "is_a_or_part_of_ancestor_projection"
                mapping_rows.append({
                    "original_uberon_id": anatomy_id,
                    "original_uberon_term": original_term,
                    "top_level_uberon_id": target,
                    "top_level_system": terms[target]["name"],
                    "display_system": display_name(terms[target]["name"]),
                    "mapping_relation": relation_type,
                    "mapping_source": f"Uberon {version}; is_a/part_of ancestor traversal",
                    "mapping_status": status,
                    "original_term_is_direct_system_candidate": original_is_candidate,
                    "top_level_candidate_was_collapsed": terminal.get(anatomy_id, anatomy_id) != anatomy_id,
                })
        else:
            mapping_rows.append({
                "original_uberon_id": anatomy_id,
                "original_uberon_term": original_term,
                "top_level_uberon_id": "",
                "top_level_system": "",
                "display_system": "",
                "mapping_relation": "",
                "mapping_source": f"Uberon {version}; is_a/part_of ancestor traversal",
                "mapping_status": status,
                "original_term_is_direct_system_candidate": original_is_candidate,
                "top_level_candidate_was_collapsed": False,
            })
    mapping = pd.DataFrame(mapping_rows).sort_values(["mapping_status", "top_level_system", "original_uberon_term"])
    mapping.to_csv(OUT / "D2_original_to_top_level_mapping.tsv", sep="\t", index=False)

    mapped = mapping[mapping["mapping_status"] == "mapped_to_terminal_top_level_system"][["original_uberon_id", "top_level_uberon_id", "top_level_system", "display_system"]]
    anatomy_projected = anatomy.merge(mapped, left_on="anatomy_id", right_on="original_uberon_id", how="inner")
    triples = relation.merge(
        anatomy_projected[["canonical_disease_id", "top_level_uberon_id", "top_level_system", "display_system"]].drop_duplicates(),
        on="canonical_disease_id", how="inner",
    ).drop_duplicates(["target_uniprot_id", "canonical_disease_id", "top_level_uberon_id"])

    collapsed_source_terms = (
        mapping[mapping["original_uberon_id"] != mapping["top_level_uberon_id"]]
        .groupby("top_level_uberon_id")["original_uberon_id"].nunique().rename("n_original_anatomical_terms_collapsed")
    )
    summary = (
        triples.groupby(["top_level_uberon_id", "top_level_system", "display_system"], as_index=False)
        .agg(
            n_unique_disease_associated_membrane_proteins=("target_uniprot_id", "nunique"),
            n_unique_diseases=("canonical_disease_id", "nunique"),
            n_unique_protein_disease_system_relations=("canonical_disease_id", "size"),
        )
        .merge(collapsed_source_terms, on="top_level_uberon_id", how="left")
    )
    summary["n_original_anatomical_terms_collapsed"] = summary["n_original_anatomical_terms_collapsed"].fillna(0).astype(int)
    summary["has_descendant_terms_collapsed"] = summary["n_original_anatomical_terms_collapsed"] > 0
    summary = summary.sort_values("n_unique_disease_associated_membrane_proteins", ascending=False).reset_index(drop=True)
    summary.to_csv(OUT / "D2_top_level_anatomical_system_summary.tsv", sep="\t", index=False)

    unresolved = mapping[mapping["mapping_status"] != "mapped_to_terminal_top_level_system"].copy()
    unresolved_counts = anatomy.merge(
        unresolved[["original_uberon_id", "mapping_status"]], left_on="anatomy_id", right_on="original_uberon_id", how="inner"
    ).merge(relation, on="canonical_disease_id", how="left")
    unresolved_summary = (
        # ``anatomy_name`` is the formal V7.2 term label; it remains the source
        # label for unresolved terms even when the downloaded ontology omits it.
        unresolved_counts.groupby(["original_uberon_id", "anatomy_name", "mapping_status"], as_index=False)
        .agg(
            n_formal_anatomy_assertions=("canonical_disease_id", "size"),
            n_unique_diseases=("canonical_disease_id", "nunique"),
            n_unique_disease_associated_membrane_proteins=("target_uniprot_id", "nunique"),
        )
        .rename(columns={"anatomy_name": "original_uberon_term"})
        .sort_values(["n_unique_disease_associated_membrane_proteins", "n_unique_diseases"], ascending=False)
    )
    unresolved_summary.to_csv(OUT / "D2_unresolved_anatomical_terms.tsv", sep="\t", index=False)

    # Multi-label terms must not be summed and presented as a disease total.
    # Keep both the union and the truly display-excluded subset explicit.
    unresolved_disease_ids = set(unresolved_counts["canonical_disease_id"])
    unresolved_protein_ids = set(unresolved_counts.loc[unresolved_counts["target_uniprot_id"].notna(), "target_uniprot_id"])
    retained_disease_ids = set(anatomy_projected["canonical_disease_id"])
    retained_protein_ids = set(
        relation.loc[relation["canonical_disease_id"].isin(retained_disease_ids), "target_uniprot_id"]
    )
    unresolved_disease_only = unresolved_disease_ids - retained_disease_ids
    unresolved_protein_only = unresolved_protein_ids - retained_protein_ids

    # Figure: all nonzero terminal systems remain visible; no display threshold is applied.
    figure_data = summary[summary["n_unique_disease_associated_membrane_proteins"] > 0].sort_values(
        "n_unique_disease_associated_membrane_proteins", ascending=True
    )
    fig_height = max(5.3, 0.44 * len(figure_data) + 1.9)
    fig, ax = plt.subplots(figsize=(7.0, fig_height))
    burgundy = "#B76E79"
    bars = ax.barh(
        figure_data["display_system"], figure_data["n_unique_disease_associated_membrane_proteins"],
        color=burgundy, edgecolor="#934F5A", linewidth=0.55, height=0.64,
    )
    xmax = max(figure_data["n_unique_disease_associated_membrane_proteins"]) * 1.36
    ax.set_xlim(0, xmax)
    ax.set_xlabel("Unique canonical disease-associated membrane proteins", fontsize=9, color="#24364B", labelpad=7)
    ax.set_ylabel("")
    ax.xaxis.grid(True, color="#E3E8EC", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#748392")
    ax.tick_params(axis="y", length=0, labelsize=8.5, colors="#24364B")
    ax.tick_params(axis="x", labelsize=8, colors="#506174")
    for bar, (_, row) in zip(bars, figure_data.iterrows()):
        label = f"{fmt(row['n_unique_disease_associated_membrane_proteins'])} proteins · {fmt(row['n_unique_diseases'])} diseases"
        ax.text(bar.get_width() + xmax * 0.012, bar.get_y() + bar.get_height() / 2, label,
                va="center", ha="left", fontsize=8.1, color="#3E5165")
    fig.suptitle("Anatomical-system coverage of disease-associated membrane proteins",
                 x=0.125, y=0.985, ha="left", fontsize=11.2, fontweight="bold", color="#24364B")
    ax.set_title("Bar length indicates unique canonical membrane proteins; labels report unique mapped diseases.",
                 loc="left", fontsize=8.2, color="#5E7182", pad=12)
    unresolved_terms = unresolved_summary["original_uberon_id"].nunique()
    fig.text(0.125, 0.012,
             f"Uberon {version} projection; parent/child system candidates are not co-displayed. "
             f"Diseases may map to multiple systems. {unresolved_terms} anatomy terms; {fmt(len(unresolved_disease_only))} diseases and {fmt(len(unresolved_protein_only))} proteins occur only in unresolved anatomy contexts and are reported separately.",
             ha="left", va="bottom", fontsize=6.8, color="#5E7182")
    fig.subplots_adjust(left=0.33, right=0.95, top=0.89, bottom=0.10)
    base = OUT / "D2_top_level_anatomical_system_coverage"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)

    # QC and field provenance.
    displayed_ids = set(figure_data["top_level_uberon_id"])
    parent_child_pairs = []
    for system_id in displayed_ids:
        for ancestor in ancestor_distances(system_id, terms):
            if ancestor in displayed_ids:
                parent_child_pairs.append((system_id, ancestor))
    source_triples = relation.merge(anatomy[["canonical_disease_id", "anatomy_id"]].drop_duplicates(), on="canonical_disease_id", how="inner")
    qc_rows = [
        ("formal_protein_disease_relations", len(relation), "unique canonical protein × disease before anatomy mapping", "PASS"),
        ("formal_disease_anatomy_assertions", len(anatomy), "read-only V7.2 source assertions", "PASS"),
        ("original_unique_anatomy_terms", anatomy["anatomy_id"].nunique(), "terms in formal disease-anatomy table", "PASS"),
        ("top_level_systems_with_nonzero_coverage", len(figure_data), "terminal direct-system candidates after ontology projection", "PASS"),
        ("direct_system_candidates_collapsed", sum(k != v for k, v in terminal.items()), "collapsed only by Uberon is_a/part_of ancestry", "PASS"),
        ("displayed_parent_child_pairs", len(parent_child_pairs), "must equal zero", "PASS" if not parent_child_pairs else "FAIL"),
        ("duplicate_protein_disease_system_triples", int(triples.duplicated(["target_uniprot_id", "canonical_disease_id", "top_level_uberon_id"]).sum()), "must equal zero", "PASS"),
        ("unresolved_original_anatomy_terms", int(unresolved_summary["original_uberon_id"].nunique()), "not forced into a system", "ACCEPTABLE_AS_UNRESOLVED"),
        ("unresolved_disease_union", len(unresolved_disease_ids), "distinct diseases with at least one unresolved term", "ACCEPTABLE_AS_UNRESOLVED"),
        ("unresolved_protein_union", len(unresolved_protein_ids), "distinct proteins linked to a disease with at least one unresolved term", "ACCEPTABLE_AS_UNRESOLVED"),
        ("unresolved_diseases_only", len(unresolved_disease_only), "distinct diseases with no retained system mapping", "ACCEPTABLE_AS_UNRESOLVED"),
        ("unresolved_proteins_only", len(unresolved_protein_only), "distinct proteins not represented by any retained system", "ACCEPTABLE_AS_UNRESOLVED"),
        ("uses_disease_name_keyword_mapping", 0, "ontology IDs and graph only", "PASS"),
        ("figure_numbers_match_summary", True, "labels created directly from summary table", "PASS"),
        ("uberon_snapshot_data_version", version, "downloaded and SHA-256 frozen for visualization", "PASS"),
        ("uberon_snapshot_sha256", snapshot_sha256, "downloaded and SHA-256 frozen for visualization", "PASS"),
    ]
    with (OUT / "D2_QC.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["check", "value", "definition", "status"])
        writer.writerows(qc_rows)

    with (OUT / "D2_fields_used.txt").open("w", encoding="utf-8") as handle:
        handle.write("D2 top-level anatomical-system coverage: fields and derivation\n\n")
        handle.write("Protein–disease relation (read-only formal V7.2): protein_disease_relation_v72.tsv.gz\n")
        handle.write("- target_uniprot_id\n- canonical_disease_id\n\n")
        handle.write("Disease–anatomy assertion (read-only formal V7.2): disease_anatomy_v72.tsv.gz\n")
        handle.write("- canonical_disease_id\n- anatomy_id\n- anatomy_name\n- assertion_status\n- source_disease_id\n- hierarchy_distance\n- mapping_predicate\n- axiom_type\n- source_ontology\n- uberon_source\n\n")
        handle.write(f"Visualization-only ontology snapshot: {SNAPSHOT.name}\n")
        handle.write(f"- data-version: {version}\n- SHA-256: {snapshot_sha256}\n")
        handle.write("Projection: anatomy_id → all terminal direct is_a children of UBERON:0000467 reachable through is_a/part_of ancestors.\n")
        handle.write("Unit: unique canonical protein × canonical disease × harmonized top-level anatomical system.\n")
        handle.write("No disease names, keywords, or manual anatomical merges were used.\n")


if __name__ == "__main__":
    main()
