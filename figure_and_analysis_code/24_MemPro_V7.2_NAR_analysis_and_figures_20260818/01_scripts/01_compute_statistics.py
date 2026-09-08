from __future__ import annotations

import itertools
import json
import math
import os
import random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, kruskal, mannwhitneyu


SOURCE = Path(os.environ["MEMPRO_DATA_ROOT"]) / "01_core_v72" / "01_release_tables"
ROOT = Path(os.environ["MEMPRO_FIGURE_OUTPUT_ROOT"])
ANALYSIS = ROOT / "02_analysis_data"
SD = ROOT / "03_source_data"
QA = ROOT / "08_QA"
SEED = 42
random.seed(SEED)
np.random.seed(SEED)

PLACEHOLDERS = {
    "unclassified",
    "no_specialist_classification",
    "family_defined_membrane_role_unresolved",
    "",
}


def load(name: str, directory: Path = ANALYSIS, **kwargs) -> pd.DataFrame:
    return pd.read_csv(directory / name, sep="\t", compression="infer", dtype=str, **kwargs)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def split_values(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    return [x.strip() for x in str(value).replace("|", ";").split(";") if x.strip()]


def bh_fdr(pvalues: list[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


def gini(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x) & (x >= 0)]
    if len(x) == 0 or x.sum() == 0:
        return float("nan")
    x = np.sort(x)
    n = len(x)
    return float((2 * np.sum((np.arange(1, n + 1)) * x) / (n * x.sum())) - (n + 1) / n)


def contingency_stats(df: pd.DataFrame, row: str, col: str, prefix: str) -> dict:
    table = pd.crosstab(df[row], df[col])
    chi2, p, dof, expected = chi2_contingency(table.values)
    n = table.values.sum()
    phi2 = chi2 / n
    r, k = table.shape
    corrected_phi2 = max(0.0, phi2 - ((k - 1) * (r - 1)) / max(n - 1, 1))
    corrected_r = r - ((r - 1) ** 2) / max(n - 1, 1)
    corrected_k = k - ((k - 1) ** 2) / max(n - 1, 1)
    cramer_v = math.sqrt(corrected_phi2 / max(min(corrected_k - 1, corrected_r - 1), 1e-12))
    residual = (table.values - expected) / np.sqrt(np.maximum(expected, 1e-12))

    cell_rows = []
    pvals = []
    for i, rv in enumerate(table.index):
        for j, cv in enumerate(table.columns):
            a = int(table.iloc[i, j])
            b = int(table.iloc[i, :].sum() - a)
            c = int(table.iloc[:, j].sum() - a)
            d = int(n - a - b - c)
            odds, fp = fisher_exact([[a, b], [c, d]])
            pvals.append(fp)
            cell_rows.append({
                row: rv, col: cv, "observed": a, "expected": expected[i, j],
                "standardized_residual": residual[i, j], "odds_ratio": odds,
                "log2_odds_ratio": math.log2(odds) if odds > 0 and math.isfinite(odds) else np.nan,
                "fisher_p": fp,
            })
    qvals = bh_fdr(pvals)
    for rec, q in zip(cell_rows, qvals):
        rec["fisher_bh_fdr"] = q
    out = pd.DataFrame(cell_rows)
    out.to_csv(SD / f"{prefix}_cell_statistics.tsv", sep="\t", index=False)
    table.to_csv(SD / f"{prefix}_contingency.tsv", sep="\t")
    return {"chi2": chi2, "p": p, "dof": int(dof), "n": int(n), "cramers_v_bias_corrected": cramer_v}


def protein_statistics(protein: pd.DataFrame) -> dict:
    pd.crosstab(protein["membrane_class_v7"], protein["membrane_evidence_level_v7"]).to_csv(SD / "F2A_membrane_class_by_evidence.tsv", sep="\t")

    base = protein.copy()
    for col in ["primary_membrane_role", "primary_molecular_function", "primary_biological_process"]:
        base[col] = base[col].fillna("unclassified")
    role_function = contingency_stats(base, "primary_membrane_role", "primary_molecular_function", "F2B_role_function")
    role_process = contingency_stats(base, "primary_membrane_role", "primary_biological_process", "F2B_role_process")

    # Compact alluvial edges use top-level primary labels only; no record is lost,
    # because rare labels are explicitly grouped as Other.
    def aggregate(values: pd.Series, n: int) -> pd.Series:
        top = set(values.value_counts().head(n).index)
        return values.where(values.isin(top), "Other")

    flow = pd.DataFrame({
        "membrane_class": base["membrane_class_v7"],
        "membrane_role": aggregate(base["primary_membrane_role"], 10),
        "molecular_function": aggregate(base["primary_molecular_function"], 8),
        "biological_process": aggregate(base["primary_biological_process"], 10),
    })
    edges = []
    for left, right in [("membrane_class", "membrane_role"), ("membrane_role", "molecular_function"), ("molecular_function", "biological_process")]:
        counts = flow.groupby([left, right]).size().reset_index(name="protein_count")
        for rec in counts.to_dict("records"):
            edges.append({"source_axis": left, "source_label": rec[left], "target_axis": right, "target_label": rec[right], "protein_count": rec["protein_count"]})
    pd.DataFrame(edges).to_csv(SD / "F2C_five_axis_alluvial_edges.tsv", sep="\t", index=False)
    return {"role_function": role_function, "role_process": role_process}


def tau_statistics(tau: pd.DataFrame) -> dict:
    tau["tau_consensus_tissue_rna"] = numeric(tau["tau_consensus_tissue_rna"])
    tau = tau.dropna(subset=["tau_consensus_tissue_rna", "primary_membrane_role"])
    counts = tau["primary_membrane_role"].value_counts()
    roles = list(counts[counts >= 30].index)
    groups = [tau.loc[tau["primary_membrane_role"].eq(role), "tau_consensus_tissue_rna"].to_numpy() for role in roles]
    kw = kruskal(*groups) if len(groups) > 1 else None
    comparisons = []
    pvals = []
    for a, b in itertools.combinations(roles, 2):
        xa = tau.loc[tau["primary_membrane_role"].eq(a), "tau_consensus_tissue_rna"].to_numpy()
        xb = tau.loc[tau["primary_membrane_role"].eq(b), "tau_consensus_tissue_rna"].to_numpy()
        stat, p = mannwhitneyu(xa, xb, alternative="two-sided")
        pvals.append(p)
        comparisons.append({"role_a": a, "role_b": b, "n_a": len(xa), "n_b": len(xb), "median_a": np.median(xa), "median_b": np.median(xb), "mannwhitney_u": stat, "p": p})
    if comparisons:
        for row, q in zip(comparisons, bh_fdr(pvals)):
            row["bh_fdr"] = q
    pd.DataFrame(comparisons).to_csv(SD / "F2D_tau_pairwise_statistics.tsv", sep="\t", index=False)
    summary = tau.groupby("primary_membrane_role")["tau_consensus_tissue_rna"].agg(["count", "median", "mean", "std", "min", "max"]).reset_index()
    summary.to_csv(SD / "F2D_tau_role_summary.tsv", sep="\t", index=False)
    return {"kruskal_h": float(kw.statistic) if kw else None, "kruskal_p": float(kw.pvalue) if kw else None, "groups_n_ge_30": len(roles)}


def compound_and_pair_statistics(protein: pd.DataFrame, pairs: pd.DataFrame) -> dict:
    role_map = protein.set_index("canonical_uniprot_accession")["primary_membrane_role"].fillna("unclassified").to_dict()
    pairs["membrane_role"] = pairs["target_uniprot_id"].map(role_map).fillna("unclassified")
    pairs["evidence_count"] = numeric(pairs["evidence_count"]).fillna(0)
    for col in ["distinct_database_count", "independent_experiment_count", "independent_structure_count", "independent_pubchem_assay_count", "independent_evidence_modality_count"]:
        pairs[col] = numeric(pairs[col]).fillna(0).astype(int)

    # Compound target breadth and target-family entropy.
    role_counts = pairs.groupby(["compound_internal_id", "membrane_role"]).size().rename("target_role_count").reset_index()
    rows = []
    for compound_id, grp in role_counts.groupby("compound_internal_id", sort=False):
        counts = grp["target_role_count"].to_numpy(float)
        p = counts / counts.sum()
        entropy = float(-(p * np.log(p)).sum())
        norm = float(entropy / np.log(len(p))) if len(p) > 1 else 0.0
        rows.append({"compound_internal_id": compound_id, "target_count": int(counts.sum()), "target_role_count": len(counts), "target_role_entropy": entropy, "target_role_entropy_normalized": norm})
    entropy_df = pd.DataFrame(rows)
    entropy_df.to_csv(ANALYSIS / "compound_target_entropy.tsv.gz", sep="\t", index=False, compression="gzip")
    entropy_df.to_csv(SD / "F3C_compound_target_entropy.tsv", sep="\t", index=False)

    # Load only interaction-linked compounds for property and preference analyses.
    chunks = []
    for chunk in pd.read_csv(ANALYSIS / "compound_features.tsv.gz", sep="\t", compression="gzip", dtype=str, chunksize=100_000):
        take = chunk[chunk["interaction_linked"].eq("1")].copy()
        if not take.empty:
            chunks.append(take)
    compounds = pd.concat(chunks, ignore_index=True)
    compounds = compounds.merge(entropy_df, on="compound_internal_id", how="left", suffixes=("", "_derived"))

    statuses = [
        ("Approved drug", "is_approved_drug"),
        ("Clinical candidate", "is_clinical_candidate"),
        ("Endogenous ligand", "is_endogenous_ligand"),
        ("Natural product", "is_natural_product"),
        ("Chemical probe", "is_chemical_probe"),
    ]
    property_rows = []
    for label, flag in statuses:
        subset = compounds[compounds[flag].eq("1")]
        for prop in ["molecular_weight", "xlogp", "tpsa", "rotatable_bond_count"]:
            vals = numeric(subset[prop]).dropna()
            for value in vals:
                property_rows.append({"biological_status": label, "property": prop, "value": value})
    # Deterministic background sample prevents an uninformative 240k-row file.
    background = compounds.sample(min(20_000, len(compounds)), random_state=SEED)
    for prop in ["molecular_weight", "xlogp", "tpsa", "rotatable_bond_count"]:
        for value in numeric(background[prop]).dropna():
            property_rows.append({"biological_status": "All interaction-linked", "property": prop, "value": value})
    pd.DataFrame(property_rows).to_csv(SD / "F3B_physicochemical_by_status.tsv", sep="\t", index=False)

    compound_class = compounds.set_index("compound_internal_id")["broad_structure_class_recomputed"].to_dict()
    pref = pairs[["compound_internal_id", "target_uniprot_id", "membrane_role"]].drop_duplicates()
    pref["broad_structure_class"] = pref["compound_internal_id"].map(compound_class)
    pref = pref.dropna(subset=["broad_structure_class"])
    class_role = contingency_stats(pref, "broad_structure_class", "membrane_role", "F3D_structure_role")

    # Evidence architecture.
    source_sets = defaultdict(set)
    pair_sources = {}
    for row in pairs[["pair_id", "source_databases"]].itertuples(index=False):
        values = split_values(row.source_databases)
        pair_sources[row.pair_id] = values
        for src in values:
            source_sets[src].add(row.pair_id)
    top_sources = [x for x, _ in sorted(((s, len(v)) for s, v in source_sets.items()), key=lambda x: x[1], reverse=True)[:6]]
    intersections = Counter()
    for pid, sources in pair_sources.items():
        key = tuple(src for src in top_sources if src in sources)
        intersections[key if key else ("Other sources only",)] += 1
    inter_rows = []
    for key, count in intersections.most_common():
        inter_rows.append({"intersection": ";".join(key), "pair_count": count, **{src: int(src in key) for src in top_sources}})
    pd.DataFrame(inter_rows).to_csv(SD / "F4A_source_upset.tsv", sep="\t", index=False)

    jac_rows = []
    for a, b in itertools.product(top_sources, repeat=2):
        union = source_sets[a] | source_sets[b]
        jac_rows.append({"source_a": a, "source_b": b, "jaccard": len(source_sets[a] & source_sets[b]) / len(union) if union else np.nan, "intersection": len(source_sets[a] & source_sets[b]), "union": len(union)})
    pd.DataFrame(jac_rows).to_csv(SD / "F4B_source_jaccard.tsv", sep="\t", index=False)

    dbexp = pairs.groupby(["distinct_database_count", "independent_experiment_count"]).size().reset_index(name="pair_count")
    dbexp.to_csv(SD / "F4C_database_vs_putative_experiment.tsv", sep="\t", index=False)

    vals = np.sort(pairs["evidence_count"].to_numpy(float))
    cumulative = np.cumsum(vals) / vals.sum()
    lorenz = pd.DataFrame({"cumulative_pair_fraction": np.arange(1, len(vals) + 1) / len(vals), "cumulative_evidence_fraction": cumulative})
    lorenz = pd.concat([pd.DataFrame({"cumulative_pair_fraction": [0.0], "cumulative_evidence_fraction": [0.0]}), lorenz], ignore_index=True)
    lorenz.to_csv(SD / "F4D_evidence_lorenz.tsv", sep="\t", index=False)
    return {"class_role": class_role, "pair_evidence_gini": gini(vals), "top_sources": top_sources, "interaction_linked_compounds": len(compounds)}


def site_statistics() -> dict:
    sites = load("site_analysis.tsv.gz")
    sites["mapping_fraction"] = numeric(sites["mapping_fraction"]).fillna(0).clip(0, 1)
    ecdf_rows = []
    for tier, grp in sites.groupby("site_tier_v72"):
        values = np.sort(grp["mapping_fraction"].to_numpy())
        for x, y in zip(values, np.arange(1, len(values) + 1) / len(values)):
            ecdf_rows.append({"site_tier": tier, "mapping_fraction": x, "ecdf": y})
    pd.DataFrame(ecdf_rows).to_csv(SD / "F4E_site_mapping_ecdf.tsv", sep="\t", index=False)
    sites["membrane_side"].fillna("unknown").value_counts().rename_axis("membrane_side").reset_index(name="site_count").to_csv(SD / "F4F_membrane_side.tsv", sep="\t", index=False)
    sites["site_tier_v72"].value_counts().rename_axis("site_tier").reset_index(name="site_count").to_csv(SD / "S_site_tier_counts.tsv", sep="\t", index=False)
    return {"site_rows": len(sites), "unknown_membrane_side": int(sites["membrane_side"].fillna("unknown").eq("unknown").sum())}


def disease_statistics(disease: pd.DataFrame) -> dict:
    disease["best_evidence_level"].value_counts().rename_axis("evidence_level").reset_index(name="relation_count").to_csv(SD / "F5D_disease_evidence_levels.tsv", sep="\t", index=False)

    disease_sets = disease[["canonical_disease_id", "anatomical_systems", "therapeutic_areas"]].drop_duplicates("canonical_disease_id")
    fractional = Counter()
    raw = Counter()
    for row in disease_sets.itertuples(index=False):
        systems = split_values(row.anatomical_systems)
        if not systems:
            systems = ["Unmapped"]
        for system in systems:
            raw[system] += 1
            fractional[system] += 1 / len(systems)
    pd.DataFrame([{"anatomical_system": s, "raw_multilabel_disease_count": raw[s], "fractional_disease_count": fractional[s]} for s in sorted(raw)]).to_csv(SD / "F5A_anatomy_fractional_counts.tsv", sep="\t", index=False)

    # Disease-level therapeutic-area x anatomy associations.
    ta_anat_rows = []
    for row in disease_sets.itertuples(index=False):
        systems = split_values(row.anatomical_systems)
        areas = split_values(row.therapeutic_areas)
        for area in areas:
            for system in systems:
                ta_anat_rows.append({"canonical_disease_id": row.canonical_disease_id, "therapeutic_area": area, "anatomical_system": system})
    ta_anat = pd.DataFrame(ta_anat_rows).drop_duplicates()
    ta_anat_stats = contingency_stats(ta_anat, "therapeutic_area", "anatomical_system", "F5B_therapeutic_anatomy") if not ta_anat.empty else {}

    channels = ["human_genetic", "clinical_genetic", "somatic_mutation", "functional", "animal_model", "expression", "literature", "uniprot_disease"]
    channel_counter = Counter(tuple(x for x in channels if x in set(split_values(v))) for v in disease["evidence_channels"])
    top = channel_counter.most_common(25)
    upset = pd.DataFrame([{"intersection": ";".join(key) if key else "none", "relation_count": count, **{c: int(c in key) for c in channels}} for key, count in top])
    upset.to_csv(SD / "F5C_disease_channel_upset.tsv", sep="\t", index=False)

    # Sparse shared-target disease network.
    target_sets = disease.groupby("canonical_disease_id")["target_uniprot_id"].agg(set).to_dict()
    inverted = defaultdict(list)
    for did, targets in target_sets.items():
        if len(targets) >= 2:
            for target in targets:
                inverted[target].append(did)
    intersections = Counter()
    for diseases in inverted.values():
        if len(diseases) > 80:
            continue
        for a, b in itertools.combinations(sorted(set(diseases)), 2):
            intersections[(a, b)] += 1
    edges = []
    for (a, b), inter in intersections.items():
        if inter < 2:
            continue
        union = len(target_sets[a] | target_sets[b])
        j = inter / union
        if j >= 0.20:
            edges.append({"disease_a": a, "disease_b": b, "shared_targets": inter, "jaccard": j, "degree_a": len(target_sets[a]), "degree_b": len(target_sets[b])})
    pd.DataFrame(sorted(edges, key=lambda x: (x["jaccard"], x["shared_targets"]), reverse=True)[:500]).to_csv(SD / "F5E_disease_shared_target_network.tsv", sep="\t", index=False)
    return {"therapeutic_anatomy": ta_anat_stats, "network_edges_retained": min(len(edges), 500)}


SYSTEM_TISSUE_KEYWORDS = {
    "nervous system": ["brain", "cerebell", "spinal cord"],
    "musculoskeletal system": ["skeletal muscle"],
    "entire sense organ system": ["retina", "eye"],
    "circulatory system": ["heart"],
    "cardiovascular system": ["heart"],
    "vascular system": ["artery", "heart"],
    "integumental system": ["skin"],
    "digestive system": ["liver", "stomach", "colon", "small intestine", "duodenum", "rectum", "pancreas", "gallbladder"],
    "alimentary part of gastrointestinal system": ["stomach", "colon", "small intestine", "duodenum", "rectum"],
    "excretory system": ["kidney", "urinary bladder"],
    "reproductive system": ["testis", "ovary", "endometrium", "cervix", "prostate", "fallopian"],
    "hematopoietic system": ["bone marrow", "spleen"],
    "lymphoid system": ["lymph", "tonsil", "spleen"],
    "respiratory system": ["lung", "bronch"],
    "exocrine system": ["salivary", "pancreas"],
    "lacrimal apparatus": ["lacrimal"],
    "neuroendocrine system": ["pituitary", "adrenal"],
}


def expression_concordance(disease: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    cols = ["target_uniprot_id", "source_dataset", "tissue", "value", "mapped_measured_denominator_eligible_v71", "measurement_status_v71"]
    chunks = []
    for chunk in pd.read_csv(SOURCE / "expression_measurement_v72.tsv.gz", sep="\t", compression="gzip", dtype=str, usecols=cols, chunksize=250_000):
        take = chunk[
            chunk["source_dataset"].eq("HPA_consensus_tissue_RNA")
            & chunk["mapped_measured_denominator_eligible_v71"].eq("1")
            & chunk["measurement_status_v71"].eq("MEASURED")
            & chunk["tissue"].notna()
        ].copy()
        take["value"] = numeric(take["value"])
        chunks.append(take.dropna(subset=["value"])[["target_uniprot_id", "tissue", "value"]])
    expr = pd.concat(chunks, ignore_index=True).groupby(["target_uniprot_id", "tissue"], as_index=False)["value"].median()
    expr_map = {pid: dict(zip(g["tissue"].str.lower(), g["value"])) for pid, g in expr.groupby("target_uniprot_id")}

    disease_systems = disease[["canonical_disease_id", "anatomical_systems"]].drop_duplicates().set_index("canonical_disease_id")["anatomical_systems"].map(split_values).to_dict()
    edges = disease[["target_uniprot_id", "canonical_disease_id"]].drop_duplicates()

    def score(protein_id: str, disease_id: str) -> float:
        values = expr_map.get(protein_id)
        systems = disease_systems.get(disease_id, [])
        if not values or not systems:
            return np.nan
        maximum = max(values.values()) if values else 0
        if maximum <= 0:
            return np.nan
        matched = []
        for system in systems:
            kws = SYSTEM_TISSUE_KEYWORDS.get(system.lower(), [])
            for tissue, value in values.items():
                if any(k in tissue for k in kws):
                    matched.append(value)
        return float(max(matched) / maximum) if matched else np.nan

    observed = edges.copy()
    observed["anatomical_expression_concordance"] = [score(p, d) for p, d in observed.itertuples(index=False)]
    observed.to_csv(ANALYSIS / "disease_expression_concordance.tsv.gz", sep="\t", index=False, compression="gzip")
    observed.to_csv(SD / "F6B_disease_expression_concordance_observed.tsv", sep="\t", index=False)

    # Degree-preserving bipartite edge swaps. Duplicate edges are rejected.
    edge_list = list(map(tuple, edges.to_numpy()))
    edge_set = set(edge_list)
    null_means = []
    rng = random.Random(SEED)
    for perm in range(200):
        work = list(edge_list)
        work_set = set(edge_set)
        attempts = min(50_000, len(work) * 8)
        for _ in range(attempts):
            i, j = rng.sample(range(len(work)), 2)
            p1, d1 = work[i]
            p2, d2 = work[j]
            if p1 == p2 or d1 == d2:
                continue
            n1, n2 = (p1, d2), (p2, d1)
            if n1 in work_set or n2 in work_set:
                continue
            work_set.remove(work[i]); work_set.remove(work[j])
            work[i], work[j] = n1, n2
            work_set.add(n1); work_set.add(n2)
        values = np.asarray([score(p, d) for p, d in work], dtype=float)
        null_means.append(float(np.nanmean(values)))
    obs_mean = float(observed["anatomical_expression_concordance"].mean())
    null = np.asarray(null_means)
    p_empirical = float((1 + np.sum(null >= obs_mean)) / (len(null) + 1))
    z = float((obs_mean - null.mean()) / null.std(ddof=1)) if null.std(ddof=1) > 0 else np.nan
    pd.DataFrame({"permutation": np.arange(1, len(null) + 1), "null_mean_concordance": null}).to_csv(SD / "F6B_concordance_permutation_null.tsv", sep="\t", index=False)
    result = {"observed_mean": obs_mean, "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)), "z_score": z, "empirical_p": p_empirical, "permutations": len(null), "observed_edges_with_score": int(observed["anatomical_expression_concordance"].notna().sum())}
    return observed, result


def integrated_statistics(protein: pd.DataFrame, pairs: pd.DataFrame, disease: pd.DataFrame, concordance: pd.DataFrame) -> dict:
    metrics = protein[["canonical_uniprot_accession", "approved_symbol", "membrane_class_v7", "membrane_evidence_level_v7", "primary_membrane_role"]].rename(columns={"canonical_uniprot_accession": "target_uniprot_id"})

    chem = pairs.groupby("target_uniprot_id").agg(
        chemical_pair_count=("pair_id", "nunique"),
        compound_count=("compound_internal_id", "nunique"),
        BE1_pair_count=("best_evidence_tier", lambda x: int((x == "BE1").sum())),
        mean_contributing_database_count=("distinct_database_count", "mean"),
    ).reset_index()
    dis = disease.groupby("target_uniprot_id").agg(
        disease_count=("canonical_disease_id", "nunique"),
        very_high_disease_relations=("best_evidence_level", lambda x: int((x == "very_high").sum())),
        high_disease_relations=("best_evidence_level", lambda x: int((x == "high").sum())),
        medium_disease_relations=("best_evidence_level", lambda x: int((x == "medium").sum())),
        max_open_targets_score=("max_ot_overall_score", lambda x: numeric(x).max()),
    ).reset_index()
    dis["transparent_disease_evidence_index"] = 3 * dis["very_high_disease_relations"] + 2 * dis["high_disease_relations"] + dis["medium_disease_relations"]
    site = load("site_analysis.tsv.gz")
    site["coordinate_docking_eligible"] = numeric(site["coordinate_docking_eligible"]).fillna(0)
    siteagg = site.groupby("target_uniprot_id").agg(site_assertion_count=("site_record_id", "nunique"), coordinate_ready_site_count=("coordinate_docking_eligible", "sum")).reset_index()
    conc = concordance.groupby("target_uniprot_id")["anatomical_expression_concordance"].mean().reset_index(name="mean_anatomical_expression_concordance")
    metrics = metrics.merge(chem, on="target_uniprot_id", how="left").merge(dis, on="target_uniprot_id", how="left").merge(siteagg, on="target_uniprot_id", how="left").merge(conc, on="target_uniprot_id", how="left")
    for col in ["chemical_pair_count", "compound_count", "BE1_pair_count", "disease_count", "very_high_disease_relations", "high_disease_relations", "medium_disease_relations", "transparent_disease_evidence_index", "site_assertion_count", "coordinate_ready_site_count"]:
        metrics[col] = numeric(metrics[col]).fillna(0)
    metrics["log10_chemical_pair_count_plus1"] = np.log10(metrics["chemical_pair_count"] + 1)
    metrics["log10_disease_evidence_index_plus1"] = np.log10(metrics["transparent_disease_evidence_index"] + 1)

    disease_positive = metrics[metrics["disease_count"] > 0]
    disease_threshold = disease_positive["transparent_disease_evidence_index"].quantile(0.75)
    # A pair-level site can only exist after at least one chemical relation is
    # known. Therefore the chemically underexplored threshold is estimated among
    # disease-linked proteins with >=1 formal pair; proteins with zero pairs remain
    # a separate structural-darkness group rather than being forced into this set.
    chemically_observed = disease_positive[disease_positive["chemical_pair_count"] > 0]
    chem_threshold = chemically_observed["chemical_pair_count"].quantile(0.25)
    metrics["hypothesis_generation_candidate"] = (
        (metrics["transparent_disease_evidence_index"] >= disease_threshold)
        & (metrics["chemical_pair_count"] > 0)
        & (metrics["chemical_pair_count"] <= chem_threshold)
        & (metrics["coordinate_ready_site_count"] > 0)
    ).astype(int)
    metrics.to_csv(ANALYSIS / "integrated_target_landscape.tsv.gz", sep="\t", index=False, compression="gzip")
    metrics.to_csv(SD / "F6A_integrated_target_landscape.tsv", sep="\t", index=False)
    metrics[metrics["hypothesis_generation_candidate"].eq(1)].sort_values(["transparent_disease_evidence_index", "chemical_pair_count"], ascending=[False, True]).head(50).to_csv(SD / "F6D_hypothesis_generation_candidates.tsv", sep="\t", index=False)
    return {"disease_evidence_q75": float(disease_threshold), "chemical_pair_q25_among_nonzero": float(chem_threshold), "candidate_count": int(metrics["hypothesis_generation_candidate"].sum())}


def main() -> None:
    protein = load("protein_analysis.tsv.gz")
    pairs = load("pair_analysis.tsv.gz")
    tau = load("expression_tau.tsv.gz")
    disease = load("disease_relation_analysis.tsv.gz")

    results = {
        "protein": protein_statistics(protein),
        "tau": tau_statistics(tau),
        "compound_pair_evidence": compound_and_pair_statistics(protein, pairs),
        "site": site_statistics(),
        "disease": disease_statistics(disease),
    }
    concordance, concordance_result = expression_concordance(disease)
    results["anatomical_expression_concordance"] = concordance_result
    results["integrated"] = integrated_statistics(protein, pairs, disease, concordance)

    (ANALYSIS / "analysis_results.json").write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False))


if __name__ == "__main__":
    main()
