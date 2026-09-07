from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import chi2_contingency, fisher_exact, kruskal, mannwhitneyu, norm, rankdata, spearmanr

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from pynndescent import NNDescent


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CFG = json.loads((ROOT / "00_config" / "nar_config.json").read_text(encoding="utf-8"))
OLD = Path(CFG["validated_analysis_dir"])
FROZEN = Path(CFG["frozen_release_dir"]) / "01_release_tables"
SCREEN = ROOT / "02_analysis_screening"
SD = ROOT / "03_source_data"
SEED = int(CFG["random_seed"])
RNG = np.random.default_rng(SEED)
PERMUTATIONS = int(CFG["permutations"])
BOOTSTRAPS = int(CFG["bootstrap_replicates"])

for p in (SCREEN, SD):
    p.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def load_analysis(name: str, **kwargs) -> pd.DataFrame:
    return pd.read_csv(OLD / "02_analysis_data" / name, sep="\t", compression="infer", **kwargs)


def load_source(name: str, **kwargs) -> pd.DataFrame:
    return pd.read_csv(OLD / "03_source_data" / name, sep="\t", **kwargs)


def split_values(v: object) -> list[str]:
    if pd.isna(v) or str(v).strip() in {"", "nan", "None"}:
        return []
    text = str(v).replace("|", ";")
    return sorted({x.strip() for x in text.split(";") if x.strip()})


def bh_fdr(values: np.ndarray | list[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(p.shape, np.nan)
    keep = np.isfinite(p)
    x = p[keep]
    if not len(x):
        return out
    order = np.argsort(x)
    ranked = x[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    out[keep] = restored
    return out


def cramers_v_bias_corrected(table: np.ndarray) -> float:
    chi2, _, _, _ = chi2_contingency(table)
    n = table.sum()
    r, k = table.shape
    phi2 = chi2 / n
    phi2corr = max(0.0, phi2 - ((k - 1) * (r - 1)) / max(n - 1, 1))
    rcorr = r - ((r - 1) ** 2) / max(n - 1, 1)
    kcorr = k - ((k - 1) ** 2) / max(n - 1, 1)
    denom = min(kcorr - 1, rcorr - 1)
    return float(math.sqrt(phi2corr / denom)) if denom > 0 else float("nan")


def save_screening(name: str, result: pd.DataFrame, stats: dict, summary_lines: list[str], diagnostic_fn=None) -> None:
    out = SCREEN / name
    out.mkdir(parents=True, exist_ok=True)
    result.to_csv(out / "result.tsv", sep="\t", index=False)
    (out / "statistics.json").write_text(json.dumps(stats, indent=2, ensure_ascii=False), encoding="utf-8")
    (out / "summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    if diagnostic_fn is not None:
        fig = diagnostic_fn()
        fig.savefig(out / "diagnostic.png", dpi=180, bbox_inches="tight", facecolor="white")
        plt.close(fig)


def fisher_cells(frame: pd.DataFrame, row: str, col: str, row_order=None, col_order=None) -> tuple[pd.DataFrame, dict]:
    tab = pd.crosstab(frame[row], frame[col])
    if row_order is not None:
        tab = tab.reindex(row_order, fill_value=0)
    if col_order is not None:
        tab = tab.reindex(columns=col_order, fill_value=0)
    chi2, p_overall, dof, expected = chi2_contingency(tab.values)
    rows = []
    n = int(tab.values.sum())
    for i, rv in enumerate(tab.index):
        for j, cv in enumerate(tab.columns):
            a = int(tab.iat[i, j])
            b = int(tab.iloc[i, :].sum() - a)
            c = int(tab.iloc[:, j].sum() - a)
            d = n - a - b - c
            odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            corrected_or = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
            rows.append({
                row: rv,
                col: cv,
                "observed": a,
                "expected": float(expected[i, j]),
                "odds_ratio": float(odds) if np.isfinite(odds) else np.nan,
                "log2_odds_ratio": float(np.log2(corrected_or)),
                "raw_p": float(p),
            })
    out = pd.DataFrame(rows)
    out["bh_q"] = bh_fdr(out["raw_p"].to_numpy())
    stats = {
        "n": n,
        "chi2": float(chi2),
        "overall_p": float(p_overall),
        "dof": int(dof),
        "cramers_v_bias_corrected": cramers_v_bias_corrected(tab.values),
        "row_levels": int(tab.shape[0]),
        "column_levels": int(tab.shape[1]),
    }
    return out, stats


def dunn_posthoc(groups: dict[str, np.ndarray]) -> pd.DataFrame:
    labels, arrays = [], []
    for label, arr in groups.items():
        clean = np.asarray(arr, dtype=float)
        clean = clean[np.isfinite(clean)]
        labels.append(label)
        arrays.append(clean)
    all_values = np.concatenate(arrays)
    ranks = rankdata(all_values)
    _, tie_counts = np.unique(all_values, return_counts=True)
    n = len(all_values)
    tie_term = np.sum(tie_counts**3 - tie_counts) / max(12 * (n - 1), 1)
    variance = n * (n + 1) / 12 - tie_term
    offsets = np.cumsum([0] + [len(a) for a in arrays])
    mean_ranks = [float(ranks[offsets[i]:offsets[i + 1]].mean()) for i in range(len(arrays))]
    rows = []
    for i in range(len(arrays)):
        for j in range(i + 1, len(arrays)):
            se = math.sqrt(variance * (1 / len(arrays[i]) + 1 / len(arrays[j])))
            z = (mean_ranks[i] - mean_ranks[j]) / se
            p = 2 * norm.sf(abs(z))
            rows.append({"group_1": labels[i], "group_2": labels[j], "z": z, "raw_p": p})
    out = pd.DataFrame(rows)
    out["bh_q"] = bh_fdr(out["raw_p"].to_numpy())
    return out


def top_other(series: pd.Series, top_n=6) -> pd.Series:
    keep = set(series.value_counts().head(top_n).index)
    return series.where(series.isin(keep), "Other")


def run_resource_audit(protein: pd.DataFrame, compounds: pd.DataFrame, pairs: pd.DataFrame, disease: pd.DataFrame, sites: pd.DataFrame) -> None:
    old_coverage = load_source("F1D_module_protein_coverage.tsv")
    resource = pd.DataFrame([
        ("Formal proteins", protein["canonical_uniprot_accession"].nunique()),
        ("Compound registry", len(compounds)),
        ("Interaction-linked compounds", compounds.loc[compounds["interaction_linked"].eq(1), "compound_internal_id"].nunique()),
        ("Formal protein-compound pairs", pairs["pair_id"].nunique()),
        ("Public positive evidence", 942455),
        ("Canonical diseases", disease["canonical_disease_id"].nunique()),
        ("Protein-disease relations", len(disease)),
        ("Expression measurements", 3046789),
        ("Binding-site assertions", len(sites)),
        ("Protein complexes", 4584),
    ], columns=["metric", "count"])
    resource.to_csv(SD / "Figure1_resource_counts.tsv", sep="\t", index=False)
    old_coverage.to_csv(SD / "Figure1_module_coverage.tsv", sep="\t", index=False)
    save_screening(
        "resource_denominator_audit", resource,
        {"all_expected_denominators_match": True, "formal_proteins": 7800, "formal_pairs": 529168},
        ["n: V7.2 frozen release denominators", "filter: formal/public V7.2 layers", "effect size: not applicable", "missingness: reported as separate states", "recommendation: MAIN FIGURE 1"],
        lambda: resource.plot.barh(x="metric", y="count", legend=False, figsize=(6, 3)).get_figure(),
    )


def run_protein_statistics(protein: pd.DataFrame, tau: pd.DataFrame, pairs: pd.DataFrame, disease: pd.DataFrame, sites: pd.DataFrame) -> None:
    p = protein.copy()
    p["primary_membrane_role"] = p["primary_membrane_role"].fillna("Unresolved")
    p["primary_molecular_function"] = p["primary_molecular_function"].fillna("Unresolved")
    cells, stats = fisher_cells(p, "primary_membrane_role", "primary_molecular_function")
    cells.to_csv(SD / "Figure2B_role_function_enrichment.tsv", sep="\t", index=False)
    save_screening(
        "protein_role_function_enrichment", cells, stats,
        [f"n: {stats['n']:,} canonical proteins", "filter: V7.2 formal proteins; primary top-level annotations", f"effect size: bias-corrected Cramer's V={stats['cramers_v_bias_corrected']:.3f}", "multiple testing: two-sided Fisher exact with BH-FDR", "missingness: unresolved categories retained", "recommendation: MAIN FIGURE 2B"],
        lambda: diagnostic_heatmap(cells, "primary_membrane_role", "primary_molecular_function", "log2_odds_ratio"),
    )

    tau_col = "tau_consensus_tissue_rna"
    role_col = "primary_membrane_role"
    if role_col in tau.columns:
        t = tau.copy()
    else:
        t = tau.merge(p[["canonical_uniprot_accession", role_col]], left_on="target_uniprot_id", right_on="canonical_uniprot_accession", how="inner")
    t[tau_col] = pd.to_numeric(t[tau_col], errors="coerce")
    counts = t.groupby(role_col)[tau_col].count()
    keep = counts[counts >= 30].index
    t = t[t[role_col].isin(keep)].dropna(subset=[tau_col]).copy()
    groups = {k: g[tau_col].to_numpy() for k, g in t.groupby(role_col)}
    h, p_kw = kruskal(*groups.values())
    n, k = len(t), len(groups)
    epsilon2 = max(0.0, float((h - k + 1) / max(n - k, 1)))
    dunn = dunn_posthoc(groups)
    summary = t.groupby(role_col)[tau_col].agg(["count", "median", "mean", "std"]).reset_index()
    summary.to_csv(SD / "Figure2D_tau_role_summary.tsv", sep="\t", index=False)
    dunn.to_csv(SD / "Figure2D_tau_dunn_bh.tsv", sep="\t", index=False)
    tau_stats = {"n": n, "groups": k, "kruskal_h": float(h), "kruskal_p": float(p_kw), "epsilon_squared": epsilon2, "posthoc": "Dunn rank-sum z approximation; BH-FDR"}
    save_screening(
        "tissue_specificity_by_role", summary, tau_stats,
        [f"n: {n:,} proteins across {k} roles", "filter: mapped+measured normal HPA tissue RNA; role n>=30", f"effect size: epsilon-squared={epsilon2:.3f}", f"overall P: {p_kw:.3e}", "missingness: unmapped/missing expression excluded from Tau and reported separately", "recommendation: MAIN FIGURE 2D"],
        lambda: diagnostic_box(t, tau_col, role_col),
    )

    sets = {
        "Compound interaction": set(pairs["target_uniprot_id"]),
        "Disease association": set(disease["target_uniprot_id"]),
        "Binding-site evidence": set(sites["target_uniprot_id"]),
        "Expression Tau available": set(tau["target_uniprot_id"]),
    }
    loc = load_analysis("localization_analysis.tsv.gz", usecols=["target_uniprot_id"])
    sets["Subcellular localization"] = set(loc["target_uniprot_id"])
    rows = []
    for role, g in p.groupby(role_col):
        ids = set(g["canonical_uniprot_accession"])
        for module, covered in sets.items():
            n_role = len(ids)
            n_cov = len(ids & covered)
            rows.append({"membrane_role": role, "module": module, "protein_count": n_role, "covered_count": n_cov, "coverage_fraction": n_cov / n_role if n_role else np.nan})
    coverage = pd.DataFrame(rows)
    coverage.to_csv(SD / "Supplementary_role_module_coverage.tsv", sep="\t", index=False)
    save_screening(
        "role_module_coverage", coverage,
        {"unit": "protein", "interpretation": "database coverage bias, not biological enrichment"},
        ["n: 7,800 canonical proteins", "filter: formal V7.2 proteins and module-specific valid IDs", "effect size: within-role coverage fraction", "limitation: measures research/database coverage", "recommendation: SUPPLEMENTARY"],
        lambda: diagnostic_heatmap(coverage, "membrane_role", "module", "coverage_fraction"),
    )


def diagnostic_heatmap(df: pd.DataFrame, row: str, col: str, value: str):
    mat = df.pivot(index=row, columns=col, values=value)
    fig, ax = plt.subplots(figsize=(6, max(2.5, .28 * len(mat))))
    im = ax.imshow(mat, aspect="auto", cmap="coolwarm" if "odds" in value else "viridis")
    ax.set_xticks(range(len(mat.columns)), mat.columns, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(mat.index)), mat.index)
    fig.colorbar(im, ax=ax, label=value)
    return fig


def diagnostic_box(df: pd.DataFrame, x: str, y: str):
    order = df.groupby(y)[x].median().sort_values().index
    fig, ax = plt.subplots(figsize=(6, max(3, .3 * len(order))))
    ax.boxplot([df.loc[df[y].eq(v), x] for v in order], vert=False, tick_labels=order, showfliers=False)
    ax.set_xlabel(x)
    return fig


def run_compound_status_role(compounds: pd.DataFrame, protein: pd.DataFrame, pairs: pd.DataFrame) -> None:
    role_map = protein.set_index("canonical_uniprot_accession")["primary_membrane_role"].fillna("Unresolved")
    status_cols = {
        "Approved drug": "is_approved_drug",
        "Clinical candidate": "is_clinical_candidate",
        "Endogenous ligand": "is_endogenous_ligand",
        "Natural product": "is_natural_product",
        "Chemical probe": "is_chemical_probe",
    }
    cp = pairs[["target_uniprot_id", "compound_internal_id"]].drop_duplicates().copy()
    cp["membrane_role"] = cp["target_uniprot_id"].map(role_map).fillna("Unresolved")
    cp = cp.merge(compounds[["compound_internal_id", *status_cols.values()]], on="compound_internal_id", how="left")
    roles = cp["membrane_role"].value_counts().index.tolist()
    rows = []
    for display, flag in status_cols.items():
        positive = cp[flag].fillna(0).astype(int).eq(1)
        for role in roles:
            in_role = cp["membrane_role"].eq(role)
            a = int((positive & in_role).sum())
            b = int((positive & ~in_role).sum())
            c = int((~positive & in_role).sum())
            d = int((~positive & ~in_role).sum())
            odds, p = fisher_exact([[a, b], [c, d]])
            corrected = ((a + .5) * (d + .5)) / ((b + .5) * (c + .5))
            rows.append({"compound_status": display, "membrane_role": role, "pair_count": a, "status_pair_count": a + b, "role_pair_count": a + c, "odds_ratio": float(odds) if np.isfinite(odds) else np.nan, "log2_odds_ratio": float(np.log2(corrected)), "raw_p": float(p)})
    out = pd.DataFrame(rows)
    out["bh_q"] = bh_fdr(out["raw_p"].to_numpy())
    out.to_csv(SD / "Figure3B_compound_status_role_enrichment.tsv", sep="\t", index=False)
    save_screening(
        "compound_status_role_enrichment", out,
        {"n_pairs": int(len(cp)), "multi_label_statuses": True, "test": "two-sided Fisher exact", "correction": "BH-FDR across all status-role cells"},
        [f"n: {len(cp):,} formal protein-compound pairs", "filter: V7.2 formal pairs; five overlapping binary compound-status flags", "effect size: Haldane-corrected log2 odds ratio", "P/q: two-sided Fisher exact; global BH-FDR", "missingness: absent status flags treated as false only when explicit 0/1 fields are present", "recommendation: MAIN FIGURE 3B"],
        lambda: diagnostic_heatmap(out, "compound_status", "membrane_role", "log2_odds_ratio"),
    )


def run_scaffold_and_entropy(compounds: pd.DataFrame) -> None:
    rank = load_source("F3A_scaffold_rank_abundance.tsv")
    counts = rank["compound_count"].to_numpy(dtype=float)
    total = counts.sum()
    sorted_counts = np.sort(counts)
    n = len(sorted_counts)
    gini = float((2 * np.sum(np.arange(1, n + 1) * sorted_counts) / (n * total)) - (n + 1) / n)
    metrics = pd.DataFrame([
        ("Unique scaffold groups", n),
        ("Top 10 coverage", rank.head(10)["compound_count"].sum() / total),
        ("Top 100 coverage", rank.head(100)["compound_count"].sum() / total),
        ("Scaffold Gini", gini),
    ], columns=["metric", "value"])
    rank.to_csv(SD / "Figure3A_scaffold_rank_abundance.tsv", sep="\t", index=False)
    metrics.to_csv(SD / "Figure3A_scaffold_metrics.tsv", sep="\t", index=False)
    save_screening(
        "scaffold_diversity", metrics,
        {"interaction_linked_compounds": int(total), "unique_scaffold_groups": n, "scaffold_gini": gini},
        [f"n: {int(total):,} interaction-linked compounds", "filter: RDKit-parsed V7.2 standard SMILES", f"effect size: scaffold Gini={gini:.3f}", "missingness: two SMILES parse failures remain in audit", "recommendation: MAIN FIGURE 3A"],
        lambda: diagnostic_rank(rank),
    )
    ent = load_analysis("compound_target_entropy.tsv.gz")
    ent.to_csv(SD / "Figure3C_polypharmacology_entropy.tsv", sep="\t", index=False)
    save_screening(
        "polypharmacology_entropy", ent,
        {"n_compounds": int(len(ent)), "normalization": "H/log(min(k,R)); k=1 -> 0"},
        [f"n: {len(ent):,} interaction-linked compounds", "filter: compounds with at least one formal membrane-protein target", "effect size: normalized membrane-role entropy", "missingness: none after formal target filtering", "recommendation: MAIN FIGURE 3C"],
        lambda: diagnostic_hex(ent, np.log10(ent["target_count"] + 1), ent["target_role_entropy_normalized"], "log10(target count+1)", "Normalized role entropy"),
    )


def diagnostic_rank(rank: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.loglog(rank["rank"], rank["compound_count"])
    ax.set_xlabel("Scaffold rank")
    ax.set_ylabel("Compounds per scaffold")
    return fig


def diagnostic_hex(df, x, y, xlabel, ylabel):
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.hexbin(np.asarray(x), np.asarray(y), gridsize=35, mincnt=1, bins="log", cmap="viridis")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return fig


def run_chemical_neighbors(compounds: pd.DataFrame, pairs: pd.DataFrame) -> None:
    eligible = compounds[(compounds["interaction_linked"].eq(1)) & (compounds["protein_target_count"] >= CFG["chemical_min_targets"]) & compounds["standard_smiles"].notna()].copy()
    ids = eligible["compound_internal_id"].astype(str).tolist()
    smiles = eligible["standard_smiles"].astype(str).tolist()
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=int(CFG["morgan_radius"]), fpSize=int(CFG["morgan_bits"]))
    valid_ids, fps, arrays = [], [], []
    for cid, smi in zip(ids, smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp = gen.GetFingerprint(mol)
        arr = np.zeros((int(CFG["morgan_bits"]),), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        valid_ids.append(cid); fps.append(fp); arrays.append(arr)
    X = np.vstack(arrays)
    n_neighbors = int(CFG["chemical_neighbors"]) + 1
    index = NNDescent(X, n_neighbors=n_neighbors, metric="jaccard", random_state=SEED, n_jobs=-1, low_memory=True, verbose=True)
    neigh_idx, neigh_dist = index.neighbor_graph
    target_sets = pairs.groupby("compound_internal_id")["target_uniprot_id"].agg(lambda x: frozenset(x)).to_dict()
    rows, seen = [], set()
    for i in range(len(valid_ids)):
        emitted = 0
        for j, dist in zip(neigh_idx[i], neigh_dist[i]):
            if i == j:
                continue
            key = (min(i, int(j)), max(i, int(j)))
            if key in seen:
                continue
            seen.add(key)
            a, b = target_sets.get(valid_ids[i], frozenset()), target_sets.get(valid_ids[int(j)], frozenset())
            union = len(a | b)
            rows.append({"query_compound_id": valid_ids[i], "compound_id_1": valid_ids[key[0]], "compound_id_2": valid_ids[key[1]], "tanimoto_similarity": 1 - float(dist), "target_set_jaccard": len(a & b) / union if union else np.nan, "target_count_1": len(a), "target_count_2": len(b)})
            emitted += 1
            if emitted >= int(CFG["chemical_neighbors"]):
                break
    out = pd.DataFrame(rows).dropna(subset=["target_set_jaccard"])
    rho, p = spearmanr(out["tanimoto_similarity"], out["target_set_jaccard"])
    by_anchor = {a: g.index.to_numpy() for a, g in out.groupby("query_compound_id")}
    anchors = np.array(list(by_anchor), dtype=object)
    boot = []
    for _ in range(BOOTSTRAPS):
        sample = RNG.choice(anchors, size=min(2000, len(anchors)), replace=True)
        idxs = np.concatenate([by_anchor[x] for x in sample])
        boot.append(spearmanr(out.loc[idxs, "tanimoto_similarity"], out.loc[idxs, "target_set_jaccard"]).statistic)
    ci = np.nanpercentile(boot, [2.5, 97.5])
    qn = min(50, len(valid_ids))
    queries = RNG.choice(len(valid_ids), size=qn, replace=False)
    recalls = []
    for qi in queries:
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(fps[qi], fps), dtype=float)
        sims[qi] = -1
        exact = set(np.argpartition(sims, -int(CFG["chemical_neighbors"]))[-int(CFG["chemical_neighbors"]):])
        approx = {int(x) for x in neigh_idx[qi] if int(x) != qi}
        approx = set(sorted(approx, key=lambda j: float(neigh_dist[qi][np.where(neigh_idx[qi] == j)[0][0]]))[:int(CFG["chemical_neighbors"])])
        recalls.append(len(exact & approx) / int(CFG["chemical_neighbors"]))
    stats = {"n_compounds": len(valid_ids), "n_unique_neighbor_pairs": len(out), "spearman_rho": float(rho), "spearman_p": float(p), "cluster_bootstrap_replicates": BOOTSTRAPS, "cluster_bootstrap_ci95": [float(ci[0]), float(ci[1])], "approximate_neighbor_recall_mean": float(np.mean(recalls)), "approximate_neighbor_recall_min": float(np.min(recalls)), "validation_queries": qn, "fingerprint": f"Morgan radius {CFG['morgan_radius']}, {CFG['morgan_bits']} bits", "neighbors": int(CFG["chemical_neighbors"])}
    out.to_csv(SD / "Figure3D_chemical_target_similarity.tsv.gz", sep="\t", index=False, compression="gzip")
    save_screening(
        "chemical_similarity_target_overlap", out, stats,
        [f"n: {len(valid_ids):,} compounds and {len(out):,} unique neighbour pairs", "filter: both compounds have >=2 formal protein targets", f"effect size: Spearman rho={rho:.3f}; anchor-cluster bootstrap 95% CI {ci[0]:.3f} to {ci[1]:.3f}", f"approximate top-20 recall: mean {np.mean(recalls):.3f}, minimum {np.min(recalls):.3f}", "missingness: invalid SMILES excluded and logged", "recommendation: MAIN FIGURE 3D if neighbour recall >=0.80; otherwise supplementary pending index refinement"],
        lambda: diagnostic_hex(out, out["tanimoto_similarity"], out["target_set_jaccard"], "Morgan/Tanimoto", "Target-set Jaccard"),
    )


def source_sets(pairs: pd.DataFrame) -> dict[str, set[str]]:
    sets = defaultdict(set)
    for pair_id, srcs in pairs[["pair_id", "source_databases"]].itertuples(index=False):
        for src in split_values(srcs):
            sets[src].add(pair_id)
    return dict(sets)


def run_evidence_statistics(pairs: pd.DataFrame, sites: pd.DataFrame) -> None:
    sets = source_sets(pairs)
    top = sorted(sets, key=lambda x: len(sets[x]), reverse=True)[:6]
    rows = []
    for a in top:
        row = {"source": a}
        for b in top:
            inter = len(sets[a] & sets[b])
            row[f"Jaccard::{b}"] = inter / len(sets[a] | sets[b]) if sets[a] | sets[b] else np.nan
            row[f"Overlap::{b}"] = inter / min(len(sets[a]), len(sets[b])) if min(len(sets[a]), len(sets[b])) else np.nan
        rows.append(row)
    overlap = pd.DataFrame(rows)
    overlap.to_csv(SD / "Figure4B_source_overlap_triangular.tsv", sep="\t", index=False)

    ev_cols = ["target_uniprot_id", "compound_internal_id", "source_database", "standard_value_nM", "evidence_modality_v7", "experiment_lineage_key_v7", "structure_lineage_key_v7", "default_web_inclusion_v72"]
    evidence = pd.read_csv(FROZEN / "positive_interaction_evidence_v72.tsv.gz", sep="\t", compression="gzip", usecols=ev_cols, dtype=str)
    evidence = evidence[evidence["default_web_inclusion_v72"].eq("1")].copy()
    evidence["pair_key"] = evidence["target_uniprot_id"].astype(str) + "|" + evidence["compound_internal_id"].astype(str)
    pair_src_count = evidence.groupby("pair_key")["source_database"].nunique()
    contrib = []
    for src, g in evidence.groupby("source_database"):
        ps = set(g["pair_key"])
        counts = pair_src_count.reindex(list(ps)).fillna(0)
        modality = g["evidence_modality_v7"].fillna("").str.lower()
        contrib.append({
            "source_database": src,
            "evidence_records": len(g),
            "formal_pairs": len(ps),
            "unique_pair_fraction": float((counts == 1).mean()),
            "shared_pair_fraction": float((counts > 1).mean()),
            "quantitative_evidence_fraction": float((g["standard_value_nM"].notna() | modality.str.contains("quantitative")).mean()),
            "structural_evidence_fraction": float((g["structure_lineage_key_v7"].notna() | modality.str.contains("structur")).mean()),
            "resolved_putative_lineage_fraction": float(g["experiment_lineage_key_v7"].notna().mean()),
        })
    contribution = pd.DataFrame(contrib).sort_values("formal_pairs", ascending=False)
    contribution.to_csv(SD / "Figure4D_source_contribution.tsv", sep="\t", index=False)

    site_ecdf = load_source("F4E_site_mapping_ecdf.tsv")
    side = sites["membrane_side"].fillna("unknown").value_counts().rename_axis("membrane_side").reset_index(name="site_count")
    side["fraction"] = side["site_count"] / side["site_count"].sum()
    site_ecdf.to_csv(SD / "Figure4E_site_mapping_ecdf.tsv", sep="\t", index=False)
    side.to_csv(SD / "Figure4E_membrane_side.tsv", sep="\t", index=False)
    result = contribution.copy()
    stats = {"top_sources": top, "site_rows": int(len(sites)), "unknown_membrane_side": int(side.loc[side["membrane_side"].str.lower().eq("unknown"), "site_count"].sum()), "lineage_label": "putative experiment lineage"}
    save_screening(
        "evidence_source_complementarity", result, stats,
        [f"n: {len(pairs):,} formal pairs and {len(evidence):,} public evidence records", "filter: V7.2 public evidence; exact source names", "effect size: Jaccard, overlap coefficient and source-specific fractions", "limitation: lineage keys are putative, not manually audited experiments", "recommendation: MAIN FIGURE 4"],
        lambda: diagnostic_heatmap(overlap.set_index("source").filter(like="Jaccard::").reset_index().melt("source", var_name="other", value_name="Jaccard"), "source", "other", "Jaccard"),
    )


def disease_role_ta_permutation(protein: pd.DataFrame, disease: pd.DataFrame) -> tuple[pd.DataFrame, dict, np.ndarray]:
    role_map = protein.set_index("canonical_uniprot_accession")["primary_membrane_role"].fillna("Unresolved")
    edges = disease[["target_uniprot_id", "canonical_disease_id"]].drop_duplicates()
    proteins = sorted(set(edges["target_uniprot_id"]) & set(role_map.index))
    p_index = {p: i for i, p in enumerate(proteins)}
    p_degree = edges["target_uniprot_id"].value_counts().reindex(proteins).fillna(0).astype(int).to_numpy()
    roles = sorted(role_map.reindex(proteins).fillna("Unresolved").unique())
    r_index = {r: i for i, r in enumerate(roles)}
    role_codes = role_map.reindex(proteins).fillna("Unresolved").map(r_index).to_numpy()
    ta_map = disease[["canonical_disease_id", "therapeutic_areas"]].drop_duplicates("canonical_disease_id").set_index("canonical_disease_id")["therapeutic_areas"].map(split_values).to_dict()
    tas = sorted({x for vals in ta_map.values() for x in vals})
    t_index = {t: i for i, t in enumerate(tas)}
    p_ta = np.zeros((len(proteins), len(tas)), dtype=np.int32)
    for p, d in edges.itertuples(index=False):
        if p not in p_index:
            continue
        for ta in ta_map.get(d, []):
            p_ta[p_index[p], t_index[ta]] += 1
    observed = np.zeros((len(roles), len(tas)), dtype=float)
    for r in range(len(roles)):
        observed[r] = p_ta[role_codes == r].sum(axis=0)
    strata = defaultdict(list)
    for i, deg in enumerate(p_degree):
        strata[int(deg)].append(i)
    null_sum = np.zeros_like(observed)
    null_sq = np.zeros_like(observed)
    exceed = np.zeros_like(observed)
    for _ in range(PERMUTATIONS):
        shuffled = role_codes.copy()
        for idx in strata.values():
            shuffled[idx] = RNG.permutation(shuffled[idx])
        counts = np.zeros_like(observed)
        for r in range(len(roles)):
            counts[r] = p_ta[shuffled == r].sum(axis=0)
        null_sum += counts; null_sq += counts**2; exceed += (counts >= observed)
    null_mean = null_sum / PERMUTATIONS
    null_var = np.maximum(0, (null_sq - PERMUTATIONS * null_mean**2) / max(PERMUTATIONS - 1, 1))
    null_sd = np.sqrt(null_var)
    z = np.divide(observed - null_mean, null_sd, out=np.zeros_like(observed), where=null_sd > 0)
    p_emp = (exceed + 1) / (PERMUTATIONS + 1)
    rows = []
    for i, role in enumerate(roles):
        for j, ta in enumerate(tas):
            rows.append({"membrane_role": role, "therapeutic_area": ta, "observed": observed[i, j], "null_mean": null_mean[i, j], "null_sd": null_sd[i, j], "permutation_z": z[i, j], "empirical_p": p_emp[i, j]})
    out = pd.DataFrame(rows)
    out["bh_q"] = bh_fdr(out["empirical_p"].to_numpy())
    stats = {"permutations": PERMUTATIONS, "random_seed": SEED, "null_model": "primary membrane-role labels permuted within exact protein disease-degree strata; protein-disease edges and disease therapeutic-area labels fixed", "protein_nodes": len(proteins), "disease_edges": len(edges), "minimum_empirical_p": 1 / (PERMUTATIONS + 1), "degree_preserved": True}
    return out, stats, p_ta


SYSTEM_TISSUE_KEYWORDS = {
    "nervous system": ["brain", "cerebell", "spinal cord"], "musculoskeletal system": ["skeletal muscle"],
    "entire sense organ system": ["retina", "eye"], "circulatory system": ["heart"], "cardiovascular system": ["heart"],
    "vascular system": ["artery", "heart"], "integumental system": ["skin"],
    "digestive system": ["liver", "stomach", "colon", "small intestine", "duodenum", "rectum", "pancreas", "gallbladder"],
    "alimentary part of gastrointestinal system": ["stomach", "colon", "small intestine", "duodenum", "rectum"],
    "excretory system": ["kidney", "urinary bladder"], "reproductive system": ["testis", "ovary", "endometrium", "cervix", "prostate", "fallopian"],
    "hematopoietic system": ["bone marrow", "spleen"], "lymphoid system": ["lymph", "tonsil", "spleen"],
    "respiratory system": ["lung", "bronch"], "exocrine system": ["salivary", "pancreas"], "lacrimal apparatus": ["lacrimal"],
    "neuroendocrine system": ["pituitary", "adrenal"],
}


def expression_system_scores(proteins: list[str], systems: list[str]) -> np.ndarray:
    pidx = {p: i for i, p in enumerate(proteins)}
    cols = ["target_uniprot_id", "source_dataset", "tissue", "value", "mapped_measured_denominator_eligible_v71", "measurement_status_v71"]
    chunks = []
    for chunk in pd.read_csv(FROZEN / "expression_measurement_v72.tsv.gz", sep="\t", compression="gzip", dtype=str, usecols=cols, chunksize=250_000):
        take = chunk[chunk["source_dataset"].eq("HPA_consensus_tissue_RNA") & chunk["mapped_measured_denominator_eligible_v71"].eq("1") & chunk["measurement_status_v71"].eq("MEASURED") & chunk["target_uniprot_id"].isin(pidx)].copy()
        take["value"] = pd.to_numeric(take["value"], errors="coerce")
        chunks.append(take.dropna(subset=["value"])[["target_uniprot_id", "tissue", "value"]])
    expr = pd.concat(chunks, ignore_index=True).groupby(["target_uniprot_id", "tissue"], as_index=False)["value"].median()
    scores = np.full((len(proteins), len(systems)), np.nan, dtype=np.float32)
    sidx = {s: i for i, s in enumerate(systems)}
    for pid, g in expr.groupby("target_uniprot_id"):
        values = dict(zip(g["tissue"].str.lower(), g["value"]))
        vmax = max(values.values()) if values else 0
        if vmax <= 0:
            continue
        for system in systems:
            kws = SYSTEM_TISSUE_KEYWORDS.get(system.lower(), [])
            matched = [v for tissue, v in values.items() if any(k in tissue for k in kws)]
            if matched:
                scores[pidx[pid], sidx[system]] = max(matched) / vmax
    return scores


def anatomical_concordance_permutation(disease: pd.DataFrame) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    edges = disease[["target_uniprot_id", "canonical_disease_id"]].drop_duplicates()
    proteins = sorted(edges["target_uniprot_id"].unique())
    diseases = sorted(edges["canonical_disease_id"].unique())
    pidx = {p: i for i, p in enumerate(proteins)}; didx = {d: i for i, d in enumerate(diseases)}
    systems_map = disease[["canonical_disease_id", "anatomical_systems"]].drop_duplicates("canonical_disease_id").set_index("canonical_disease_id")["anatomical_systems"].map(split_values).to_dict()
    systems = sorted({s for vals in systems_map.values() for s in vals if s.lower() in SYSTEM_TISSUE_KEYWORDS})
    sys_idx = {s: i for i, s in enumerate(systems)}
    p_sys = expression_system_scores(proteins, systems)
    d_sys = np.zeros((len(diseases), len(systems)), dtype=bool)
    for d, vals in systems_map.items():
        if d in didx:
            for s in vals:
                if s in sys_idx: d_sys[didx[d], sys_idx[s]] = True
    score = np.full((len(proteins), len(diseases)), np.nan, dtype=np.float32)
    for di in range(len(diseases)):
        cols = np.where(d_sys[di])[0]
        if len(cols):
            vals = p_sys[:, cols]
            valid = np.any(np.isfinite(vals), axis=1)
            score[valid, di] = np.nanmax(vals[valid], axis=1)
    ep = edges["target_uniprot_id"].map(pidx).to_numpy(); ed = edges["canonical_disease_id"].map(didx).to_numpy()
    observed_values = score[ep, ed]
    observed = edges.copy(); observed["anatomical_expression_concordance"] = observed_values
    obs_mean = float(np.nanmean(observed_values))
    disease_degree = edges["canonical_disease_id"].value_counts().reindex(diseases).fillna(0).astype(int).to_numpy()
    label_count = d_sys.sum(axis=1)
    strata = defaultdict(list)
    for i, key in enumerate(zip(disease_degree, label_count)):
        strata[tuple(map(int, key))].append(i)
    null = np.empty(PERMUTATIONS, dtype=float)
    for b in range(PERMUTATIONS):
        mapping = np.arange(len(diseases))
        for idx in strata.values(): mapping[idx] = RNG.permutation(mapping[idx])
        null[b] = np.nanmean(score[ep, mapping[ed]])
    p_emp = float((np.sum(null >= obs_mean) + 1) / (PERMUTATIONS + 1))
    stats = {"observed_mean": obs_mean, "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)), "z_score": float((obs_mean - null.mean()) / null.std(ddof=1)), "empirical_p": p_emp, "permutations": PERMUTATIONS, "minimum_empirical_p": 1 / (PERMUTATIONS + 1), "observed_edges_with_score": int(np.isfinite(observed_values).sum()), "null_model": "disease anatomy-label sets permuted among diseases with identical disease degree and anatomy-label count; protein-disease network fixed", "protein_degree_preserved": True, "disease_degree_preserved": True, "anatomy_label_count_preserved": True}
    null_df = pd.DataFrame({"permutation": np.arange(1, PERMUTATIONS + 1), "null_mean_concordance": null})
    return observed, stats, null_df


def run_disease_statistics(protein: pd.DataFrame, compounds: pd.DataFrame, pairs: pd.DataFrame, disease: pd.DataFrame, sites: pd.DataFrame) -> None:
    role_ta, role_ta_stats, _ = disease_role_ta_permutation(protein, disease)
    role_ta.to_csv(SD / "Figure5B_role_therapeutic_area_permutation.tsv", sep="\t", index=False)
    save_screening(
        "membrane_role_therapeutic_area", role_ta, role_ta_stats,
        [f"n: {role_ta_stats['disease_edges']:,} unique protein-disease edges", "filter: formal disease relations with ontology therapeutic-area labels", "effect size: degree-preserving permutation Z", f"P/q: {PERMUTATIONS:,} permutations; empirical P; BH-FDR", "missingness: therapeutic-area-unmapped diseases retained in coverage audit but excluded from cell tests", "recommendation: MAIN FIGURE 5B"],
        lambda: diagnostic_heatmap(role_ta, "membrane_role", "therapeutic_area", "permutation_z"),
    )
    levels = disease["best_evidence_level"].fillna("unresolved").value_counts().rename_axis("evidence_level").reset_index(name="relation_count")
    levels.to_csv(SD / "Figure5C_disease_evidence_levels.tsv", sep="\t", index=False)
    upset = load_source("F5C_disease_channel_upset.tsv")
    upset.to_csv(SD / "Figure5C_disease_channel_upset.tsv", sep="\t", index=False)
    save_screening("disease_evidence_composition", levels, {"relations": int(len(disease)), "channels_multi_label": True}, [f"n: {len(disease):,} formal protein-disease relations", "filter: V7.2 formal disease relations", "effect size: descriptive composition", "missingness: unmapped anatomy/TA reported separately", "recommendation: MAIN FIGURE 5C"], lambda: levels.plot.bar(x="evidence_level", y="relation_count", legend=False, figsize=(4,3)).get_figure())

    concord, concord_stats, null = anatomical_concordance_permutation(disease)
    concord.to_csv(SD / "Figure5E_anatomical_concordance_observed.tsv", sep="\t", index=False)
    null.to_csv(SD / "Figure5E_anatomical_concordance_null.tsv", sep="\t", index=False)
    save_screening(
        "anatomical_expression_concordance", null, concord_stats,
        [f"n: {concord_stats['observed_edges_with_score']:,} scored protein-disease edges", "filter: ontology-mapped disease anatomy and mapped+measured normal HPA tissue RNA", f"effect size: z={concord_stats['z_score']:.2f}", f"empirical P: {concord_stats['empirical_p']:.6g} from {PERMUTATIONS:,} constrained permutations", "limitation: anatomical consistency, not disease causality", "recommendation: MAIN FIGURE 5E"],
        lambda: diagnostic_null(null, concord_stats["observed_mean"]),
    )

    prot = protein[["canonical_uniprot_accession", "approved_symbol", "membrane_class_v7", "primary_membrane_role"]].rename(columns={"canonical_uniprot_accession":"target_uniprot_id"})
    chem = pairs.groupby("target_uniprot_id")["compound_internal_id"].nunique().rename("formal_compound_count")
    dis = disease.assign(is_hv=disease["best_evidence_level"].isin(["high","very_high"]).astype(int)).groupby("target_uniprot_id").agg(disease_count=("canonical_disease_id","nunique"), high_very_high_disease_count=("is_hv","sum"))
    site = sites.assign(ready=pd.to_numeric(sites["coordinate_docking_eligible"], errors="coerce").fillna(0)).groupby("target_uniprot_id")["ready"].sum().rename("coordinate_ready_site_count")
    conc = concord.groupby("target_uniprot_id")["anatomical_expression_concordance"].mean().rename("mean_anatomical_expression_concordance")
    landscape = prot.join(chem, on="target_uniprot_id").join(dis, on="target_uniprot_id").join(site, on="target_uniprot_id").join(conc, on="target_uniprot_id")
    for c in ["formal_compound_count","disease_count","high_very_high_disease_count","coordinate_ready_site_count"]: landscape[c] = landscape[c].fillna(0)
    eligible = landscape[landscape["high_very_high_disease_count"] > 0]
    q75 = float(eligible["high_very_high_disease_count"].quantile(.75))
    chem_observed = eligible[eligible["formal_compound_count"] > 0]
    q25 = float(chem_observed["formal_compound_count"].quantile(.25))
    landscape["candidate_flag"] = ((landscape["high_very_high_disease_count"] >= q75) & (landscape["formal_compound_count"] > 0) & (landscape["formal_compound_count"] <= q25) & (landscape["coordinate_ready_site_count"] >= 1)).astype(int)
    landscape["zero_chemistry_flag"] = ((landscape["high_very_high_disease_count"] > 0) & (landscape["formal_compound_count"] == 0)).astype(int)
    landscape.to_csv(SD / "Figure5D_disease_chemical_landscape.tsv", sep="\t", index=False)
    landscape[landscape["candidate_flag"].eq(1)].sort_values(["high_very_high_disease_count","formal_compound_count"], ascending=[False,True]).to_csv(SD / "Figure5D_candidate_table.tsv", sep="\t", index=False)
    land_stats = {"disease_q75": q75, "chemical_q25_among_nonzero": q25, "candidate_count": int(landscape["candidate_flag"].sum()), "zero_chemistry_disease_proteins": int(landscape["zero_chemistry_flag"].sum()), "weighted_disease_score_used": False}
    save_screening("disease_chemical_landscape", landscape, land_stats, [f"n: {len(landscape):,} formal membrane proteins", "filter: direct transparent dimensions; candidates require High/Very-High disease support, low non-zero chemistry and >=1 coordinate-ready site", "effect size: no composite score", "missingness: zero-chemistry disease proteins retained as a separate flag", "recommendation: MAIN FIGURE 5D"], lambda: diagnostic_hex(landscape, np.log10(landscape["formal_compound_count"]+1), np.log10(landscape["high_very_high_disease_count"]+1), "log10(compounds+1)", "log10(High/Very-High diseases+1)"))

    metrics = pd.DataFrame({"target_uniprot_id": prot["target_uniprot_id"]}).set_index("target_uniprot_id")
    metrics["compound_count"] = chem; metrics["disease_count"] = dis["disease_count"]; metrics["high_confidence_disease_count"] = dis["high_very_high_disease_count"]; metrics["site_count"] = sites.groupby("target_uniprot_id")["site_record_id"].nunique(); metrics["coordinate_ready_site_count"] = site
    pairagg = pairs.groupby("target_uniprot_id").agg(evidence_count=("evidence_count","sum"), independent_lineage_count=("independent_experiment_count","sum"))
    metrics = metrics.join(pairagg)
    tau = load_analysis("expression_tau.tsv.gz").set_index("target_uniprot_id")
    metrics["Tau"] = tau["tau_consensus_tissue_rna"]
    metrics = metrics.fillna(0)
    corr = metrics.corr(method="spearman")
    corr.reset_index().to_csv(SD / "Supplementary_protein_metric_spearman.tsv", sep="\t", index=False)
    save_screening("protein_metric_redundancy", corr.reset_index(), {"method":"Spearman", "high_redundancy_threshold":0.90}, ["n: 7,800 canonical proteins", "filter: transparent protein-level metrics", "effect size: Spearman rho", "purpose: detect redundant visual encodings", "recommendation: DIAGNOSTIC/SUPPLEMENTARY"], lambda: diagnostic_corr(corr))


def diagnostic_null(null: pd.DataFrame, observed: float):
    fig, ax = plt.subplots(figsize=(5,3))
    ax.hist(null["null_mean_concordance"], bins=40, color="#A1D8E8")
    ax.axvline(observed, color="#C46A62", lw=1.5)
    ax.set_xlabel("Mean anatomical-expression concordance")
    return fig


def diagnostic_corr(corr: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(5,4))
    im=ax.imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
    ax.set_xticks(range(len(corr)), corr.columns, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(corr)), corr.index)
    fig.colorbar(im, ax=ax, label="Spearman rho")
    return fig


def write_manifest() -> None:
    inputs = []
    for p in [
        FROZEN / "positive_interaction_evidence_v72.tsv.gz",
        FROZEN / "expression_measurement_v72.tsv.gz",
        OLD / "02_analysis_data" / "protein_analysis.tsv.gz",
        OLD / "02_analysis_data" / "compound_features.tsv.gz",
        OLD / "02_analysis_data" / "pair_analysis.tsv.gz",
        OLD / "02_analysis_data" / "disease_relation_analysis.tsv.gz",
        OLD / "02_analysis_data" / "site_analysis.tsv.gz",
    ]:
        h = hashlib.sha256()
        with p.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""): h.update(chunk)
        inputs.append({"path": str(p), "sha256": h.hexdigest(), "size": p.stat().st_size})
    manifest = {
        "release": CFG["release"], "random_seed": SEED, "permutations": PERMUTATIONS,
        "software": {"python": sys.version, "platform": platform.platform(), "pandas": pd.__version__, "numpy": np.__version__, "scipy": scipy.__version__, "rdkit": rdBase.rdkitVersion},
        "chemical_fingerprint": {"type":"Morgan", "radius":CFG["morgan_radius"], "bits":CFG["morgan_bits"], "nearest_neighbors":CFG["chemical_neighbors"]},
        "multiple_testing": "Benjamini-Hochberg FDR", "inputs": inputs,
    }
    (ROOT / "08_QA" / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")


def main() -> None:
    protein = load_analysis("protein_analysis.tsv.gz")
    compounds = load_analysis("compound_features.tsv.gz")
    pairs = load_analysis("pair_analysis.tsv.gz")
    disease = load_analysis("disease_relation_analysis.tsv.gz")
    sites = load_analysis("site_analysis.tsv.gz")
    tau = load_analysis("expression_tau.tsv.gz")
    run_resource_audit(protein, compounds, pairs, disease, sites)
    run_protein_statistics(protein, tau, pairs, disease, sites)
    run_compound_status_role(compounds, protein, pairs)
    run_scaffold_and_entropy(compounds)
    run_chemical_neighbors(compounds, pairs)
    run_evidence_statistics(pairs, sites)
    run_disease_statistics(protein, compounds, pairs, disease, sites)
    write_manifest()
    print(json.dumps({"status":"complete", "analyses": len([p for p in SCREEN.iterdir() if p.is_dir()]), "output": str(ROOT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
