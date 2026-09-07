"""
Add scoring/ranking fields to docking manifest.
Computes a composite priority score for ordering docking tasks.
"""
import csv, gzip, os
from collections import defaultdict

MANIFEST = r'D:\finale\negative_upgrade\docking_package\docking_manifest.tsv'
EVIDENCE = r'D:\finale\01_正式数据_V6.2\binding_evidence_master_v6_2.tsv.gz'
PROTEIN = r'D:\finale\01_正式数据_V6.2\human_membrane_protein_master_v6_2.tsv'
OUT_MANIFEST = r'D:\finale\negative_upgrade\docking_package\docking_manifest_ranked.tsv'

print('Loading protein E tiers...')
prot_e = {}
with open(PROTEIN, encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        prot_e[row['target_uniprot_id']] = row.get('evidence_level_v52', '')

print('Loading manifest...')
tasks = []
with open(MANIFEST, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    manifest_cols = reader.fieldnames
    for row in reader:
        tasks.append(row)
print(f'  {len(tasks):,} tasks')

# ============================================================
# Collect best evidence-level data per pair
# ============================================================
print('Collecting evidence-level scoring fields...')
pair_evidence = {}  # (cpd, up) -> best evidence row data

with gzip.open(EVIDENCE, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        cpd = row.get('compound_internal_id', '')
        up = row.get('target_uniprot_id', '')
        key = (cpd, up)

        be = row.get('evidence_tier', '')
        snm_str = row.get('standard_value_nM', '')

        # Parse activity
        try:
            snm = float(snm_str) if snm_str else None
        except ValueError:
            snm = None

        # Determine best evidence per pair: prefer BE1, then lower nM
        if key not in pair_evidence:
            pair_evidence[key] = row
        else:
            existing = pair_evidence[key]
            existing_be = existing.get('evidence_tier', '')
            tier_order = {'BE1': 3, 'BE2': 2, 'BE3': 1}

            # Better BE tier wins
            if tier_order.get(be, 0) > tier_order.get(existing_be, 0):
                pair_evidence[key] = row
            elif tier_order.get(be, 0) == tier_order.get(existing_be, 0) and snm:
                # Same tier: prefer lower (more potent) nM
                existing_snm = None
                try:
                    existing_snm = float(existing.get('standard_value_nM', '') or '')
                except ValueError:
                    pass
                if existing_snm is None or (snm and snm < existing_snm):
                    pair_evidence[key] = row

print(f'  Best evidence collected for {len(pair_evidence):,} pairs')

# ============================================================
# Compute priority score
# ============================================================
def compute_priority_score(task, ev_row):
    """Composite score 0-100, higher = higher docking priority."""
    score = 0.0

    # 1. BE Tier (0-35 points)
    be = ev_row.get('evidence_tier', '')
    be_scores = {'BE1': 35, 'BE2': 20, 'BE3': 5}
    score += be_scores.get(be, 0)

    # 2. Activity potency (0-20 points)
    # pAct = -log10(M), cap at 10 (1nM) ~ 3 (1mM)
    snm_str = ev_row.get('standard_value_nM', '')
    try:
        snm = float(snm_str) if snm_str else None
    except ValueError:
        snm = None
    if snm and snm > 0:
        import math
        p_act = 9 - math.log10(snm)  # 1nM=9, 1000nM=6, 1000000nM=3
        score += max(0, min(20, p_act * 2))  # scale: pAct 10->20pts, pAct 5->10pts
    score += 0  # no nM = no activity points

    # 3. Bio-status flags (0-15 points)
    bio = task.get('bio_flags', '')
    flag_scores = {
        'is_approved_drug': 15,
        'is_clinical_candidate': 12,
        'is_endogenous_ligand': 8,
        'is_natural_product': 5,
        'is_chemical_probe': 10,
    }
    flag_pts = []
    for flag, pts in flag_scores.items():
        if flag in bio:
            flag_pts.append(pts)
    if flag_pts:
        score += max(flag_pts)  # best flag wins

    # 4. PDB source quality (0-15 points)
    pdb_src = task.get('pdb_source', '')
    src_scores = {'direct': 15, 'fallback': 8, 'none': 0}
    score += src_scores.get(pdb_src, 0)

    # 5. Has binding residues (0-10 points)
    has_res = task.get('has_binding_residues', '')
    if has_res == 'True':
        score += 10

    # 6. Evidence type quality (0-5 points)
    ev_type = ev_row.get('evidence_type', '')
    if 'structure' in ev_type.lower():
        score += 5
    elif 'quantitative' in ev_type.lower():
        score += 3
    elif 'curated' in ev_type.lower():
        score += 4

    return round(score, 1)


# ============================================================
# Add scoring fields to all tasks
# ============================================================
print('Computing priority scores...')

new_fields = [
    'priority_score',
    'best_evidence_tier',
    'best_activity_nM',
    'activity_relation',
    'activity_outcome',
    'evidence_directness',
    'pubmed_ids',
    'doi',
    'source_database',
    'protein_e_tier',
]

all_cols = list(manifest_cols) + new_fields

scored_tasks = []
for task in tasks:
    cpd = task['compound_internal_id']
    up = task['target_uniprot']
    key = (cpd, up)

    ev = pair_evidence.get(key, {})

    if not ev:
        # No evidence found for this pair - should not happen
        score = 0
        task.update({
            'priority_score': '0',
            'best_evidence_tier': '',
            'best_activity_nM': '',
            'activity_relation': '',
            'activity_outcome': '',
            'evidence_directness': '',
            'pubmed_ids': '',
            'doi': '',
            'source_database': '',
            'protein_e_tier': prot_e.get(up, ''),
        })
    else:
        score = compute_priority_score(task, ev)
        task.update({
            'priority_score': str(score),
            'best_evidence_tier': ev.get('evidence_tier', ''),
            'best_activity_nM': ev.get('standard_value_nM', ''),
            'activity_relation': ev.get('activity_relation', ''),
            'activity_outcome': ev.get('activity_outcome_v60', ''),
            'evidence_directness': ev.get('evidence_directness', ''),
            'pubmed_ids': ev.get('pubmed_ids', ''),
            'doi': ev.get('doi', ''),
            'source_database': ev.get('source_database', ''),
            'protein_e_tier': prot_e.get(up, ''),
        })

    scored_tasks.append(task)

# Sort by priority score (highest first)
scored_tasks.sort(key=lambda t: float(t.get('priority_score', 0)), reverse=True)

# Re-assign task_ids in sorted order
for i, task in enumerate(scored_tasks):
    task['task_id'] = f'DOCK-{i+1:06d}'

# ============================================================
# Write output
# ============================================================
print(f'Writing ranked manifest: {OUT_MANIFEST}')
with open(OUT_MANIFEST, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=all_cols, delimiter='\t', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(scored_tasks)

# ============================================================
# Score distribution
# ============================================================
print(f'\n{"="*60}')
print(f'PRIORITY SCORE DISTRIBUTION')
print(f'{"="*60}')

from collections import Counter
score_bins = Counter()
for t in scored_tasks:
    s = float(t.get('priority_score', 0))
    if s >= 80:
        score_bins['80-100 (excellent)'] += 1
    elif s >= 60:
        score_bins['60-80 (very good)'] += 1
    elif s >= 40:
        score_bins['40-60 (good)'] += 1
    elif s >= 20:
        score_bins['20-40 (fair)'] += 1
    else:
        score_bins['0-20 (low)'] += 1

for bin_label in ['80-100 (excellent)', '60-80 (very good)', '40-60 (good)', '20-40 (fair)', '0-20 (low)']:
    cnt = score_bins.get(bin_label, 0)
    bar = '#' * (cnt // 5000)
    print(f'  {bin_label:25s}: {cnt:>8,}  {bar}')

# Top-10 examples
print(f'\nTop 10 priority tasks:')
for t in scored_tasks[:10]:
    print(f"  {t['task_id']}: score={t['priority_score']}  "
          f"{t['target_symbol']:10s} x {t['preferred_name'][:40]:40s}  "
          f"BE={t['best_evidence_tier']}  nM={t['best_activity_nM']}  "
          f"PDB={t['pdb_source']}  bio={t['bio_flags'][:30]}")

# Bottom 5
print(f'\nBottom 5 priority tasks:')
for t in scored_tasks[-5:]:
    print(f"  {t['task_id']}: score={t['priority_score']}  "
          f"{t['target_symbol']:10s} x {t['preferred_name'][:40]:40s}  "
          f"BE={t['best_evidence_tier']}  nM={t.get('best_activity_nM','')}  "
          f"PDB={t['pdb_source']}")

# Score components for top vs bottom
print(f'\n{"="*60}')
print(f'SCORE COMPONENTS: Top-100 vs Bottom-100')
print(f'{"="*60}')
for label, subset in [('Top 100', scored_tasks[:100]), ('Bottom 100', scored_tasks[-100:])]:
    avg_be = sum(1 for t in subset if t['best_evidence_tier'] == 'BE1') / len(subset) * 100
    avg_pdb = sum(1 for t in subset if t['pdb_source'] == 'direct') / len(subset) * 100
    avg_res = sum(1 for t in subset if t['has_binding_residues'] == 'True') / len(subset) * 100
    has_activity = sum(1 for t in subset if t.get('best_activity_nM', '')) / len(subset) * 100
    print(f'  {label}: BE1={avg_be:.0f}%  direct_PDB={avg_pdb:.0f}%  has_res={avg_res:.0f}%  has_nM={has_activity:.0f}%')

print(f'\nDone. Output: {OUT_MANIFEST}')
