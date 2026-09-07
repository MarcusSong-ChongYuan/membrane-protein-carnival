"""Fast tier breakdown: single pass, minimal memory."""
import csv, gzip
from collections import Counter, defaultdict

EVIDENCE = r'D:\finale\01_正式数据_V6.2\binding_evidence_master_v6_2.tsv.gz'
PROTEIN = r'D:\finale\01_正式数据_V6.2\human_membrane_protein_master_v6_2.tsv'

# Load membrane protein set
all_mp = set()
with open(PROTEIN, encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        all_mp.add(row['target_uniprot_id'])
print(f'Membrane proteins: {len(all_mp):,}')

# Single pass: count pairs, best tier per pair
print('Scanning evidence (single pass)...')
pair_best_be = {}  # (cpd,up) -> best BE tier
pair_has_pdb = set()
pair_has_res = set()
prot_best_pdb = {}  # up -> pdb (prefer with residues)
pair_best_e = {}  # (cpd,up) -> best E tier
tier_order = {'BE1': 3, 'BE2': 2, 'BE3': 1}
row_count = 0

with gzip.open(EVIDENCE, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        row_count += 1
        if row_count % 1000000 == 0:
            print(f'  {row_count:,} rows...')

        cpd = row.get('compound_internal_id', '')
        up = row.get('target_uniprot_id', '')
        if up not in all_mp:
            continue

        key = (cpd, up)
        be = row.get('evidence_tier', '')
        e_tier = row.get('protein_release_tier', '')

        # Track best BE tier for this pair
        be_val = tier_order.get(be, 0)
        if key not in pair_best_be or be_val > tier_order.get(pair_best_be.get(key, ''), 0):
            pair_best_be[key] = be
            if e_tier:
                pair_best_e[key] = e_tier
        elif not pair_best_e.get(key) and e_tier:
            pair_best_e[key] = e_tier

        # Track PDB access
        pdb_str = row.get('pdb_ids', '').strip()
        if pdb_str:
            pair_has_pdb.add(key)
            # Protein-level PDB fallback
            for p in pdb_str.split(';'):
                p = p.strip().lower()
                if len(p) == 4 and p[0].isdigit():
                    has_res = bool(row.get('binding_site_residues', '').strip())
                    if up not in prot_best_pdb or has_res:
                        prot_best_pdb[up] = p

        if row.get('binding_site_residues', '').strip():
            pair_has_res.add(key)

print(f'  Done: {row_count:,} rows scanned')
total_pairs = len(pair_best_be)
print(f'  Unique pairs: {total_pairs:,}')
print(f'  Unique compounds: {len(set(k[0] for k in pair_best_be)):,}')
print(f'  Unique proteins: {len(set(k[1] for k in pair_best_be)):,}')
print(f'  Proteins with any PDB: {len(prot_best_pdb):,}')

# ============================================================
# Statistics
# ============================================================
print(f'\n{"="*60}')
print(f'By Evidence Tier (BE)')
print(f'{"="*60}')
be_counts = Counter(pair_best_be.values())
for tier in ['BE1', 'BE2', 'BE3']:
    cnt = be_counts.get(tier, 0)
    print(f'  {tier}: {cnt:,} ({cnt/total_pairs*100:.1f}%)')
empty_be = sum(1 for v in pair_best_be.values() if not v)
print(f'  unclassified: {empty_be:,} ({empty_be/total_pairs*100:.1f}%)')

print(f'\n{"="*60}')
print(f'By Protein Release Tier (E)')
print(f'{"="*60}')
e_counts = Counter(pair_best_e.values())
for tier in ['E0', 'E1', 'E2', 'E3']:
    cnt = e_counts.get(tier, 0)
    print(f'  {tier}: {cnt:,} ({cnt/total_pairs*100:.1f}%)')
empty_e = sum(1 for v in pair_best_e.values() if not v)
print(f'  unclassified: {empty_e:,} ({empty_e/total_pairs*100:.1f}%)')

# Cross-tabulation
print(f'\n{"="*60}')
print(f'Cross: BE x E')
print(f'{"="*60}')
cross = Counter()
for key in pair_best_be:
    be = pair_best_be[key] or 'BE-'
    e = pair_best_e.get(key) or 'E-'
    cross[(be, e)] += 1

print(f'{"BE\\E":>8} {"E0":>8} {"E1":>8} {"E2":>8} {"E3":>8} {"E-":>8} {"TOTAL":>8}')
for be in ['BE1', 'BE2', 'BE3', 'BE-']:
    vals = [cross.get((be, e), 0) for e in ['E0', 'E1', 'E2', 'E3', 'E-']]
    row_total = sum(vals)
    row = ' '.join(f'{v:>8,}' for v in vals)
    print(f'{be:>8} {row} {row_total:>8,}')

# ============================================================
# PDB Coverage by Scenario
# ============================================================
print(f'\n{"="*60}')
print(f'SCENARIO ANALYSIS')
print(f'{"="*60}')

scenarios = [
    ('ALL pairs (maximum)',
     lambda k: True),
    ('BE1 only (structural evidence)',
     lambda k: pair_best_be.get(k) == 'BE1'),
    ('BE1+BE2 (drop PubChem HTS)',
     lambda k: pair_best_be.get(k) in ('BE1', 'BE2')),
    ('BE1+BE2 + E1+E2+E3 (with protein validation)',
     lambda k: pair_best_be.get(k) in ('BE1', 'BE2') and pair_best_e.get(k) in ('E1', 'E2', 'E3')),
    ('BE1+BE2 + E1+E2 (strictest)',
     lambda k: pair_best_be.get(k) in ('BE1', 'BE2') and pair_best_e.get(k) in ('E1', 'E2')),
]

for label, fn in scenarios:
    keys = [k for k in pair_best_be if fn(k)]
    n = len(keys)
    if n == 0:
        print(f'\n  {label}: 0 pairs')
        continue
    direct_pdb = sum(1 for k in keys if k in pair_has_pdb)
    fallback_pdb = sum(1 for k in keys if k in pair_has_pdb or k[1] in prot_best_pdb)
    with_res = sum(1 for k in keys if k in pair_has_res)
    prots = len(set(k[1] for k in keys))
    cpds = len(set(k[0] for k in keys))

    print(f'\n  [{label}]')
    print(f'    Pairs: {n:,}')
    print(f'    Proteins: {prots:,} | Compounds: {cpds:,}')
    print(f'    Direct PDB: {direct_pdb:,} ({direct_pdb/n*100:.1f}%)')
    print(f'    PDB available (incl fallback): {fallback_pdb:,} ({fallback_pdb/n*100:.1f}%)')
    print(f'    Has binding residues: {with_res:,} ({with_res/n*100:.1f}%)')

# ============================================================
# Recommendation
# ============================================================
print(f'\n{"="*60}')
print(f'RECOMMENDATION')
print(f'{"="*60}')
print(f'''
Evidence Tier (BE) 判断的是证据类型:
  BE1 = PDB 共晶结构 → 蛋白结构 + 结合位点明确 → 最适合 docking
  BE2 = ChEMBL/DrugCentral 文献 curated 数据 → 有 Ki/Kd/IC50 值 → 可以 docking（需 PDB）
  BE3 = PubChem HTS 高通量筛选 → 假阳性率较高 → 不建议作为 docking 依据

Protein Release Tier (E) 判断的是蛋白层面的证据丰富度:
  E1 = 多来源交叉验证 → 最高置信
  E2 = 较好验证
  E3 = 证据有限
  E0 = 最低
  unclassified = 未进入分级体系

Docking 推荐策略:
  BE1+BE2: 去掉 PubChem HTS 假阳性，保留有结构和文献支撑的
  E 等级可以放宽 (E1+E2+E3+unclassified): 因为只要有 PDB 就能做 docking
''')

# Final recommendation
for be_filt, e_filt, label in [
    (('BE1','BE2'), None, 'BE1+BE2, all E (推荐)'),
    (('BE1','BE2'), ('E1','E2','E3'), 'BE1+BE2 + E1+E2+E3'),
    (('BE1',), None, 'BE1 only (最保守)'),
    (None, None, 'ALL (最大范围)'),
]:
    if be_filt:
        keys = [k for k in pair_best_be if pair_best_be.get(k) in be_filt]
    else:
        keys = list(pair_best_be.keys())
    if e_filt:
        keys = [k for k in keys if pair_best_e.get(k) in e_filt]
    n = len(keys)
    pdb_ok = sum(1 for k in keys if k in pair_has_pdb or k[1] in prot_best_pdb)
    print(f'  {label}: {n:,} pairs ({pdb_ok:,} with PDB = {pdb_ok/n*100:.1f}%)')

print()
