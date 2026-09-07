"""
Correct tier analysis: BE from evidence table + E (evidence_level_v52) from protein master.
Single pass, accurate pair-level stats.
"""
import csv, gzip
from collections import Counter, defaultdict

EVIDENCE = r'D:\finale\01_正式数据_V6.2\binding_evidence_master_v6_2.tsv.gz'
PROTEIN = r'D:\finale\01_正式数据_V6.2\human_membrane_protein_master_v6_2.tsv'

# ============================================================
# 1. Load protein master: E tier + BE flags
# ============================================================
print('Loading protein master...')
prot_info = {}  # uniprot -> {e_tier, has_be1, has_be2, has_be3, symbol, mp_class}
with open(PROTEIN, encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        up = row['target_uniprot_id']
        prot_info[up] = {
            'e_tier': row.get('evidence_level_v52', ''),
            'has_be1': row.get('has_BE1_structural_binding', '0') == '1',
            'has_be2': row.get('has_BE2_direct_binding', '0') == '1',
            'has_be3': row.get('has_BE3_pharmacology', '0') == '1',
            'symbol': row.get('approved_symbol', ''),
            'scope': row.get('membrane_scope', ''),
        }

print(f'  Loaded {len(prot_info):,} membrane proteins')

# E tier distribution
e_dist = Counter(p['e_tier'] for p in prot_info.values())
for tier in ['E0', 'E1', 'E2', 'E3']:
    print(f'    {tier}: {e_dist.get(tier, 0):,}')

# ============================================================
# 2. Single pass over evidence: build pair-level stats
# ============================================================
print('\nScanning binding evidence (single pass)...')
pair_be = {}  # (cpd, up) -> best BE tier
pair_has_pdb = set()
pair_has_residues = set()
prot_pdb = {}  # up -> best PDB (prefer with residues)

# for PDB fallback collection
prot_pdb_with_res = {}

row_count = 0
tier_order = {'BE1': 3, 'BE2': 2, 'BE3': 1}

with gzip.open(EVIDENCE, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        row_count += 1
        if row_count % 1000000 == 0:
            print(f'  {row_count:,} rows...')

        cpd = row.get('compound_internal_id', '')
        up = row.get('target_uniprot_id', '')
        if up not in prot_info:
            continue

        key = (cpd, up)
        be = row.get('evidence_tier', '')
        be_val = tier_order.get(be, 0)

        # Track best BE per pair
        if key not in pair_be or be_val > tier_order.get(pair_be[key], 0):
            pair_be[key] = be

        # Track PDB availability
        pdb_str = row.get('pdb_ids', '').strip()
        has_res = bool(row.get('binding_site_residues', '').strip())
        if pdb_str:
            pair_has_pdb.add(key)
            for p in pdb_str.split(';'):
                p = p.strip().lower()
                if len(p) == 4 and p[0].isdigit():
                    if up not in prot_pdb_with_res or has_res:
                        prot_pdb_with_res[up] = p
                        if has_res:
                            break  # prefer residue-bearing PDB
        if has_res:
            pair_has_residues.add(key)

# Fill prot_pdb with fallback (any PDB)
for up in prot_info:
    if up in prot_pdb_with_res:
        prot_pdb[up] = prot_pdb_with_res[up]

# Second pass for any PDB (not just residue-bearing)
with gzip.open(EVIDENCE, 'rt', encoding='utf-8') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        up = row.get('target_uniprot_id', '')
        pdb_str = row.get('pdb_ids', '').strip()
        if up not in prot_info or up in prot_pdb:
            continue
        if pdb_str:
            for p in pdb_str.split(';'):
                p = p.strip().lower()
                if len(p) == 4 and p[0].isdigit():
                    prot_pdb[up] = p
                    break

print(f'  Done: {row_count:,} rows scanned')
print(f'  Unique pairs: {len(pair_be):,}')
print(f'  Unique compounds: {len(set(k[0] for k in pair_be)):,}')
print(f'  Unique proteins: {len(set(k[1] for k in pair_be)):,}')
print(f'  Proteins with any PDB: {len(prot_pdb):,}')

# ============================================================
# 3. Cross-tab: E tier (protein) x BE tier (evidence)
# ============================================================
print(f'\n{"="*60}')
print(f'CROSS: Protein E-tier x Pair BE-tier')
print(f'{"="*60}')

cross = Counter()
e_tier_be = defaultdict(Counter)

for (cpd, up), be in pair_be.items():
    e = prot_info[up]['e_tier']
    cross[(e, be)] += 1
    e_tier_be[e][be] += 1

print(f'\n{"E\\BE":>6} {"BE1":>10} {"BE2":>10} {"BE3":>10} {"TOTAL":>10}')
print(f'{"-"*46}')
for e in ['E0', 'E1', 'E2', 'E3']:
    vals = [e_tier_be[e].get(be, 0) for be in ['BE1', 'BE2', 'BE3']]
    row_total = sum(vals)
    row = ' '.join(f'{v:>10,}' for v in vals)
    print(f'{e:>6} {row} {row_total:>10,}')

# Column totals
col_totals = [sum(e_tier_be[e].get(be, 0) for e in ['E0', 'E1', 'E2', 'E3']) for be in ['BE1', 'BE2', 'BE3']]
grand_total = sum(col_totals)
row_str = ' '.join(f'{v:>10,}' for v in col_totals)
print(f'{"-"*46}')
print(f'{"TOTAL":>6} {row_str} {grand_total:>10,}')

# ============================================================
# 4. Scenario analysis: E + BE combinations
# ============================================================
print(f'\n{"="*60}')
print(f'DOCKING SCENARIOS (with PDB availability)')
print(f'{"="*60}')

scenarios = [
    # (label, BE filter, E filter)
    ('ALL pairs',
     None, None),
    ('BE1+BE2, ALL E tiers',
     ['BE1', 'BE2'], None),
    ('BE1+BE2, E1 only (strictest)',
     ['BE1', 'BE2'], ['E1']),
    ('BE1+BE2, E1+E2 (recommended)',
     ['BE1', 'BE2'], ['E1', 'E2']),
    ('BE1+BE2, E1+E2+E3',
     ['BE1', 'BE2'], ['E1', 'E2', 'E3']),
    ('BE1 only, ALL E',
     ['BE1'], None),
    ('BE1 only, E1+E2',
     ['BE1'], ['E1', 'E2']),
    ('BE1 only, E1+E2+E3',
     ['BE1'], ['E1', 'E2', 'E3']),
]

for label, be_filt, e_filt in scenarios:
    keys = []
    for (cpd, up), be in pair_be.items():
        be_ok = (be_filt is None) or (be in be_filt)
        e_ok = (e_filt is None) or (prot_info[up]['e_tier'] in e_filt)
        if be_ok and e_ok:
            keys.append((cpd, up))

    n = len(keys)
    if n == 0:
        print(f'\n  [{label}]: 0 pairs')
        continue

    n_direct_pdb = sum(1 for k in keys if k in pair_has_pdb)
    n_any_pdb = sum(1 for k in keys if k in pair_has_pdb or k[1] in prot_pdb)
    n_res = sum(1 for k in keys if k in pair_has_residues)
    n_prots = len(set(k[1] for k in keys))
    n_cpds = len(set(k[0] for k in keys))

    # E tier breakdown of these pairs
    e_breakdown = Counter(prot_info[k[1]]['e_tier'] for k in keys)
    e_str = ', '.join(f'{t}:{e_breakdown.get(t,0):,}' for t in ['E0','E1','E2','E3'])

    # BE breakdown
    be_breakdown = Counter(pair_be[k] for k in keys)
    be_str = ', '.join(f'{t}:{be_breakdown.get(t,0):,}' for t in ['BE1','BE2','BE3'])

    print(f'\n  [{label}]')
    print(f'    Pairs: {n:,}')
    print(f'    Proteins: {n_prots:,} | Compounds: {n_cpds:,}')
    print(f'    PDB direct: {n_direct_pdb:,} ({n_direct_pdb/n*100:.1f}%)')
    print(f'    PDB any: {n_any_pdb:,} ({n_any_pdb/n*100:.1f}%)')
    print(f'    Has residues: {n_res:,} ({n_res/n*100:.1f}%)')
    print(f'    E breakdown: {e_str}')
    print(f'    BE breakdown: {be_str}')

# ============================================================
# 5. What do E0/E1/E2/E3 mean here?
# ============================================================
print(f'\n{"="*60}')
print(f'E TIER MEANING (evidence_level_v52)')
print(f'{"="*60}')

# Show what proportion of proteins at each E tier have BE evidence
for e in ['E0', 'E1', 'E2', 'E3']:
    prots_in_tier = [up for up, p in prot_info.items() if p['e_tier'] == e]
    n_total = len(prots_in_tier)
    n_with_be1 = sum(1 for up in prots_in_tier if prot_info[up]['has_be1'])
    n_with_be2 = sum(1 for up in prots_in_tier if prot_info[up]['has_be2'])
    n_with_be3 = sum(1 for up in prots_in_tier if prot_info[up]['has_be3'])
    n_with_any = sum(1 for up in prots_in_tier if any([prot_info[up]['has_be1'], prot_info[up]['has_be2'], prot_info[up]['has_be3']]))
    n_pairs = e_tier_be[e].total() if e in e_tier_be else 0

    print(f'\n  {e}: {n_total:,} proteins')
    print(f'    With BE1: {n_with_be1:,} | BE2: {n_with_be2:,} | BE3: {n_with_be3:,} | Any: {n_with_any:,}')
    print(f'    Pairs (cpd x protein): {n_pairs:,}')
    # PDB coverage
    prots_with_pdb = sum(1 for up in prots_in_tier if up in prot_pdb)
    print(f'    Proteins with PDB: {prots_with_pdb:,}/{n_total:,} ({prots_with_pdb/n_total*100:.1f}%)')

# ============================================================
# 6. Recommendation
# ============================================================
print(f'\n{"="*60}')
print(f'RECOMMENDATION')
print(f'{"="*60}')

recs = [
    ('BE1+BE2, E1+E2 (RECOMMENDED)',
     lambda k: pair_be[k] in ('BE1','BE2') and prot_info[k[1]]['e_tier'] in ('E1','E2'),
     '高质量证据 + 高质量蛋白，平衡数量和质量'),
    ('BE1+BE2, E1+E2+E3 (extended)',
     lambda k: pair_be[k] in ('BE1','BE2') and prot_info[k[1]]['e_tier'] in ('E1','E2','E3'),
     '含 E3 蛋白，证据弱但有binding记录'),
    ('BE1+BE2, ALL E (broad)',
     lambda k: pair_be[k] in ('BE1','BE2'),
     '最大范围的高质量证据，不限制蛋白等级'),
]

for label, fn, desc in recs:
    keys = [(cpd, up) for (cpd, up) in pair_be if fn((cpd, up))]
    n = len(keys)
    n_any_pdb = sum(1 for k in keys if k in pair_has_pdb or k[1] in prot_pdb)
    n_prots = len(set(k[1] for k in keys))
    n_cpds = len(set(k[0] for k in keys))
    print(f'\n  {label}')
    print(f'    {n:,} pairs | {n_prots:,} proteins | {n_cpds:,} compounds')
    print(f'    PDB: {n_any_pdb:,} ({n_any_pdb/n*100:.1f}%)')
    print(f'    {desc}')

print()
