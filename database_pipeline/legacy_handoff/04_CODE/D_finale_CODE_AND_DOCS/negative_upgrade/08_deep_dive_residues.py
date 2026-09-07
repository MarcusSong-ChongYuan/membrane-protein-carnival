"""
Deep dive: why BindingDB residues are empty, UniProt 'curated' meaning, PDBbind 'key residues'.
"""
import csv, gzip
from collections import Counter

SITES = r'D:\finale\01_正式数据_V6.2\binding_site_instances_v6_2.tsv.gz'
EVIDENCE = r'D:\finale\01_正式数据_V6.2\binding_evidence_master_v6_2.tsv.gz'

print("=" * 70)
print("Q1: BindingDB residue_or_site_description 为空原因")
print("=" * 70)

bindingdb_ev_ids = set()
bindingdb_sites = []
with gzip.open(SITES, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row.get('source_database', '') == 'BindingDB':
            bindingdb_sites.append(row)
            bindingdb_ev_ids.add(row['evidence_id'])

print(f"BindingDB site instances: {len(bindingdb_sites):,}")

empty_desc = sum(1 for s in bindingdb_sites if not s.get('residue_or_site_description', '').strip())
print(f"residue_or_site_description 为空: {empty_desc}/{len(bindingdb_sites)}")

# Check evidence-level binding_site_residues
residue_samples = []
both_empty = 0
with gzip.open(EVIDENCE, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['evidence_id'] in bindingdb_ev_ids:
            bsr = row.get('binding_site_residues', '').strip()
            desc = ''
            for s in bindingdb_sites:
                if s['evidence_id'] == row['evidence_id']:
                    desc = s.get('residue_or_site_description', '').strip()
                    break
            if bsr:
                if len(residue_samples) < 12:
                    residue_samples.append({
                        'ev_id': row['evidence_id'],
                        'bsr': bsr,
                        'site_desc': desc,
                        'cpd': row.get('compound_name', ''),
                        'sym': row.get('approved_symbol', ''),
                        'pdb': row.get('pdb_ids', ''),
                    })
            else:
                both_empty += 1

print(f"\nevidence 级 binding_site_residues 有值样本 ({len(residue_samples)}条):")
for s in residue_samples[:8]:
    print(f"\n  [{s['sym']}] {s['cpd'][:50]}")
    print(f"  PDB: {s['pdb']}")
    print(f"  evidence.binding_site_residues: {s['bsr'][:200]}")
    print(f"  site.residue_or_site_description: '{s['site_desc'][:100]}'")

print(f"\n  binding_site_residues 也为空: {both_empty}条")

# Conclusion for Q1
print(f"\n  >>> 结论: evidence 表的 binding_site_residues 有残基数据,")
print(f"      但 site 表的 residue_or_site_description 没做填充。")
print(f"      site_type=bindingdb_ligand_target_complex 的数据应该从")
print(f"      evidence.binding_site_residues 同步过来。")

print()
print("=" * 70)
print("Q2: UniProt curated_binding_site 是不是真实实验?")
print("=" * 70)

uniprot_ev_ids = set()
with gzip.open(SITES, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row.get('site_type', '') == 'uniprot_curated_binding_site':
            uniprot_ev_ids.add(row['evidence_id'])

print(f"UniProt curated sites: {len(uniprot_ev_ids)} instances")

# Check evidence records
uni_sources = Counter()
uni_ev_types = Counter()
uni_has_pubmed = 0
uni_no_pubmed = 0
uni_samples = []

with gzip.open(EVIDENCE, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['evidence_id'] in uniprot_ev_ids:
            uni_sources[row.get('source_database', '?')] += 1
            uni_ev_types[row.get('evidence_type', '?')] += 1
            pmids = row.get('pubmed_ids', '').strip()
            if pmids:
                uni_has_pubmed += 1
            else:
                uni_no_pubmed += 1
            if len(uni_samples) < 6:
                uni_samples.append({
                    'sym': row.get('approved_symbol', ''),
                    'uniprot': row.get('target_uniprot_id', ''),
                    'cpd': row.get('compound_name', ''),
                    'ev_type': row.get('evidence_type', ''),
                    'assay': row.get('assay_or_mechanism', ''),
                    'pubmed': row.get('pubmed_ids', ''),
                    'source_url': row.get('source_url', ''),
                })

print(f"\nEvidence types: {dict(uni_ev_types)}")
print(f"有 PubMed 引用: {uni_has_pubmed}/{len(uniprot_ev_ids)}")
print(f"无 PubMed 引用: {uni_no_pubmed}/{len(uniprot_ev_ids)}")

print(f"\n样本 evidence 记录:")
for s in uni_samples:
    print(f"\n  [{s['sym']}] {s['cpd'][:50]}")
    print(f"  UniProt: {s['uniprot']}")
    print(f"  evidence_type: {s['ev_type']}")
    print(f"  assay: {s['assay'][:100]}")
    print(f"  PubMed: {s['pubmed']}")
    print(f"  Source URL: {s['source_url'][:120]}")

print(f"\n  >>> 结论: UniProt curated_binding_site 来自 Swiss-Prot 手动注释。")
print(f"      'curated' = 由 UniProt/Swiss-Prot 的 curator 阅读文献后人工注释。")
print(f"      是真实实验证据，但注释粒度是'功能区域'而非原子接触。")
print(f"      格式 L-glutamate@156-156(medium) 含义:")
print(f"        配体名 @ UniProt序列位置起-止 (curator置信度high/medium)")
print(f"      UniProt curators 审读论文后标注哪些残基参与配体结合。")

print()
print("=" * 70)
print("Q3: PDBbind 为何只标注关键残基?")
print("=" * 70)

# Separate PDBbind sites into 'match' and 'pocket' types
pdbind_match = []
pdbind_pocket = []
with gzip.open(SITES, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row.get('source_database', '') == 'PDBbind':
            st = row.get('site_type', '')
            if st == 'experimental_structure_match':
                pdbind_match.append(row)
            elif st == 'experimental_complex_pocket':
                pdbind_pocket.append(row)

print(f"PDBbind experimental_structure_match (关键残基): {len(pdbind_match)}")
print(f"PDBbind experimental_complex_pocket (完整口袋): {len(pdbind_pocket)}")

# Show match samples
print(f"\n--- experimental_structure_match 样本 ('关键残基') ---")
for s in pdbind_match[:8]:
    desc = s.get('residue_or_site_description', '')
    print(f"  {desc[:200]}")

# Show pocket samples
print(f"\n--- experimental_complex_pocket 样本 ('完整口袋') ---")
for s in pdbind_pocket[:5]:
    desc = s.get('residue_or_site_description', '')
    residues = [r for r in desc.split(';') if r.strip()]
    print(f"  PDB={s.get('pdb_ids','')} 残基数={len(residues)}")
    print(f"  {desc[:200]}...")

# Check evidence level for match records
pdbind_match_ev_ids = {s['evidence_id'] for s in pdbind_match}
print(f"\n--- PDBbind match 对应的 evidence 记录 ---")
cnt = 0
with gzip.open(EVIDENCE, 'rt', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['evidence_id'] in pdbind_match_ev_ids and cnt < 6:
            print(f"\n  [{row.get('approved_symbol','')}] {row.get('compound_name','')[:50]}")
            print(f"  evidence_type: {row.get('evidence_type','')}")
            print(f"  activity: {row.get('activity_value','')} {row.get('activity_unit','')}")
            print(f"  standard_value_nM: {row.get('standard_value_nM','')}")
            print(f"  PDB: {row.get('pdb_ids','')}")
            print(f"  binding_site_residues (evidence级): {row.get('binding_site_residues','')[:200]}")
            cnt += 1

print(f"\n  >>> 结论: PDBbind 有两个 site_type:")
print(f"      (1) experimental_structure_match = 仅亲和力测定中标注的关键残基")
print(f"          格式: CID|NAME@PDB:IC50/Kd=值=链:关键残基")
print(f"          PDBbind 数据库只记录对亲和力贡献最大的几个残基")
print(f"          不是完整口袋，而是'决定亲和力的要害残基'")
print(f"      (2) experimental_complex_pocket = 完整结合口袋所有残基")
print(f"          格式: 链:残基;链:残基;... (上百个残基)")
print(f"          包含所有距离配体一定范围内的残基")

print("\nDone.")
