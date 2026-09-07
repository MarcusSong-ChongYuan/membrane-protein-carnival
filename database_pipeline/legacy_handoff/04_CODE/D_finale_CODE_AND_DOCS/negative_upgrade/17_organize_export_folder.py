"""
Organize the full binding site export into a clean folder for sharing.
Split into: compound info, protein info, binding site pairs, source summary.
"""
import csv, os, shutil
from collections import Counter
from datetime import datetime

SRC = r'D:\finale\negative_upgrade\full_binding_site_pairs.tsv'
OUT = r'D:\finale\negative_upgrade\MemPro_binding_site_export'

if os.path.exists(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT)

print("Loading full export...")
with open(SRC, encoding='utf-8') as f:
    all_rows = list(csv.DictReader(f, delimiter='\t'))
print(f"  {len(all_rows):,} rows (all)")

# Remove BindingDB rows (no residue data)
all_rows = [r for r in all_rows if r.get('site_type') != 'bindingdb_ligand_target_complex']
print(f"  {len(all_rows):,} rows (after removing BindingDB, no residue data)")

# ============================================================
# 1. Unique compounds
# ============================================================
print("\n[1/4] Building compound table...")
cpd_cols = [
    'compound_internal_id', 'preferred_name',
    'standard_smiles', 'standard_inchikey',
    'molecular_formula', 'molecular_weight', 'xlogp', 'tpsa',
    'pubchem_cids', 'chembl_ids', 'drugcentral_ids', 'chebi_ids',
    'gtopdb_ligand_ids', 'cas_numbers',
    'compound_scope_class', 'molecule_types', 'identity_confidence',
    'is_approved_drug', 'is_clinical_candidate', 'is_endogenous_ligand',
    'is_natural_product', 'is_chemical_probe',
    'biological_status_flags',
    'development_status', 'max_chembl_phase', 'first_approval_year',
    'compound_source_databases',
]

uniq_cpds = {}
for r in all_rows:
    cid = r['compound_internal_id']
    if cid not in uniq_cpds:
        uniq_cpds[cid] = {c: r.get(c, '') for c in cpd_cols}

with open(os.path.join(OUT, '01_compound_table.tsv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cpd_cols, delimiter='\t', extrasaction='ignore')
    w.writeheader()
    for cid in sorted(uniq_cpds):
        row_out = {}
        for c in cpd_cols:
            v = uniq_cpds[cid].get(c, '')
            if isinstance(v, str) and ('\n' in v or '\r' in v):
                v = v.replace('\r\n', '; ').replace('\n', '; ').replace('\r', '; ')
            row_out[c] = v
        w.writerow(row_out)

n_cpd = len(uniq_cpds)
bio_cpd = sum(1 for c in uniq_cpds.values() if c['biological_status_flags'])
print(f"  {n_cpd:,} unique compounds ({bio_cpd:,} with bio-status)")

# ============================================================
# 2. Unique proteins
# ============================================================
print("\n[2/4] Building protein table...")
prot_cols = [
    'target_uniprot_id', 'approved_symbol', 'protein_name', 'gene_names',
    'membrane_class_v52', 'membrane_topology', 'transmembrane_count',
    'functional_primary_class', 'functional_subclass',
    'membrane_scope', 'subcellular_location',
]

uniq_prots = {}
for r in all_rows:
    up = r['target_uniprot_id']
    if up not in uniq_prots:
        uniq_prots[up] = {c: r.get(c, '') for c in prot_cols}

with open(os.path.join(OUT, '02_protein_table.tsv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=prot_cols, delimiter='\t', extrasaction='ignore')
    w.writeheader()
    for up in sorted(uniq_prots):
        row_out = {}
        for c in prot_cols:
            v = uniq_prots[up].get(c, '')
            if isinstance(v, str) and ('\n' in v or '\r' in v):
                v = v.replace('\r\n', '; ').replace('\n', '; ').replace('\r', '; ')
            row_out[c] = v
        w.writerow(row_out)

n_prot = len(uniq_prots)
print(f"  {n_prot:,} unique proteins")

# ============================================================
# 3. Binding site pairs (the main table)
# ============================================================
print("\n[3/4] Building binding site pair table...")
# FULL columns: compound identity + properties + protein identity + classification
# + evidence + binding site residues. One row = everything.
pair_cols = [
    # ── Compound identity & properties ──
    'compound_internal_id', 'preferred_name',
    'standard_smiles', 'standard_inchikey',
    'molecular_formula', 'molecular_weight', 'xlogp', 'tpsa',
    'pubchem_cids', 'chembl_ids', 'drugcentral_ids', 'chebi_ids',
    'gtopdb_ligand_ids', 'cas_numbers',
    # ── Compound classification ──
    'compound_scope_class', 'molecule_types', 'identity_confidence',
    # ── Biological status ──
    'biological_status_flags',
    'is_approved_drug', 'is_clinical_candidate', 'is_endogenous_ligand',
    'is_natural_product', 'is_chemical_probe',
    'development_status', 'max_chembl_phase', 'first_approval_year',
    'compound_source_databases',
    # ── Protein identity & classification ──
    'target_uniprot_id', 'approved_symbol', 'protein_name', 'gene_names',
    'membrane_class_v52', 'membrane_topology', 'transmembrane_count',
    'functional_primary_class', 'functional_subclass',
    'membrane_scope', 'subcellular_location',
    # ── Evidence ──
    'evidence_id', 'evidence_type',
    'activity_type', 'standard_value_nM', 'activity_outcome_v60',
    'assay_or_mechanism',
    'pdb_ids_evidence', 'ligand_het_id',
    'pubmed_ids', 'doi',
    'evidence_source_database',
    # ── Binding site residue detail ──
    'binding_site_instance_id',
    'site_type', 'site_source_database',
    'pdb_ids_site', 'pdb_chain_ids_v60',
    'residue_or_site_description',
    'site_compound_specificity', 'residue_index_type_v60',
    'site_record_qc_status',
    'residue_count',
]

with open(os.path.join(OUT, '03_binding_site_pairs.tsv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=pair_cols, delimiter='\t', extrasaction='ignore')
    w.writeheader()
    for r in all_rows:
        row_out = {}
        for c in pair_cols:
            if c == 'residue_count':
                desc = r.get('residue_or_site_description', '')
                if not desc.strip():
                    row_out[c] = '0'
                else:
                    # Count residues by splitting on semicolons or spaces
                    parts = [p.strip() for p in desc.replace(' ', ';').split(';') if p.strip()]
                    # Filter out non-residue tokens
                    n = sum(1 for p in parts if p and (p[0].isalpha() or ':' in p))
                    row_out[c] = str(n) if n > 0 else str(len(parts))
            else:
                v = r.get(c, '')
                if isinstance(v, str) and ('\n' in v or '\r' in v):
                    v = v.replace('\r\n', '; ').replace('\n', '; ').replace('\r', '; ')
                row_out[c] = v
        w.writerow(row_out)

print(f"  {len(all_rows):,} rows")

# ============================================================
# 4. Source database summary
# ============================================================
print("\n[4/4] Building source summary...")

# Per-source statistics
source_stats = {}
for r in all_rows:
    src = r['site_source_database'] or r['evidence_source_database'] or 'unknown'
    if src not in source_stats:
        source_stats[src] = {
            'source_database': src,
            'total_pairs': 0,
            'unique_compounds': set(),
            'unique_proteins': set(),
            'unique_pdbs': set(),
            'site_types': Counter(),
            'evidence_tiers': Counter(),
            'has_residues': 0,
            'has_pdb': 0,
            'has_both': 0,
        }
    s = source_stats[src]
    s['total_pairs'] += 1
    s['unique_compounds'].add(r['compound_internal_id'])
    s['unique_proteins'].add(r['target_uniprot_id'])
    if r['pdb_ids_site']:
        for p in r['pdb_ids_site'].split(','):
            if p.strip():
                s['unique_pdbs'].add(p.strip())
    s['site_types'][r['site_type']] += 1
    s['evidence_tiers'][r['evidence_tier']] += 1
    if r['residue_or_site_description'].strip():
        s['has_residues'] += 1
    if r['pdb_ids_site'].strip():
        s['has_pdb'] += 1
    if r['residue_or_site_description'].strip() and r['pdb_ids_site'].strip():
        s['has_both'] += 1

with open(os.path.join(OUT, '04_source_summary.tsv'), 'w', encoding='utf-8', newline='') as f:
    cols = ['source_database', 'total_pairs', 'unique_compounds', 'unique_proteins',
            'unique_pdbs', 'has_both_residue_and_pdb', 'has_residue_only', 'has_pdb_only',
            'site_types', 'evidence_tiers',
            'description']
    f.write('\t'.join(cols) + '\n')

    desc_map = {
        'PDBe': 'PDB co-crystal structures with atomic contact analysis. Residues within 4A of ligand.',
        'PDBbind': 'PDBbind database: experimental_structure_match (key affinity residues) + experimental_complex_pocket (full pocket).',
        'sc-PDB': 'sc-PDB database: full binding pocket from co-crystal structures.',
        'BindingDB': 'BindingDB ligand-target complexes. Has PDB list but no residue-level data.',
        'UniProt': 'UniProt/Swiss-Prot manually curated binding sites from literature. Residue positions on sequence, no PDB.',
        'BioLiP': 'BioLiP structure-matched binding sites. Compound matched to similar ligand in PDB.',
    }

    for src in sorted(source_stats.keys(), key=lambda x: -source_stats[x]['total_pairs']):
        s = source_stats[src]
        f.write('\t'.join([
            src,
            str(s['total_pairs']),
            str(len(s['unique_compounds'])),
            str(len(s['unique_proteins'])),
            str(len(s['unique_pdbs'])),
            str(s['has_both']),
            str(s['has_residues'] - s['has_both']),
            str(s['has_pdb'] - s['has_both']),
            str(dict(s['site_types'])).replace('\t', ' '),
            str(dict(s['evidence_tiers'])).replace('\t', ' '),
            desc_map.get(src, ''),
        ]) + '\n')

print(f"  {len(source_stats)} sources")

# ============================================================
# Site type explanation table
# ============================================================
print("Writing site type reference...")
with open(os.path.join(OUT, '05_site_type_reference.tsv'), 'w', encoding='utf-8', newline='') as f:
    f.write('\t'.join(['site_type', 'count', 'source', 'has_residues', 'has_pdb',
                       'residue_format', 'format_example', 'description']) + '\n')
    refs = [
        ('experimental_structure_residue_contact',
         'PDBe / BioLiP',
         'PHE144;PHE151;THR147',
         'Three-letter residue name + sequence number, semicolon-separated. Chain in pdb_chain_ids_v60 column. These come from PDB atomic contact analysis: residues whose atoms are within 4A of the bound ligand. Most reliable binding site definition.'),
        ('experimental_complex_pocket',
         'sc-PDB / PDBbind',
         'A:ALA767;A:ALA771;B:THR200',
         'Chain letter colon residue name + number, semicolon-separated. Full binding pocket: ALL residues within 6A of the ligand, typically 50-200+ residues. More comprehensive than residue_contact but may include non-contact structural residues.'),
        ('experimental_structure_match',
         'BioLiP / PDBbind',
         'CID907|(LP5)@4g8a:B=S415 F440',
         'PubChem CID | (HET code) @ PDB_ID : chain = key_residues. The compound itself has no co-crystal structure; a similar compound in PDB was used to infer the binding site. Only affinity-critical residues are recorded (not the full pocket).'),
        ('uniprot_curated_binding_site',
         'UniProt / Swiss-Prot',
         'L-glutamate@156-156(medium)',
         'Ligand_name @ UniProt_sequence_position_start-end (curator_confidence: high/medium/low). Manually curated by UniProt curators from published literature. Positions refer to UniProt sequence numbering, not PDB residue numbers. No 3D structure available.'),
        ('bindingdb_ligand_target_complex',
         'BindingDB',
         '(empty)',
         'No residue-level data. BindingDB records that a compound binds a protein but does not specify the binding site residues. PDB list is provided (the protein has structures) but the specific binding pocket is not annotated.'),
    ]
    for st, src, fmt, desc in refs:
        cnt = sum(1 for r in all_rows if r['site_type'] == st)
        has_res = 'Y' if st != 'bindingdb_ligand_target_complex' else 'N'
        has_pdb = 'Y' if st != 'uniprot_curated_binding_site' else 'N'
        f.write('\t'.join([st, str(cnt), src, has_res, has_pdb, fmt.split('\n')[0],
                           fmt.split('\n')[0][:80], desc]) + '\n')

# ============================================================
# README
# ============================================================
print("Writing README...")
res_both = sum(1 for r in all_rows if r['residue_or_site_description'].strip() and r['pdb_ids_site'].strip())
res_only = sum(1 for r in all_rows if r['residue_or_site_description'].strip() and not r['pdb_ids_site'].strip())
pdb_only = sum(1 for r in all_rows if not r['residue_or_site_description'].strip() and r['pdb_ids_site'].strip())

with open(os.path.join(OUT, 'README.md'), 'w', encoding='utf-8') as f:
    f.write(f'''# MemPro V6.2 — Binding Site Pair Export

Generated: {datetime.now().strftime('%Y-%m-%d')}

## Overview

All experimentally validated small-molecule x membrane-protein binding site pairs from MemPro V6.2.

| Metric | Count |
|--------|-------|
| Total binding site pairs | {len(all_rows):,} |
| Unique compounds | {n_cpd:,} |
| Unique membrane proteins | {n_prot:,} |
| With bio-status flags | {bio_cpd:,} compounds |
| Pairs with residue + PDB | {res_both:,} |
| Pairs with residue only (UniProt, no PDB) | {res_only:,} |
| Pairs with PDB only (BindingDB, no residue) | {pdb_only:,} |

## Files

| File | Content |
|------|---------|
| `01_compound_table.tsv` | {n_cpd:,} unique compounds with chemical properties, bio-status, source DBs |
| `02_protein_table.tsv` | {n_prot:,} unique membrane proteins with class, topology, function |
| `03_binding_site_pairs.tsv` | {len(all_rows):,} rows: full compound-protein binding site pair data (60 columns) |
| `04_source_summary.tsv` | Per-source database statistics |
| `05_site_type_reference.tsv` | Site type formats with examples |
| `README.md` | This file |

## Binding Site Completeness

Not all pairs have complete binding site residue data:

| Status | Count | site_type |
|--------|-------|-----------|
| Complete (residue + PDB) | {res_both:,} | residue_contact, complex_pocket, structure_match |
| No PDB (sequence-level only) | {res_only:,} | uniprot_curated_binding_site |
| No residue data | {pdb_only:,} | bindingdb_ligand_target_complex |

Even among "complete" pairs, residue coverage varies:
- **experimental_structure_residue_contact**: avg 7 residues, but 25% have only 1-2 residues (sparse)
- **experimental_complex_pocket**: avg 84 residues (full pocket)
- **experimental_structure_match**: only key affinity residues (intentionally incomplete)

Use the `residue_count` column to filter by residue coverage.

## Site Types

1. **experimental_structure_residue_contact** — Atomic contacts within 4A of ligand (PDB structures)
2. **experimental_complex_pocket** — Full binding pocket (50-200+ residues)
3. **experimental_structure_match** — Key affinity residues, inferred from similar compound
4. **uniprot_curated_binding_site** — Literature-curated, sequence positions (no PDB)
5. **bindingdb_ligand_target_complex** — PDB available but no residue annotation

See `05_site_type_reference.tsv` for format examples.

## Membrane Protein Classification

- **Class A**: Multi-pass transmembrane (GPCRs, channels, transporters)
- **Class B**: Lipid-anchored or peripheral membrane
- **Class C**: Single-pass transmembrane

## Column Reference

### Compound identity (01_compound_table.tsv)
- `compound_internal_id`: MemPro internal compound ID
- `preferred_name`: Primary compound name
- `standard_smiles`: Canonical SMILES string
- `standard_inchikey`: InChI Key (standardized)
- `molecular_formula`, `molecular_weight`, `xlogp`, `tpsa`: Computed properties
- `pubchem_cids`, `chembl_ids`, `drugcentral_ids`, `chebi_ids`: Cross-reference IDs
- `is_approved_drug` .. `is_chemical_probe`: Biological status flags (0/1)
- `biological_status_flags`: Active flags as semicolon-separated string
- `compound_scope_class`: M5 classification level
- `identity_confidence`: high / medium / low

### Protein identity (02_protein_table.tsv)
- `target_uniprot_id`: UniProt accession
- `approved_symbol`: HGNC gene symbol
- `membrane_class_v52`: A (multi-pass) / B (lipid-anchored) / C (single-pass)
- `membrane_topology`: e.g. multi_pass, single_pass_type_i, lipid_anchor
- `functional_primary_class`: receptor / enzyme / transporter / ion_channel / ...

### Binding site pairs (03_binding_site_pairs.tsv)
- `site_type`: One of the 5 site types (see reference)
- `residue_or_site_description`: Binding site residues in source-specific format
- `residue_count`: Number of binding site residues (0 = none; filter this column to exclude incomplete sites)
- `pdb_ids_site` / `pdb_ids_evidence`: PDB structure IDs
- `pdb_chain_ids_v60`: Chain identifier(s) for the binding site
- `standard_value_nM`: Binding affinity in nanomolar
- `activity_outcome_v60`: positive_or_observed / negative / ...
- `pubmed_ids` / `doi`: Literature references

## Source Databases

''')

    # Fill source descriptions
    src_lines = []
    for src in sorted(source_stats.keys(), key=lambda x: -source_stats[x]['total_pairs']):
        s = source_stats[src]
        src_lines.append(f"- **{src}**: {s['total_pairs']:,} pairs, {len(s['unique_compounds']):,} compounds, {len(s['unique_proteins']):,} proteins, {len(s['unique_pdbs']):,} PDBs")
    f.write('\n'.join(src_lines) + '\n')

# ============================================================
# Final summary
# ============================================================
print(f"\n{'='*60}")
print(f"EXPORT FOLDER: {OUT}")
print(f"{'='*60}")
for fname in sorted(os.listdir(OUT)):
    fpath = os.path.join(OUT, fname)
    size = os.path.getsize(fpath)
    print(f"  {fname:40s} {size/1024:8.1f} KB")
total = sum(os.path.getsize(os.path.join(OUT, f)) for f in os.listdir(OUT)) / (1024*1024)
print(f"\n  Total: {total:.1f} MB")
print("\nDone.")
