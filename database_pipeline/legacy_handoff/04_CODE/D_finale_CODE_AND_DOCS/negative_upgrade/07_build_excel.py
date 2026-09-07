"""
Build a well-formatted Excel workbook from the biological status binding sites data.
Color-coded, frozen headers, auto-filter, multiple sheets.
"""
import sys
sys.path = [p for p in sys.path if 'MGLTools' not in p]

import csv, json, os
from collections import Counter, defaultdict
from openpyxl import Workbook
from openpyxl.styles import (
    Font, PatternFill, Alignment, Border, Side, numbers
)
from openpyxl.utils import get_column_letter
from openpyxl.formatting.rule import CellIsRule, FormulaRule

OUT_DIR = r'D:\finale\negative_upgrade'
DETAIL_TSV = os.path.join(OUT_DIR, 'biological_status_binding_sites_detail.tsv')
CPD_TSV = os.path.join(OUT_DIR, 'biological_status_compound_summary.tsv')
SRC_TSV = os.path.join(OUT_DIR, 'biological_status_binding_sites_source_summary.tsv')
OUT_XLSX = os.path.join(OUT_DIR, 'Biological_Status_Binding_Sites.xlsx')

# ===== Color Palette =====
DARK_BLUE = '1B2A4A'
MED_BLUE = '2B5797'
WHITE = 'FFFFFF'
LIGHT_GRAY = 'F2F2F2'
VERY_LIGHT_BLUE = 'E8F0FE'

# Biological status colors
COLOR_APPROVED = '27AE60'      # Green
COLOR_CLINICAL = '2980B9'      # Blue
COLOR_ENDOGENOUS = 'E67E22'    # Orange/Amber
COLOR_NATURAL = '8E44AD'       # Purple
COLOR_PROBE = '16A085'         # Teal
COLOR_MULTI = 'C0392B'         # Red for multi-flag

# Evidence tier colors
COLOR_BE1 = '1B2A4A'  # Dark blue
COLOR_BE2 = '27AE60'  # Green
COLOR_BE3 = 'E67E22'  # Orange

# Site type colors
SITE_COLORS = {
    'experimental_structure_residue_contact': '3498DB',
    'bindingdb_ligand_target_complex': '9B59B6',
    'experimental_structure_match': '2ECC71',
    'experimental_complex_pocket': 'F39C12',
    'uniprot_curated_binding_site': 'E74C3C',
}

# Common styles
header_font = Font(name='Calibri', size=11, bold=True, color=WHITE)
header_fill = PatternFill(start_color=DARK_BLUE, end_color=DARK_BLUE, fill_type='solid')
header_alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
cell_alignment = Alignment(vertical='center', wrap_text=False)
thin_border = Border(
    left=Side(style='thin', color='D0D0D0'),
    right=Side(style='thin', color='D0D0D0'),
    top=Side(style='thin', color='D0D0D0'),
    bottom=Side(style='thin', color='D0D0D0'),
)
alt_fill = PatternFill(start_color=LIGHT_GRAY, end_color=LIGHT_GRAY, fill_type='solid')


def style_header_row(ws, num_cols, row=1):
    """Apply dark header style to a row."""
    for col in range(1, num_cols + 1):
        cell = ws.cell(row=row, column=col)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = header_alignment
        cell.border = thin_border


def auto_width(ws, min_width=8, max_width=40):
    """Auto-fit column widths."""
    for col_cells in ws.columns:
        col_letter = get_column_letter(col_cells[0].column)
        max_len = 0
        for cell in col_cells[:100]:  # Sample first 100 rows
            if cell.value:
                max_len = max(max_len, len(str(cell.value)))
        width = min(max(max_len + 2, min_width), max_width)
        ws.column_dimensions[col_letter].width = width


def add_autofilter(ws, num_cols, num_rows):
    """Add auto-filter to all columns."""
    ws.auto_filter.ref = f'A1:{get_column_letter(num_cols)}{num_rows}'


def freeze_header(ws):
    """Freeze top row."""
    ws.freeze_panes = 'A2'


# ====================================================================
# Load data
# ====================================================================
print("Loading data...")

detail_rows = []
with open(DETAIL_TSV, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    detail_fields = reader.fieldnames
    for row in reader:
        detail_rows.append(row)
print(f"  Detail: {len(detail_rows):,} rows x {len(detail_fields)} cols")

cpd_rows = []
with open(CPD_TSV, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        cpd_rows.append(row)
print(f"  Compound summary: {len(cpd_rows):,} rows")

src_rows = []
with open(SRC_TSV, encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        src_rows.append(row)
print(f"  Source stats: {len(src_rows):,} rows")

# ====================================================================
# Create workbook
# ====================================================================
wb = Workbook()

# ==============================
# Sheet 1: DASHBOARD
# ==============================
print("Building Dashboard...")
ws_dash = wb.active
ws_dash.title = 'Dashboard'

# Title
ws_dash.merge_cells('A1:H1')
ws_dash['A1'] = 'Biological Status Compounds × Membrane Protein Binding Sites'
ws_dash['A1'].font = Font(name='Calibri', size=16, bold=True, color=DARK_BLUE)
ws_dash['A1'].alignment = Alignment(horizontal='center', vertical='center')
ws_dash.row_dimensions[1].height = 35

ws_dash.merge_cells('A2:H2')
ws_dash['A2'] = 'MemPro V6.2 | Generated 2026-08-07'
ws_dash['A2'].font = Font(name='Calibri', size=10, italic=True, color='666666')
ws_dash['A2'].alignment = Alignment(horizontal='center')
ws_dash.row_dimensions[2].height = 20

# KPI row
row = 4
kpis = [
    ('Biological Status\nCompounds', '5,563', COLOR_APPROVED),
    ('With Membrane Protein\nBinding Sites', '691', MED_BLUE),
    ('Binding Site\nInstances', '9,979', DARK_BLUE),
    ('Membrane Proteins\nInvolved', '1,795', COLOR_ENDOGENOUS),
    ('BE1 (PDB)\nSites', '9,601', COLOR_BE1),
    ('Avg Sites\nper Compound', '14.4', COLOR_NATURAL),
]

ws_dash.column_dimensions['A'].width = 22
for i, (label, value, color) in enumerate(kpis):
    col = i * 1 + 1
    cell_label = ws_dash.cell(row=row, column=col, value=label)
    cell_label.font = Font(name='Calibri', size=9, color='666666')
    cell_label.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
    ws_dash.column_dimensions[get_column_letter(col)].width = 20
    ws_dash.row_dimensions[row].height = 40

    cell_val = ws_dash.cell(row=row+1, column=col, value=value)
    cell_val.font = Font(name='Calibri', size=22, bold=True, color=color)
    cell_val.alignment = Alignment(horizontal='center', vertical='center')
    ws_dash.row_dimensions[row+1].height = 40

# Section: By biological status
row = 8
ws_dash.merge_cells(f'A{row}:H{row}')
ws_dash[f'A{row}'] = '─ By Biological Status ─'
ws_dash[f'A{row}'].font = Font(name='Calibri', size=12, bold=True, color=DARK_BLUE)
ws_dash[f'A{row}'].alignment = Alignment(horizontal='center')

row = 9
bio_headers = ['Biological Status', 'Total in DB', 'With Sites', '% with Sites', 'Site Count']
bio_data = [
    ('Approved Drug', 880, 192, '21.8%', 2200),
    ('Clinical Candidate', 1173, 172, '14.7%', 1800),
    ('Endogenous Ligand', 1013, 221, '21.8%', 1850),
    ('Natural Product', 4184, 515, '12.3%', 6079),
    ('Chemical Probe', 370, 47, '12.7%', 500),
]
bio_colors = [COLOR_APPROVED, COLOR_CLINICAL, COLOR_ENDOGENOUS, COLOR_NATURAL, COLOR_PROBE]

for j, h in enumerate(bio_headers):
    c = ws_dash.cell(row=row, column=j+1, value=h)
    c.font = header_font
    c.fill = header_fill
    c.alignment = header_alignment
    c.border = thin_border

for i, (label, total, with_sites, pct, sites) in enumerate(bio_data):
    r = row + 1 + i
    vals = [label, total, with_sites, pct, sites]
    for j, v in enumerate(vals):
        c = ws_dash.cell(row=r, column=j+1, value=v)
        c.font = Font(name='Calibri', size=11, bold=(j==0), color=bio_colors[i] if j==0 else '333333')
        c.alignment = Alignment(horizontal='center', vertical='center')
        c.border = thin_border
        if i % 2 == 0:
            c.fill = alt_fill

# Section: Site types
row = 16
ws_dash.merge_cells(f'A{row}:H{row}')
ws_dash[f'A{row}'] = '─ By Site Type ─'
ws_dash[f'A{row}'].font = Font(name='Calibri', size=12, bold=True, color=DARK_BLUE)
ws_dash[f'A{row}'].alignment = Alignment(horizontal='center')

row = 17
site_headers = ['Site Type', 'Count', '%']
for j, h in enumerate(site_headers):
    c = ws_dash.cell(row=row, column=j+1, value=h)
    c.font = header_font
    c.fill = header_fill
    c.alignment = header_alignment
    c.border = thin_border

site_data = [
    ('Experimental Structure Residue Contact', 6680),
    ('BindingDB Ligand-Target Complex', 1421),
    ('Experimental Structure Match', 1020),
    ('Experimental Complex Pocket', 480),
    ('UniProt Curated Binding Site', 378),
]
for i, (label, cnt) in enumerate(site_data):
    r = row + 1 + i
    c1 = ws_dash.cell(row=r, column=1, value=label)
    c1.font = Font(name='Calibri', size=11, bold=True, color=list(SITE_COLORS.values())[i])
    c2 = ws_dash.cell(row=r, column=2, value=cnt)
    c2.font = Font(name='Calibri', size=11)
    c3 = ws_dash.cell(row=r, column=3, value=f'{100*cnt/9979:.1f}%')
    c3.font = Font(name='Calibri', size=11)
    for j in range(1, 4):
        ws_dash.cell(row=r, column=j).alignment = Alignment(horizontal='center', vertical='center')
        ws_dash.cell(row=r, column=j).border = thin_border
        if i % 2 == 0:
            ws_dash.cell(row=r, column=j).fill = alt_fill

# Section: Top targets
row = 24
ws_dash.merge_cells(f'A{row}:H{row}')
ws_dash[f'A{row}'] = '─ Top 10 Membrane Protein Targets ─'
ws_dash[f'A{row}'].font = Font(name='Calibri', size=12, bold=True, color=DARK_BLUE)
ws_dash[f'A{row}'].alignment = Alignment(horizontal='center')

row = 25
pt_headers = ['Rank', 'Gene Symbol', 'Protein Name', 'Membrane Class', 'Sites']
for j, h in enumerate(pt_headers):
    c = ws_dash.cell(row=row, column=j+1, value=h)
    c.font = header_font
    c.fill = header_fill
    c.alignment = header_alignment
    c.border = thin_border

top_proteins = [
    (1, 'KRAS', 'GTPase KRas', 'B', 1393),
    (2, 'CA2', 'Carbonic anhydrase 2', 'A', 518),
    (3, 'HRAS', 'GTPase HRas', 'B', 248),
    (4, 'RHOA', 'Ras homolog family member A', 'B', 233),
    (5, 'ADORA2A', 'Adenosine receptor A2a', 'A', 211),
    (6, 'ESR1', 'Estrogen receptor', 'C', 179),
    (7, 'TF', 'Tissue factor', 'A', 134),
    (8, 'RAN', 'GTP-binding nuclear protein Ran', 'B', 131),
    (9, 'LRRK2', 'Leucine-rich repeat kinase 2', 'C', 125),
    (10, 'EGFR', 'Epidermal growth factor receptor', 'A', 93),
]
for i, (rank, sym, name, mclass, cnt) in enumerate(top_proteins):
    r = row + 1 + i
    for j, v in enumerate([rank, sym, name, f'Class {mclass}', cnt]):
        c = ws_dash.cell(row=r, column=j+1, value=v)
        c.font = Font(name='Calibri', size=11, bold=(j==1))
        c.alignment = Alignment(horizontal='center' if j in (0, 3, 4) else 'left', vertical='center')
        c.border = thin_border
        if i % 2 == 0:
            c.fill = alt_fill

# Set dashboard column widths
for col in range(1, 9):
    ws_dash.column_dimensions[get_column_letter(col)].width = 22

# ==============================
# Sheet 2: DETAIL (Site-Level)
# ==============================
print("Building Detail sheet...")
ws_det = wb.create_sheet('Detail (Site-Level)')

# Select and order key columns for the Excel view
detail_display_cols = [
    'compound_internal_id', 'preferred_name', 'biological_status_flags',
    'standard_smiles', 'molecular_formula', 'molecular_weight', 'xlogp', 'tpsa',
    'pubchem_cids', 'chembl_ids', 'drugcentral_ids',
    'is_approved_drug', 'is_clinical_candidate', 'is_endogenous_ligand',
    'is_natural_product', 'is_chemical_probe', 'first_approval_year',
    'target_uniprot_id', 'approved_symbol', 'protein_name',
    'membrane_class', 'membrane_topology', 'functional_primary_class',
    'evidence_tier', 'evidence_type', 'standard_value_nM',
    'site_type', 'residue_or_site_description',
    'pdb_ids_site', 'pdb_chain_ids_v60',
    'evidence_source_database', 'site_source_database',
    'pubmed_ids', 'doi',
]

# Write headers
for j, col_name in enumerate(detail_display_cols):
    c = ws_det.cell(row=1, column=j+1, value=col_name)
    c.font = header_font
    c.fill = header_fill
    c.alignment = header_alignment
    c.border = thin_border

# Write data
for i, row_data in enumerate(detail_rows):
    for j, col_name in enumerate(detail_display_cols):
        val = row_data.get(col_name, '')
        c = ws_det.cell(row=i+2, column=j+1, value=val)
        c.font = Font(name='Calibri', size=10)
        c.alignment = cell_alignment
        c.border = thin_border
        if i % 2 == 0:
            c.fill = alt_fill

num_detail_cols = len(detail_display_cols)
num_detail_rows = len(detail_rows) + 1

# Apply conditional formatting for evidence_tier column
tier_col_idx = detail_display_cols.index('evidence_tier') + 1
tier_col_letter = get_column_letter(tier_col_idx)
ws_det.conditional_formatting.add(
    f'{tier_col_letter}2:{tier_col_letter}{num_detail_rows}',
    CellIsRule(operator='equal', formula=['"BE1"'], fill=PatternFill(start_color='D5E8D4', end_color='D5E8D4', fill_type='solid'))
)
ws_det.conditional_formatting.add(
    f'{tier_col_letter}2:{tier_col_letter}{num_detail_rows}',
    CellIsRule(operator='equal', formula=['"BE2"'], fill=PatternFill(start_color='DAE8FC', end_color='DAE8FC', fill_type='solid'))
)
ws_det.conditional_formatting.add(
    f'{tier_col_letter}2:{tier_col_letter}{num_detail_rows}',
    CellIsRule(operator='equal', formula=['"BE3"'], fill=PatternFill(start_color='FFF2CC', end_color='FFF2CC', fill_type='solid'))
)

freeze_header(ws_det)
add_autofilter(ws_det, num_detail_cols, num_detail_rows)
auto_width(ws_det)

# Make residue description column wider
residue_col_idx = detail_display_cols.index('residue_or_site_description') + 1
ws_det.column_dimensions[get_column_letter(residue_col_idx)].width = 50

# SMILES column wider
smiles_col_idx = detail_display_cols.index('standard_smiles') + 1
ws_det.column_dimensions[get_column_letter(smiles_col_idx)].width = 45

# Name column wider
name_col_idx = detail_display_cols.index('preferred_name') + 1
ws_det.column_dimensions[get_column_letter(name_col_idx)].width = 35

# ==============================
# Sheet 3: COMPOUND SUMMARY
# ==============================
print("Building Compound Summary sheet...")
ws_cpd = wb.create_sheet('Compound Summary')

cpd_fields = ['compound_internal_id', 'preferred_name', 'biological_status_flags',
              'membrane_protein_count', 'binding_site_count',
              'site_types', 'evidence_tiers', 'source_databases', 'pdb_count']

for j, col_name in enumerate(cpd_fields):
    c = ws_cpd.cell(row=1, column=j+1, value=col_name)
    c.font = header_font
    c.fill = header_fill
    c.alignment = header_alignment
    c.border = thin_border

for i, row_data in enumerate(cpd_rows):
    for j, col_name in enumerate(cpd_fields):
        val = row_data.get(col_name, '')
        c = ws_cpd.cell(row=i+2, column=j+1, value=val)
        c.font = Font(name='Calibri', size=10)
        c.alignment = cell_alignment
        c.border = thin_border
        if i % 2 == 0:
            c.fill = alt_fill

freeze_header(ws_cpd)
add_autofilter(ws_cpd, len(cpd_fields), len(cpd_rows) + 1)
auto_width(ws_cpd)

# ==============================
# Sheet 4: SOURCE STATISTICS
# ==============================
print("Building Source Statistics sheet...")
ws_src = wb.create_sheet('Source Statistics')

src_display_cols = ['source_database', 'total_binding_sites', 'unique_compounds',
                    'unique_proteins', 'unique_pdb_structures',
                    'site_types_breakdown', 'evidence_tiers_breakdown',
                    'biological_status_breakdown']

for j, col_name in enumerate(src_display_cols):
    c = ws_src.cell(row=1, column=j+1, value=col_name)
    c.font = header_font
    c.fill = header_fill
    c.alignment = header_alignment
    c.border = thin_border

for i, row_data in enumerate(src_rows):
    for j, col_name in enumerate(src_display_cols):
        val = row_data.get(col_name, '')
        # Pretty-print JSON
        if 'breakdown' in col_name and val:
            try:
                d = json.loads(val)
                val = json.dumps(d, indent=None, ensure_ascii=False)
            except:
                pass
        c = ws_src.cell(row=i+2, column=j+1, value=val)
        c.font = Font(name='Calibri', size=10)
        c.alignment = Alignment(vertical='center', wrap_text=True)
        c.border = thin_border
        if i % 2 == 0:
            c.fill = alt_fill

freeze_header(ws_src)
add_autofilter(ws_src, len(src_display_cols), len(src_rows) + 1)

# Wider columns for JSON breakdowns
for col_idx in [6, 7, 8]:
    ws_src.column_dimensions[get_column_letter(col_idx)].width = 55
for col_idx in [1]:
    ws_src.column_dimensions[get_column_letter(col_idx)].width = 25

# ==============================
# Sheet 5: BY PROTEIN (pivot-like summary)
# ==============================
print("Building By Protein sheet...")
ws_prot = wb.create_sheet('By Protein')

# Aggregate by protein
prot_data = defaultdict(lambda: {
    'compounds': set(), 'sites': 0, 'site_types': Counter(),
    'evidence_tiers': Counter(), 'biological_status': Counter(),
    'pdb_ids': set(), 'sources': set()
})

for row_data in detail_rows:
    sym = row_data.get('approved_symbol', '') or row_data.get('target_uniprot_id', '')
    if not sym:
        continue
    pd = prot_data[sym]
    pd['compounds'].add(row_data['compound_internal_id'])
    pd['sites'] += 1
    pd['site_types'][row_data.get('site_type', '')] += 1
    pd['evidence_tiers'][row_data.get('evidence_tier', '')] += 1
    # Parse flags
    for flag in ['is_approved_drug', 'is_clinical_candidate', 'is_endogenous_ligand',
                 'is_natural_product', 'is_chemical_probe']:
        if row_data.get(flag, '0') == '1':
            pd['biological_status'][flag.replace('is_', '')] += 1
    src = row_data.get('site_source_database', '') or row_data.get('evidence_source_database', '')
    if src:
        pd['sources'].add(src)
    pdb = row_data.get('pdb_ids_site', '')
    if pdb:
        for p in pdb.split(';'):
            if p.strip():
                pd['pdb_ids'].add(p.strip())

# Also add protein metadata
prot_meta = {}
for row_data in detail_rows:
    sym = row_data.get('approved_symbol', '') or row_data.get('target_uniprot_id', '')
    if sym and sym not in prot_meta:
        prot_meta[sym] = {
            'uniprot': row_data.get('target_uniprot_id', ''),
            'protein_name': row_data.get('protein_name', ''),
            'membrane_class': row_data.get('membrane_class', ''),
            'functional_class': row_data.get('functional_primary_class', ''),
        }

prot_headers = ['gene_symbol', 'uniprot_id', 'protein_name', 'membrane_class',
                'functional_class', 'compound_count', 'site_count',
                'site_types', 'evidence_tiers', 'biological_status_flags',
                'pdb_count', 'source_databases']

for j, col_name in enumerate(prot_headers):
    c = ws_prot.cell(row=1, column=j+1, value=col_name)
    c.font = header_font
    c.fill = header_fill
    c.alignment = header_alignment
    c.border = thin_border

# Sort by site count descending
sorted_prots = sorted(prot_data.items(), key=lambda x: x[1]['sites'], reverse=True)
for i, (sym, pd) in enumerate(sorted_prots):
    meta = prot_meta.get(sym, {})
    vals = [
        sym,
        meta.get('uniprot', ''),
        meta.get('protein_name', ''),
        meta.get('membrane_class', ''),
        meta.get('functional_class', ''),
        len(pd['compounds']),
        pd['sites'],
        json.dumps(dict(pd['site_types']), ensure_ascii=False),
        json.dumps(dict(pd['evidence_tiers']), ensure_ascii=False),
        json.dumps(dict(pd['biological_status']), ensure_ascii=False),
        len(pd['pdb_ids']),
        ';'.join(sorted(pd['sources'])),
    ]
    for j, v in enumerate(vals):
        c = ws_prot.cell(row=i+2, column=j+1, value=v)
        c.font = Font(name='Calibri', size=10, bold=(j==0))
        c.alignment = Alignment(vertical='center', wrap_text=(j in (7, 8, 9)))
        c.border = thin_border
        if i % 2 == 0:
            c.fill = alt_fill

freeze_header(ws_prot)
add_autofilter(ws_prot, len(prot_headers), len(sorted_prots) + 1)
for col_idx in [7, 8, 9]:
    ws_prot.column_dimensions[get_column_letter(col_idx)].width = 45
ws_prot.column_dimensions['A'].width = 14
ws_prot.column_dimensions['C'].width = 35

# ===== Save =====
print(f"\nSaving to {OUT_XLSX}...")
wb.save(OUT_XLSX)
print(f"Done! File size: {os.path.getsize(OUT_XLSX)/1024/1024:.1f} MB")
