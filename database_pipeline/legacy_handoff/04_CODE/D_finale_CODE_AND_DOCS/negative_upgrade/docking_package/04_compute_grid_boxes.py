#!/usr/bin/env python3
"""
Compute Vina grid boxes from binding site residues.  V4.2

Mapping priority (tried in order; stop when observed >= 50% of requested):
  1. Direct author-number match (manifest residue numbers already PDB numbering)
  2. SIFTS per-residue mapping (from sifts_mappings.tsv, if present)
  3. DBREF segment-based mapping (equal-length, no-insertion segments only)
  4. Same-UniProt homologous chain candidate → REVIEW, no grid written

Hard rules:
  - REVIEW / FAIL / NO_SITE_DATA tasks → grid fields LEFT BLANK
  - Cross-chain is NEVER auto-trusted (always REVIEW, never MEDIUM)
  - DBREF non-equal-length segments are NEVER extrapolated
  - DBREF cross-chain borrowing is forbidden
  - Direct match validates chain accession against DBREF when available
  - mapped_residue_count = residues that have a valid mapping (not same as observed)
"""
import csv, os, re, glob, gzip, json
from collections import defaultdict
from xml.etree import ElementTree as ET

PDB_DIR    = "pdb_structures"
MANIFEST   = "docking_manifest.tsv"
SIFTS_TSV  = "sifts_mappings.tsv"         # optional – only used if present
MARGIN     = 10.0
MIN_SIZE   = 20.0

# 1-letter → 3-letter amino acid codes
AA1TO3 = { 'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY',
           'H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN',
           'P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL',
           'W':'TRP','Y':'TYR' }

# ═══════════════════════════════════════════════════════════════
#  SIFTS loader
# ═══════════════════════════════════════════════════════════════

def load_sifts(tsv_path):
    """Parse sifts_mappings.tsv into a dict.
       Key:   (pdb_id_upper, author_chain, uniprot_acc, uniprot_resnum)
       Value: (pdb_resnum, insertion_code, pdb_resname, observed_flag)

       Skip comment/empty lines.  Robust to missing columns."""
    sifts = {}
    if not os.path.exists(tsv_path):
        return sifts

    with open(tsv_path, encoding='utf-8', newline='') as f:
        header = None
        for line in f:
            line = line.rstrip('\r\n')
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if header is None:
                header = cols  # first data-like row is header
                continue
            if len(cols) < 6:
                continue

            pdb_id  = cols[0].strip().upper()
            chain   = cols[1].strip().upper()
            # cols[2] = struct_asym_id (skip)
            uniprot = cols[3].strip()
            uniprot_num_str = cols[4].strip()
            pdb_num_str     = cols[5].strip()
            icode   = cols[6].strip() if len(cols) > 6 else ''
            resname = cols[7].strip() if len(cols) > 7 else ''
            obs     = cols[8].strip() if len(cols) > 8 else ''

            if not uniprot or not pdb_id or not chain:
                continue

            try:
                uniprot_num = int(uniprot_num_str)
                pdb_resnum  = int(pdb_num_str)
            except (ValueError, AttributeError):
                continue

            key = (pdb_id, chain, uniprot, uniprot_num)
            sifts[key] = (pdb_resnum, icode, resname.upper(), obs)
    return sifts


# ═══════════════════════════════════════════════════════════════
#  DBREF loader (handles DBREF / DBREF1+DBREF2 pairs)
# ═══════════════════════════════════════════════════════════════
#
#  PDB format v3.30 column positions (1-indexed → 0-indexed):
#
#            DBREF / DBREF1               DBREF2
#  chain     col 13  →  line[12]         col 13  →  line[12]
#  pdb_b     cols 15-18 → line[14:18]    cols 15-18 → line[14:18]
#  pdb_bi    col 19 → line[18]           col 19 → line[18]
#  pdb_e     cols 21-24 → line[20:24]    cols 21-24 → line[20:24]
#  pdb_ei    col 25 → line[24]           col 25 → line[24]
#  db_acc    cols 34-41 → line[33:41]    -- not present --
#  db_b      cols 56-60 → line[55:60]    cols 56-60 → line[55:60]
#  db_bi     col 61 → line[60]           col 61 → line[60]
#  db_e      cols 63-67 → line[62:67]    cols 63-67 → line[62:67]
#  db_ei     col 68 → line[67]           col 68 → line[67]

def build_dbref_map(pdb_path):
    """Return {chain: {'accessions': set, 'segments': [seg, ...]}}.

       DBREF1 and DBREF2 are complementary paired records: DBREF1 carries
       the database accession, DBREF2 carries the database residue range.
       They are paired by matching (chain, pdb_begin, pdb_end).

       Standalone DBREF records carry both — used directly."""
    raw_dbref  = []   # standalone DBREF
    raw_dbref1 = []   # pairs with DBREF2
    raw_dbref2 = []   # pairs with DBREF1

    with open(pdb_path, encoding='utf-8', errors='ignore') as f:
        for line in f:
            rectype = line[0:6].strip()
            if rectype not in ('DBREF', 'DBREF1', 'DBREF2'):
                continue
            if len(line) < 68:
                continue

            try:
                ch = line[12].strip().upper()
                if not ch:
                    continue

                rec = {
                    'type':    rectype,
                    'chain':   ch,
                    'pdb_b':   int(line[14:18].strip()),
                    'pdb_bi':  line[18].strip() if len(line) > 18 else '',
                    'pdb_e':   int(line[20:24].strip()),
                    'pdb_ei':  line[24].strip() if len(line) > 24 else '',
                    'db_b':    int(line[55:60].strip()),
                    'db_bi':   line[60].strip() if len(line) > 60 else '',
                    'db_e':    int(line[62:67].strip()),
                    'db_ei':   line[67].strip() if len(line) > 67 else '',
                }

                if rectype == 'DBREF':
                    rec['db_acc'] = line[33:41].strip()
                elif rectype == 'DBREF1':
                    rec['db_acc'] = line[33:41].strip()
                    # DBREF1 carries accession; DBREF2 carries db range
                    rec['db_b'] = -999999   # placeholder – real values from DBREF2
                    rec['db_e'] = -999999
                else:  # DBREF2
                    rec['db_acc'] = ''      # DBREF2 has no accession
            except (ValueError, IndexError):
                continue

            if rectype == 'DBREF':
                raw_dbref.append(rec)
            elif rectype == 'DBREF1':
                raw_dbref1.append(rec)
            else:
                raw_dbref2.append(rec)

    mapping = {}

    # ── Standalone DBREF: use directly ──
    for rec in raw_dbref:
        ch = rec['chain']
        mapping.setdefault(ch, {'accessions': set(), 'segments': []})
        mapping[ch]['accessions'].add(rec['db_acc'])
        mapping[ch]['segments'].append(_make_seg(rec))

    # ── Pair DBREF1 + DBREF2 by (chain, pdb_begin, pdb_end) ──
    # Index DBREF2 by (chain, pdb_begin, pdb_end)
    dbref2_by_key = {}
    for rec in raw_dbref2:
        key = (rec['chain'], rec['pdb_b'], rec['pdb_e'])
        dbref2_by_key[key] = rec

    for rec1 in raw_dbref1:
        key = (rec1['chain'], rec1['pdb_b'], rec1['pdb_e'])
        rec2 = dbref2_by_key.get(key)
        if not rec2:
            continue  # unpaired DBREF1 — skip, can't resolve without DBREF2

        ch = rec1['chain']
        mapping.setdefault(ch, {'accessions': set(), 'segments': []})
        mapping[ch]['accessions'].add(rec1['db_acc'])

        # Merge: DBREF1 has accession, DBREF2 has database residue range
        merged = dict(rec1)
        merged['db_b']  = rec2['db_b']
        merged['db_bi'] = rec2['db_bi']
        merged['db_e']  = rec2['db_e']
        merged['db_ei'] = rec2['db_ei']
        mapping[ch]['segments'].append(_make_seg(merged))

    return mapping


def _make_seg(rec):
    """Convert a DBREF record dict into a segment dict."""
    return {
        'pdb_begin':     rec['pdb_b'],
        'pdb_begin_ins': rec['pdb_bi'],
        'pdb_end':       rec['pdb_e'],
        'pdb_end_ins':   rec['pdb_ei'],
        'db_begin':      rec['db_b'],
        'db_begin_ins':  rec['db_bi'],
        'db_end':        rec['db_e'],
        'db_end_ins':    rec['db_ei'],
        'accession':     rec.get('db_acc', ''),
    }


# ═══════════════════════════════════════════════════════════════
#  Residue parser  (handles 3-letter and 1-letter codes)
# ═══════════════════════════════════════════════════════════════

def parse_residue_parts(desc):
    """Parse residue descriptor strings like 'TYR598;TYR604;VAL654'.
       Returns list of (aa3, resnum, icode) tuples."""
    results = []
    if not desc or not desc.strip():
        return results
    desc = desc.strip()

    # Check for PDBbind range notation: "456-470@..." -- skip range
    if '@' in desc and '-' in desc.split('@')[-1].split('(')[0]:
        return results

    if '@' in desc:
        desc = desc.split('@')[-1]
    if '=' in desc:
        _, desc = desc.split('=', 1)
    elif ':' in desc:
        # Two distinct formats use colons:
        #   A) Chain-prefixed:  A:GLN71;A:TYR72;...  (':' per residue token)
        #   B) PDB/chain header: 6o4x:A=GLN71 TYR72   (one ':' in prefix)
        tokens_test = [t for t in desc.replace(';', ' ').split() if t.strip()]
        n_colon_tokens = sum(1 for t in tokens_test if ':' in t)
        # When nearly every token contains a colon, it's chain-prefixed
        if n_colon_tokens >= len(tokens_test) * 0.5 and len(tokens_test) >= 2:
            cleaned = []
            for t in tokens_test:
                if ':' in t:
                    _, residue = t.split(':', 1)
                    cleaned.append(residue.strip())
                else:
                    cleaned.append(t.strip())
            desc = ' '.join(cleaned)
        else:
            parts = desc.split(':')
            if len(parts) >= 2:
                desc = parts[-1]

    tokens = [t.strip() for t in desc.replace(';', ' ').split() if t.strip()]
    for token in tokens:
        m = re.match(r'([A-Za-z]{1,3})(-?\d+)([A-Za-z]?)', token)
        if not m:
            continue
        aa_code = m.group(1).upper()
        resnum  = int(m.group(2))
        icode   = m.group(3).strip()

        if len(aa_code) == 3:
            aa3 = aa_code
        elif aa_code in AA1TO3:
            aa3 = AA1TO3[aa_code]
        else:
            continue

        results.append((aa3, resnum, icode))
    return results


# ═══════════════════════════════════════════════════════════════
#  Atom index builder + residue matcher
# ═══════════════════════════════════════════════════════════════

def build_atom_index(pdb_path):
    """Build index of all heavy/backbone atoms in a PDB file.
       Returns: {(chain, resname, resnum, icode): [(x, y, z, is_heavy)]}"""
    idx = defaultdict(list)
    with open(pdb_path, encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not (line.startswith('ATOM  ') or line.startswith('HETATM')):
                continue
            if len(line) < 55:
                continue
            try:
                ch      = line[21:22].strip().upper()
                rn      = line[17:20].strip().upper()
                rn_num  = int(line[22:26].strip())
                rn_ico  = line[26:27].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                aname   = line[12:16].strip()
                heavy   = aname not in ('N', 'CA', 'C', 'O', 'H', 'HA', 'HN')
            except (ValueError, IndexError):
                continue
            idx[(ch, rn, rn_num, rn_ico)].append((x, y, z, heavy))
    return idx


def match_residues(atom_index, residue_list, chain_filter):
    """Match a list of (aa3, resnum, icode) against the atom index.

       Returns (coords, mapped, observed):
         coords:   list of (x,y,z) for matched atoms (heavy preferred, then CA)
         mapped:   number of input residues that could be mapped (i.e. had a
                   valid residue number in the mapping scheme)
         observed: number of mapped residues that actually have atom coordinates

       Note: mapped ≤ len(residue_list), observed ≤ mapped."""
    coords   = []
    mapped   = 0
    observed = 0

    for aa3, resnum, icode in residue_list:
        found = False
        # Search for exact match within chain filter
        for (ch, rn, rn_num, rn_ico), atoms in atom_index.items():
            if chain_filter and ch != chain_filter.upper():
                continue
            if rn == aa3 and rn_num == resnum and rn_ico == icode:
                mapped   += 1
                observed += 1
                found     = True
                # Prefer heavy atoms for centroid
                heavy = [(x, y, z) for x, y, z, h in atoms if h]
                ca    = [(x, y, z) for x, y, z, h in atoms if not h]
                coords.extend(heavy if heavy else ca)
                break

    return coords, mapped, observed


# ═══════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════

print("=== Loading manifest ===")
tasks = []
with open(MANIFEST, encoding='utf-8', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        tasks.append(row)
n_tasks = len(tasks)
print(f"Tasks: {n_tasks}")

# ── Load SIFTS (optional) ──
print("\n=== Loading SIFTS ===")
sifts = load_sifts(SIFTS_TSV)
print(f"SIFTS residue-level mappings loaded: {len(sifts)}")

# ── Build DBREF maps ──
print("\n=== Building DBREF maps ===")
dbref_cache = {}
pdb_files = sorted(glob.glob(os.path.join(PDB_DIR, '*.pdb')))
for fp in pdb_files:
    pdb_id = os.path.basename(fp).replace('.pdb', '')
    dbref_cache[pdb_id] = build_dbref_map(fp)
n_dbref = sum(1 for v in dbref_cache.values() if v)
print(f"PDB files: {len(pdb_files)}  with DBREF: {n_dbref}")

# ── Ensure output columns exist ──
fieldnames = list(tasks[0].keys())
for col in ['grid_method', 'mapped_residue_count', 'observed_residue_count',
            'requested_residue_count', 'grid_qc_status', 'grid_confidence',
            'grid_fallback_reason', 'mapping_source', 'actual_chain_used']:
    if col not in fieldnames:
        fieldnames.append(col)

grid_file = "docking_manifest_with_grid.tsv"
audit_file = "grid_audit_report.tsv"

qc_stats   = defaultdict(int)
written_ok = 0

with open(grid_file, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                            extrasaction='ignore')
    writer.writeheader()

    for task in tasks:
        task_id        = task.get('task_id', '').strip()
        pdb_id         = task.get('pdb_id', '').strip()
        chain_manifest = task.get('chain', '').strip().upper()
        target_uniprot = task.get('target_uniprot', '').strip()

        # ── Clear old grid fields (zero-init every run) ──
        for fld in ['grid_center_x', 'grid_center_y', 'grid_center_z',
                    'grid_size_x', 'grid_size_y', 'grid_size_z']:
            task[fld] = ''
        task['grid_method']            = ''
        task['mapped_residue_count']   = '0'
        task['observed_residue_count'] = '0'
        task['requested_residue_count'] = '0'
        task['grid_qc_status']         = 'unchecked'
        task['grid_confidence']        = ''
        task['grid_fallback_reason']   = ''
        task['mapping_source']         = ''
        task['actual_chain_used']      = ''

        # ── No PDB ID ──
        if not pdb_id:
            task['grid_qc_status']       = 'no_pdb'
            task['grid_confidence']      = 'FAIL'
            task['grid_fallback_reason'] = 'No PDB ID in manifest'
            qc_stats['no_pdb'] += 1
            writer.writerow(task)
            continue

        # ── Residue descriptor ──
        residue_desc = (
            task.get('binding_residues', '').strip() or
            task.get('site_residues', '').strip() or
            task.get('evidence_residues', '').strip()
        )
        if not residue_desc:
            task['grid_qc_status']       = 'NO_SITE_DATA'
            task['grid_confidence']      = 'NO_SITE_DATA'
            task['grid_fallback_reason'] = 'No binding residue data in manifest'
            qc_stats['NO_SITE_DATA'] += 1
            writer.writerow(task)
            continue

        # Determine grid method priority
        if task.get('binding_residues', '').strip():
            task['grid_method'] = 'binding_residues'
        elif task.get('site_residues', '').strip():
            task['grid_method'] = 'site_residues'
        else:
            task['grid_method'] = 'evidence_residues'

        # ── Parse residues ──
        requested = parse_residue_parts(residue_desc)
        total_req = len(requested)
        task['requested_residue_count'] = str(total_req)

        if not requested:
            task['grid_qc_status']       = 'no_parseable_residues'
            task['grid_confidence']      = 'FAIL'
            task['grid_fallback_reason'] = 'Cannot parse residue descriptor'
            qc_stats['no_parseable'] += 1
            writer.writerow(task)
            continue

        # ── Locate PDB file ──
        pdb_path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
        if not os.path.exists(pdb_path):
            task['grid_qc_status']       = 'pdb_missing'
            task['grid_confidence']      = 'FAIL'
            task['grid_fallback_reason'] = f'PDB file missing: {pdb_path}'
            qc_stats['pdb_missing'] += 1
            writer.writerow(task)
            continue

        # ── Build atom index (once per task — cached per PDB?) ──
        # For simplicity we rebuild each time; DBREF lookups are already cached.
        atom_index = build_atom_index(pdb_path)

        # ── Which chains exist in the PDB? ──
        all_pdb_chains = sorted(set(ch for (ch, _, _, _) in atom_index))
        chain_exists   = (not chain_manifest) or (chain_manifest in all_pdb_chains)

        # ── Get DBREF for this PDB ──
        dbref          = dbref_cache.get(pdb_id, {})
        chain_dbref    = dbref.get(chain_manifest) if chain_exists else None
        dbref_accs     = chain_dbref.get('accessions', set()) if chain_dbref else set()

        # ── Try mappings in priority order ──
        coords         = []
        mapped         = 0     # residues that were successfully mapped to PDB numbering
        observed       = 0     # mapped residues with actual atom coordinates
        obs_frac       = 0.0
        mapping_source = 'none'
        actual_chain   = chain_manifest  # which chain was actually used for matching
        best_attempt   = ''
        fallback_detail = ''

        # 1) Direct match
        c_direct, m_direct, o_direct = match_residues(
            atom_index, requested, chain_manifest)
        obs_frac_direct = o_direct / total_req if total_req > 0 else 0

        # Accession validation for direct matches:
        # If the manifest specifies a UniProt accession and the PDB chain has
        # DBREF accessions that don't include it, flag as accession_mismatch.
        accession_ok = True
        if target_uniprot and dbref_accs and obs_frac_direct > 0:
            if target_uniprot not in dbref_accs:
                accession_ok = False

        if obs_frac_direct >= 0.5 and accession_ok:
            coords         = c_direct
            mapped         = total_req   # direct: all residues have PDB numbering
            observed       = o_direct
            obs_frac       = obs_frac_direct
            mapping_source = 'direct'
            best_attempt   = 'direct'
        else:
            best_attempt = f'direct({o_direct}/{total_req})'
            if not accession_ok:
                best_attempt += f'_acc_mismatch(dbref={dbref_accs}, manifest={target_uniprot})'

            # 2) SIFTS per-residue mapping
            sifts_mapped_list = []
            sifts_remap_count = 0
            for aa3, resnum, icode in requested:
                sifts_key = (pdb_id.upper(), chain_manifest,
                             target_uniprot, resnum)
                if sifts_key in sifts:
                    pdb_num, sifts_icode, sifts_resname, _ = sifts[sifts_key]
                    sifts_mapped_list.append((aa3, pdb_num, sifts_icode))
                    sifts_remap_count += 1
                elif chain_manifest:
                    # Try chain-insensitive lookup (same PDB, same UniProt, same UniProt num)
                    alt_found = False
                    for (spdb, sch, suni, snum), val in sifts.items():
                        if (spdb == pdb_id.upper() and suni == target_uniprot
                                and snum == resnum):
                            sifts_mapped_list.append((aa3, val[0], val[1]))
                            sifts_remap_count += 1
                            alt_found = True
                            break
                    if not alt_found:
                        # Keep original numbering for this residue
                        sifts_mapped_list.append((aa3, resnum, icode))

            if sifts_remap_count > 0:
                c_sifts, m_sifts, o_sifts = match_residues(
                    atom_index, sifts_mapped_list, chain_manifest)
                obs_frac_sifts = o_sifts / total_req if total_req > 0 else 0
                if obs_frac_sifts > obs_frac:  # better than previous best
                    coords         = c_sifts
                    mapped         = sifts_remap_count
                    observed       = o_sifts
                    obs_frac       = obs_frac_sifts
                    mapping_source = 'sifts'
                    best_attempt   = f'sifts({o_sifts}/{total_req})'

            # 3) DBREF segment-based (equal-length, no-insertion only)
            if obs_frac < 0.5 and chain_dbref:
                dbref_residues = []
                dbref_remap_count = 0
                for aa3, resnum, icode in requested:
                    mapped_num = None
                    for seg in chain_dbref.get('segments', []):
                        pdb_len = seg['pdb_end'] - seg['pdb_begin']
                        db_len  = seg['db_end'] - seg['db_begin']
                        no_ins  = not (seg.get('pdb_begin_ins', '') or
                                       seg.get('pdb_end_ins', '') or
                                       seg.get('db_begin_ins', '') or
                                       seg.get('db_end_ins', ''))
                        # Only equal-length, no-insertion segments
                        if no_ins and pdb_len == db_len and seg['db_begin'] <= resnum <= seg['db_end']:
                            offset     = resnum - seg['db_begin']
                            mapped_num = seg['pdb_begin'] + offset
                            dbref_remap_count += 1
                            break
                    dbref_residues.append(
                        (aa3, mapped_num if mapped_num is not None else resnum, icode))

                if dbref_remap_count > 0:
                    c_db, m_db, o_db = match_residues(
                        atom_index, dbref_residues, chain_manifest)
                    obs_frac_db = o_db / total_req if total_req > 0 else 0
                    if obs_frac_db > obs_frac:
                        coords         = c_db
                        mapped         = dbref_remap_count
                        observed       = o_db
                        obs_frac       = obs_frac_db
                        mapping_source = 'dbref_segment'
                        best_attempt   = f'dbref({o_db}/{total_req})'

            # 4) Same-UniProt homologous chain (REVIEW only, no grid written)
            #     Runs even when requested chain doesn't exist — rescue chain_mismatch cases
            if obs_frac < 0.5 and chain_manifest:
                candidate_chains = []
                for ch, cd in dbref.items():
                    if ch == chain_manifest:
                        continue
                    if target_uniprot and target_uniprot in cd.get('accessions', set()):
                        if ch in all_pdb_chains:
                            candidate_chains.append(ch)

                if candidate_chains:
                    for cand_ch in candidate_chains:
                        c_cross, m_cross, o_cross = match_residues(
                            atom_index, requested, cand_ch)
                        obs_frac_cross = o_cross / total_req if total_req > 0 else 0
                        if obs_frac_cross > obs_frac:
                            coords         = c_cross
                            mapped         = total_req
                            observed       = o_cross
                            obs_frac       = obs_frac_cross
                            mapping_source = 'homologous_chain'
                            actual_chain   = cand_ch
                            best_attempt   = f'homologous({cand_ch}:{o_cross}/{total_req})'
                            break  # take first candidate that improves

        # ═══════════════════════════════════════════════════════
        #  Write results
        # ═══════════════════════════════════════════════════════

        task['mapped_residue_count']   = str(mapped)
        task['observed_residue_count'] = str(observed)
        task['mapping_source']         = mapping_source
        task['actual_chain_used']      = actual_chain

        # ── Failure cases (observed == 0) ──
        if observed == 0:
            if not chain_exists:
                task['grid_qc_status']       = 'chain_mismatch'
                task['grid_confidence']      = 'FAIL'
                task['grid_fallback_reason'] = (
                    f'Manifest chain "{chain_manifest}" not in PDB. '
                    f'PDB chains: {all_pdb_chains}. '
                    f'Attempted: {best_attempt}')
                qc_stats['chain_mismatch'] += 1
            else:
                task['grid_qc_status']       = 'no_residue_observed'
                task['grid_confidence']      = 'FAIL'
                task['grid_fallback_reason'] = (
                    f'0/{total_req} residues observed. '
                    f'Attempted: {best_attempt}. '
                    f'Chain: {chain_manifest}. '
                    f'DBREF accessions: {dbref_accs}.')
                qc_stats['no_residue_observed'] += 1
            writer.writerow(task)
            continue

        # ── Compute grid box ──
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        sx = max(c[0] for c in coords) - min(c[0] for c in coords) + 2 * MARGIN
        sy = max(c[1] for c in coords) - min(c[1] for c in coords) + 2 * MARGIN
        sz = max(c[2] for c in coords) - min(c[2] for c in coords) + 2 * MARGIN
        sx = max(sx, MIN_SIZE)
        sy = max(sy, MIN_SIZE)
        sz = max(sz, MIN_SIZE)

        # ── Confidence & grid writing ──
        if mapping_source == 'homologous_chain':
            # Cross-chain: REVIEW only, NO grid written
            conf   = 'REVIEW'
            status = 'homologous_chain_candidate'
            task['grid_fallback_reason'] = (
                f'Matched on homologous chain {actual_chain} '
                f'(requested chain: {chain_manifest}). '
                f'Both map to UniProt {target_uniprot}. '
                f'{observed}/{total_req} residues have coordinates. '
                f'REVIEW required before use.')
            # Grid fields remain empty

        elif mapping_source in ('direct', 'sifts', 'dbref_segment'):
            if obs_frac >= 0.8:
                conf   = 'HIGH'
                status = f'mapped_{mapping_source}'
            elif obs_frac >= 0.5:
                conf   = 'MEDIUM'
                status = f'mapped_{mapping_source}_partial'
            else:
                conf   = 'REVIEW'
                status = 'low_coverage'

            # Write grid coordinates for HIGH/MEDIUM only
            if conf in ('HIGH', 'MEDIUM'):
                task['grid_center_x'] = f'{cx:.3f}'
                task['grid_center_y'] = f'{cy:.3f}'
                task['grid_center_z'] = f'{cz:.3f}'
                task['grid_size_x']   = f'{sx:.1f}'
                task['grid_size_y']   = f'{sy:.1f}'
                task['grid_size_z']   = f'{sz:.1f}'
                written_ok += 1
            else:
                task['grid_fallback_reason'] = (
                    f'Only {observed}/{total_req} residues observed '
                    f'(source: {mapping_source}). '
                    f'Below 50% threshold for grid writing.')
        else:
            conf   = 'REVIEW'
            status = mapping_source if mapping_source != 'none' else 'no_mapping'
            task['grid_fallback_reason'] = (
                f'No usable mapping found. '
                f'{observed}/{total_req} residues observed. '
                f'Attempts: {best_attempt}')

        task['grid_confidence'] = conf
        task['grid_qc_status']  = status
        qc_stats[f'confidence_{conf}'] += 1
        writer.writerow(task)

# ── Summary ──
print(f"\n{'='*60}")
print(f"Grid boxes with HIGH/MEDIUM: {written_ok}/{n_tasks}")
print(f"\nQC breakdown:")
for k, v in sorted(qc_stats.items()):
    print(f"  {k:30s}: {v}")

# ── Audit report ──
print(f"\nWriting audit report ...")
with open(audit_file, 'w', encoding='utf-8', newline='') as af:
    af_fields = ['task_id', 'pdb_id', 'chain', 'target_uniprot',
                 'grid_qc_status', 'grid_confidence', 'mapping_source',
                 'actual_chain_used',
                 'requested_residue_count', 'mapped_residue_count',
                 'observed_residue_count',
                 'grid_fallback_reason']
    writer = csv.DictWriter(af, fieldnames=af_fields, delimiter='\t',
                            extrasaction='ignore')
    writer.writeheader()
    for task in tasks:
        writer.writerow(task)

print(f"Output: {grid_file}")
print(f"Audit:  {audit_file}")
