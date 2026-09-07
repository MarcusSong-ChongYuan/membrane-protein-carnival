#!/usr/bin/env python3
"""Conservative target-chain receptor preparation for the MemPro D1/G1 pilot.

One receptor is keyed by PDB + actual author chain.  Only polymer ATOM records are
converted.  Non-water HETATM records in or near the grid are reported and block
automatic configuration; they are never silently deleted from a production task.
"""
import csv, hashlib, json, os, subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OBABEL = "/home/csong/miniconda3/envs/docking/bin/obabel"
SRC = ROOT / "docking_manifest_D1_G1_prefilter_pass.tsv"
PDB_CACHE = Path("/home/csong/docking/handoff/MemPro_Docking_Dscope_RCSBclean_14333_20260811/pdb_structures")
LIG_CACHE = Path("/home/csong/docking/handoff/MemPro_Docking_Dscope_RCSBclean_14333_20260811/ligands_pdbqt")
CHAIN_PDB = ROOT / "receptors_chain_pdb_v2"
RAW = ROOT / "receptors_pdbqt_raw_chain_v2"
CLEAN = ROOT / "receptors_pdbqt_clean_chain_v2"
LIG = ROOT / "ligands_pdbqt_chain_v2"
CONF = ROOT / "vina_configs_chain_v2"
META = ROOT / "metadata_chain_v2"
for d in (CHAIN_PDB, RAW, CLEAN, LIG, CONF, META): d.mkdir(exist_ok=True)

def fnum(line, a, b):
    try: return float(line[a:b])
    except ValueError: return None

def in_shell(x, y, z, row, margin=5.0):
    if None in (x, y, z): return True
    for axis, value in zip("xyz", (x, y, z)):
        c = float(row[f"grid_center_{axis}"])
        half = float(row[f"grid_size_{axis}"]) / 2.0 + margin
        if abs(value - c) > half: return False
    return True

def blank_type(line):
    return line.startswith(("ATOM  ", "HETATM")) and (len(line) < 79 or not line[77:79].strip())

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""): h.update(block)
    return h.hexdigest()

with SRC.open(encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    base_fields, rows = reader.fieldnames, list(reader)

# Pre-extract each PDB/author-chain once. PDB legacy files have one-character chain IDs.
keys = sorted({(r["pdb_id"].lower(), r["actual_chain_used"].strip()) for r in rows})
prep = {}
for pdb, chain in keys:
    key = f"{pdb}__{chain or 'blank'}"
    src, chain_pdb = PDB_CACHE / f"{pdb}.pdb", CHAIN_PDB / f"{key}.pdb"
    raw, clean = RAW / f"{key}.pdbqt", CLEAN / f"{key}.pdbqt"
    if not src.exists(): prep[(pdb, chain)] = {"status": "PDB_MISSING"}; continue
    atom_lines, het = [], []
    with src.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            rec = line[:6]
            if rec in ("ATOM  ", "HETATM") and line[21:22].strip() == chain:
                if rec == "ATOM  ": atom_lines.append(line.rstrip("\r\n"))
                elif line[17:20].strip() not in {"HOH", "WAT", "DOD"}:
                    het.append((line[17:20].strip(), line[12:16].strip(), fnum(line,30,38), fnum(line,38,46), fnum(line,46,54)))
    if not atom_lines: prep[(pdb, chain)] = {"status": "TARGET_CHAIN_NO_ATOM", "het": het}; continue
    chain_pdb.write_text("\n".join(atom_lines) + "\nTER\nEND\n", encoding="ascii", errors="ignore")
    tmp = RAW / f"{key}.tmp.pdbqt"
    cp = subprocess.run(["timeout", "300", OBABEL, str(chain_pdb), "-O", str(tmp), "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if cp.returncode != 0 or not tmp.exists():
        prep[(pdb, chain)] = {"status": "OBABEL_FAILED", "het": het, "stderr": cp.stderr[-300:]}; continue
    kept = []
    for line in tmp.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM", "TER")): kept.append(line)
    raw.write_text("\n".join(kept) + "\n", encoding="utf-8")
    tmp.unlink(missing_ok=True)
    blanks = []
    for line in kept:
        if blank_type(line): blanks.append((fnum(line,30,38), fnum(line,38,46), fnum(line,46,54), line[17:20].strip(), line[12:16].strip()))
    # A clean file may omit unsupported atoms only after every task using this chain
    # confirms those atoms are outside its grid shell. Creation is deferred below.
    prep[(pdb, chain)] = {"status": "CONVERTED_AWAITING_TASK_QC", "raw": raw, "clean": clean, "het": het, "blanks": blanks}

extra = ["chain_receptor_key","chain_prep_status","nonwater_hetatm_total_same_chain","nonwater_hetatm_in_grid_or_5A","blank_autodock_type_total","blank_autodock_type_in_grid_or_5A","ligand_pdbqt_status","config_status"]
out = []
key_all_safe = {k: True for k in prep}
task_qc = []
for r in rows:
    pdb, chain = r["pdb_id"].lower(), r["actual_chain_used"].strip()
    key = f"{pdb}__{chain or 'blank'}"; p = prep[(pdb, chain)]
    het_near = [a for a in p.get("het", []) if in_shell(a[2],a[3],a[4],r)]
    blank_near = [a for a in p.get("blanks", []) if in_shell(a[0],a[1],a[2],r)]
    status = p["status"]
    if status == "CONVERTED_AWAITING_TASK_QC":
        if blank_near: status = "PREP_FAILED_UNSUPPORTED_ATOM_NEAR_POCKET"
        elif het_near: status = "HOLD_NONWATER_HETATM_NEAR_POCKET"
        else: status = "PASS_CHAIN_CHEMISTRY_QC"
    if status != "PASS_CHAIN_CHEMISTRY_QC": key_all_safe[(pdb,chain)] = False
    ligsrc = LIG_CACHE / f"{r['compound_internal_id']}.pdbqt"
    ligstatus = "PASS" if ligsrc.exists() and ligsrc.stat().st_size else "LIGAND_PDBQT_MISSING"
    x = dict(r); x.update({"chain_receptor_key":key,"chain_prep_status":status,"nonwater_hetatm_total_same_chain":str(len(p.get('het',[]))),"nonwater_hetatm_in_grid_or_5A":str(len(het_near)),"blank_autodock_type_total":str(len(p.get('blanks',[]))),"blank_autodock_type_in_grid_or_5A":str(len(blank_near)),"ligand_pdbqt_status":ligstatus,"config_status":"HOLD_PENDING_FINAL_CHAIN_QC"})
    task_qc.append((x, ligsrc)); out.append(x)

# Materialize receptor/config only where every task using the chain receptor is safe.
for x, ligsrc in task_qc:
    pdb, chain = x["pdb_id"].lower(), x["actual_chain_used"].strip(); p = prep[(pdb,chain)]
    if x["chain_prep_status"] != "PASS_CHAIN_CHEMISTRY_QC" or x["ligand_pdbqt_status"] != "PASS" or not key_all_safe[(pdb,chain)]: continue
    clean = p["clean"]
    if not clean.exists():
        lines = [line for line in p["raw"].read_text(encoding="utf-8", errors="ignore").splitlines() if not blank_type(line)]
        clean.write_text("\n".join(lines)+"\n", encoding="utf-8")
    dst = LIG / ligsrc.name
    if not dst.exists(): os.link(ligsrc, dst)
    tid = x["task_id"]
    conf = (f"receptor = receptors_pdbqt_clean_chain_v2/{x['chain_receptor_key']}.pdbqt\n"
            f"ligand = ligands_pdbqt_chain_v2/{x['compound_internal_id']}.pdbqt\n"
            f"center_x = {x['grid_center_x']}\ncenter_y = {x['grid_center_y']}\ncenter_z = {x['grid_center_z']}\n"
            f"size_x = {x['grid_size_x']}\nsize_y = {x['grid_size_y']}\nsize_z = {x['grid_size_z']}\n"
            "thread = 8000\nnum_modes = 9\nenergy_range = 3\n"
            f"out = results_chain_v2/{tid}_out.pdbqt\n")
    (CONF / f"{tid}.conf").write_text(conf, encoding="ascii")
    x["config_status"] = "READY_FOR_PARAMETER_BENCHMARK"

with (ROOT / "docking_manifest_D1_G1_chain_receptor_qc_v2.tsv").open("w", encoding="utf-8", newline="") as f:
    w=csv.DictWriter(f, fieldnames=base_fields+extra, delimiter="\t"); w.writeheader(); w.writerows(out)
summary={"input_tasks":len(out),"unique_pdb_chain":len(keys),"ready_configs":sum(x["config_status"].startswith("READY") for x in out),"chain_prep_status":dict(Counter(x["chain_prep_status"] for x in out)),"ligand_status":dict(Counter(x["ligand_pdbqt_status"] for x in out)),"note":"No docking launched; READY means eligible for parameter benchmark only."}
(META / "D1_G1_CHAIN_RECEPTOR_QC_V2_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))


