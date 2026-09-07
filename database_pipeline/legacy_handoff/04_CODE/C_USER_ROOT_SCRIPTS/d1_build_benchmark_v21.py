#!/usr/bin/env python3
from pathlib import Path
src=(Path(__file__).resolve().parent/'d1_build_benchmark_v2.py').read_text(encoding='utf-8')
src=src.replace("lig=ROOT/'ligands_pdbqt_chain_v2'/f\"{r['compound_internal_id']}.pdbqt\";q=ligand_qc(lig);reason=[]","repaired=ROOT/'ligands_pdbqt_repaired_v2'/f\"{r['compound_internal_id']}.pdbqt\";lig=repaired if repaired.exists() else ROOT/'ligands_pdbqt_chain_v2'/f\"{r['compound_internal_id']}.pdbqt\";q=ligand_qc(lig);reason=[]")
src=src.replace("if key in {'size_x','size_y','size_z','thread','num_modes','energy_range','out'}:continue\n   lines.append(line)","if key in {'size_x','size_y','size_z','thread','num_modes','energy_range','out'}:continue\n   if key=='ligand':line=f'ligand = {lig}'\n   lines.append(line)")
exec(compile(src,'d1_build_benchmark_v21.generated.py','exec'))
