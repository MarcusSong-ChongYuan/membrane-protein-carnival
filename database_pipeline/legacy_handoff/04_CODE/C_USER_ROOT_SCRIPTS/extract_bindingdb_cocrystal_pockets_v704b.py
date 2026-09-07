import csv,gzip
from pathlib import Path
evidence_ids=set();rem=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites\binding_site_coordinate_remaining_v703.tsv.gz")
with gzip.open(rem,'rt',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['uniprot_rescue_status_v703']=='REMAIN_NO_RESIDUE_ANNOTATION':evidence_ids.add(r['evidence_id'])
evidence_het={};evid=Path(r"D:\finale\13_MemPro_V7_final_data_release_20260813\01_data\protein_compound_evidence_v7.tsv.gz")
with gzip.open(evid,'rt',encoding='utf-8',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['evidence_id'] in evidence_ids:evidence_het[r['evidence_id']]=r.get('ligand_het_id','').strip().upper()
src=Path(r"C:\Users\Administrator\extract_bindingdb_cocrystal_pockets_v704.py").read_text(encoding='utf-8')
src=src.replace("r.get('ligand_het_id','').strip().upper()","evidence_het.get(r['evidence_id'],'')")
src=src.replace("v704_cocrystal","v704b_cocrystal")
src=src.replace("POCKET_V704_REPORT","POCKET_V704B_REPORT")
src=src.replace("pocket_v704.tsv.gz","pocket_v704b.tsv.gz")
src=src.replace("remaining_v704.tsv.gz","remaining_v704b.tsv.gz")
src=src.replace("pocket_rescue_status_v704","pocket_rescue_status_v704b")
src=src.replace("pocket_rescue_rule_v704","pocket_rescue_rule_v704b")
exec(compile(src,'extract_bindingdb_cocrystal_pockets_v704b.generated.py','exec'))
