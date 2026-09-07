from pathlib import Path

p=Path(r"C:\Users\Administrator\build_v7_isoform_processed_exact.py")
s=p.read_text(encoding="utf-8")
old="""for acc,p in protein.items():
 rec=records.get(acc,{});seq=rec.get('sequence',{}).get('value','');segments=[]
 for ft in rec.get('features',[]):
  typ=ft.get('type','');loc=ft.get('location',{});start,sm=pos(loc,'start');end,em=pos(loc,'end')
  if typ in {'Transmembrane','Intramembrane'} and start and end and sm=='EXACT' and em=='EXACT':segments.append((typ,start,end,seq[start-1:end]))
  if typ not in {'Chain','Peptide','Propeptide','Signal','Transit peptide'}:continue
"""
new="""for acc,p in protein.items():
 rec=records.get(acc,{});seq=rec.get('sequence',{}).get('value','');segments=[]
 for ft in rec.get('features',[]):
  typ=ft.get('type','');loc=ft.get('location',{});start,sm=pos(loc,'start');end,em=pos(loc,'end')
  if typ in {'Transmembrane','Intramembrane'} and start and end and sm=='EXACT' and em=='EXACT':segments.append((typ,start,end,seq[start-1:end]))
 for ft in rec.get('features',[]):
  typ=ft.get('type','');loc=ft.get('location',{});start,sm=pos(loc,'start');end,em=pos(loc,'end')
  if typ not in {'Chain','Peptide','Propeptide','Signal','Transit peptide'}:continue
"""
if old not in s:
    raise RuntimeError("Expected parser block not found")
exec(compile(s.replace(old,new),"build_v7_isoform_processed_exact_v2.py","exec"))
