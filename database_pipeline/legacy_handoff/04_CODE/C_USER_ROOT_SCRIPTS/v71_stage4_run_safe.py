from pathlib import Path

p=Path(r'C:\Users\Administrator\v71_stage4_expression_negative.py')+s=p.read_text(encoding='utf-8')
s=s.replace('import csv,gzip,json,hashlib','import csv,gzip,json,hashlib,itertools')
old=""" eg=expression_rows(estats);firste=next(eg);ef=list(firste.keys());ec=write(OUT/'expression_measurement_v71.tsv.gz',ef,(x for x in [firste]));
 # append remaining rows without recompressing through a second file
 # Re-run once to preserve simple deterministic streaming and overwrite complete table.
 estats={'mapping':Counter(),'measurement':Counter(),'detection':Counter(),'assay':Counter(),'denom':Counter()};eg=expression_rows(estats);firste=next(eg);ef=list(firste.keys());ec=write(OUT/'expression_measurement_v71.tsv.gz',ef,(x for x in [firste]+list(eg)))
 lstats=Counter();lg=localization_rows(lstats);firstl=next(lg);lf=list(firstl.keys());lc=write(OUT/'subcellular_localization_v71.tsv.gz',lf,(x for x in [firstl]+list(lg)))"""
new=""" eg=expression_rows(estats);firste=next(eg);ef=list(firste.keys());ec=write(OUT/'expression_measurement_v71.tsv.gz',ef,itertools.chain([firste],eg))
 lstats=Counter();lg=localization_rows(lstats);firstl=next(lg);lf=list(firstl.keys());lc=write(OUT/'subcellular_localization_v71.tsv.gz',lf,itertools.chain([firstl],lg))"""
if old not in s: raise SystemExit('safe replacement target not found')
s=s.replace(old,new)
exec(compile(s,str(p),'exec'),{'__name__':'__main__','__file__':str(p)})
