from pathlib import Path

p=Path(r"C:\Users\Administrator\MemPro_narrative_v3\build_stage4_m7_m8.py")
s=p.read_text(encoding="utf-8")
s=s.replace('sl=sl.loc[sl.sum(1).nlargest(10).index];sl.to_csv', 'sl=sl.loc[sl.sum(1).nlargest(10).index];sl=sl.rename(columns={"database_record_lineage_not_independent_experiment":"DB record","literature_experiment_proxy_key":"Literature key","pubchem_assay_proxy_key":"PubChem key","structure_key":"Structure key"});sl.to_csv')
s=s.replace('mt=mt.loc[mt.sum(1).nlargest(10).index];mt.to_csv', 'mt=mt.loc[mt.sum(1).nlargest(10).index];mt.index=[textwrap.shorten(str(x),24,placeholder="...") for x in mt.index];mt.to_csv')
s=s.replace('left=.075,right=.985,top=.87,bottom=.09', 'left=.075,right=.985,top=.87,bottom=.14')
s=s.replace('"M7_evidence_profiles_and_lineage"', '"M7_evidence_profiles_and_lineage_final"')
s=s.replace('if __name__=="__main__":draw_m7();draw_m8();', 'if __name__=="__main__":draw_m7();')
exec(compile(s,str(p),"exec"),{"__name__":"__main__"})
