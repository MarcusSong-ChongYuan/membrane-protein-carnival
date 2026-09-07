from pathlib import Path

p=Path(r"C:\Users\Administrator\MemPro_narrative_v3\build_final_m1_m2.py")
s=p.read_text(encoding="utf-8")
s=s.replace('share=tab.div(tab.sum(1),axis=0)', 'share=tab.div(tab.sum(1),axis=0);share=share.rename(columns={"structural_family":"SF","molecular_function":"MF","biological_process":"BP","membrane_role":"MR","specialist_classification":"SC"})')
s=s.replace('"M2_multi_axis_annotation_space"', '"M2_multi_axis_annotation_space_final2"')
s=s.replace('if __name__=="__main__":m1();m2();print("final M1/M2 complete")','if __name__=="__main__":m2();print("clean M2 complete")')
exec(compile(s,str(p),"exec"),{"__name__":"__main__"})
