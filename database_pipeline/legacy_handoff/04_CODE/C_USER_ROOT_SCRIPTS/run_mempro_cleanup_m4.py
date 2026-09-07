from pathlib import Path

p=Path(r"C:\Users\Administrator\MemPro_narrative_v3\build_stage2_m3_m4.py")
s=p.read_text(encoding="utf-8")
s=s.replace('share=tab.div(tab.sum(1),axis=0);share.to_csv', 'share=tab.div(tab.sum(1),axis=0);share=share.rename(columns={"receptor":"Receptor","membrane_associated_enzyme":"Enzyme","ion_channel":"Channel","transporter":"Transporter","family_defined_membrane_role_unresolved":"Role unresolved","membrane_scaffold_or_linker":"Scaffold","signaling_regulator":"Signaling","adhesion_recognition":"Adhesion","immune_or_cell_recognition":"Immune"});share.index=[textwrap.shorten(str(x),22,placeholder="...") for x in share.index];share.to_csv')
s=s.replace('ev=ev.loc[ev.sum(1).nlargest(10).index];sns.heatmap', 'ev=ev.loc[ev.sum(1).nlargest(10).index];ev.index=[textwrap.shorten(str(x),22,placeholder="...") for x in ev.index];sns.heatmap')
s=s.replace('"M4_disease_expression_concordance"', '"M4_disease_expression_concordance_final"')
s=s.replace('if __name__=="__main__":draw_m3();draw_m4();', 'if __name__=="__main__":draw_m4();')
exec(compile(s,str(p),"exec"),{"__name__":"__main__"})
