from pathlib import Path

p=Path(r"C:\Users\Administrator\MemPro_narrative_v3\build_stage3_m5_m6.py")
s=p.read_text(encoding="utf-8")
s=s.replace('["kNN retention","cluster purity","assigned fraction"]*2', '["kNN","purity","assigned"]*2')
s=s.replace('"M5_chemical_space_embedding_comparison"', '"M5_chemical_space_embedding_comparison_final"')
old='sc=sc[~sc.membrane_side.fillna("unknown").str.lower().eq("unknown")];sc.to_csv'
new='sc=sc[~sc.membrane_side.fillna("unknown").str.lower().eq("unknown")];sc["membrane_side"]=sc["membrane_side"].replace({"extramembrane_side_unresolved":"Extra, unresolved","membrane_interface_or_mixed":"Interface/mixed","intramembrane":"Intramembrane","cytoplasmic_side":"Cytoplasmic","non_cytoplasmic_side":"Non-cytoplasmic","both_sides_or_mixed":"Both/mixed"});sc["xlogp"]=sc["xlogp"].clip(sc["xlogp"].quantile(.01),sc["xlogp"].quantile(.99));sc["tpsa"]=sc["tpsa"].clip(sc["tpsa"].quantile(.01),sc["tpsa"].quantile(.99));sc.to_csv'
s=s.replace(old,new)
s=s.replace('"M6_polypharmacology_and_membrane_side"', '"M6_polypharmacology_and_membrane_side_final"')
exec(compile(s,str(p),"exec"),{"__name__":"__main__"})
