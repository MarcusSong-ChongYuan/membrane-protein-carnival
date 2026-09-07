from pathlib import Path

p = Path(r"C:\Users\Administrator\MemPro_redraw_v2\build_clean_m1_m8.py")
s = p.read_text(encoding="utf-8")
s = s.replace(
    'anatomy(ax,organ,"detected_proteins","A  Broad tissue RNA coverage spans major organs","YlGnBu","proteins",8)',
    'anatomy(ax,organ,"detected_proteins","YlGnBu","A  Broad tissue RNA coverage spans major organs","proteins",8)',
)
s = s.replace(
    'anatomy(ax,organ,"protein_disease_pairs","A  Disease burden is distributed across organ systems","YlOrRd","pairs",8)',
    'anatomy(ax,organ,"protein_disease_pairs","YlOrRd","A  Disease burden is distributed across organ systems","pairs",8)',
)
exec(compile(s, str(p), "exec"), {"__name__": "__main__"})
