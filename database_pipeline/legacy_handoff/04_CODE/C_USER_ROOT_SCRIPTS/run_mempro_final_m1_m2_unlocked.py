from pathlib import Path

p = Path(r"C:\Users\Administrator\MemPro_narrative_v3\build_final_m1_m2.py")
s = p.read_text(encoding="utf-8")
s = s.replace('"M2_multi_axis_annotation_space"', '"M2_multi_axis_annotation_space_final"')
exec(compile(s, str(p), "exec"), {"__name__": "__main__"})
