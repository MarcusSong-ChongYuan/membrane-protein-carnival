from pathlib import Path

p = Path(r"C:\Users\Administrator\MemPro_narrative_v3\build_stage1_m1_m2.py")
s = p.read_text(encoding="utf-8")
s = s.replace("[:10]]", "[:6]]", 1)
s = s.replace('figsize=(13.6,9.2)', 'figsize=(13.6,10.4)')
s = s.replace('hspace=.43', 'hspace=.62', 1)
s = s.replace('ia=ax.inset_axes([0,-.42,1,.34])', 'ia=ax.inset_axes([0,-.30,1,.22])')
s = s.replace('fontsize=6.8)', 'fontsize=7.2)', 1)
s = s.replace('ax.set_xlabel("");ax.set_ylabel("community")', 'ax.set_xlabel("");ax.set_ylabel("community");ax.tick_params(axis="x",rotation=45)')
old = 'fig.text(.56,.105,f"Embedding input: {shape[0]:,} proteins x {shape[1]:,} recurrent annotation labels; seed=42.",fontsize=7,color=C["muted"])'
s = s.replace(old, '')
exec(compile(s, str(p), "exec"), {"__name__": "__main__"})
