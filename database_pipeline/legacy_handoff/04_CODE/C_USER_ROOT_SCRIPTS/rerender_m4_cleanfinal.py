from pathlib import Path

source = Path(r"C:\Users\Administrator\MemPro_narrative_v3\render_clean_final_problem_panels.py")
text = source.read_text(encoding="utf-8")
text = text.replace(
    'panel(ax, "B", "Observed disease–expression overlap exceeds matched random expectation")',
    'panel(ax, "B", "Disease–expression overlap exceeds random expectation")',
)
namespace = {"__name__": "mempro_clean_renderer"}
exec(compile(text, str(source), "exec"), namespace)
namespace["draw_m4"]()
