from pathlib import Path

source = Path(r"C:\Users\Administrator\MemPro_narrative_v3\qa_full_suite.py")
text = source.read_text(encoding="utf-8")
text = text.replace(
    '"M3": "M3_expression_space_and_RNA_IHC",',
    '"M3": "M3_expression_space_and_RNA_IHC_cleanfinal",',
)
namespace = {"__name__": "mempro_qa_v2"}
exec(compile(text, str(source), "exec"), namespace)
namespace["main"]()
