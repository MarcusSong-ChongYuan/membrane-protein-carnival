from pathlib import Path

source_path = Path(__file__).with_name("build_revised_figures.py")
source = source_path.read_text(encoding="utf-8")
source = source.replace(
    "import matplotlib as mpl\n",
    "import matplotlib as mpl\nfrom matplotlib import font_manager as fm\nfm.fontManager.addfont(r'C:\\Windows\\Fonts\\msyh.ttc')\n",
)
source = source.replace(
    '"font.family": ["Microsoft YaHei", "Arial", "DejaVu Sans"],',
    '"font.family": fm.FontProperties(fname=r"C:\\Windows\\Fonts\\msyh.ttc").get_name(),',
)
source = source.replace(
    'sns.set_theme(style="whitegrid", rc=mpl.rcParams)',
    'sns.set_theme(style="whitegrid", rc=mpl.rcParams)\nmpl.rcParams["font.family"] = fm.FontProperties(fname=r"C:\\Windows\\Fonts\\msyh.ttc").get_name()',
)
source = source.replace(
    '"无正–负冲突": ~dock["positive_negative_conflict_flag_v62"].isin(["1", "true", "True"]),',
    '"无正–负冲突": int((~dock["positive_negative_conflict_flag_v62"].isin(["1", "true", "True"])).sum()),',
)
source = source.replace(
    '    readiness.iloc[-1] = int(readiness.iloc[-1].sum()) if hasattr(readiness.iloc[-1], "sum") else int(readiness.iloc[-1])',
    '',
)
source = source.replace('joined["hpa_mapping_status_v62_expr"]', 'joined["hpa_mapping_status_v62"]')
source = source.replace('joined["hpa_tissue_detected_count_v62_expr"]', 'joined["hpa_tissue_detected_count_v62"]')
source = source.replace('labels=[role_cn[r] for r in wanted]', 'tick_labels=[role_cn[r] for r in wanted]')
source = source.replace('labels=order, patch_artist=True', 'tick_labels=order, patch_artist=True')
exec(compile(source, str(source_path), "exec"), {"__name__": "__main__", "__file__": str(source_path)})
