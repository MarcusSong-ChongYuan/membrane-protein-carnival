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
source = source.replace('ta = ta_counts.sort_values(ta_counts.columns[-1]).tail(10)', 'ta = ta_counts.sort_values(ta_counts.columns[-1]).tail(8)')
source = source.replace('an = an_counts.sort_values(an_counts.columns[-1]).tail(10)', 'an = an_counts.sort_values(an_counts.columns[-1]).tail(8)')
source = source.replace(
    'ax.set_yticks(np.arange(len(ta)), [wrap(x, 24) for x in ta.iloc[:, 0]])',
    '''ta_labels = [str(x).replace("genetic, familial or congenital disease", "遗传/先天").replace("cancer or benign tumor", "肿瘤").replace("nervous system disorder", "神经系统").replace("musculoskeletal or connective tissue disease", "肌肉骨骼/结缔组织").replace("nutritional or metabolic disease", "营养/代谢").replace("gastrointestinal disease", "胃肠道").replace("respiratory or thoracic disease", "呼吸/胸部").replace("hematologic disease", "血液") for x in ta.iloc[:, 0]]\n    ax.set_yticks(np.arange(len(ta)), ta_labels, fontsize=9.0)''',
)
source = source.replace(
    'ax.set_yticks(np.arange(len(an)), [wrap(x, 24) for x in an.iloc[:, 0]])',
    '''an_labels = [str(x).replace("nervous system", "神经系统").replace("musculoskeletal system", "肌肉骨骼").replace("digestive system", "消化系统").replace("alimentary part of gastrointestinal system", "消化道").replace("gastrointestinal system", "胃肠系统").replace("respiratory system", "呼吸系统").replace("cardiovascular system", "心血管").replace("genitourinary system", "泌尿生殖").replace("integumental system", "皮肤/体被") for x in an.iloc[:, 0]]\n    ax.set_yticks(np.arange(len(an)), an_labels, fontsize=9.0)''',
)
namespace = {"__name__": "mempro_figure_module", "__file__": str(source_path)}
exec(compile(source, str(source_path), "exec"), namespace)
namespace["make_disease"]()
