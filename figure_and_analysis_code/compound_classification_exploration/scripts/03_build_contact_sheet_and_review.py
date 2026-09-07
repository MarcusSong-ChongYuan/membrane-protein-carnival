from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

ROOT = Path(r"D:\finale\compound_classification_exploration")
OUT = ROOT / "08_contact_sheet"
REVIEW = ROOT / "09_review"
OUT.mkdir(exist_ok=True)
REVIEW.mkdir(exist_ok=True)

ITEMS = [
    ("V1", ROOT/"02_classification_tables/V1B_chemical_taxonomy_hierarchy.png",
     "How are interaction-linked compounds distributed across broad local chemical regimes?",
     "240,054 structure-valid interaction-linked compounds", "SUPPLEMENTARY",
     "Useful overview, but the regime is a deterministic local taxonomy rather than ClassyFire/ChEBI ontology."),
    ("V2", ROOT/"03_scaffold/V2A_exact_scaffold_rank_abundance.png",
     "Is scaffold use concentrated or long-tailed?",
     "240,054 structure-valid interaction-linked compounds", "MAIN CANDIDATE",
     "Directly quantifies scaffold diversity and concentration without treating clusters as classes."),
    ("V3", ROOT/"03_scaffold/V3_top_scaffold_gallery.png",
     "Which recurrent Bemis-Murcko scaffolds dominate the linked chemical collection?",
     "Top recurrent scaffolds among 240,054 compounds", "MAIN CANDIDATE",
     "Chemically intuitive; use with rank-abundance rather than as a stand-alone frequency gallery."),
    ("V4", ROOT/"04_heatmaps/V4B_chemical_class_role_profile_heatmap.png",
     "Do broad chemical regimes show different membrane-role profiles?",
     "529,168 formal pairs; row-normalized descriptive profiles", "SUPPLEMENTARY",
     "Shows preference patterns, but broad local classes may hide scaffold-level specificity."),
    ("V5", ROOT/"04_heatmaps/V5C_scaffold_role_heatmap_biclustered.png",
     "Which abundant scaffolds share membrane-target role profiles?",
     "Top 50 generic scaffolds across 529,168 formal pairs", "MAIN CANDIDATE",
     "Strong cross-module view; retain support counts and avoid interpreting clustering as causal taxonomy."),
    ("V6", ROOT/"04_heatmaps/V6_physchem_class_heatmap.png",
     "How do physicochemical properties vary across local chemical regimes?",
     "240,054 structure-valid interaction-linked compounds", "SUPPLEMENTARY",
     "Compact descriptor comparison; medians are descriptive and do not establish drug-likeness."),
    ("V7", ROOT/"05_chemical_space/V7A_PCA_chemical_superclass.png",
     "How much variance in standardized descriptors is captured by a linear projection?",
     "Deterministic stratified sample n=50,000", "SUPPLEMENTARY",
     "Axes are interpretable through loadings, but overlap indicates coarse regimes are not cleanly separable."),
    ("V8", ROOT/"05_chemical_space/V8B_UMAP_scaffold_family.png",
     "Do Morgan-fingerprint neighborhoods form local scaffold islands?",
     "Deterministic stratified sample n=20,000", "EXPLORATORY ONLY",
     "Useful for local-neighborhood inspection; UMAP axes and global distances have no direct chemical meaning."),
    ("V9", ROOT/"06_target_breadth/V9A_target_breadth_distribution.png",
     "How selective or polypharmacological are linked compounds by observed target breadth?",
     "240,054 interaction-linked compounds; formal-pair target counts", "MAIN CANDIDATE",
     "Directly supports target-breadth narrative; observed breadth remains coverage-dependent."),
    ("V10", ROOT/"07_optional_3d/V10_3D_physchem_space.png",
     "Can users interactively inspect MW-LogP-TPSA space and target breadth?",
     "Deterministic stratified sample n=20,000", "EXPLORATORY ONLY",
     "Interactive exploration aid; perspective and overplotting make it unsuitable as primary evidence."),
]

plt.rcParams.update({"font.family":"Arial", "font.size":8, "pdf.fonttype":42, "svg.fonttype":"none"})

def make_page(items, page_no):
    fig = plt.figure(figsize=(11.69, 8.27), facecolor="white")
    gs = fig.add_gridspec(2, 3, left=.035, right=.98, top=.92, bottom=.055, wspace=.18, hspace=.28)
    fig.suptitle(f"MemPro compound-classification exploration | candidate views {page_no}", x=.035, ha="left", fontsize=14, fontweight="bold")
    for i, (code, path, question, denom, rec, note) in enumerate(items):
        ax = fig.add_subplot(gs[i//3, i%3]); ax.axis("off")
        im = Image.open(path).convert("RGB")
        ax.imshow(im)
        ax.set_title(f"{code}  {rec}", loc="left", fontsize=9, fontweight="bold", color="#31445a", pad=4)
        ax.text(0, -0.06, question, transform=ax.transAxes, fontsize=7.2, va="top", wrap=True)
        ax.text(0, -0.16, f"Denominator: {denom}", transform=ax.transAxes, fontsize=6.7, va="top", color="#555555", wrap=True)
        ax.text(0, -0.26, note, transform=ax.transAxes, fontsize=6.7, va="top", color="#555555", wrap=True)
    return fig

pdf_path = OUT / "COMPOUND_VISUAL_EXPLORATION_CONTACT_SHEET.pdf"
with PdfPages(pdf_path) as pdf:
    pdf.savefig(make_page(ITEMS[:6], "1/2"), bbox_inches="tight")
    plt.close()
    pdf.savefig(make_page(ITEMS[6:], "2/2"), bbox_inches="tight")
    plt.close()

# One-page raster overview for fast review.
fig = plt.figure(figsize=(13, 8), facecolor="white")
gs = fig.add_gridspec(2,5,left=.02,right=.99,top=.93,bottom=.04,wspace=.08,hspace=.16)
fig.suptitle("MemPro compound-classification exploration | V1-V10", x=.02, ha="left", fontsize=16, fontweight="bold")
for i,(code,path,question,denom,rec,note) in enumerate(ITEMS):
    ax=fig.add_subplot(gs[i//5,i%5]); ax.axis("off"); ax.imshow(Image.open(path).convert("RGB")); ax.set_title(f"{code} | {rec}",loc="left",fontsize=8,fontweight="bold")
fig.savefig(OUT/"COMPOUND_VISUAL_EXPLORATION_CONTACT_SHEET.png",dpi=220,bbox_inches="tight")
plt.close(fig)

review_lines = [
    "# MemPro compound-classification visualization review",
    "",
    "This is an exploratory selection report, not the final manuscript Figure 3. All structural analyses use 240,054 structure-valid interaction-linked compounds; the quarantined identity-conflict record remains in the registry but is excluded from structure-derived analyses.",
    "",
    "## Denominator guardrails",
    "",
    "- Compound registry: **646,670** identities.",
    "- Interaction-linked compounds: **240,055** identities.",
    "- Structure-valid interaction-linked compounds: **240,054**.",
    "- Formal protein-compound pairs: **529,168**.",
    "- These denominators are not interchangeable.",
    "",
]
for code,path,question,denom,rec,note in ITEMS:
    inferential = "No; descriptive/exploratory" if code not in {"V4","V5"} else "Descriptive association view; no causal inference"
    review_lines += [
        f"## {code} | {rec}", "",
        f"- **Question:** {question}",
        f"- **Denominator:** {denom}",
        f"- **Analysis type:** {inferential}.",
        f"- **Strength:** {note.split(';')[0]}.",
        f"- **Limitation:** {note}",
        f"- **Recommended placement:** {rec}.",
        "",
    ]
review_lines += [
    "## Recommended Figure 3 core", "",
    "The strongest nonredundant manuscript candidates are V2 (scaffold concentration), V3 (representative recurrent scaffolds), V5 (scaffold-by-membrane-role profiles), and V9 (observed target breadth). V1, V4, V6, and V7 are useful supporting context. V8 and V10 should remain exploratory because their geometry is projection- or perspective-dependent.",
    "",
    "## Classification boundary", "",
    "`chemical_regime_local` is a transparent deterministic classification derived locally from structure/descriptor rules. It must not be described as ClassyFire, ChEBI, or an experimentally validated chemical ontology. Bemis-Murcko scaffold identity and Morgan/Butina neighborhoods answer different questions and must not be merged into one 'class' field.",
]
(REVIEW/"COMPOUND_VISUALIZATION_REVIEW.md").write_text("\n".join(review_lines)+"\n", encoding="utf-8")

summary = {
    "status":"PASS", "candidate_views":10,
    "recommendations":{"MAIN CANDIDATE":4,"SUPPLEMENTARY":4,"EXPLORATORY ONLY":2},
    "contact_sheet_pdf":str(pdf_path),
}
(REVIEW/"COMPOUND_VISUALIZATION_REVIEW.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
