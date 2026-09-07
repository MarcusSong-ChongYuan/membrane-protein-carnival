import sys, os
# Remove MGLTools contamination
sys.path = [p for p in sys.path if 'MGLTools' not in p]
# Ensure venv site-packages is first
venv_sp = r'C:\Users\Administrator\plotenv\Lib\site-packages'
if venv_sp not in sys.path:
    sys.path.insert(0, venv_sp)

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import matplotlib.patches as mpatches

# ---- Style ----
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Microsoft YaHei', 'Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# ---- Load data ----
with open('C:/Users/Administrator/residue_recall_v2.jsonl') as f:
    data = [json.loads(l) for l in f if l.strip()]

direct = [r for r in data if r['mapping'] == 'direct']
dbref = [r for r in data if r['mapping'] == 'dbref_segment']

dr = [r['recall']*100 for r in direct]
dp = [r['precision']*100 for r in direct]
d2r = [r['recall']*100 for r in dbref]
d2p = [r['precision']*100 for r in dbref]

C1 = '#2563EB'  # G1 blue
C1L = '#93C5FD'
C2 = '#EA580C'  # G2 orange
C2L = '#FDBA74'
CG = '#16A34A'  # green for precision
CR = '#DC2626'  # red accent

# ---- Figure ----
fig = plt.figure(figsize=(16, 12))
fig.suptitle('Docking Residue Accuracy  —  Recall · Precision · Hit Rate',
             fontsize=18, fontweight='bold', y=0.98)

# ===== Panel 1: G1 Recall histogram =====
ax1 = fig.add_subplot(2, 2, 1)
bins = np.linspace(0, 100, 51)
ax1.hist(dr, bins=bins, color=C1, edgecolor='white', alpha=0.88, linewidth=0.3)
ax1.axvline(np.mean(dr), color='#991B1B', linestyle='--', linewidth=2,
            label=f'Mean = {np.mean(dr):.1f}%')
ax1.axvline(np.median(dr), color='#7F1D1D', linestyle='-', linewidth=1.8,
            label=f'Median = {np.median(dr):.1f}%')
ax1.set_xlabel('Recall (%)')
ax1.set_ylabel('Task Count')
ax1.set_title(f'G1 Recall Distribution  (N={len(dr):,})', fontweight='bold')
ax1.legend(loc='upper right', frameon=True, fancybox=True)
ax1.set_xlim(0, 100)
ax1.set_ylim(bottom=0)
# Inset: zoomed 0-50
axins = ax1.inset_axes([0.55, 0.45, 0.42, 0.42])
axins.hist(dr, bins=np.linspace(0, 100, 101), color=C1, edgecolor='white', alpha=0.9, linewidth=0.2)
axins.set_xlim(0, 50)
axins.set_title('Zoom 0-50%', fontsize=8)
axins.tick_params(labelsize=7)
ax1.indicate_inset_zoom(axins, edgecolor='gray', alpha=0.5)

# ===== Panel 2: G1 Precision histogram =====
ax2 = fig.add_subplot(2, 2, 2)
ax2.hist(dp, bins=bins, color=CG, edgecolor='white', alpha=0.88, linewidth=0.3)
ax2.axvline(np.mean(dp), color='#991B1B', linestyle='--', linewidth=2,
            label=f'Mean = {np.mean(dp):.1f}%')
ax2.axvline(np.median(dp), color='#7F1D1D', linestyle='-', linewidth=1.8,
            label=f'Median = {np.median(dp):.1f}%')
ax2.set_xlabel('Precision (%)')
ax2.set_ylabel('Task Count')
ax2.set_title(f'G1 Precision Distribution  (N={len(dp):,})', fontweight='bold')
ax2.legend(loc='upper right', frameon=True, fancybox=True)
ax2.set_xlim(0, 100)
ax2.set_ylim(bottom=0)

# ===== Panel 3: Recall + Precision + Hit Rate by exp_n bins =====
ax3 = fig.add_subplot(2, 2, 3)
exp_bins = [(1,5), (5,10), (10,20), (20,50), (50,150)]
x_labels = ['1–5', '5–10', '10–20', '20–50', '50+']
x = np.arange(len(x_labels))
w = 0.22

g1_rec, g1_prec, g1_hit = [], [], []
for lo, hi in exp_bins:
    sub = [r for r in direct if lo <= r['exp_n'] < hi]
    g1_rec.append(np.mean([r['recall']*100 for r in sub]))
    g1_prec.append(np.mean([r['precision']*100 for r in sub]))
    g1_hit.append(sum(1 for r in sub if r['intersect_n']>0)/len(sub)*100 if sub else 0)

b1 = ax3.bar(x - w, g1_rec, w, color=C1, alpha=0.9, label='Recall', zorder=3)
b2 = ax3.bar(x, g1_prec, w, color=CG, alpha=0.9, label='Precision', zorder=3)
b3 = ax3.bar(x + w, g1_hit, w, color='#F97316', alpha=0.9, label='Hit Rate', zorder=3)

# Value labels
for bars, vals in [(b1, g1_rec), (b2, g1_prec), (b3, g1_hit)]:
    for bar, v in zip(bars, vals):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.6,
                 f'{v:.1f}', ha='center', va='bottom', fontsize=7, fontweight='bold')

ax3.set_xticks(x)
ax3.set_xticklabels(x_labels)
ax3.set_xlabel('Experimental Residue Count (exp_n)')
ax3.set_ylabel('%')
ax3.set_title('G1 Metrics by Residue Set Size', fontweight='bold')
ax3.legend(loc='upper left', frameon=True, fancybox=True)
ax3.set_ylim(0, 105)
ax3.grid(axis='y', alpha=0.3, linewidth=0.5)

# ===== Panel 4: G1 vs G2 Recall CDF overlay =====
ax4 = fig.add_subplot(2, 2, 4)

for label, vals, color, ls in [('G1 Recall', dr, C1, '-'),
                                ('G1 Precision', dp, CG, '--'),
                                ('G2 Recall', d2r, C2, '-'),
                                ('G2 Precision', d2p, '#F97316', '--')]:
    sv = np.sort(vals)
    cdf = np.arange(1, len(sv)+1) / len(sv) * 100
    ax4.plot(sv, cdf, color=color, linestyle=ls, linewidth=2.0, label=label, alpha=0.9)

# Add percentile markers
for v, pct, c in [(np.median(dr), 50, C1), (np.median(dp), 50, CG)]:
    y = np.searchsorted(np.sort(dr), v) / len(dr) * 100 if c == C1 else np.searchsorted(np.sort(dp), v) / len(dp) * 100
    ax4.axvline(v, ymax=y/100, color=c, linestyle=':', alpha=0.5, linewidth=1)

ax4.set_xlabel('%')
ax4.set_ylabel('Cumulative % of Tasks')
ax4.set_title('CDF: G1 vs G2  —  Recall & Precision', fontweight='bold')
ax4.legend(loc='lower right', frameon=True, fancybox=True)
ax4.set_xlim(0, 100)
ax4.set_ylim(0, 105)
ax4.grid(True, alpha=0.3, linewidth=0.5)

# ---- Summary text box ----
summary_text = (
    f'G1 (N={len(dr):,}):  '
    f'Recall μ={np.mean(dr):.1f}%  M={np.median(dr):.1f}%  |  '
    f'Precision μ={np.mean(dp):.1f}%  M={np.median(dp):.1f}%  |  '
    f'Hit Rate {sum(1 for r in direct if r["intersect_n"]>0)/len(direct)*100:.1f}%\n'
    f'G2 (N={len(d2r):,}):  '
    f'Recall μ={np.mean(d2r):.1f}%  M={np.median(d2r):.1f}%  |  '
    f'Precision μ={np.mean(d2p):.1f}%  M={np.median(d2p):.1f}%  |  '
    f'Hit Rate {sum(1 for r in dbref if r["intersect_n"]>0)/len(dbref)*100:.1f}%'
)
fig.text(0.5, 0.01, summary_text, ha='center', va='bottom', fontsize=8.5,
         fontfamily='monospace', bbox=dict(boxstyle='round,pad=0.5', facecolor='#F8FAFC', edgecolor='#CBD5E1'))

plt.tight_layout(rect=[0, 0.04, 1, 0.95])
plt.savefig('C:/Users/Administrator/docking_accuracy_report.png')
plt.savefig('D:/finale/negative_upgrade/docking_package/docking_accuracy_report.png')
print('Done! Saved docking_accuracy_report.png')
