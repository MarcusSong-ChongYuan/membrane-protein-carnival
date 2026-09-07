import sys, os
sys.path = [p for p in sys.path if 'MGLTools' not in p]
venv_sp = r'C:\Users\Administrator\plotenv\Lib\site-packages'
if venv_sp not in sys.path: sys.path.insert(0, venv_sp)

import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Microsoft YaHei', 'Arial'],
    'font.size': 10, 'axes.titlesize': 13, 'axes.labelsize': 11,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 150, 'savefig.dpi': 200, 'savefig.bbox': 'tight',
    'axes.spines.top': False, 'axes.spines.right': False,
})

with open('C:/Users/Administrator/residue_recall_55.jsonl') as f:
    data = [json.loads(l) for l in f if l.strip()]

direct = [r for r in data if r['mapping'] == 'direct']
dbref = [r for r in data if r['mapping'] == 'dbref_segment']

g1_r = [r['recall']*100 for r in direct]
g1_p = [r['precision']*100 for r in direct]
g2_r = [r['recall']*100 for r in dbref]
g2_p = [r['precision']*100 for r in dbref]
g1_hit = sum(1 for r in direct if r['intersect_n']>0)/len(direct)*100
g2_hit = sum(1 for r in dbref if r['intersect_n']>0)/len(dbref)*100

C1, C1L = '#2563EB', '#93C5FD'
C2, C2L = '#EA580C', '#FDBA74'
CG = '#16A34A'
CO = '#F97316'

fig = plt.figure(figsize=(18, 12))
fig.suptitle('Docking Residue Accuracy — 5.5Å Contact Cutoff', fontsize=18, fontweight='bold', y=0.98)

# ===== Panel 1: G1 Recall =====
ax1 = fig.add_subplot(2, 3, 1)
bins = np.linspace(0, 100, 51)
ax1.hist(g1_r, bins=bins, color=C1, edgecolor='white', alpha=0.88, linewidth=0.3)
ax1.axvline(np.mean(g1_r), color='#991B1B', ls='--', lw=2, label=f'Mean={np.mean(g1_r):.1f}%')
ax1.axvline(np.median(g1_r), color='#7F1D1D', ls='-', lw=1.8, label=f'Median={np.median(g1_r):.1f}%')
ax1.set_xlabel('Recall (%)'); ax1.set_ylabel('Task Count')
ax1.set_title(f'G1 Recall  (N={len(g1_r):,})', fontweight='bold')
ax1.legend(loc='upper right', frameon=True); ax1.set_xlim(0, 100)

# ===== Panel 2: G1 Precision =====
ax2 = fig.add_subplot(2, 3, 2)
ax2.hist(g1_p, bins=bins, color=CG, edgecolor='white', alpha=0.88, linewidth=0.3)
ax2.axvline(np.mean(g1_p), color='#991B1B', ls='--', lw=2, label=f'Mean={np.mean(g1_p):.1f}%')
ax2.axvline(np.median(g1_p), color='#7F1D1D', ls='-', lw=1.8, label=f'Median={np.median(g1_p):.1f}%')
ax2.set_xlabel('Precision (%)'); ax2.set_ylabel('Task Count')
ax2.set_title(f'G1 Precision  (N={len(g1_p):,})', fontweight='bold')
ax2.legend(loc='upper right', frameon=True); ax2.set_xlim(0, 100)

# ===== Panel 3: Hit Rate pie =====
ax3 = fig.add_subplot(2, 3, 3)
sizes = [g1_hit, 100-g1_hit]
colors = [C1, '#E5E7EB']
wedges, texts, autotexts = ax3.pie(sizes, labels=['Hit', 'Miss'], colors=colors,
    autopct='%1.1f%%', startangle=90, explode=(0.02, 0),
    textprops={'fontsize': 11})
autotexts[0].set_fontweight('bold'); autotexts[0].set_fontsize(14)
ax3.set_title(f'G1 Hit Rate\n(any exp residue contacted)', fontweight='bold')

# ===== Panel 4: CDF G1 vs G2 =====
ax4 = fig.add_subplot(2, 3, 4)
for label, vals, color, ls in [('G1 Recall', g1_r, C1, '-'),
                                ('G1 Precision', g1_p, CG, '--'),
                                ('G2 Recall', g2_r, C2, '-'),
                                ('G2 Precision', g2_p, '#F97316', '--')]:
    sv = np.sort(vals); cdf = np.arange(1, len(sv)+1)/len(sv)*100
    ax4.plot(sv, cdf, color=color, ls=ls, lw=2, label=label, alpha=0.9)
ax4.set_xlabel('%'); ax4.set_ylabel('Cumulative % of Tasks')
ax4.set_title('CDF: G1 vs G2', fontweight='bold')
ax4.legend(loc='lower right', frameon=True); ax4.set_xlim(0, 100); ax4.set_ylim(0, 105)
ax4.grid(True, alpha=0.3, lw=0.5)

# ===== Panel 5: G1 metrics by exp_n =====
ax5 = fig.add_subplot(2, 3, 5)
exp_bins = [(1,5),(5,10),(10,20),(20,50),(50,150)]
xl = ['1–5','5–10','10–20','20–50','50+']
x = np.arange(len(xl)); w = 0.22
g1r_b, g1p_b, g1h_b = [], [], []
for lo, hi in exp_bins:
    sub = [r for r in direct if lo <= r['exp_n'] < hi]
    g1r_b.append(np.mean([r['recall']*100 for r in sub]))
    g1p_b.append(np.mean([r['precision']*100 for r in sub]))
    g1h_b.append(sum(1 for r in sub if r['intersect_n']>0)/len(sub)*100 if sub else 0)
b1 = ax5.bar(x-w, g1r_b, w, color=C1, alpha=0.9, label='Recall')
b2 = ax5.bar(x, g1p_b, w, color=CG, alpha=0.9, label='Precision')
b3 = ax5.bar(x+w, g1h_b, w, color=CO, alpha=0.9, label='Hit Rate')
for bars, vals in [(b1,g1r_b),(b2,g1p_b),(b3,g1h_b)]:
    for bar, v in zip(bars, vals):
        ax5.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.6,
                 f'{v:.1f}', ha='center', va='bottom', fontsize=7, fontweight='bold')
ax5.set_xticks(x); ax5.set_xticklabels(xl)
ax5.set_xlabel('exp_n'); ax5.set_ylabel('%')
ax5.set_title('G1 by Residue Set Size', fontweight='bold')
ax5.legend(loc='upper left', frameon=True); ax5.set_ylim(0, 105); ax5.grid(axis='y', alpha=0.3, lw=0.5)

# ===== Panel 6: 4Å vs 5.5Å comparison =====
ax6 = fig.add_subplot(2, 3, 6)
# Load 4A data for comparison
with open('C:/Users/Administrator/residue_recall_v2.jsonl') as f:
    d4 = [json.loads(l) for l in f if l.strip()]
d4r = [r['recall']*100 for r in d4 if r['mapping']=='direct']

# Paired comparison: mean recall
cutoffs = ['4Å', '5.5Å']
g1_means = [np.mean(d4r), np.mean(g1_r)]
g1_meds = [np.median(d4r), np.median(g1_r)]
g1_hits = [sum(1 for r in d4 if r['mapping']=='direct' and r['intersect_n']>0)/sum(1 for r in d4 if r['mapping']=='direct')*100, g1_hit]
g2_means = [np.mean([r['recall']*100 for r in d4 if r['mapping']=='dbref_segment']), np.mean(g2_r)]
g2_hits = [sum(1 for r in d4 if r['mapping']=='dbref_segment' and r['intersect_n']>0)/sum(1 for r in d4 if r['mapping']=='dbref_segment')*100, g2_hit]

x2 = np.arange(len(cutoffs)); w2 = 0.2
ax6.bar(x2-1.5*w2, g1_means, w2, color=C1, alpha=0.9, label='G1 Recall')
ax6.bar(x2-0.5*w2, g1_hits, w2, color=C1L, alpha=0.9, label='G1 Hit Rate')
ax6.bar(x2+0.5*w2, g2_means, w2, color=C2, alpha=0.9, label='G2 Recall')
ax6.bar(x2+1.5*w2, g2_hits, w2, color=C2L, alpha=0.9, label='G2 Hit Rate')
for i in range(2):
    for j, vals in enumerate([g1_means, g1_hits, g2_means, g2_hits]):
        ax6.text(x2[i]+(j-1.5)*w2, vals[i]+0.6, f'{vals[i]:.1f}', ha='center', fontsize=7, fontweight='bold')
ax6.set_xticks(x2); ax6.set_xticklabels(cutoffs)
ax6.set_ylabel('%'); ax6.set_ylim(0, 75)
ax6.set_title('4Å vs 5.5Å Cutoff Comparison', fontweight='bold')
ax6.legend(loc='upper right', frameon=True, fontsize=8)
ax6.grid(axis='y', alpha=0.3, lw=0.5)

# Summary
summary = (
    f'Cutoff: 5.5Å  |  '
    f'G1 (N={len(g1_r):,}): Recall μ={np.mean(g1_r):.1f}% M={np.median(g1_r):.1f}%  '
    f'Precision μ={np.mean(g1_p):.1f}% M={np.median(g1_p):.1f}%  Hit={g1_hit:.1f}%  |  '
    f'G2 (N={len(g2_r):,}): Recall μ={np.mean(g2_r):.1f}% M={np.median(g2_r):.1f}%  Hit={g2_hit:.1f}%'
)
fig.text(0.5, 0.01, summary, ha='center', va='bottom', fontsize=8,
         fontfamily='monospace', bbox=dict(boxstyle='round,pad=0.5', facecolor='#F8FAFC', edgecolor='#CBD5E1'))

plt.tight_layout(rect=[0, 0.04, 1, 0.95])
plt.savefig('C:/Users/Administrator/docking_55A_report.png')
plt.savefig('D:/finale/negative_upgrade/docking_package/docking_55A_report.png')
print('Done!')
print(f'G1 5.5A: Recall mean={np.mean(g1_r):.1f}% med={np.median(g1_r):.1f}%  Prec mean={np.mean(g1_p):.1f}%  Hit={g1_hit:.1f}%')
print(f'G2 5.5A: Recall mean={np.mean(g2_r):.1f}% med={np.median(g2_r):.1f}%  Hit={g2_hit:.1f}%')
