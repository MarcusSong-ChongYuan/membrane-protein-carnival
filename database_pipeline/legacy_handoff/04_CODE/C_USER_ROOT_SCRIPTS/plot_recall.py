import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

with open('/home/csong/docking/docking_package_trial/residue_recall_v2.jsonl') as f:
    data = [json.loads(l) for l in f if l.strip()]

direct = [r['recall']*100 for r in data if r['mapping'] == 'direct']
dbref = [r['recall']*100 for r in data if r['mapping'] == 'dbref_segment']

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('Binding Residue Recall: Docked vs Experimental', fontsize=16, fontweight='bold')

# ---- Panel 1: Histogram G1 direct ----
ax = axes[0, 0]
bins = np.linspace(0, 100, 51)
ax.hist(direct, bins=bins, color='#2196F3', edgecolor='white', alpha=0.85)
ax.axvline(np.mean(direct), color='red', linestyle='--', linewidth=2, label=f'Mean={np.mean(direct):.1f}%')
ax.axvline(np.median(direct), color='darkred', linestyle='-', linewidth=2, label=f'Median={np.median(direct):.1f}%')
ax.set_xlabel('Recall (%)')
ax.set_ylabel('Task Count')
ax.set_title(f'G1: direct mapping (N={len(direct):,})')
ax.legend(fontsize=9)
ax.set_xlim(0, 100)

# ---- Panel 2: Histogram G2 dbref_segment ----
ax = axes[0, 1]
ax.hist(dbref, bins=bins, color='#FF9800', edgecolor='white', alpha=0.85)
ax.axvline(np.mean(dbref), color='red', linestyle='--', linewidth=2, label=f'Mean={np.mean(dbref):.1f}%')
ax.axvline(np.median(dbref), color='darkred', linestyle='-', linewidth=2, label=f'Median={np.median(dbref):.1f}%')
ax.set_xlabel('Recall (%)')
ax.set_ylabel('Task Count')
ax.set_title(f'G2: dbref_segment (N={len(dbref):,})')
ax.legend(fontsize=9)
ax.set_xlim(0, 100)

# ---- Panel 3: Overlaid CDF ----
ax = axes[1, 0]
for label, vals, color in [('G1: direct', direct, '#2196F3'), ('G2: dbref_segment', dbref, '#FF9800')]:
    sorted_vals = np.sort(vals)
    cdf = np.arange(1, len(sorted_vals)+1) / len(sorted_vals) * 100
    ax.plot(sorted_vals, cdf, color=color, linewidth=2, label=f'{label} (N={len(vals):,})')
ax.set_xlabel('Recall (%)')
ax.set_ylabel('Cumulative % of Tasks')
ax.set_title('Cumulative Distribution (CDF)')
ax.legend(fontsize=9)
ax.set_xlim(0, 100)
ax.grid(True, alpha=0.3)
# Add reference lines
for thresh, ls in [(25, '--'), (50, '-.'), (75, ':')]:
    ax.axvline(thresh, color='gray', linestyle=ls, alpha=0.4)
    ax.text(thresh+0.5, 95, f'{thresh}%', fontsize=7, color='gray')

# ---- Panel 4: Recall by exp_n bin (boxplot style as bar) ----
ax = axes[1, 1]
exp_bins = [(1,5), (5,10), (10,20), (20,50), (50,150)]
x_labels = []
g1_means, g1_meds = [], []
g2_means, g2_meds = [], []
for lo, hi in exp_bins:
    g1_sub = [r['recall']*100 for r in data if r['mapping']=='direct' and lo <= r['exp_n'] < hi]
    g2_sub = [r['recall']*100 for r in data if r['mapping']=='dbref_segment' and lo <= r['exp_n'] < hi]
    x_labels.append(f'{lo}-{hi}')
    g1_means.append(np.mean(g1_sub) if g1_sub else 0)
    g1_meds.append(np.median(g1_sub) if g1_sub else 0)
    g2_means.append(np.mean(g2_sub) if g2_sub else 0)
    g2_meds.append(np.median(g2_sub) if g2_sub else 0)

x = np.arange(len(x_labels))
w = 0.35
bars1 = ax.bar(x - w/2, g1_means, w, color='#2196F3', alpha=0.85, label='G1 mean')
bars2 = ax.bar(x + w/2, g2_means, w, color='#FF9800', alpha=0.85, label='G2 mean')
# Overlay median markers
ax.scatter(x - w/2, g1_meds, marker='_', color='darkblue', s=200, linewidth=2, zorder=5)
ax.scatter(x + w/2, g2_meds, marker='_', color='darkred', s=200, linewidth=2, zorder=5)
ax.set_xticks(x)
ax.set_xticklabels(x_labels)
ax.set_xlabel('Experimental Residue Count (exp_n)')
ax.set_ylabel('Recall (%)')
ax.set_title('Recall by Experimental Residue Set Size')
ax.legend(fontsize=9)

# Add value labels on bars
for bar, val in zip(bars1, g1_means):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, f'{val:.1f}%',
            ha='center', va='bottom', fontsize=8, fontweight='bold')
for bar, val in zip(bars2, g2_means):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, f'{val:.1f}%',
            ha='center', va='bottom', fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig('/home/csong/docking/docking_package_trial/recall_plot.png', dpi=150, bbox_inches='tight')
print('Saved to recall_plot.png')
print()
print('G1 direct stats:')
print('  N=%d  Mean=%.1f%%  Median=%.1f%%  Q1=%.1f%%  Q3=%.1f%%' % (
    len(direct), np.mean(direct), np.median(direct),
    np.percentile(direct, 25), np.percentile(direct, 75)))
print('  0%% recall: %d (%.1f%%)' % (sum(1 for d in direct if d == 0), 100*sum(1 for d in direct if d == 0)/len(direct)))
print('  >50%% recall: %d (%.1f%%)' % (sum(1 for d in direct if d > 50), 100*sum(1 for d in direct if d > 50)/len(direct)))
print()
print('G2 dbref_segment stats:')
print('  N=%d  Mean=%.1f%%  Median=%.1f%%' % (len(dbref), np.mean(dbref), np.median(dbref)))
print('  0%% recall: %d (%.1f%%)' % (sum(1 for d in dbref if d == 0), 100*sum(1 for d in dbref if d == 0)/len(dbref)))
