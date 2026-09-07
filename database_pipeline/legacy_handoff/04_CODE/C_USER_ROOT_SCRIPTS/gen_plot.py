import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

with open('C:/Users/Administrator/residue_recall_v2.jsonl') as f:
    data = [json.loads(l) for l in f if l.strip()]

direct = [r['recall']*100 for r in data if r['mapping'] == 'direct']
dbref = [r['recall']*100 for r in data if r['mapping'] == 'dbref_segment']

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('Binding Residue Recall: Docked vs Experimental (N=7,839)', fontsize=16, fontweight='bold')

# Panel 1: G1 direct histogram
ax = axes[0, 0]
bins = np.linspace(0, 100, 51)
ax.hist(direct, bins=bins, color='#2196F3', edgecolor='white', alpha=0.85)
ax.axvline(np.mean(direct), color='red', linestyle='--', linewidth=2,
           label='Mean = {:.1f}%'.format(np.mean(direct)))
ax.axvline(np.median(direct), color='darkred', linestyle='-', linewidth=2,
           label='Median = {:.1f}%'.format(np.median(direct)))
ax.set_xlabel('Recall (%)')
ax.set_ylabel('Task Count')
ax.set_title('G1: direct mapping (N={:,})'.format(len(direct)))
ax.legend(fontsize=9)
ax.set_xlim(0, 100)

# Panel 2: G2 dbref histogram
ax = axes[0, 1]
ax.hist(dbref, bins=bins, color='#FF9800', edgecolor='white', alpha=0.85)
ax.axvline(np.mean(dbref), color='red', linestyle='--', linewidth=2,
           label='Mean = {:.1f}%'.format(np.mean(dbref)))
ax.axvline(np.median(dbref), color='darkred', linestyle='-', linewidth=2,
           label='Median = {:.1f}%'.format(np.median(dbref)))
ax.set_xlabel('Recall (%)')
ax.set_ylabel('Task Count')
ax.set_title('G2: dbref_segment (N={:,})'.format(len(dbref)))
ax.legend(fontsize=9)
ax.set_xlim(0, 100)

# Panel 3: CDF overlay
ax = axes[1, 0]
for label, vals, color in [('G1: direct', direct, '#2196F3'),
                            ('G2: dbref_segment', dbref, '#FF9800')]:
    sv = np.sort(vals)
    cdf = np.arange(1, len(sv)+1) / len(sv) * 100
    ax.plot(sv, cdf, color=color, linewidth=2.5, label='{} (N={:,})'.format(label, len(vals)))
ax.set_xlabel('Recall (%)')
ax.set_ylabel('Cumulative % of Tasks')
ax.set_title('Cumulative Distribution (CDF)')
ax.legend(fontsize=9)
ax.set_xlim(0, 100)
ax.grid(True, alpha=0.3)
for thresh in [25, 50, 75]:
    ax.axvline(thresh, color='gray', linestyle='--', alpha=0.3)

# Panel 4: Recall & hit rate by exp_n bins
ax = axes[1, 1]
exp_bins = [(1,5), (5,10), (10,20), (20,50), (50,150)]
x_labels = ['{}-{}'.format(lo, hi) for lo, hi in exp_bins]
g1_means, g1_pcts = [], []
g2_means, g2_pcts = [], []
g1_ns, g2_ns = [], []
for lo, hi in exp_bins:
    g1s = [r['recall']*100 for r in data if r['mapping']=='direct' and lo <= r['exp_n'] < hi]
    g2s = [r['recall']*100 for r in data if r['mapping']=='dbref_segment' and lo <= r['exp_n'] < hi]
    g1_means.append(np.mean(g1s) if g1s else 0)
    g2_means.append(np.mean(g2s) if g2s else 0)
    g1_pcts.append(sum(1 for v in g1s if v > 0)/len(g1s)*100 if g1s else 0)
    g2_pcts.append(sum(1 for v in g2s if v > 0)/len(g2s)*100 if g2s else 0)
    g1_ns.append(len(g1s))
    g2_ns.append(len(g2s))

x = np.arange(len(x_labels))
w = 0.30

bars1 = ax.bar(x - w/2, g1_means, w, color='#2196F3', alpha=0.9, label='G1: mean recall')
bars3 = ax.bar(x - w/2, g1_pcts, w*0.5, color='#1565C0', alpha=0.6, label='G1: hit rate (any contact)')
bars2 = ax.bar(x + w/2, g2_means, w, color='#FF9800', alpha=0.9, label='G2: mean recall')
ax.set_xticks(x)
ax.set_xticklabels(x_labels)
ax.set_xlabel('Experimental Residue Count (exp_n)')
ax.set_ylabel('%')
ax.set_title('Recall & Hit Rate by Exp Residue Set Size')
ax.legend(fontsize=7)
ax.set_ylim(0, max(max(g1_pcts)+10, 80))

# Add N labels
for i, (n1, n2) in enumerate(zip(g1_ns, g2_ns)):
    ax.text(i - w/2, 3, 'n={}'.format(n1), ha='center', fontsize=6, color='darkblue')
    ax.text(i + w/2, 3, 'n={}'.format(n2), ha='center', fontsize=6, color='darkred')

plt.tight_layout()
plt.savefig('C:/Users/Administrator/recall_plot.png', dpi=150, bbox_inches='tight')
print('Saved to recall_plot.png')
print()
print('G1 direct: N={} Mean={:.1f}% Median={:.1f}%'.format(len(direct), np.mean(direct), np.median(direct)))
print('  Zero recall: {} ({:.1f}%)'.format(sum(1 for d in direct if d==0), 100*sum(1 for d in direct if d==0)/len(direct)))
print('  Recall >50%: {} ({:.1f}%)'.format(sum(1 for d in direct if d>50), 100*sum(1 for d in direct if d>50)/len(direct)))
print()
print('G2 dbref_segment: N={} Mean={:.1f}% Median={:.1f}%'.format(len(dbref), np.mean(dbref), np.median(dbref)))
print('  Zero recall: {} ({:.1f}%)'.format(sum(1 for d in dbref if d==0), 100*sum(1 for d in dbref if d==0)/len(dbref)))
print('  Recall >50%: {} ({:.1f}%)'.format(sum(1 for d in dbref if d>50), 100*sum(1 for d in dbref if d>50)/len(dbref)))
