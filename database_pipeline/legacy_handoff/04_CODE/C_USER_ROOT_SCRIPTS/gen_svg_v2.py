import json, math

with open('C:/Users/Administrator/residue_recall_v2.jsonl') as f:
    data = [json.loads(l) for l in f if l.strip()]

direct = [r for r in data if r['mapping'] == 'direct']
dbref = [r for r in data if r['mapping'] == 'dbref_segment']

def mean(arr): return sum(arr)/len(arr) if arr else 0
def median(arr):
    s = sorted(arr); return s[len(s)//2] if s else 0
def pc(arr, p):
    s = sorted(arr); return s[int(len(s)*p/100)] if s else 0

SW, SH = 500, 380
ML, MR, MT, MB = 60, 20, 55, 60
PW = SW - ML - MR
PH = SH - MT - MB

def histogram(data, color, title, label_x='%'):
    nbins = 40
    max_cnt = 0
    counts = [0]*nbins
    for v in data:
        idx = min(int(v * nbins / 100), nbins-1)
        counts[idx] += 1
        max_cnt = max(max_cnt, counts[idx])
    bar_w = PW / nbins
    scale = PH / max(max_cnt, 1)
    svg = ''
    for i in range(nbins):
        h = counts[i] * scale
        if h > 0:
            x = ML + i * bar_w
            y = MT + PH - h
            svg += f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w-0.5:.1f}" height="{h:.1f}" fill="{color}" opacity="0.85"/>\n'
    m = mean(data); med = median(data)
    mx = ML + m/100 * PW
    medx = ML + med/100 * PW
    svg += f'<line x1="{mx:.1f}" y1="{MT}" x2="{mx:.1f}" y2="{MT+PH}" stroke="red" stroke-width="1.5" stroke-dasharray="5,3"/>\n'
    svg += f'<text x="{mx:.1f}" y="{MT-12}" text-anchor="middle" fill="red" font-size="10">Mean={m:.1f}{label_x}</text>\n'
    svg += f'<line x1="{medx:.1f}" y1="{MT}" x2="{medx:.1f}" y2="{MT+PH}" stroke="darkred" stroke-width="1"/>\n'
    svg += f'<text x="{medx:.1f}" y="{MT-24}" text-anchor="middle" fill="darkred" font-size="9">Med={med:.1f}{label_x}</text>\n'
    # axes
    svg += f'<line x1="{ML}" y1="{MT+PH}" x2="{ML+PW}" y2="{MT+PH}" stroke="black" stroke-width="1"/>\n'
    svg += f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+PH}" stroke="black" stroke-width="1"/>\n'
    for pct in [0, 25, 50, 75, 100]:
        x = ML + pct/100 * PW
        svg += f'<text x="{x:.1f}" y="{MT+PH+16}" text-anchor="middle" font-size="9">{pct}{label_x}</text>\n'
        svg += f'<line x1="{x:.1f}" y1="{MT+PH}" x2="{x:.1f}" y2="{MT+PH+4}" stroke="black" stroke-width="0.5"/>\n'
    svg += f'<text x="{SW/2:.1f}" y="{MT+PH+38}" text-anchor="middle" font-size="11">{label_x}</text>\n'
    svg += f'<text x="{ML-38:.1f}" y="{MT+PH/2:.1f}" text-anchor="middle" font-size="10" transform="rotate(-90,{ML-38},{MT+PH/2})">Tasks</text>\n'
    svg += f'<text x="{SW/2:.1f}" y="{MT-35}" text-anchor="middle" font-size="13" font-weight="bold">{title}</text>\n'
    stats = f'N={len(data)}  Mean={m:.1f}{label_x}  Med={med:.1f}{label_x}  P25={pc(data,25):.0f}{label_x}  P75={pc(data,75):.0f}{label_x}'
    svg += f'<text x="{SW/2:.1f}" y="{MT-10}" text-anchor="middle" font-size="9" fill="#555">{stats}</text>\n'
    return svg

def grouped_bar(x_labels, groups, colors, legend_names, title, ymax=100):
    """groups = [[vals_per_x], ...] for each series"""
    n_groups = len(x_labels)
    n_bars = len(groups)
    bar_w = PW / (n_groups * n_bars + n_groups) * 0.65
    gap = bar_w * 0.5
    max_val = ymax

    svg = ''
    for gi in range(n_groups):
        for bi in range(n_bars):
            v = groups[bi][gi]
            h = v / max_val * PH
            x = ML + gi * (PW / n_groups) + bi * (bar_w + gap)
            y = MT + PH - h
            svg += f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{colors[bi]}" opacity="0.88" rx="1"/>\n'
            svg += f'<text x="{x+bar_w/2:.1f}" y="{y-3:.1f}" text-anchor="middle" font-size="8" font-weight="bold">{v:.1f}</text>\n'

    # legend
    for bi in range(n_bars):
        svg += f'<rect x="{ML+4}" y="{MT-40+bi*13}" width="9" height="9" fill="{colors[bi]}" rx="1"/>\n'
        svg += f'<text x="{ML+16}" y="{MT-31+bi*13}" font-size="9">{legend_names[bi]}</text>\n'

    svg += f'<line x1="{ML}" y1="{MT+PH}" x2="{ML+PW}" y2="{MT+PH}" stroke="black" stroke-width="1"/>\n'
    svg += f'<line x1="{ML}" y1="{MT}" x2="{ML}" y2="{MT+PH}" stroke="black" stroke-width="1"/>\n'
    for gi in range(n_groups):
        x = ML + gi * (PW / n_groups) + (n_bars-1) * (bar_w + gap) / 2
        svg += f'<text x="{x:.1f}" y="{MT+PH+15}" text-anchor="middle" font-size="9">{x_labels[gi]}</text>\n'
    svg += f'<text x="{SW/2:.1f}" y="{MT+PH+35}" text-anchor="middle" font-size="11">exp_n (Experimental Residue Count)</text>\n'
    svg += f'<text x="15" y="{MT+PH/2:.1f}" text-anchor="middle" font-size="10" transform="rotate(-90,15,{MT+PH/2})">%</text>\n'
    svg += f'<text x="{SW/2:.1f}" y="{MT-32}" text-anchor="middle" font-size="13" font-weight="bold">{title}</text>\n'
    return svg

# ---- Compute data ----
g1_recall = [r['recall']*100 for r in direct]
g1_precision = [r['precision']*100 for r in direct]
g2_recall = [r['recall']*100 for r in dbref]
g2_precision = [r['precision']*100 for r in dbref]

# By exp_n bins
exp_bins = [(1,5), (5,10), (10,20), (20,50), (50,150)]
x_labels = ['1-5', '5-10', '10-20', '20-50', '50+']
g1_recall_bin, g1_prec_bin, g1_hit_bin = [], [], []
g2_recall_bin, g2_prec_bin, g2_hit_bin = [], [], []
g1_count_bin, g2_count_bin = [], []

for lo, hi in exp_bins:
    g1s = [r for r in direct if lo <= r['exp_n'] < hi]
    g2s = [r for r in dbref if lo <= r['exp_n'] < hi]
    g1_count_bin.append(len(g1s))
    g2_count_bin.append(len(g2s))
    g1_recall_bin.append(mean([r['recall']*100 for r in g1s]))
    g1_prec_bin.append(mean([r['precision']*100 for r in g1s]))
    g1_hit_bin.append(sum(1 for r in g1s if r['intersect_n']>0)/len(g1s)*100 if g1s else 0)
    g2_recall_bin.append(mean([r['recall']*100 for r in g2s]))
    g2_prec_bin.append(mean([r['precision']*100 for r in g2s]))
    g2_hit_bin.append(sum(1 for r in g2s if r['intersect_n']>0)/len(g2s)*100 if g2s else 0)

# ---- Build SVG ----
TW = SW * 2
TH = SH * 2 + 30

def panel(svg, col, row):
    return f'<g transform="translate({col*SW},{row*SH+30})">\n{svg}\n</g>\n'

# Panel 1: G1 Recall
p1 = histogram(g1_recall, '#2196F3', 'G1 Recall (direct mapping)')
# Panel 2: G1 Precision
p2 = histogram(g1_precision, '#4CAF50', 'G1 Precision (direct mapping)')
# Panel 3: G1 metrics by exp_n
p3 = grouped_bar(x_labels,
    [g1_recall_bin, g1_prec_bin, g1_hit_bin],
    ['#2196F3', '#4CAF50', '#FF5722'],
    ['Recall', 'Precision', 'Hit Rate'],
    'G1: Recall + Precision + Hit Rate by exp_n')
# Panel 4: G1 vs G2 recall comparison by exp_n
p4 = grouped_bar(x_labels,
    [g1_recall_bin, g2_recall_bin],
    ['#2196F3', '#FF9800'],
    ['G1 Recall', 'G2 Recall'],
    'G1 vs G2: Recall by exp_n')

full = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{TW}" height="{TH}" viewBox="0 0 {TW} {TH}">
<rect width="{TW}" height="{TH}" fill="white"/>
<text x="{TW/2}" y="22" text-anchor="middle" font-size="16" font-weight="bold">Docking Residue Accuracy: Recall · Precision · Hit Rate (N=7,839)</text>
{panel(p1, 0, 0)}
{panel(p2, 1, 0)}
{panel(p3, 0, 1)}
{panel(p4, 1, 1)}
</svg>'''

with open('C:/Users/Administrator/recall_precision_hit.svg', 'w', encoding='utf-8') as f:
    f.write(full)
with open('D:/finale/negative_upgrade/docking_package/recall_precision_hit.svg', 'w', encoding='utf-8') as f:
    f.write(full)

print(f'Done: {len(full):,} bytes')
print()
print('=== Summary Stats ===')
print(f'G1 (N={len(direct)}):')
print(f'  Recall:    Mean={mean(g1_recall):.1f}%  Med={median(g1_recall):.1f}%')
print(f'  Precision: Mean={mean(g1_precision):.1f}%  Med={median(g1_precision):.1f}%')
print(f'  Hit Rate:  {sum(1 for r in direct if r["intersect_n"]>0)/len(direct)*100:.1f}%')
print(f'G2 (N={len(dbref)}):')
print(f'  Recall:    Mean={mean(g2_recall):.1f}%  Med={median(g2_recall):.1f}%')
print(f'  Precision: Mean={mean(g2_precision):.1f}%  Med={median(g2_precision):.1f}%')
print(f'  Hit Rate:  {sum(1 for r in dbref if r["intersect_n"]>0)/len(dbref)*100:.1f}%')
print()
print('By exp_n bins (G1):')
for i, (lo, hi) in enumerate(exp_bins):
    print(f'  exp_n {lo:3d}-{hi:3d}: N={g1_count_bin[i]:5d}  Recall={g1_recall_bin[i]:5.1f}%  Prec={g1_prec_bin[i]:5.1f}%  Hit={g1_hit_bin[i]:5.1f}%')
