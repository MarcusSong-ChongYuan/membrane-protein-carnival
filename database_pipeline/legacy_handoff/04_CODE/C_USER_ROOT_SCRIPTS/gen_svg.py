import json, math
from collections import defaultdict

with open('C:/Users/Administrator/residue_recall_v2.jsonl') as f:
    data = [json.loads(l) for l in f if l.strip()]

direct = [r['recall']*100 for r in data if r['mapping'] == 'direct']
dbref = [r['recall']*100 for r in data if r['mapping'] == 'dbref_segment']

def mean(arr): return sum(arr)/len(arr) if arr else 0
def median(arr):
    s = sorted(arr)
    return s[len(s)//2] if s else 0
def percentile(arr, p):
    s = sorted(arr)
    return s[int(len(s)*p/100)] if s else 0

W, H = 1000, 700
MARGIN_L, MARGIN_R, MARGIN_T, MARGIN_B = 80, 30, 50, 80
PLOT_W = W - MARGIN_L - MARGIN_R
PLOT_H = H - MARGIN_T - MARGIN_B

def histogram_svg(data, title, color, nbins=40):
    max_count = 0
    counts = [0]*nbins
    for v in data:
        idx = min(int(v * nbins / 100), nbins-1)
        counts[idx] += 1
        max_count = max(max_count, counts[idx])

    bar_w = PLOT_W / nbins
    scale_y = PLOT_H / max(max_count, 1)

    svg = ''
    for i in range(nbins):
        h = counts[i] * scale_y
        x = MARGIN_L + i * bar_w
        y = MARGIN_T + PLOT_H - h
        svg += f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w-1:.1f}" height="{h:.1f}" fill="{color}" opacity="0.85"/>\n'

    # Mean and median lines
    m = mean(data)
    med = median(data)
    mx = MARGIN_L + m/100 * PLOT_W
    medx = MARGIN_L + med/100 * PLOT_W
    svg += f'<line x1="{mx:.1f}" y1="{MARGIN_T}" x2="{mx:.1f}" y2="{MARGIN_T+PLOT_H}" stroke="red" stroke-width="2" stroke-dasharray="6,3"/>\n'
    svg += f'<text x="{mx:.1f}" y="{MARGIN_T-10}" text-anchor="middle" fill="red" font-size="11">Mean={m:.1f}%</text>\n'
    svg += f'<line x1="{medx:.1f}" y1="{MARGIN_T}" x2="{medx:.1f}" y2="{MARGIN_T+PLOT_H}" stroke="darkred" stroke-width="1.5"/>\n'
    svg += f'<text x="{medx:.1f}" y="{MARGIN_T-22}" text-anchor="middle" fill="darkred" font-size="10">Med={med:.1f}%</text>\n'

    # Axes
    svg += f'<line x1="{MARGIN_L}" y1="{MARGIN_T+PLOT_H}" x2="{MARGIN_L+PLOT_W}" y2="{MARGIN_T+PLOT_H}" stroke="black" stroke-width="1"/>\n'
    svg += f'<line x1="{MARGIN_L}" y1="{MARGIN_T}" x2="{MARGIN_L}" y2="{MARGIN_T+PLOT_H}" stroke="black" stroke-width="1"/>\n'

    # X labels
    for pct in [0, 25, 50, 75, 100]:
        x = MARGIN_L + pct/100 * PLOT_W
        svg += f'<text x="{x:.1f}" y="{MARGIN_T+PLOT_H+20}" text-anchor="middle" font-size="11">{pct}%</text>\n'
        svg += f'<line x1="{x:.1f}" y1="{MARGIN_T+PLOT_H}" x2="{x:.1f}" y2="{MARGIN_T+PLOT_H+5}" stroke="black" stroke-width="1"/>\n'

    svg += f'<text x="{W/2:.1f}" y="{MARGIN_T+PLOT_H+45}" text-anchor="middle" font-size="13">Recall (%)</text>\n'
    svg += f'<text x="{MARGIN_L-50:.1f}" y="{MARGIN_T+PLOT_H/2:.1f}" text-anchor="middle" font-size="13" transform="rotate(-90,{MARGIN_L-50},{MARGIN_T+PLOT_H/2})">Task Count</text>\n'
    svg += f'<text x="{W/2:.1f}" y="{MARGIN_T-30}" text-anchor="middle" font-size="15" font-weight="bold">{title}</text>\n'

    # Stats box
    stats = f'N={len(data)}  Mean={m:.1f}%  Median={med:.1f}%  P25={percentile(data,25):.0f}%  P75={percentile(data,75):.0f}%'
    svg += f'<text x="{W/2:.1f}" y="{MARGIN_T-10}" text-anchor="middle" font-size="11" fill="#555">{stats}</text>\n'

    return svg, max_count

def cdf_svg(data_sets, title):
    svg = ''
    for data, label, color in data_sets:
        s = sorted(data)
        n = len(s)
        points = []
        for i in range(0, n, max(1, n//200)):
            x = s[i]
            y = (i+1) / n * 100
            points.append((x, y))
        # Ensure endpoints
        if points[0][0] > 0:
            points.insert(0, (0, 0))
        points.append((100, 100))

        d = f'M{MARGIN_L + points[0][0]/100*PLOT_W:.1f},{MARGIN_T + PLOT_H - points[0][1]/100*PLOT_H:.1f}'
        for x, y in points[1:]:
            d += f' L{MARGIN_L + x/100*PLOT_W:.1f},{MARGIN_T + PLOT_H - y/100*PLOT_H:.1f}'

        svg += f'<path d="{d}" fill="none" stroke="{color}" stroke-width="2.5"/>\n'
        svg += f'<text x="{MARGIN_L+80}" y="{MARGIN_T+30+40*data_sets.index((data,label,color))}" fill="{color}" font-size="12">{label} (N={len(data):,})</text>\n'

    # Reference lines
    for pct in [25, 50, 75]:
        y = MARGIN_T + PLOT_H - pct/100 * PLOT_H
        svg += f'<line x1="{MARGIN_L}" y1="{y:.1f}" x2="{MARGIN_L+PLOT_W}" y2="{y:.1f}" stroke="#ddd" stroke-width="0.5" stroke-dasharray="4,4"/>\n'
        svg += f'<text x="{MARGIN_L-5}" y="{y+4:.1f}" text-anchor="end" font-size="9" fill="#999">{pct}%</text>\n'

    # Axes
    svg += f'<line x1="{MARGIN_L}" y1="{MARGIN_T+PLOT_H}" x2="{MARGIN_L+PLOT_W}" y2="{MARGIN_T+PLOT_H}" stroke="black" stroke-width="1"/>\n'
    svg += f'<line x1="{MARGIN_L}" y1="{MARGIN_T}" x2="{MARGIN_L}" y2="{MARGIN_T+PLOT_H}" stroke="black" stroke-width="1"/>\n'
    for pct in [0, 25, 50, 75, 100]:
        x = MARGIN_L + pct/100 * PLOT_W
        svg += f'<text x="{x:.1f}" y="{MARGIN_T+PLOT_H+20}" text-anchor="middle" font-size="11">{pct}%</text>\n'
    svg += f'<text x="{W/2:.1f}" y="{MARGIN_T+PLOT_H+45}" text-anchor="middle" font-size="13">Recall (%)</text>\n'
    svg += f'<text x="{MARGIN_L-55:.1f}" y="{MARGIN_T+PLOT_H/2:.1f}" text-anchor="middle" font-size="13" transform="rotate(-90,{MARGIN_L-55},{MARGIN_T+PLOT_H/2})">Cumulative % of Tasks</text>\n'
    svg += f'<text x="{W/2:.1f}" y="{MARGIN_T-30}" text-anchor="middle" font-size="15" font-weight="bold">{title}</text>\n'
    return svg

def bar_chart_svg(groups_data, title):
    """groups_data = [(label, [values], color), ...]"""
    n_groups = len(groups_data[0][1])  # number of bar groups
    n_bars = len(groups_data)  # bars per group
    bar_w = PLOT_W / (n_groups * n_bars + n_groups) * 0.7
    gap = bar_w * 0.4

    max_val = max(max(vals) for _, vals, _ in groups_data) * 1.15

    svg = ''
    for gi in range(n_groups):
        for bi, (label, vals, color) in enumerate(groups_data):
            v = vals[gi]
            h = v / max_val * PLOT_H
            x = MARGIN_L + gi * (PLOT_W / n_groups) + bi * (bar_w + gap)
            y = MARGIN_T + PLOT_H - h
            svg += f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{color}" opacity="0.85"/>\n'
            if v > 0:
                svg += f'<text x="{x+bar_w/2:.1f}" y="{y-3:.1f}" text-anchor="middle" font-size="9" font-weight="bold">{v:.1f}%</text>\n'

    # Axes
    svg += f'<line x1="{MARGIN_L}" y1="{MARGIN_T+PLOT_H}" x2="{MARGIN_L+PLOT_W}" y2="{MARGIN_T+PLOT_H}" stroke="black" stroke-width="1"/>\n'
    svg += f'<line x1="{MARGIN_L}" y1="{MARGIN_T}" x2="{MARGIN_L}" y2="{MARGIN_T+PLOT_H}" stroke="black" stroke-width="1"/>\n'

    # Legend
    for bi, (label, vals, color) in enumerate(groups_data):
        svg += f'<rect x="{MARGIN_L+10}" y="{MARGIN_T-40+bi*18}" width="12" height="12" fill="{color}"/>\n'
        svg += f'<text x="{MARGIN_L+26}" y="{MARGIN_T-30+bi*18}" font-size="10">{label}</text>\n'

    svg += f'<text x="{W/2:.1f}" y="{MARGIN_T+PLOT_H+45}" text-anchor="middle" font-size="13">Experimental Residue Count (exp_n)</text>\n'
    svg += f'<text x="{MARGIN_L-55:.1f}" y="{MARGIN_T+PLOT_H/2:.1f}" text-anchor="middle" font-size="13" transform="rotate(-90,{MARGIN_L-55},{MARGIN_T+PLOT_H/2})">%</text>\n'
    svg += f'<text x="{W/2:.1f}" y="{MARGIN_T-30}" text-anchor="middle" font-size="15" font-weight="bold">{title}</text>\n'

    # X labels
    x_labels = ['1-5', '5-10', '10-20', '20-50', '50+']
    for gi in range(n_groups):
        x = MARGIN_L + gi * (PLOT_W / n_groups) + (n_bars-1) * (bar_w + gap) / 2
        svg += f'<text x="{x:.1f}" y="{MARGIN_T+PLOT_H+18}" text-anchor="middle" font-size="10">{x_labels[gi]}</text>\n'

    return svg

# Build composite SVG: 2 rows x 2 cols
SW, SH = 500, 380  # per panel
COLS, ROWS = 2, 2
TW = SW * COLS
TH = SH * ROWS + 30

def make_panel(svg_content, col, row):
    x = col * SW
    y = row * SH + 30
    return f'<g transform="translate({x},{y})">\n{svg_content}\n</g>\n'

# We need to adjust the global W/H for each panel
W_orig, H_orig = W, H
PLOT_W_orig, PLOT_H_orig = PLOT_W, PLOT_H
MARGIN_L_orig = MARGIN_L

# Panel settings
W, H = SW, SH
MARGIN_L, MARGIN_R, MARGIN_T, MARGIN_B = 65, 20, 55, 60
PLOT_W = W - MARGIN_L - MARGIN_R
PLOT_H = H - MARGIN_T - MARGIN_B

# Panel 1: G1 histogram
svg1, _ = histogram_svg(direct, 'G1: direct mapping', '#2196F3')

# Panel 2: G2 histogram
svg2, _ = histogram_svg(dbref, 'G2: dbref_segment', '#FF9800')

# Panel 3: CDF overlay (need to adjust for smaller panel)
# Redo CDF with new dimensions
W_use, H_use = W, H
MARGIN_L_use, MARGIN_R_use = 65, 20
MARGIN_T_use, MARGIN_B_use = 55, 60
PLOT_W_use = W_use - MARGIN_L_use - MARGIN_R_use
PLOT_H_use = H_use - MARGIN_T_use - MARGIN_B_use

def make_cdf_panel(data_sets, title):
    svg = ''
    for di, (data, label, color) in enumerate(data_sets):
        s = sorted(data)
        n = len(s)
        points = [(0, 0)]
        for i in range(max(1, n//150), n, max(1, n//150)):
            x = s[i]
            y = i / n * 100
            if points and x != points[-1][0]:
                points.append((x, y))
        points.append((100, 100))

        d = f'M{MARGIN_L_use + points[0][0]/100*PLOT_W_use:.1f},{MARGIN_T_use + PLOT_H_use - points[0][1]/100*PLOT_H_use:.1f}'
        for x, y in points[1:]:
            d += f' L{MARGIN_L_use + x/100*PLOT_W_use:.1f},{MARGIN_T_use + PLOT_H_use - y/100*PLOT_H_use:.1f}'
        svg += f'<path d="{d}" fill="none" stroke="{color}" stroke-width="2"/>\n'
        svg += f'<text x="{MARGIN_L_use+15}" y="{MARGIN_T_use+15+di*16}" fill="{color}" font-size="10">{label}</text>\n'

    for pct in [25, 50, 75]:
        y = MARGIN_T_use + PLOT_H_use - pct/100 * PLOT_H_use
        svg += f'<line x1="{MARGIN_L_use}" y1="{y:.1f}" x2="{MARGIN_L_use+PLOT_W_use}" y2="{y:.1f}" stroke="#eee" stroke-width="0.5"/>\n'
    svg += f'<line x1="{MARGIN_L_use}" y1="{MARGIN_T_use+PLOT_H_use}" x2="{MARGIN_L_use+PLOT_W_use}" y2="{MARGIN_T_use+PLOT_H_use}" stroke="black" stroke-width="1"/>\n'
    svg += f'<line x1="{MARGIN_L_use}" y1="{MARGIN_T_use}" x2="{MARGIN_L_use}" y2="{MARGIN_T_use+PLOT_H_use}" stroke="black" stroke-width="1"/>\n'
    for pct in [0, 25, 50, 75, 100]:
        x = MARGIN_L_use + pct/100 * PLOT_W_use
        svg += f'<text x="{x:.1f}" y="{MARGIN_T_use+PLOT_H_use+15}" text-anchor="middle" font-size="9">{pct}%</text>\n'
    svg += f'<text x="{W_use/2:.1f}" y="{MARGIN_T_use+PLOT_H_use+35}" text-anchor="middle" font-size="11">Recall (%)</text>\n'
    svg += f'<text x="15" y="{MARGIN_T_use+PLOT_H_use/2:.1f}" text-anchor="middle" font-size="11" transform="rotate(-90,15,{MARGIN_T_use+PLOT_H_use/2})">Cumulative %</text>\n'
    svg += f'<text x="{W_use/2:.1f}" y="{MARGIN_T_use-25}" text-anchor="middle" font-size="13" font-weight="bold">{title}</text>\n'
    return svg

svg3 = make_cdf_panel([(direct, 'G1: direct', '#2196F3'), (dbref, 'G2: dbref_segment', '#FF9800')], 'Cumulative Distribution (CDF)')

# Panel 4: Recall by exp_n
def make_bar_panel():
    exp_bins = [(1,5), (5,10), (10,20), (20,50), (50,150)]
    x_labels = ['1-5', '5-10', '10-20', '20-50', '50+']

    groups = []
    for label_key in ['G1 recall', 'G1 hit%', 'G2 recall']:
        groups.append([])

    for lo, hi in exp_bins:
        g1s = [r['recall']*100 for r in data if r['mapping']=='direct' and lo <= r['exp_n'] < hi]
        g2s = [r['recall']*100 for r in data if r['mapping']=='dbref_segment' and lo <= r['exp_n'] < hi]
        groups[0].append(mean(g1s))
        groups[1].append(sum(1 for v in g1s if v > 0)/len(g1s)*100 if g1s else 0)
        groups[2].append(mean(g2s))

    colors = ['#2196F3', '#64B5F6', '#FF9800']
    labels = ['G1: mean recall', 'G1: hit rate', 'G2: mean recall']

    n_groups = len(x_labels)
    n_bars = 3
    bar_w = PLOT_W_use / (n_groups * n_bars + n_groups) * 0.65
    gap = bar_w * 0.5
    max_val = max(max(g) for g in groups) * 1.2

    svg = ''
    for gi in range(n_groups):
        for bi in range(n_bars):
            v = groups[bi][gi]
            h = v / max_val * PLOT_H_use
            x = MARGIN_L_use + gi * (PLOT_W_use / n_groups) + bi * (bar_w + gap)
            y = MARGIN_T_use + PLOT_H_use - h
            svg += f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{colors[bi]}" opacity="0.85"/>\n'
            if v > 0:
                svg += f'<text x="{x+bar_w/2:.1f}" y="{y-2:.1f}" text-anchor="middle" font-size="8" font-weight="bold">{v:.1f}</text>\n'

    for bi in range(n_bars):
        svg += f'<rect x="{MARGIN_L_use+5}" y="{MARGIN_T_use-40+bi*14}" width="10" height="10" fill="{colors[bi]}"/>\n'
        svg += f'<text x="{MARGIN_L_use+18}" y="{MARGIN_T_use-31+bi*14}" font-size="9">{labels[bi]}</text>\n'

    svg += f'<line x1="{MARGIN_L_use}" y1="{MARGIN_T_use+PLOT_H_use}" x2="{MARGIN_L_use+PLOT_W_use}" y2="{MARGIN_T_use+PLOT_H_use}" stroke="black" stroke-width="1"/>\n'
    svg += f'<line x1="{MARGIN_L_use}" y1="{MARGIN_T_use}" x2="{MARGIN_L_use}" y2="{MARGIN_T_use+PLOT_H_use}" stroke="black" stroke-width="1"/>\n'
    for gi in range(n_groups):
        x = MARGIN_L_use + gi * (PLOT_W_use / n_groups) + (n_bars-1) * (bar_w + gap) / 2
        svg += f'<text x="{x:.1f}" y="{MARGIN_T_use+PLOT_H_use+15}" text-anchor="middle" font-size="10">{x_labels[gi]}</text>\n'
    svg += f'<text x="{W_use/2:.1f}" y="{MARGIN_T_use+PLOT_H_use+35}" text-anchor="middle" font-size="11">Experimental Residue Count (exp_n)</text>\n'
    svg += f'<text x="15" y="{MARGIN_T_use+PLOT_H_use/2:.1f}" text-anchor="middle" font-size="11" transform="rotate(-90,15,{MARGIN_T_use+PLOT_H_use/2})">%</text>\n'
    svg += f'<text x="{W_use/2:.1f}" y="{MARGIN_T_use-25}" text-anchor="middle" font-size="13" font-weight="bold">Recall & Hit Rate by Exp Residue Count</text>\n'
    return svg

svg4 = make_bar_panel()

# Assemble full page
full_svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{TW}" height="{TH}" viewBox="0 0 {TW} {TH}">
<rect width="{TW}" height="{TH}" fill="white"/>
<text x="{TW/2}" y="22" text-anchor="middle" font-size="16" font-weight="bold">Binding Residue Recall: Docked vs Experimental (N=7,839)</text>
{make_panel(svg1, 0, 0)}
{make_panel(svg2, 1, 0)}
{make_panel(svg3, 0, 1)}
{make_panel(svg4, 1, 1)}
</svg>'''

with open('C:/Users/Administrator/recall_plot.svg', 'w', encoding='utf-8') as f:
    f.write(full_svg)

with open('D:/finale/negative_upgrade/docking_package/recall_plot.svg', 'w', encoding='utf-8') as f:
    f.write(full_svg)

print(f'Done! SVG size: {len(full_svg):,} bytes')
print(f'Saved to recall_plot.svg')
