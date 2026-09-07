import json

# Load data
with open('C:/Users/Administrator/residue_recall_v2.jsonl') as f:
    data = [json.loads(l) for l in f if l.strip()]

# Read template
with open('C:/Users/Administrator/recall_chart.html', encoding='utf-8') as f:
    html = f.read()

# Inject data
json_str = json.dumps(data)
html = html.replace('__DATA_PLACEHOLDER__', json_str)

# Write final
with open('C:/Users/Administrator/recall_chart_final.html', 'w', encoding='utf-8') as f:
    f.write(html)

print('Done: ' + str(len(html)) + ' bytes')
