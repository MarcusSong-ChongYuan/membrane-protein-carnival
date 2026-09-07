import csv,json,hashlib,shutil
from pathlib import Path
from PIL import Image,ImageOps,ImageDraw

ROOT=Path(r'D:\finale\21_MemPro_V7.1.1_figures_20260814');FIG=ROOT/'03_figures';QA=ROOT/'05_QA';REP=ROOT/'06_reproducibility';QA.mkdir(exist_ok=True);REP.mkdir(exist_ok=True)
names=['M1_V711_architecture_scale_sources','M2_V711_membrane_proteome_five_axis','M3_V711_expression_localization_corrected','M4_V711_disease_ontology_multilabel','M5_V711_compound_upset_structure_properties','M6_V711_binding_sites_alluvial_docking','M7_V711_lineage_negative_QA','M8_V711_release_docking_readiness']
checks=[]
for n in names:
 for ext in ['png','svg','pdf']:
  p=FIG/ext/f'{n}.{ext}';checks.append({'figure':n,'format':ext,'exists':int(p.exists()),'bytes':p.stat().st_size if p.exists() else 0})
# Montage for rapid review.
thumbs=[]
for n in names:
 im=Image.open(FIG/'png'/f'{n}.png').convert('RGB');im.thumbnail((760,500),Image.Resampling.LANCZOS);canvas=Image.new('RGB',(780,540),'white');canvas.paste(im,((780-im.width)//2,30));ImageDraw.Draw(canvas).text((12,8),n,fill='#30343b');thumbs.append(canvas)
mont=Image.new('RGB',(1560,2160),'white')
for i,im in enumerate(thumbs):mont.paste(im,((i%2)*780,(i//2)*540))
mont.save(QA/'V711_M1_M8_contact_sheet.png',quality=92)
# Reproducibility scripts.
for name in ['build_v711_figure_stats.py','draw_v711_m1_m8.py','finalize_v711_figures.py']:
 p=Path(r'C:\Users\Administrator')/name
 if p.exists():shutil.copy2(p,REP/name)
with (QA/'FIGURE_FILE_QA.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=checks[0].keys(),delimiter='\t');w.writeheader();w.writerows(checks)
report={'release':'MemPro V7.1.1 figures','figure_count':8,'formats_per_figure':3,'all_expected_files_exist':all(x['exists'] for x in checks),'zero_byte_files':sum(x['bytes']==0 for x in checks),'statistics_qa':json.loads((QA/'V711_FIGURE_STATS_QA.json').read_text(encoding='utf-8')),'status':'PASS' if all(x['exists'] and x['bytes']>0 for x in checks) else 'FAIL'}
(QA/'V711_FIGURE_FINAL_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
(ROOT/'README.md').write_text('''# MemPro V7.1.1 figure suite\n\nEight manuscript-oriented composite figures were rebuilt from the frozen V7.1.1 release. Each figure is available as 320-dpi PNG, editable SVG and PDF. Statistics are frozen in `02_statistics`, captions in `04_captions`, source/statistical-unit provenance in `PANEL_SOURCE_VERSION_UNIT_V711.tsv`, and executable scripts in `06_reproducibility`.\n\nVisual policy: Arial; restrained Seaborn Paired/rocket palettes; transparent overlays; directly labelled statistical units; source footer on every figure. M1 uses architecture and scale; M2 composition/lollipop; M3 stacked composition and raincloud-style distributions; M4 ontology heatmap; M5 UpSet and chemical distributions; M6 alluvial and funnel; M7 lineage/QA; M8 release and docking-readiness funnel.\n''',encoding='utf-8')
def sha(p):
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
records=[]
for p in sorted(ROOT.rglob('*')):
 if p.is_file() and p.name!='MANIFEST.tsv':records.append({'relative_path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
with (ROOT/'MANIFEST.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['relative_path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(records)
print(json.dumps({'status':report['status'],'files':len(records),'contact_sheet':str(QA/'V711_M1_M8_contact_sheet.png')},ensure_ascii=False))
