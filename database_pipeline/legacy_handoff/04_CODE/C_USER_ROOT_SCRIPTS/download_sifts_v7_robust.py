import csv,gzip,hashlib,json,time,urllib.request,urllib.error
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT=Path(r"D:\finale\12_V7_final_completion_20260813\03_binding_sites")
IDS=[x.strip().lower() for x in (ROOT/'pdb_ids_requiring_sifts_v7.txt').read_text().splitlines() if x.strip()]
CP=ROOT/'sifts_residue_mapping_v7.checkpoint.tsv';STATUS=ROOT/'sifts_download_status_v7.tsv';SUMMARY=ROOT/'SIFTS_DOWNLOAD_SUMMARY.json'
header=['pdb_id','chain','struct_asym_id','uniprot_acc','uniprot_resnum','pdb_resnum','pdb_icode','pdb_resname','observed','mapping_status','xml_sha256']
if not CP.exists():
 with CP.open('w',encoding='utf-8',newline='') as f:csv.writer(f,delimiter='\t').writerow(header)
if not STATUS.exists():
 with STATUS.open('w',encoding='utf-8',newline='') as f:csv.writer(f,delimiter='\t').writerow(['pdb_id','status','mapping_rows','xml_sha256','attempts','message'])
done=set()
with STATUS.open(encoding='utf-8') as f:
 for r in csv.DictReader(f,delimiter='\t'):done.add(r['pdb_id'].lower())
pending=[x for x in IDS if x not in done]

def lname(tag):return tag.split('}')[-1]
def fetch(pid):
 url=f'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{pid}.xml.gz';last=''
 for attempt,delay in enumerate([0,1,5],1):
  if delay:time.sleep(delay)
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'MemPro-V7-SIFTS/1.0'})
   with urllib.request.urlopen(req,timeout=60) as resp:b=resp.read()
   raw=gzip.decompress(b) if b[:2]==b'\x1f\x8b' else b;sha=hashlib.sha256(raw).hexdigest();root=ET.fromstring(raw);rows=[]
   for entity in root.iter():
    if lname(entity.tag)!='entity':continue
    asym=entity.get('entityId','')
    for residue in entity.iter():
     if lname(residue.tag)!='residue':continue
     pdb={};uni={};details=[]
     for child in residue.iter():
      nm=lname(child.tag)
      if nm=='crossRefDb':
       src=child.get('dbSource','').lower()
       if src in {'pdb','pdbe'} and not pdb:pdb=dict(child.attrib)
       elif src=='uniprot' and not uni:uni=dict(child.attrib)
      elif nm=='residueDetail':details.append((child.text or '')+' '+child.get('property',''))
     if not uni:continue
     chain=pdb.get('dbChainId') or residue.get('dbChainId','');rn=pdb.get('dbResNum') or residue.get('dbResNum','');icode=''
     # PDB residue numbers may contain an insertion-code suffix.
     if rn and rn[-1:].isalpha():icode=rn[-1];rn=rn[:-1]
     detail=' '.join(details).lower();obs='unobserved' if 'not_observed' in detail or 'unobserved' in detail else 'observed'
     rows.append([pid.upper(),chain,asym,uni.get('dbAccessionId',''),uni.get('dbResNum',''),rn,icode,pdb.get('dbResName') or residue.get('dbResName',''),obs,'mapped' if uni.get('dbResNum') else 'no_uniprot_number',sha])
   return pid,'ok' if rows else 'no_uniprot_mapping',rows,sha,attempt,''
  except urllib.error.HTTPError as e:
   last=f'HTTP {e.code}'
   if e.code==404:return pid,'not_found',[],'',attempt,last
  except Exception as e:last=repr(e)[:300]
 return pid,'failed',[],'',3,last

results=[]
with ThreadPoolExecutor(max_workers=6) as ex:
 futs={ex.submit(fetch,p):p for p in pending}
 for i,fut in enumerate(as_completed(futs),1):
  result=fut.result();results.append(result)
  if i%50==0:print(f'{i}/{len(pending)}',flush=True)
# Single atomic append phase avoids interleaved checkpoint writes.
with CP.open('a',encoding='utf-8',newline='') as fm,STATUS.open('a',encoding='utf-8',newline='') as fs:
 wm=csv.writer(fm,delimiter='\t');ws=csv.writer(fs,delimiter='\t')
 for pid,st,rows,sha,attempt,msg in sorted(results):
  wm.writerows(rows);ws.writerow([pid.upper(),st,len(rows),sha,attempt,msg])
status=list(csv.DictReader(STATUS.open(encoding='utf-8'),delimiter='\t'));counts={}
for r in status:counts[r['status']]=counts.get(r['status'],0)+1
report={'requested_pdb':len(IDS),'completed_status_rows':len(status),'pending':len(IDS)-len({r['pdb_id'].lower() for r in status}),'status_counts':counts,'mapping_rows':sum(int(r['mapping_rows']) for r in status),'source':'SIFTS XML snapshot fetched from EBI','parser':'namespace-safe crossRefDb PDB/UniProt parser v1.0'}
SUMMARY.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
