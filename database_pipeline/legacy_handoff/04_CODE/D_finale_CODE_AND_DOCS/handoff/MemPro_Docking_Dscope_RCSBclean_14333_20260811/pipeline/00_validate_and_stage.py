#!/usr/bin/env python3
from pathlib import Path
import csv,hashlib,json,sys
r=Path(__file__).resolve().parents[1];p=r/'docking_manifest.tsv';x=list(csv.DictReader(open(p,encoding="utf-8"),delimiter="\t"));k=[(a['compound_internal_id'],a['target_uniprot']) for a in x];bad=[]
if len(x)!=14333:bad.append(f"rows={len(x)}")
if len(set(k))!=len(k):bad.append(f"duplicate_pairs={len(k)-len(set(k))}")
print(json.dumps({'rows':len(x),'unique_pairs':len(set(k)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'problems':bad},indent=2));sys.exit(2 if bad else 0)
