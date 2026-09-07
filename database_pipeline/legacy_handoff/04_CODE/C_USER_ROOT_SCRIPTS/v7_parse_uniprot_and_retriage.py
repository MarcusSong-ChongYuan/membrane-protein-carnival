#!/usr/bin/env python3
import csv,json
from collections import Counter
from pathlib import Path

ROOT=Path(r"D:\finale\09_V7_data_freeze_working_20260813")
INIT=ROOT/"02_protein_audit"/"protein_identity_audit_all_v7_initial.tsv"
RAW=ROOT/"02_protein_audit"/"uniprot_current_jsonl"
OUT=ROOT/"02_protein_audit"/"protein_identity_audit_all_v7_uniprot_parsed.tsv"

def feature_location(x):
    loc=x.get('location',{})
    def pos(k):
        v=loc.get(k,{})
        return str(v.get('value',''))+((':'+v.get('modifier','')) if v.get('modifier') else '')
    return pos('start')+'..'+pos('end')

def ev_codes(x):
    vals=[]
    for e in x.get('evidences',[]) or []:
        vals.append('|'.join(str(e.get(k,'')) for k in ('evidenceCode','source','id') if e.get(k)))
    return ';'.join(vals)

def locations(entry):
    out=[]
    for c in entry.get('comments',[]):
        if c.get('commentType')!='SUBCELLULAR LOCATION':continue
        for s in c.get('subcellularLocations',[]):
            bits=[]
            for k in ('location','topology','orientation'):
                v=(s.get(k) or {}).get('value','')
                if v:bits.append(v)
            ev=[]
            for e in s.get('evidences',[]) or []:
                ev.append('|'.join(str(e.get(k,'')) for k in ('evidenceCode','source','id') if e.get(k)))
            out.append(' | '.join(bits)+((' ['+';'.join(ev)+']') if ev else ''))
    return '; '.join(out)

entries={}
for p in sorted(RAW.glob('batch_*.json')):
    d=json.loads(p.read_text(encoding='utf-8'))
    for e in d.get('results',[]):entries[e['primaryAccession']]=e

with INIT.open(encoding='utf-8',newline='') as f:
    rd=csv.DictReader(f,delimiter='\t'); base=rd.fieldnames; rows=list(rd)
extra=['uniprot_entry_type_current','uniprot_sequence_length_current','uniprot_transmembrane_count_current','uniprot_intramembrane_count_current','uniprot_lipidation_count_current','uniprot_signal_count_current','uniprot_chain_count_current','uniprot_peptide_count_current','uniprot_transmembrane_features_current','uniprot_intramembrane_features_current','uniprot_lipidation_features_current','uniprot_signal_features_current','uniprot_subcellular_locations_current','uniprot_isoform_ids_current','proposed_primary_class_v7','proposed_primary_mode_v7','proposed_secondary_modes_v7','proposed_evidence_level_v7','retriage_priority_v7','retriage_flags_v7','decision_status_v7']
out=[]
for r in rows:
    e=entries.get(r['target_uniprot'])
    if not e:
        x=dict(r);x.update({k:'' for k in extra});x.update({'retriage_priority_v7':'P0','retriage_flags_v7':'UNIPROT_ENTRY_MISSING','decision_status_v7':'REVIEW_REQUIRED'});out.append(x);continue
    fs=e.get('features',[])
    by={t:[x for x in fs if x.get('type')==t] for t in ['Transmembrane','Intramembrane','Lipidation','Signal','Chain','Peptide']}
    loc=locations(e); low=loc.lower()
    modes=[]
    if by['Transmembrane']:modes.append('transmembrane_segment')
    if by['Intramembrane']:modes.append('intramembrane_segment')
    if by['Lipidation']:modes.append('covalent_lipid_anchor')
    periph='peripheral membrane protein' in low or 'peripheral membrane' in low
    if periph:modes.append('peripheral_membrane_association')
    if by['Transmembrane'] or by['Intramembrane']:
        cls='A';primary='transmembrane_or_intramembrane';secondary=[m for m in modes if m not in {'transmembrane_segment','intramembrane_segment'}]
    elif by['Lipidation']:
        cls='B';primary='covalent_lipid_anchor';secondary=['peripheral_membrane_association'] if periph else []
    elif periph:
        cls='C';primary='peripheral_membrane_association';secondary=[]
    elif 'membrane' in low or 'cell surface' in low:
        cls='UNRESOLVED';primary='location_without_binding_mode';secondary=[]
    else:
        cls='EXCLUDE_CANDIDATE';primary='no_current_membrane_binding_mode';secondary=[]
    all_ev=' '.join(ev_codes(x) for x in fs)+' '+loc
    direct=('ECO:0000269' in all_ev or 'ECO:0000305' in all_ev)
    curated=bool(modes) and e.get('entryType','').lower().startswith('uniprotkb reviewed')
    if direct and modes:lev='E1'
    elif curated:lev='E2'
    else:lev='E3'
    flags=[]
    if cls in {'UNRESOLVED','EXCLUDE_CANDIDATE'}:flags.append('NO_EXPLICIT_MEMBRANE_BINDING_MODE_CURRENT_UNIPROT')
    if r['current_class']!=cls and cls in {'A','B','C'}:flags.append('ABC_RECLASSIFICATION_PROPOSED')
    if r['current_class'] in {'A','B','C'} and cls=='EXCLUDE_CANDIDATE':flags.append('CURRENT_INCLUDED_BUT_NO_MODE')
    if r['current_class']=='A' and cls in {'B','C'}:flags.append('A_DOWNCLASS_TO_CORRECT_MECHANISM')
    if r['cross_source_conflict']=='1':flags.append('CROSS_SOURCE_CONFLICT')
    if cls=='EXCLUDE_CANDIDATE' or 'CURRENT_INCLUDED_BUT_NO_MODE' in flags:prio='P0'
    elif cls=='UNRESOLVED' or 'CROSS_SOURCE_CONFLICT' in flags:prio='P1'
    elif 'ABC_RECLASSIFICATION_PROPOSED' in flags or lev=='E2':prio='P2'
    elif lev=='E3':prio='P3'
    else:prio='P4_CONTROL'
    def desc(t):return ';'.join(f"{feature_location(x)}|{x.get('description','')}|{ev_codes(x)}" for x in by[t])
    iso=[]
    for c in e.get('comments',[]):
        if c.get('commentType')=='ALTERNATIVE PRODUCTS':
            for isoform in c.get('isoforms',[]):iso.extend(isoform.get('isoformIds',[]) or [])
    x=dict(r);x.update({'uniprot_entry_type_current':e.get('entryType',''),'uniprot_sequence_length_current':len((e.get('sequence') or {}).get('value','')),'uniprot_transmembrane_count_current':len(by['Transmembrane']),'uniprot_intramembrane_count_current':len(by['Intramembrane']),'uniprot_lipidation_count_current':len(by['Lipidation']),'uniprot_signal_count_current':len(by['Signal']),'uniprot_chain_count_current':len(by['Chain']),'uniprot_peptide_count_current':len(by['Peptide']),'uniprot_transmembrane_features_current':desc('Transmembrane'),'uniprot_intramembrane_features_current':desc('Intramembrane'),'uniprot_lipidation_features_current':desc('Lipidation'),'uniprot_signal_features_current':desc('Signal'),'uniprot_subcellular_locations_current':loc,'uniprot_isoform_ids_current':';'.join(sorted(set(iso))),'proposed_primary_class_v7':cls,'proposed_primary_mode_v7':primary,'proposed_secondary_modes_v7':';'.join(secondary),'proposed_evidence_level_v7':lev,'retriage_priority_v7':prio,'retriage_flags_v7':';'.join(flags),'decision_status_v7':'AUTO_PROPOSAL_PENDING_CROSS_SOURCE_REVIEW'})
    out.append(x)
with OUT.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=base+extra,delimiter='\t');w.writeheader();w.writerows(out)
summary={'input':len(rows),'uniprot_entries_loaded':len(entries),'proposed_class':dict(Counter(x['proposed_primary_class_v7'] for x in out)),'proposed_evidence':dict(Counter(x['proposed_evidence_level_v7'] for x in out)),'priority':dict(Counter(x['retriage_priority_v7'] for x in out)),'flags':dict(Counter(y for x in out for y in x['retriage_flags_v7'].split(';') if y)),'status':'UNIPROT_PARSED_CROSS_SOURCE_REVIEW_PENDING'}
(ROOT/'02_protein_audit'/'UNIPROT_PARSED_RETRIAGE_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
