import csv,gzip,json,re,os
from collections import Counter
from pathlib import Path

BASE=Path(r'D:\finale\16_MemPro_V7.0.2_final_20260814')
ROOT=Path(r'D:\finale\17_MemPro_V7.1_candidate_20260814')
TABLE=ROOT/'01_release_tables'/'interaction_site_v71.tsv.gz'

def readgz(p):
    with gzip.open(p,'rt',encoding='utf-8-sig',newline='') as f:yield from csv.DictReader(f,delimiter='\t')
def residues(s):
    if not s:return []
    out=[]
    for x in re.findall(r'(?i)(?:\b[A-Z]{3}\s*\d+[A-Z]?\b|\b[A-Z]\d+[A-Z]?\b|(?<![A-Za-z])\d+[A-Z]?(?![A-Za-z]))',s):
        x=re.sub(r'\s+','',x.upper())
        if x not in out:out.append(x)
    return out
def raw_index():
    d={}
    for p in [BASE/'01_data'/'binding_site_coordinate_verified_v702.tsv.gz',BASE/'02_nonpublic_resolved_status'/'binding_site_coordinate_nonpublic_v702.tsv.gz']:
        for r in readgz(p):
            d[r['binding_site_instance_id']]={
              'source_res':r.get('requested_uniprot_residue_tokens_v7') or r.get('residue_or_site_description',''),
              'pocket_status':r.get('pocket_rescue_status_v704b',''),
              'pocket_contacts':r.get('pocket_contact_uniprot_residues_v706') or r.get('pocket_contact_residues_v704','')}
    return d
def main():
    raw=raw_index();tmp=TABLE.with_name(TABLE.name+'.tmp');cs=Counter();lig=Counter();n=0
    with gzip.open(TABLE,'rt',encoding='utf-8',newline='') as fi,gzip.open(tmp,'wt',encoding='utf-8',newline='') as fo:
        rd=csv.DictReader(fi,delimiter='\t');w=csv.DictWriter(fo,fieldnames=rd.fieldnames,delimiter='\t');w.writeheader()
        for r in rd:
            n+=1;x=raw[r['site_record_id']];src=residues(x['source_res']);contacts=residues(x['pocket_contacts']);exact=x['pocket_status']=='PUBLIC_RESCUED_EXACT_COCRYSTAL_POCKET'
            r['source_reported_residues']=';'.join(src)
            if exact:
                r['site_residue_status']='RESIDUE_COMPLETE';r['structure_ligand_relationship']='EXACT_COCRYSTAL_LIGAND';r['structure_ligand_basis']='exact_HET_on_SIFTS_target_chain_with_protein_contacts'
                if contacts:r['mapped_structure_residues']=';'.join(contacts)
            elif src:
                m=int(r.get('mapped_residue_count') or 0);q=int(r.get('requested_residue_count') or len(src));r['requested_residue_count']=q or len(src)
                r['site_residue_status']='RESIDUE_COMPLETE' if m and m>=int(r['requested_residue_count']) else ('RESIDUE_PARTIAL' if m else 'RESIDUE_SOURCE_ONLY')
                r['structure_ligand_relationship']='LIGAND_IDENTITY_UNRESOLVED' if r.get('source_pdb_ids') else 'NOT_APPLICABLE'
            else:
                r['site_residue_status']='NO_RESIDUE_REPORTED';r['requested_residue_count']='0';r['mapped_residue_count']=str(len(contacts) if exact else 0);r['observed_residue_count']=r['mapped_residue_count'];r['mapping_fraction']='1.000000' if exact else '0.000000'
                if not exact:r['mapped_structure_residues']='';r['coordinate_docking_eligible']='0';r['coordinate_mapping_status']='SOURCE_STRUCTURE_ONLY' if r.get('source_pdb_ids') else 'NO_COORDINATE';r['coordinate_confidence']='NONE';r['structure_ligand_relationship']='STRUCTURE_ONLY_NO_QUERY_LIGAND' if r.get('source_pdb_ids') else 'NOT_APPLICABLE';r['structure_ligand_basis']='PDB_association_does_not_by_itself_prove_query_ligand_is_present'
            lim=[z for z in r.get('limitations','').split(';') if z and z!='source_did_not_report_residues']
            if r['site_residue_status']=='NO_RESIDUE_REPORTED':lim.append('source_did_not_report_residues')
            r['limitations']=';'.join(dict.fromkeys(lim));w.writerow(r);cs[r['site_residue_status']]+=1;lig[r['structure_ligand_relationship']]+=1
    os.replace(tmp,TABLE)
    report={'rows':n,'residue_status':dict(cs),'ligand_relationship':dict(lig),'exact_cocrystal_count':lig['EXACT_COCRYSTAL_LIGAND'],'pass':n==55610 and lig['EXACT_COCRYSTAL_LIGAND']==663}
    (ROOT/'04_QA'/'V71_STAGE2_SITE_CORRECTED_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False))
if __name__=='__main__':main()
