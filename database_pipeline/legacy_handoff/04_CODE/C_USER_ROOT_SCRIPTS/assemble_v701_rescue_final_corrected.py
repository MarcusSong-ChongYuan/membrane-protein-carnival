from pathlib import Path
src=Path(r"C:\Users\Administrator\assemble_v701_rescue_final.py").read_text(encoding='utf-8')
needle="def copy_remaining(src,out):\n shutil.copy2(src,out)\n with op(src) as f:return sum(1 for _ in csv.DictReader(f,delimiter='\\t'))"
replacement=needle+"\ndef filter_remaining(src,out,status):\n with op(src) as f:\n  rd=csv.DictReader(f,delimiter='\\t');fields=rd.fieldnames;rows=[r for r in rd if not r.get(status,'').startswith('PUBLIC')]\n with op(out,'wt') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\\t');w.writeheader();w.writerows(rows)\n return len(rows)"
src=src.replace(needle,replacement)
src=src.replace("iso_rem=copy_remaining(W/'06_remaining_modules'/'protein_isoform_v701.tsv.gz',REV/'isoform_remaining_v701.tsv.gz')","iso_rem=filter_remaining(W/'06_remaining_modules'/'protein_isoform_v701.tsv.gz',REV/'isoform_remaining_v701.tsv.gz','isoform_disposition_v701')")
src=src.replace("pro_rem=copy_remaining(W/'06_remaining_modules'/'protein_processed_form_v701.tsv.gz',REV/'processed_form_remaining_v701.tsv.gz')","pro_rem=filter_remaining(W/'06_remaining_modules'/'protein_processed_form_v701.tsv.gz',REV/'processed_form_remaining_v701.tsv.gz','assignment_status_v701')")
exec(compile(src,'assemble_v701_rescue_final.corrected.py','exec'))
