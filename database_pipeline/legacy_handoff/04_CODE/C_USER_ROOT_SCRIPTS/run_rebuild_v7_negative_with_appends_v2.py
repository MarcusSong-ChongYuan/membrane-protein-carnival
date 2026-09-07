from pathlib import Path

p=Path(r"C:\Users\Administrator\rebuild_v7_negative_with_appends.py")
s=p.read_text(encoding="utf-8")
s=s.replace(
"form_rows={};form_conflict=[];fform=None\nfor p in [FORMBASE,FORMAPP]:\n with op(p) as f:\n  rd=csv.DictReader(f,delimiter='\\t');fform=rd.fieldnames",
"form_rows={};form_conflict=[];fform=[]\nfor p in [FORMBASE,FORMAPP]:\n with op(p) as f:\n  rd=csv.DictReader(f,delimiter='\\t');fform += [x for x in rd.fieldnames if x not in fform]"
)
if "fform += [x for x in rd.fieldnames if x not in fform]" not in s:
    raise RuntimeError("schema-union patch did not apply")
exec(compile(s,"rebuild_v7_negative_with_appends_v2.py","exec"))
