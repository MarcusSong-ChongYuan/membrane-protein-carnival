import json,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(r"D:\finale\12_V7_final_completion_20260813");SUMMARY=ROOT/'03_binding_sites'/'SIFTS_FROZEN_RESCUE_DOWNLOAD_SUMMARY.json';LOG=ROOT/'logs'/'v701_automatic_finalize_v2.log';STATE=ROOT/'logs'/'v701_automatic_finalize_state.json';PY=r"C:\Users\Administrator\AppData\Local\Programs\Python\Launcher\py.exe"
def write(x):
 with LOG.open('a',encoding='utf-8') as f:f.write(datetime.now(timezone.utc).isoformat()+'\t'+x+'\n')
def run(x):
 write('START '+x);r=subprocess.run([PY,'-3',str(Path(r"C:\Users\Administrator")/x)],text=True,capture_output=True)
 with LOG.open('a',encoding='utf-8') as f:f.write(r.stdout+r.stderr)
 if r.returncode:raise RuntimeError(f'{x} exited {r.returncode}')
 write('DONE '+x)
write('WATCHER_V2_STARTED')
while True:
 if SUMMARY.exists() and int(json.loads(SUMMARY.read_text(encoding='utf-8')).get('pending',-1))==0:break
 time.sleep(60)
try:
 run('rescue_v7_frozen_sites_v701.py');run('integrate_v701_rescues.py');state={'status':'V701_RESCUE_INTEGRATED_AND_QA_COMPLETE','completed_at_utc':datetime.now(timezone.utc).isoformat()}
except Exception as e:state={'status':'FAILED','error':repr(e),'failed_at_utc':datetime.now(timezone.utc).isoformat()}
STATE.write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8');write(state['status'])
