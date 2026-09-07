import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(r"D:\finale\12_V7_final_completion_20260813")
SUMMARY = ROOT / "03_binding_sites" / "SIFTS_FROZEN_RESCUE_DOWNLOAD_SUMMARY.json"
LOG = ROOT / "logs" / "v701_automatic_finalize.log"
STATE = ROOT / "logs" / "v701_automatic_finalize_state.json"
PY = r"C:\Users\Administrator\AppData\Local\Programs\Python\Launcher\py.exe"

def write(message):
    stamp = datetime.now(timezone.utc).isoformat()
    with LOG.open("a", encoding="utf-8") as handle: handle.write(f"{stamp}\t{message}\n")

def run(script):
    write(f"START {script}")
    result = subprocess.run([PY, "-3", str(Path(r"C:\Users\Administrator") / script)], text=True, capture_output=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(result.stdout); handle.write(result.stderr)
    if result.returncode:
        raise RuntimeError(f"{script} exited {result.returncode}")
    write(f"DONE {script}")

write("WATCHER_STARTED")
while True:
    if SUMMARY.exists():
        data = json.loads(SUMMARY.read_text(encoding="utf-8"))
        if int(data.get("pending", -1)) == 0:
            break
    time.sleep(60)

try:
    run("rescue_v7_frozen_sites_v701.py")
    state = {"status":"SIFTS_COMPLETE_AND_RESCUE_RECALCULATED","completed_at_utc":datetime.now(timezone.utc).isoformat(),"next":"V7.0.1 integration and QA"}
except Exception as exc:
    state = {"status":"FAILED","error":repr(exc),"failed_at_utc":datetime.now(timezone.utc).isoformat()}
STATE.write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8")
write(state["status"])
