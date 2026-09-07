"""
Batch Vina-GPU runner. Processes configs sequentially on one A100 GPU.
Resumes: skips existing outputs. Logs progress.
"""
import subprocess, os, time, glob
from pathlib import Path

DOCK_DIR = r'D:\finale\negative_upgrade\docking_package'
CONF_DIR = os.path.join(DOCK_DIR, 'vina_configs')
RESULT_DIR = os.path.join(DOCK_DIR, 'docking_results')
LOG_DIR = os.path.join(DOCK_DIR, 'logs')
os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

VINA_GPU_BIN = os.environ.get('VINA_GPU_BIN', 'vina_gpu')
GPU_ID = int(os.environ.get('GPU_ID', '0'))

configs = sorted(glob.glob(os.path.join(CONF_DIR, '*.conf')))
total = len(configs)
print(f"Found {total:,} config files")
print(f"GPU ID: {GPU_ID}")
print(f"Vina binary: {VINA_GPU_BIN}")

done = 0; skipped = 0; failed = 0
start_time = time.time()
fail_log = open(os.path.join(LOG_DIR, 'failed_tasks.txt'), 'w')

for i, conf in enumerate(configs):
    task_id = Path(conf).stem
    out_file = os.path.join(RESULT_DIR, f'{task_id}_out.pdbqt')

    if os.path.exists(out_file) and os.path.getsize(out_file) > 100:
        skipped += 1
        continue

    print(f"[{i+1}/{total}] {task_id} ... ", end='', flush=True)

    try:
        result = subprocess.run(
            [VINA_GPU_BIN, '--config', conf, '--gpu_id', str(GPU_ID)],
            capture_output=True, text=True, timeout=600,
            cwd=DOCK_DIR
        )
        if result.returncode == 0 and os.path.exists(out_file):
            print("OK")
            done += 1
        else:
            print(f"FAILED (rc={result.returncode})")
            fail_log.write(f"{task_id}\t{result.returncode}\n")
            failed += 1
    except subprocess.TimeoutExpired:
        print("TIMEOUT")
        fail_log.write(f"{task_id}\tTIMEOUT\n")
        failed += 1
    except Exception as e:
        print(f"ERROR: {e}")
        fail_log.write(f"{task_id}\t{e}\n")
        failed += 1

    if (i + 1) % 500 == 0:
        elapsed = time.time() - start_time
        rate = (i + 1) / elapsed * 3600
        eta = (total - i - 1) / rate if rate > 0 else 0
        print(f"  [PROGRESS] {i+1:,}/{total:,} | done={done} fail={failed} skip={skipped}")
        print(f"  Rate: {rate:.0f}/hr | ETA: {eta:.1f}hr")
        fail_log.flush()

elapsed = time.time() - start_time
fail_log.close()
print(f"\n=== FINAL ===")
print(f"Total: {total:,}  Done: {done:,}  Failed: {failed:,}  Skipped: {skipped:,}")
print(f"Time: {elapsed/3600:.1f} hr  |  Rate: {total/elapsed*3600:.0f}/hr")
