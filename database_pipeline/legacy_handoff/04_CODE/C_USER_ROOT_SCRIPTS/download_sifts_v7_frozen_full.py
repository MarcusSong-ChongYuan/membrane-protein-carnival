from pathlib import Path

source = Path(r"C:\Users\Administrator\download_sifts_v7_robust.py").read_text(encoding="utf-8")
source = source.replace("ROOT/'pdb_ids_requiring_sifts_v7.txt'", "ROOT/'pdb_ids_for_frozen_rescue.txt'")
source = source.replace("SIFTS_DOWNLOAD_SUMMARY.json", "SIFTS_FROZEN_RESCUE_DOWNLOAD_SUMMARY.json")
source = source.replace("ThreadPoolExecutor(max_workers=6)", "ThreadPoolExecutor(max_workers=12)")
source = source.replace("parser v1.0", "parser v1.0; expanded frozen-site coverage")
exec(compile(source, "download_sifts_v7_frozen_full.generated.py", "exec"))
