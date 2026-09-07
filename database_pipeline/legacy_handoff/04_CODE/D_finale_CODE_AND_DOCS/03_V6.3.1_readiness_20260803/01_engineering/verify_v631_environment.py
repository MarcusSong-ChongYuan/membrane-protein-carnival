import importlib
import json
import platform
import sys

mods = ["pandas", "numpy", "pyarrow", "rdkit", "openpyxl", "sklearn", "scipy"]
result = {"python": sys.version, "platform": platform.platform(), "modules": {}, "status": "PASS"}
for name in mods:
    try:
        m = importlib.import_module(name)
        result["modules"][name] = getattr(m, "__version__", "installed")
    except Exception as exc:
        result["modules"][name] = f"ERROR: {exc}"
        result["status"] = "FAIL"
print(json.dumps(result, indent=2, ensure_ascii=False))
raise SystemExit(0 if result["status"] == "PASS" else 1)
