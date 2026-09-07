from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

WORKFLOW = Path(__file__).resolve().parent
if str(WORKFLOW) not in sys.path:
    sys.path.insert(0, str(WORKFLOW))

from src.export_abce_release_v1_1 import run_export


def main() -> int:
    parser = argparse.ArgumentParser(description="Build A/B/C/E scored archive and Rule B v1.1 release")
    parser.add_argument("--config", required=True)
    parser.add_argument("--source-policy", required=True)
    parser.add_argument("--scoring-policy", required=True)
    args = parser.parse_args()
    print(json.dumps(run_export(args.config, args.source_policy, args.scoring_policy), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
