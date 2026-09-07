from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

WORKFLOW = Path(__file__).resolve().parent
if str(WORKFLOW) not in sys.path:
    sys.path.insert(0, str(WORKFLOW))

from src.validate_abce_release_v1_1 import validate_release


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the A/B/C/E archive or Rule B-only v1.1 release")
    parser.add_argument("--config", required=True)
    parser.add_argument("--scope", choices=["archive", "rule-b"], default="archive")
    args = parser.parse_args()
    print(json.dumps(validate_release(args.config, args.scope), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
