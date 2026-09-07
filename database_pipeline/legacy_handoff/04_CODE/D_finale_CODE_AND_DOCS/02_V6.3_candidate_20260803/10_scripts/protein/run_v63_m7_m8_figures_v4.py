#!/usr/bin/env python3
"""Final ASCII-safe runner for M7/M8."""

from __future__ import annotations

import re
import sys
from pathlib import Path


source_path = Path(__file__).with_name("build_v63_m7_m8_figures.py")
source = source_path.read_text(encoding="utf-8")
source = source.replace("metadata are sparse.\\n\",encoding='utf-8')", "metadata are sparse.\\n\"\"\",encoding='utf-8')")
source = source.replace("predicted binding affinity.\\n\",encoding='utf-8')", "predicted binding affinity.\\n\"\"\",encoding='utf-8')")
source = re.sub(r'"[^"\r\n]*2 modalities":', '"At least 2 modalities":', source, count=1)
sys.argv[0] = str(source_path)
exec(compile(source, str(source_path), "exec"), {"__name__": "__main__", "__file__": str(source_path)})
