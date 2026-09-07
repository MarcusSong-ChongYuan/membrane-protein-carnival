#!/usr/bin/env python3
"""Compatibility runner correcting two unterminated caption literals in v0.1."""

from __future__ import annotations

import sys
from pathlib import Path


source_path = Path(__file__).with_name("build_v63_m7_m8_figures.py")
source = source_path.read_text(encoding="utf-8")
source = source.replace("metadata are sparse.\\n\",encoding='utf-8')", "metadata are sparse.\\n\"\"\",encoding='utf-8')")
source = source.replace("predicted binding affinity.\\n\",encoding='utf-8')", "predicted binding affinity.\\n\"\"\",encoding='utf-8')")
sys.argv[0] = str(source_path)
exec(compile(source, str(source_path), "exec"), {"__name__": "__main__", "__file__": str(source_path)})
