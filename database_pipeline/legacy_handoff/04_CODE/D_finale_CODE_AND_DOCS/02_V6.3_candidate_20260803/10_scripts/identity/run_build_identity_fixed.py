#!/usr/bin/env python3
"""Compatibility runner that counts quoted TSV records instead of physical lines."""

from __future__ import annotations

import csv
import build_identity_layer as module


def record_count(path):
    with module.open_text(path) as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


module.count_rows = record_count
module.main()
