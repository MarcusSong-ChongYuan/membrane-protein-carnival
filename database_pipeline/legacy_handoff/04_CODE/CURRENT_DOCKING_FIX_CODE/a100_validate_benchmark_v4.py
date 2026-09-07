#!/usr/bin/env python3
"""Validate paired Vina-GPU thread=5000/8000 benchmark outputs."""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata"
RESULTS = ROOT / "benchmark_results"
RUN = META / "thread_benchmark_run.tsv"
OUT = META / "thread_benchmark_validation_v4.tsv"
SUMMARY = META / "thread_benchmark_validation_v4.txt"


def parse_pose(path: Path):
    score = None
    coords = []
    in_first = False
    text = path.read_text(encoding="utf-8", errors="ignore")
    has_models = "MODEL" in text
    for line in text.splitlines():
        if line.startswith("MODEL"):
            if in_first:
                break
            in_first = True
            continue
        if line.startswith("ENDMDL") and in_first:
            break
        if "REMARK VINA RESULT:" in line and score is None:
            match = re.search(r"RESULT:\s+(-?\d+(?:\.\d+)?)", line)
            if match:
                score = float(match.group(1))
        if line.startswith(("ATOM  ", "HETATM")):
            if has_models and not in_first:
                continue
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                pass
    return score, coords


def rmsd(a, b):
    if len(a) != len(b) or not a:
        return None
    return math.sqrt(sum((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
                         for (x1, y1, z1), (x2, y2, z2) in zip(a, b)) / len(a))


def main() -> int:
    if not RUN.exists():
        SUMMARY.write_text("FAIL: benchmark run table missing\n", encoding="utf-8")
        return 2
    with RUN.open(encoding="utf-8", newline="") as handle:
        runs = list(csv.DictReader(handle, delimiter="\t"))
    exit_codes = {row["config"]: int(row["exit_code"]) for row in runs}
    task_ids = sorted({re.sub(r"_t(?:5000|8000)$", "", key) for key in exit_codes})
    rows = []
    for task_id in task_ids:
        key5, key8 = f"{task_id}_t5000", f"{task_id}_t8000"
        p5 = RESULTS / f"{key5}_out.pdbqt"
        p8 = RESULTS / f"{key8}_out.pdbqt"
        score5, xyz5 = parse_pose(p5) if p5.exists() else (None, [])
        score8, xyz8 = parse_pose(p8) if p8.exists() else (None, [])
        pose_rmsd = rmsd(xyz5, xyz8)
        valid = (exit_codes.get(key5) == 0 and exit_codes.get(key8) == 0
                 and score5 is not None and score8 is not None
                 and len(xyz5) >= 5 and len(xyz8) >= 5 and pose_rmsd is not None)
        rows.append({
            "task_id": task_id, "exit_5000": exit_codes.get(key5, -1),
            "exit_8000": exit_codes.get(key8, -1), "score_5000": score5,
            "score_8000": score8,
            "abs_score_delta": None if score5 is None or score8 is None else abs(score8-score5),
            "top_pose_direct_rmsd": pose_rmsd, "valid_pair": int(valid),
        })
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        fields = list(rows[0]) if rows else ["task_id"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    valid = [row for row in rows if row["valid_pair"]]
    score_deltas = sorted(row["abs_score_delta"] for row in valid)
    rmsds = sorted(row["top_pose_direct_rmsd"] for row in valid)
    median_score = score_deltas[len(score_deltas)//2] if score_deltas else float("inf")
    median_rmsd = rmsds[len(rmsds)//2] if rmsds else float("inf")
    stable = [row for row in valid if row["abs_score_delta"] <= 1.0 and row["top_pose_direct_rmsd"] <= 3.0]
    pass_gate = (len(rows) >= 5 and len(valid) == len(rows)
                 and median_score <= 1.0 and len(stable) / len(valid) >= .60)
    text = (
        f"status={'PASS' if pass_gate else 'FAIL'}\n"
        f"paired_tasks={len(rows)}\nvalid_pairs={len(valid)}\n"
        f"median_abs_score_delta={median_score:.4f}\n"
        f"median_top_pose_direct_rmsd={median_rmsd:.4f}\n"
        f"stable_fraction={(len(stable)/len(valid) if valid else 0):.4f}\n"
        "criteria=all outputs valid; median |score delta| <= 1.0 kcal/mol; >=60% pairs with |score delta| <=1.0 and direct RMSD <=3.0 A\n"
        "note=direct RMSD is diagnostic and not symmetry-corrected\n"
    )
    SUMMARY.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if pass_gate else 5


if __name__ == "__main__":
    sys.exit(main())
