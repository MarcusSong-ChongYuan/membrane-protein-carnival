from __future__ import annotations

import json
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(os.environ.get("MEMPRO_REVISION_ROOT", r"D:\finale\07_M1-M8图表修订_20260806"))
V62 = Path(r"D:\7.22\evidence_expansion_v2_working\releases\release_mempro_v6_2_20260730")
OUT = ROOT / "05_alphafold"
OUT.mkdir(parents=True, exist_ok=True)
CHECKPOINT = OUT / "alphafold_api_checkpoint.jsonl"
API = "https://alphafold.ebi.ac.uk/api/prediction/{}"


def fetch_json(url: str, retries: int = 4):
    for i in range(retries):
        try:
            r = requests.get(url, timeout=45, headers={"User-Agent": "MemProDB/figure-revision-20260806"})
            if r.status_code == 404:
                return None, 404
            r.raise_for_status()
            return r.json(), r.status_code
        except Exception:
            if i == retries - 1:
                return None, 599
            time.sleep(1.5 * (i + 1))


def fetch_meta(acc: str):
    data, status = fetch_json(API.format(acc))
    if not data:
        return {"target_uniprot_id": acc, "request_status": status}
    item = data[0] if isinstance(data, list) else data
    return {
        "target_uniprot_id": acc,
        "request_status": status,
        "entry_id": item.get("entryId") or item.get("modelEntityId"),
        "latest_version": item.get("latestVersion"),
        "model_created_date": item.get("modelCreatedDate"),
        "sequence_start": item.get("sequenceStart"),
        "sequence_end": item.get("sequenceEnd"),
        "global_mean_plddt": item.get("globalMetricValue"),
        "fraction_plddt_very_low": item.get("fractionPlddtVeryLow"),
        "fraction_plddt_low": item.get("fractionPlddtLow"),
        "fraction_plddt_confident": item.get("fractionPlddtConfident"),
        "fraction_plddt_very_high": item.get("fractionPlddtVeryHigh"),
        "plddt_url": item.get("plddtDocUrl"),
        "cif_url": item.get("cifUrl"),
        "license": "CC BY 4.0",
    }


def main():
    protein = pd.read_csv(V62 / "human_membrane_protein_master_v6_2.tsv", sep="\t",
                          usecols=["target_uniprot_id", "alphafolddb_ids"], dtype=str, keep_default_na=False)
    targets = sorted(protein.loc[protein.alphafolddb_ids.ne(""), "target_uniprot_id"].unique())
    done = {}
    if CHECKPOINT.exists():
        with CHECKPOINT.open("r", encoding="utf-8") as fh:
            for line in fh:
                try:
                    row = json.loads(line); done[row["target_uniprot_id"]] = row
                except Exception:
                    pass
    remaining = [x for x in targets if x not in done]
    with CHECKPOINT.open("a", encoding="utf-8", buffering=1) as out:
        with ThreadPoolExecutor(max_workers=12) as ex:
            futures = {ex.submit(fetch_meta, acc): acc for acc in remaining}
            for i, fut in enumerate(as_completed(futures), 1):
                row = fut.result(); done[row["target_uniprot_id"]] = row
                out.write(json.dumps(row, ensure_ascii=False) + "\n")
                if i % 250 == 0:
                    print(json.dumps({"metadata_completed": len(done), "total": len(targets)}), flush=True)
    meta = pd.DataFrame([done[x] for x in targets])
    meta.to_csv(OUT / "alphafold_model_metadata.tsv", sep="\t", index=False)

    side_path = ROOT / "06_binding_site_side/binding_site_membrane_side_instances.tsv.gz"
    local_rows = []
    if side_path.exists():
        sites = pd.read_csv(side_path, sep="\t", compression="gzip", dtype=str, keep_default_na=False)
        positions = {}
        for acc, grp in sites.groupby("target_uniprot_id"):
            pos = sorted({int(v) for cell in grp.residue_positions_uniprot for v in cell.split(";") if v.isdigit()})
            positions[acc] = pos
        meta_idx = meta.set_index("target_uniprot_id")

        def local_one(acc):
            if acc not in meta_idx.index or not str(meta_idx.loc[acc].get("plddt_url", "")):
                return {"target_uniprot_id": acc, "local_status": "no_model_or_url"}
            data, status = fetch_json(str(meta_idx.loc[acc]["plddt_url"]))
            if not data:
                return {"target_uniprot_id": acc, "local_status": f"http_{status}"}
            item = data[0] if isinstance(data, list) else data
            nums = item.get("residueNumber", []); scores = item.get("confidenceScore", [])
            lookup = {int(n): float(s) for n, s in zip(nums, scores)}
            vals = [lookup[p] for p in positions[acc] if p in lookup]
            return {"target_uniprot_id": acc, "local_status": "ok" if vals else "positions_outside_model",
                    "site_residue_count": len(positions[acc]), "mapped_site_residue_count": len(vals),
                    "site_mean_plddt": sum(vals)/len(vals) if vals else None,
                    "site_min_plddt": min(vals) if vals else None,
                    "site_fraction_plddt_ge_70": sum(x >= 70 for x in vals)/len(vals) if vals else None}

        with ThreadPoolExecutor(max_workers=10) as ex:
            futures = {ex.submit(local_one, acc): acc for acc in positions}
            for i, fut in enumerate(as_completed(futures), 1):
                local_rows.append(fut.result())
                if i % 100 == 0:
                    print(json.dumps({"local_plddt_completed": i, "total": len(positions)}), flush=True)
    local = pd.DataFrame(local_rows)
    local.to_csv(OUT / "alphafold_binding_site_local_plddt.tsv", sep="\t", index=False)
    qa = {
        "status": "PASS", "targets_with_alphafold_id": len(targets),
        "metadata_http_200": int((meta.request_status == 200).sum()),
        "metadata_missing_or_error": int((meta.request_status != 200).sum()),
        "local_site_targets_attempted": len(local),
        "local_site_targets_mapped": int((local.get("local_status", pd.Series(dtype=str)) == "ok").sum()),
        "source": "AlphaFold DB predictions API and per-residue confidence JSON",
        "source_version_policy": "API latestVersion captured per record; expected AFDB v6",
        "interpretation": "pLDDT is local model confidence, not experimental validation or biological correctness",
    }
    (OUT / "ALPHAFOLD_CONFIDENCE_QA.json").write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(qa, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
