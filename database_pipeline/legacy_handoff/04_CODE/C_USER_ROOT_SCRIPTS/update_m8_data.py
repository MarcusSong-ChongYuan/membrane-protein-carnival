import csv, json
from collections import Counter
from pathlib import Path

plot = Path(r"C:\Users\Administrator\mempro_ppt_final_redraw\data\plot_data.json")
pilot = Path(r"D:\finale\03_V6.3.1_readiness_20260803\04_docking_hpc\docking_pilot_manifest_v631.tsv")
data = json.loads(plot.read_text(encoding="utf-8"))
tiers=Counter(); structures=Counter(); receptor=Counter(); ligand=Counter(); boxes=Counter(); targets=set(); gates=Counter()
with pilot.open(encoding="utf-8", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        tiers[r["docking_tier"]] += 1
        structures[r["structure_grade"]] += 1
        receptor[r["receptor_preparation_status_v2"]] += 1
        ligand[r["ligand_preparation_status_v2"]] += 1
        boxes[r["box_definition_status_v2"]] += 1
        gates[r["docking_execution_status_v631"]] += 1
        targets.add(r["target_uniprot_id"])
data["M8"] = {
    "pairs": sum(tiers.values()), "targets": len(targets), "tiers": dict(tiers),
    "structure_grade": structures.most_common(), "receptor_prep": receptor.most_common(),
    "ligand_prep": ligand.most_common(), "box_status": boxes.most_common(), "execution": dict(gates),
}
plot.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
print(data["M8"])
