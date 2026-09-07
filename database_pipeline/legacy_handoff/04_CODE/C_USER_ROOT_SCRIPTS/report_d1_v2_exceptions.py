import csv,json
p='/home/csong/docking/handoff/MemPro_Docking_D1_AE1_20260813/benchmark_v2_adjusted/D1_ADJUSTED_BENCHMARK_TASKS.tsv'
with open(p,encoding='utf-8') as f:
 rows=list(csv.DictReader(f,delimiter='\t'))
keys=['task_id','preferred_name','target_symbol','ligand_atom_count_v2','ligand_rotor_count_v2','ligand_net_charge_v2','ligand_max_span_A_v2','adaptive_box_x_v2','adaptive_box_y_v2','adaptive_box_z_v2','method_class_v2','benchmark_v2_disposition','benchmark_v2_reason','ligand_qc_flags_v2']
print(json.dumps([{k:r.get(k,'') for k in keys} for r in rows if r['benchmark_v2_disposition']!='READY_FOR_ADJUSTED_BENCHMARK'],indent=2))
