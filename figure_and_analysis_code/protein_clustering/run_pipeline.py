"""One-command, restart-safe MemPro protein clustering pipeline."""
from __future__ import annotations
import subprocess,sys
from pathlib import Path
from pipeline_lib import ROOT,load_config,software_versions
STEPS=['prepare_proteins.py','extract_embeddings.py','embedding_qc.py','run_pca.py','run_clustering.py','run_umap.py','evaluate_clusters.py','plot_results.py','generate_report.py']
def main():
 cfg=load_config();(cfg['output_dir']/'software_versions.txt').write_text(software_versions(),encoding='utf-8')
 for step in STEPS:
  print(f'=== {step} ===',flush=True);subprocess.run([sys.executable,str(ROOT/step)],check=True,cwd=ROOT)
if __name__=='__main__':main()
