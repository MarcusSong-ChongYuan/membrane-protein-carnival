# MemPro protein-language-model clustering

This restart-safe Python workflow constructs one record per canonical UniProt accession from the frozen MemPro formal protein table, extracts an ESM-2 protein representation, performs PCA and HDBSCAN in representation space, and generates independent 2D/3D UMAP visualizations.

## Scientific boundary

HDBSCAN is run on standardized PCA embeddings, never on UMAP coordinates. Existing membrane-role labels are joined only after clustering for composition, ARI/NMI, and display. UMAP is a visualization layer, not evidence of a structural distance or a clustering input.

## Run

```powershell
D:\finale\_envs\protein_clustering\Scripts\python.exe run_pipeline.py
```

The embedding stage caches one vector per protein under `results/embedding_chunks/`, so rerunning resumes safely after interruption. All output is confined to `results/`; frozen MemPro tables are read only.
