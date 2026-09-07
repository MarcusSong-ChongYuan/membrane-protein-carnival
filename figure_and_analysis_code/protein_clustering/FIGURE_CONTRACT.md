Core conclusion: Protein language-model embeddings may reveal family-level organization of the 7,800 formal membrane proteins; known classes are used only after clustering to evaluate it.
Figure archetype: quantitative grid with 3D-landscape hero panel.
Target output: publication-ready PDF/SVG/600-dpi PNG plus interactive HTML.
Backend: Python only.
Clustering: Standardized ESM-2 embeddings -> PCA; HDBSCAN never receives UMAP coordinates.
Visualization: independent 2D/3D UMAP projections of the PCA representation.
Reviewer risks: model-derived representation is not structural distance; UMAP does not prove cluster separation; target class is strictly post hoc; noise remains visible.