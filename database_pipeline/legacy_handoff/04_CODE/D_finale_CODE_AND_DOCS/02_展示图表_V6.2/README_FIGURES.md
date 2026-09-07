# MemPro V6.2 scientific figure package

This package contains 14 concentrated 2×3 composite figures: eight main figures (M1–M8) and six supplementary figures (S1–S6).

## Formats

- `svg/`: editable vector artwork with live text.
- `pdf/`: publication-ready vector output.
- `png/`: 600 dpi raster output for review and submission systems.
- `data/`: panel-level plotted statistics and deterministic embedding coordinates.
- `captions/`: detailed captions, including interpretation limits.
- `qa/`: input profiles, validation reports and contact sheets.

## Visual semantics

Membrane classes A/B/C use blue/green/amber; evidence tiers E1/BE1, E2/BE2 and E3/BE3 use dark blue, green and orange; positive, negative, conflict and review states use green, grey, vermilion and magenta. Docking tiers R0/D1/D2/D3 use purple, green, ochre and slate blue. The palette is color-vision-aware and remains interpretable through labels in grayscale.

## Reproducibility

All figures were generated from the frozen V6.2 release. UMAP and rendered compound-space samples use random seed 42. Chemical-space coordinates are saved in `data/`. The full categorical totals use complete source tables; large scatter/hexbin panels use deterministic sampling solely for rendering.

Package QA status: **PASS**.
