# MemPro protein-module figure package

## Figure contract

These three standalone quantitative figures support one claim: MemPro's formal protein layer contains **7,800 canonical membrane proteins** with traceable membrane evidence, multi-source annotation coverage, and functional organization beyond a single family label.

## Output figures

1. `Figure1A_membrane_class_evidence_radial`: A/B/C membrane classes, each split into E1/E2/E3 evidence tiers. Each annular arc uses a common count scale; individual annuli are not mutually stacked.
2. `Figure1B_source_by_module_coverage_heatmap`: source-specific coverage measured as unique formal canonical proteins, rather than evidence rows. Blank cells are no direct protein-level coverage counted in the relevant input tables.
3. `Figure1C_membrane_role_molecular_function_enrichment`: Haldane–Anscombe-corrected Fisher-test log2 odds ratios. A black dot requires observed support >=20 and BH-FDR q<0.05. `family_defined_membrane_role_unresolved` is intentionally visible.

## Data integrity

- Frozen core protein denominator: 7,800.
- No rows were changed, inferred, or reclassified.
- Source identifiers in the heatmap are derived from source-linked protein fields; source rows and evidence records were not used as the statistical unit.
- The heatmap does not claim that contributing database count equals independent experimental count.
