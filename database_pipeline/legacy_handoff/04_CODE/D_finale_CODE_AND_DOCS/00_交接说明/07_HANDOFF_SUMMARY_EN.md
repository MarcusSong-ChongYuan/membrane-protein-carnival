# MemPro V6.2 handoff summary

This archive contains the frozen MemPro V6.2 derived release, publication-ready figures and source data, docking-prioritization lists, reproducibility scripts/configuration snapshots, and QA reports.

The release validation status is PASS with zero blocking errors. Use `01_正式数据_V6.2` as the immutable analysis input. Review, conflict, E3 and unmapped tables must not be mixed into the default high-confidence layer without an explicit policy.

Third-party raw archives are not redistributed. A clean rebuild requires reacquiring the source versions listed in `05_SOURCE_VERSIONS.tsv` and adapting the original Windows absolute paths in the scripts.

The docking lists are prioritization outputs only. Receptor state/chain/assembly, ligand microstates, docking boxes, redocking and protocol validation remain required before HPC production runs.
