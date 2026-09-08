# MemPro reproducibility status

## What this repository can reproduce now

With the separately distributed, checksum-verified `MemPro_Complete_Migration_20260908.zip` data package:

1. Validate the extracted FORMAL release against its release manifest.
2. Recount the main frozen tables and compare them with the FORMAL QA report.
3. Audit the 27-source registry and hash-check every historical raw snapshot
   that is actually retained in the migration package.
4. Start the website source after its static data bundle is present in `website/public/data`.
5. Run the preserved, named analysis workflows after installing the figure dependency profile.

## What is deliberately not claimed yet

- Rebuilding FORMAL from every live upstream API is not deterministic: upstream releases, licenses, and API responses change.
- Docking is not portable to a normal desktop. The HPC/Slurm workflow is separately documented and must be run on a compatible Linux cluster.
- Historical `legacy_handoff/` scripts are preserved for audit only and are **not** part of the official entry point.

## Current official commands

```powershell
# New Windows machine
git clone https://github.com/MarcusSong-ChongYuan/membrane-protein-carnival.git mempro
cd mempro
.\setup.ps1

# Point to the extracted migration-package database.
$env:MEMPRO_DATA_ROOT = 'E:\MemPro\MemPro_Complete_Migration_20260908\01_database_FORMAL'
.\run.ps1 doctor
.\run.ps1 verify
.\run.ps1 summary
.\run.ps1 snapshot-audit
```

Use `.\run.ps1 verify-full` for a complete file-by-file SHA-256 verification. It can take several minutes for the multi-gigabyte release.

`snapshot-audit` does not contact upstream APIs. It reports which of the 27
historical inputs are physically retained and validates only those declared
snapshots. See [RAW_SNAPSHOT_REPRODUCIBILITY.md](RAW_SNAPSHOT_REPRODUCIBILITY.md)
for the distinction between an integrity error and an unavailable historical
input.

## NAR figure workflow

The official figure workflow is now path-portable, but it is intentionally a
separate optional profile because RDKit, UMAP and the plotting stack are not
needed for release validation.

```powershell
.\setup.ps1 -Figures
.\run.ps1 figures -DataRoot $env:MEMPRO_DATA_ROOT `
  -FigureOutputRoot 'E:\MemPro\derived\nar_figures' `
  -FigureAssetRoot 'E:\MemPro\MemPro_Complete_Migration_20260908\04_figures_and_analysis\01_nar_figures_v72\05_graphics'
```

The output directory must be new.  The workflow does not overwrite historical
figures or frozen tables.

## Acceptance target for the next validation pass

The next pass must use a clean computer or clean virtual environment and record:

- operating system and Python/Node versions;
- result of `doctor`, sampled `verify`, and `summary`;
- whether website build succeeds after re-staging its data bundle;
- whether selected figure workflows reproduce their published source-data counts.
