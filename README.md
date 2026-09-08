# MemPro

## Official reproducible entry point

The repository now distinguishes the **FORMAL frozen release workflow** from
historical development scripts.  The separately distributed migration archive
contains the release data; GitHub intentionally contains code and lightweight
documentation only.

```powershell
git clone https://github.com/MarcusSong-ChongYuan/membrane-protein-carnival.git mempro
cd mempro
.\setup.ps1
$env:MEMPRO_DATA_ROOT = 'E:\MemPro\MemPro_Complete_Migration_20260908\01_database_FORMAL'
.\run.ps1 doctor
.\run.ps1 verify
.\run.ps1 summary
.\run.ps1 snapshot-audit
```

If an older local pipeline backup exists, its surviving source exports can be
recovered into a separate manifest-verified context package with
`./run.ps1 recover-legacy-context`; this never modifies FORMAL tables.

Read [reproducibility status](docs/REPRODUCIBILITY_STATUS.md) before using a
figure workflow, website deployment or the Linux/HPC docking material.
Read [frozen upstream-input reproducibility](docs/RAW_SNAPSHOT_REPRODUCIBILITY.md)
before claiming an exact rebuild from raw source snapshots.

`database_pipeline/legacy_handoff/` is retained as a historical audit archive;
it is not the default execution path.

MemPro is a membrane-protein–small-molecule knowledgebase project. This repository contains reproducible **code, configuration and documentation** for maintaining the data pipeline, reproducing analyses/figures, and developing the web application.

## Repository layout

| Directory | Purpose |
|---|---|
| `database_pipeline/` | Source collection, normalization, QC, entity-resolution, release and Docking handoff scripts. Historical subdirectories are retained because they document pipeline provenance. |
| `website/` | Next.js/Vite website source and lightweight assets. |
| `figure_and_analysis_code/` | Figure-generation, clustering, chemical-space and exploratory analysis code. |
| `docs/` | Project workflow, design decisions, data summary, migration and web-developer documentation. |

## Data access

The complete formal database is deliberately not committed to GitHub. It includes large evidence, expression and negative-evidence tables and is distributed as a separately checksum-verified migration archive:

```text
MemPro_Complete_Migration_20260908.zip
```

After extracting the archive, set:

```powershell
$env:MEMPRO_DATA_ROOT = 'D:\MemPro\01_database_FORMAL'
```

Then use the pipeline and website code against that read-only data root. See [`docs/README_MIGRATION.md`](docs/README_MIGRATION.md) and [`docs/README_WEB_DEVELOPER.md`](docs/README_WEB_DEVELOPER.md).

## Quick start: web application

```powershell
cd website
npm ci
npm run dev
```

## Data/release principles

- The `FORMAL` release is read-only; generate application caches and derived databases elsewhere.
- `E1/E2/E3` are membrane-evidence levels; `BE1/BE2/BE3` are small-molecule evidence levels.
- Contributing source databases are not equivalent to independent experiments.
- RNA, IHC, localization and missing/unmapped data remain separately represented.
- External-source audit candidates are not automatically formal release data.

## Large-file policy

Do not commit raw source downloads, generated database tables, build directories, `node_modules`, temporary files, keys, or Docking outputs. Use releases/object storage with a SHA-256 manifest for formal data distribution.
