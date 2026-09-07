# MemPro — complete project migration package

**Prepared:** 2026-09-08  
**Source baseline:** `D:\finale\FORMAL`  
**Purpose:** move the MemPro database project to a new computer without relying on historical working directories.

## What this package contains

| Folder | Contents |
|---|---|
| `01_database_FORMAL` | Full formal MemPro database: V7.2 core, V7.2.1 companion annotations, publication metadata/QA, frozen negative-evidence archive, source data dictionary, and the current external-source audit workspace. |
| `02_pipeline_and_api_code` | Database collection, normalization, entity-resolution, QA, website/API planning, Docking handoff and figure-generation code from the project handoff archive. |
| `03_web_application_source` | Website source code and static assets, with `package.json` / lockfile. Generated build folders, node_modules and Git history are intentionally excluded. |
| `04_figures_and_analysis` | NAR figure sources/exports, editable versions, protein/compound modules, PyMOL example, graphical-abstract drafts, clustering, EDA and network analysis. |
| `05_docs_and_handoff` | Project decisions, workflow, state, new-agent prompt, database summary and website-developer handoff. |
| `06_environment` | New-computer setup guidance and portable environment manifests. |

## Deliberate exclusions

- Superseded V3–V6 and candidate-release directories; they are historical duplicates, not needed to reproduce the current formal release.
- Temporary caches, virtual environments, `node_modules`, compiled `.next`/`dist` folders and Git object history.
- A100-only live Docking jobs/results not present on this PC.

## First steps on the new computer

1. Extract the archive to a short local path, preferably `D:\MemPro`.
2. Verify `PACKAGE_SHA256.tsv` before use:

   ```powershell
   .\VERIFY_MIGRATION_PACKAGE.ps1
   ```

3. Read `05_docs_and_handoff\README_FIRST.md`, then `02_END_TO_END_WORKFLOW.md`.
4. For the website:

   ```powershell
   cd .\03_web_application_source
   npm ci
   npm run dev
   ```

5. Treat `01_database_FORMAL` as read-only source data. Build any PostgreSQL, DuckDB, Parquet, search-index or API cache under a separate application `workspace/` folder.

## Important scientific/data semantics

- `E1/E2/E3` are membrane-evidence levels; `BE1/BE2/BE3` are small-molecule evidence levels. Never conflate them.
- Database contribution count does not mean independently replicated experiments.
- RNA, IHC, localization and missing/unmapped status must remain distinct.
- Binding-site absence does not imply absence of an interaction.
- The current DrugCentral/TTD/ChEBI audit content is an audit workspace, not automatically part of the formal public-release pair count.

## GitHub boundary

The GitHub repository contains code, documentation, configuration and tiny example schemas—not 5+ GB of data. Use this package or separate release attachments/object storage for the formal data; retain the SHA-256 manifests.
