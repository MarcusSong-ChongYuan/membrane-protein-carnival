# Frozen upstream-input reproducibility

## Scope

MemPro FORMAL is a reproducible **frozen release**, not an assertion that a
later live API response can recreate the release byte-for-byte.  Upstream
resources change versions, endpoint schemas, curation, licensing and access
rules.  Replacing an archived 2026 input with a current download would create
a different dataset and must never be described as an exact rebuild.

This repository therefore has two separate reproducibility layers:

| Layer | What it proves | Official command |
|---|---|---|
| FORMAL release integrity | The distributed frozen tables have not changed. | `./run.ps1 verify-full` |
| Upstream snapshot inventory | Which historical inputs are actually retained and SHA-256 verified. | `./run.ps1 snapshot-audit` |

## Inventory files

- `config/sources/source_registry_v1.tsv` is the 27-source publication
  registry.  It records the source/version used by FORMAL, scope, historical
  snapshot state and the retained legacy parser/workflow where one exists.
- `config/sources/snapshot_manifest_v1.tsv` lists only raw inputs physically
  present in the migration package.  Every listed `PRESENT` input has a
  relative path, byte count and SHA-256 checksum.

The source registry deliberately lists resources for which the historical raw
input is no longer in the migration package.  These rows receive
`HISTORICAL_SNAPSHOT_UNAVAILABLE`, not a fabricated live-download substitute.

## Run the audit

```powershell
$env:MEMPRO_DATA_ROOT = 'E:\MemPro\MemPro_Complete_Migration_20260908\01_database_FORMAL'
./run.ps1 snapshot-audit
```

The default raw-input root is:

```text
<MEMPRO_DATA_ROOT>/06_scripts/legacy_handoff
```

Set `MEMPRO_RAW_SNAPSHOT_ROOT` if it was restored elsewhere:

```powershell
$env:MEMPRO_RAW_SNAPSHOT_ROOT = 'F:\MemPro_recovered_raw_snapshots'
./run.ps1 snapshot-audit
```

Outputs are written outside frozen tables:

```text
outputs/reproducibility/raw_snapshot_audit.json
outputs/reproducibility/raw_snapshot_audit.tsv
```

## Interpreting the result

- `SNAPSHOT_VERIFIED`: a retained raw snapshot passes its expected byte count
  and SHA-256 check.  It may be replayed only through a schema-specific
  ingestion workflow.
- `DERIVED_FROM_FORMAL`: MemPro-derived annotations are reproducible from the
  frozen FORMAL tables and documented workflow; they are not a new upstream
  source download.
- `HISTORICAL_SNAPSHOT_UNAVAILABLE`: the release keeps source provenance and
  legacy parser references, but the exact historical input was not retained in
  the current migration package.  This is a documented limitation, not an
  integrity failure.
- `SNAPSHOT_ERROR` or audit `FAIL`: a source declared `PRESENT` is missing or
  differs from its manifest.  This is an engineering error and must be fixed
  before claiming that source snapshot is available.

`PASS_WITH_LIMITATIONS` is the expected result for the current package: it
means all declared retained snapshots validate, while unavailable historical
sources remain openly reported.  It does **not** mean an exact all-source ETL
rebuild has been achieved.

## Recovery protocol for an original snapshot

When an original raw input is recovered, do not overwrite a release table.

1. Place the original file under a versioned path beneath the snapshot root.
2. Add its metadata, byte count and SHA-256 to `snapshot_manifest_v1.tsv`.
3. Change only the matching source registry row from
   `HISTORICAL_SNAPSHOT_NOT_IN_CURRENT_PACKAGE` to `SNAPSHOT_PRESENT`.
4. Add a deterministic parser with a documented input schema and expected
   output contract.
5. Run `snapshot-audit`; then compare parser output to the frozen FORMAL
   contribution tables in a new derived output directory.
6. Commit the registry, manifest, parser and comparison report together.

Do not use a current live source to fill an unrecovered historic snapshot
without creating a separately versioned future release.
