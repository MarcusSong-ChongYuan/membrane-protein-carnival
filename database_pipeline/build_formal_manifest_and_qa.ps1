$ErrorActionPreference = 'Stop'

$formal = 'D:\finale\FORMAL'
$metadata = Join-Path $formal '05_metadata'
$qaDir = Join-Path $formal '07_QA'
$manifestPath = Join-Path $metadata 'FORMAL_MANIFEST_SHA256.tsv'
$rolePath = Join-Path $formal '03_publication_repairs\protein_membrane_role_FORMAL.tsv'
$sourcePath = Join-Path $metadata 'SOURCE_REGISTRY_FORMAL.tsv'
$negativeCatalog = Join-Path $metadata 'NEGATIVE_EVIDENCE_ARCHIVE_CATALOG.tsv'

$roles = Import-Csv -LiteralPath $rolePath -Delimiter "`t"
$sourceRegistry = Import-Csv -LiteralPath $sourcePath -Delimiter "`t"
$negativeRows = Import-Csv -LiteralPath $negativeCatalog -Delimiter "`t"

$upgraded = @($roles | Where-Object { $_.role_resolution_status -like 'UPGRADED*' }).Count
$unresolved = @($roles | Where-Object { $_.formal_primary_membrane_role -eq 'family_defined_membrane_role_unresolved' }).Count
$negativeArchiveFiles = Get-ChildItem -LiteralPath (Join-Path $formal '04_archive_negative_evidence') -File
$archiveHashSet = @{}
foreach ($f in $negativeArchiveFiles) {
    $archiveHashSet[(Get-FileHash -LiteralPath $f.FullName -Algorithm SHA256).Hash.ToLowerInvariant()] = $true
}
$missingNegativeHashes = @($negativeRows | Where-Object { -not $archiveHashSet.ContainsKey($_.sha256) })

$checks = [ordered]@{
    formal_protein_count_7800 = ($roles.Count -eq 7800)
    strict_role_upgrades_17 = ($upgraded -eq 17)
    formal_role_unresolved_1705 = ($unresolved -eq 1705)
    source_registry_27 = ($sourceRegistry.Count -eq 27)
    source_registry_unique = (($sourceRegistry.source_database | Sort-Object -Unique).Count -eq 27)
    negative_catalog_present = (Test-Path -LiteralPath $negativeCatalog)
    negative_archive_hashes_complete = ($missingNegativeHashes.Count -eq 0)
    v72_validation_pass = ((Get-Content -LiteralPath (Join-Path $formal '01_core_v72\03_QA\V72_VALIDATION_REPORT.json') -Raw | ConvertFrom-Json).status -eq 'PASS')
    v721_validation_pass = ((Get-Content -LiteralPath (Join-Path $formal '02_companion_v721\05_QA\V721_VALIDATION_REPORT.json') -Raw | ConvertFrom-Json).status -eq 'PASS')
}

$failed = @($checks.GetEnumerator() | Where-Object { -not $_.Value })
$manifestFileCount = @(Get-ChildItem -LiteralPath $formal -Recurse -File |
    Where-Object { $_.FullName -ne $manifestPath }).Count
$report = [ordered]@{
    release = 'MemPro FORMAL'
    generated_at = (Get-Date).ToString('o')
    status = if ($failed.Count -eq 0) { 'PASS' } else { 'FAIL' }
    fixed_counts = [ordered]@{
        formal_proteins = 7800
        positive_evidence = 942455
        formal_pairs = 529168
        source_registry_entries = $sourceRegistry.Count
        role_upgrades = $upgraded
        unresolved_formal_roles = $unresolved
        negative_original_paths_cataloged = $negativeRows.Count
        negative_unique_content_files = $negativeArchiveFiles.Count
    }
    checks = $checks
    failed_checks = @($failed | ForEach-Object Key)
    manifest_files = $manifestFileCount
}
$report | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath (Join-Path $qaDir 'FORMAL_VALIDATION_REPORT.json') -Encoding utf8

# Write all derived metadata before hashing it.  The earlier ordering wrote the
# manifest first, then changed FORMAL_VALIDATION_REPORT.json, which made the
# release appear corrupt on a full checksum verification.
$files = Get-ChildItem -LiteralPath $formal -Recurse -File |
    Where-Object { $_.FullName -ne $manifestPath } |
    Sort-Object FullName
$manifest = foreach ($f in $files) {
    [pscustomobject]@{
        relative_path = $f.FullName.Substring($formal.Length + 1)
        bytes = $f.Length
        sha256 = (Get-FileHash -LiteralPath $f.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    }
}
if ($manifest.Count -ne $manifestFileCount) { throw 'Manifest file count changed during generation.' }
$manifest | Export-Csv -LiteralPath $manifestPath -Delimiter "`t" -NoTypeInformation -Encoding utf8
$report | ConvertTo-Json -Depth 6
if ($failed.Count -gt 0) { exit 2 }
