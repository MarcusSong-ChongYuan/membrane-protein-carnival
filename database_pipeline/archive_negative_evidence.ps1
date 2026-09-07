$ErrorActionPreference = 'Stop'

$formal = 'D:\finale\FORMAL'
$archive = Join-Path $formal '04_archive_negative_evidence'
$catalogPath = Join-Path $formal '05_metadata\NEGATIVE_EVIDENCE_ARCHIVE_CATALOG.tsv'
$roots = @('D:\finale', 'D:\7.22')

$allowedSuffixes = @(
    '.tsv.gz', '.csv.gz', '.txt.gz', '.jsonl.gz', '.json.gz',
    '.parquet', '.feather', '.arrow', '.tsv', '.csv'
)

$files = foreach ($root in $roots) {
    if (-not (Test-Path -LiteralPath $root)) { continue }
    Get-ChildItem -LiteralPath $root -Recurse -File -ErrorAction SilentlyContinue |
        Where-Object {
            $candidateName = $_.Name.ToLowerInvariant()
            $_.FullName -notlike "$formal*" -and
            $_.Name -match '(?i)negative' -and
            ($allowedSuffixes | Where-Object { $candidateName.EndsWith($_) })
        }
}

$files = $files | Sort-Object FullName -Unique
$byHash = @{}
$rows = [System.Collections.Generic.List[object]]::new()

foreach ($file in $files) {
    $hash = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    if (-not $byHash.ContainsKey($hash)) {
        $safeName = ($file.Name -replace '[^A-Za-z0-9._-]', '_')
        $archiveName = $hash.Substring(0, 12) + '__' + $safeName
        $archivePath = Join-Path $archive $archiveName
        Copy-Item -LiteralPath $file.FullName -Destination $archivePath -Force
        $byHash[$hash] = $archiveName
    }
    $rows.Add([pscustomobject]@{
        original_path = $file.FullName
        bytes = $file.Length
        sha256 = $hash
        archive_file = $byHash[$hash]
        duplicate_content = if (($rows | Where-Object sha256 -eq $hash).Count -gt 0) { 1 } else { 0 }
    })
}

$rows | Export-Csv -LiteralPath $catalogPath -Delimiter "`t" -NoTypeInformation -Encoding utf8

$summary = [ordered]@{
    scanned_files = $rows.Count
    unique_content_files = $byHash.Count
    original_bytes = ($rows | Measure-Object bytes -Sum).Sum
    archived_bytes = (Get-ChildItem -LiteralPath $archive -File | Measure-Object Length -Sum).Sum
    catalog = $catalogPath
}
$summary | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $formal '07_QA\NEGATIVE_EVIDENCE_ARCHIVE_QA.json') -Encoding utf8
$summary | ConvertTo-Json
