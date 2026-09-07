$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
$manifest = Join-Path $root '08_MANIFESTS\PORTABLE_FILES_SHA256.tsv'

if (-not (Test-Path -LiteralPath $manifest)) {
    throw "Missing manifest: $manifest"
}

$rows = Import-Csv -LiteralPath $manifest -Delimiter "`t"
$bad = [System.Collections.Generic.List[object]]::new()
foreach ($row in $rows) {
    $path = Join-Path $root $row.relative_path
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        $bad.Add([pscustomobject]@{relative_path=$row.relative_path; status='MISSING'})
        continue
    }
    $actual = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash
    if ($actual -ne $row.sha256) {
        $bad.Add([pscustomobject]@{relative_path=$row.relative_path; status='HASH_MISMATCH'})
    }
}

if ($bad.Count -gt 0) {
    $bad | Format-Table -AutoSize
    throw "Handoff verification failed: $($bad.Count) problem(s)."
}

Write-Host "PASS: $($rows.Count) portable files match SHA-256 manifest." -ForegroundColor Green
Write-Host "Note: 07_LOCAL_FULL_WORKSPACE_LINKS contains same-machine links and is intentionally excluded." -ForegroundColor Yellow
