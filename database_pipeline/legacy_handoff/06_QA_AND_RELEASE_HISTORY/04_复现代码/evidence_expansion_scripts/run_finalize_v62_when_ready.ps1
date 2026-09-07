$ErrorActionPreference = 'Stop'
$run = 'D:\7.22\evidence_expansion_v2_working\runs\v62_completion_20260729'
$requirements = @(
    (Join-Path $run 'qa\EXPRESSION_LOCATION_V2_VALIDATION.json'),
    (Join-Path $run 'qa\SUBUNIT_ASSEMBLY_V2_VALIDATION.json'),
    (Join-Path $run 'qa\NEGATIVE_EVIDENCE_RESOLUTION_V62_VALIDATION.json')
)
while ($true) {
    $ready = $true
    foreach ($path in $requirements) {
        if (-not (Test-Path -LiteralPath $path)) {
            $ready = $false
            break
        }
        try {
            $state = Get-Content -LiteralPath $path -Raw | ConvertFrom-Json
            if (-not ([string]$state.status).StartsWith('passed')) {
                $ready = $false
                break
            }
        } catch {
            $ready = $false
            break
        }
    }
    if ($ready) { break }
    Start-Sleep -Seconds 30
}
$python = 'C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $python 'D:\7.22\evidence_expansion_v2_working\scripts\preflight_v62_qa.py'
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& $python 'D:\7.22\evidence_expansion_v2_working\scripts\finalize_v62_release.py'
exit $LASTEXITCODE
