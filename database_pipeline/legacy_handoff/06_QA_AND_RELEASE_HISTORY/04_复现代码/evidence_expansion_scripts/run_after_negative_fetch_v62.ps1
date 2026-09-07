$ErrorActionPreference = 'Stop'
$run = 'D:\7.22\evidence_expansion_v2_working\runs\v62_completion_20260729'
$progress = Join-Path $run 'qa\NEGATIVE_CID_FETCH_V62_PROGRESS.json'
$python = 'C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$env:PYTHONPATH = 'D:\7.22\evidence_expansion_v2_working\tools\python_v61'
while ($true) {
    if (Test-Path -LiteralPath $progress) {
        try {
            $state = Get-Content -LiteralPath $progress -Raw | ConvertFrom-Json
            if ($state.status -eq 'complete') { break }
            if ($state.status -eq 'incomplete') {
                throw 'Negative CID fetch ended incomplete.'
            }
        } catch {
            if ($_.Exception.Message -eq 'Negative CID fetch ended incomplete.') {
                throw
            }
        }
    }
    Start-Sleep -Seconds 30
}
& $python 'D:\7.22\evidence_expansion_v2_working\scripts\resolve_negative_evidence_v62.py'
exit $LASTEXITCODE
