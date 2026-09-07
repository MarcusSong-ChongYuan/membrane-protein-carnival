$ErrorActionPreference = 'Stop'
$python = 'C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$scripts = 'D:\7.22\evidence_expansion_v2_working\scripts'
& $python (Join-Path $scripts 'prepare_negative_identity_v62.py')
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
& $python (Join-Path $scripts 'fetch_negative_pubchem_properties_v62.py')
exit $LASTEXITCODE
