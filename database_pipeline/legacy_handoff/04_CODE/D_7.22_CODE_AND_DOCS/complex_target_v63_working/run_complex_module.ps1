param(
    [Parameter(Mandatory = $true)]
    [string]$Python,
    [Parameter(Mandatory = $true)]
    [string]$ProteinMaster,
    [Parameter(Mandatory = $true)]
    [string]$CompoundMaster,
    [Parameter(Mandatory = $true)]
    [string]$PdbeAssembly
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Raw = Join-Path $Root "raw"
$V01 = Join-Path $Root "output_v0_1"
$V02 = Join-Path $Root "output_v0_2"

& $Python (Join-Path $Root "scripts\build_complex_target_module_v0_1.py") `
    --protein-master $ProteinMaster `
    --complexportal-curated (Join-Path $Raw "complexportal_9606_20260114.tsv") `
    --complexportal-predicted (Join-Path $Raw "complexportal_9606_predicted_20260114.tsv") `
    --corum-human (Join-Path $Raw "corum_5_3_humanComplexes.txt") `
    --pdbe-assembly $PdbeAssembly `
    --output $V01

if ($LASTEXITCODE -ne 0) {
    throw "Complex module V0.1 failed."
}

& $Python (Join-Path $Root "scripts\finalize_complex_target_module_v0_2.py") `
    --protein-master $ProteinMaster `
    --compound-master $CompoundMaster `
    --v01 $V01 `
    --output $V02

if ($LASTEXITCODE -ne 0) {
    throw "Complex module V0.2 failed."
}

Write-Host "Complex-target V0.2 completed: $V02"
