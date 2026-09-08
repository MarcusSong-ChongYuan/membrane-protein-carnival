[CmdletBinding()]
param(
    [string]$Python = "py -3.13",
    [switch]$Figures
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$Venv = Join-Path $RepoRoot '.venv'

if (-not (Test-Path (Join-Path $Venv 'Scripts\\python.exe'))) {
    Invoke-Expression "$Python -m venv `"$Venv`""
}

$VenvPython = Join-Path $Venv 'Scripts\\python.exe'
# The core package has no third-party runtime dependencies.  Do not upgrade
# pip here: an offline lab workstation should still be able to verify a
# previously delivered FORMAL release.
& $VenvPython -m pip install --no-deps -e $RepoRoot
if ($Figures) {
    & $VenvPython -m pip install -r (Join-Path $RepoRoot 'requirements\\requirements-figures.txt')
}

Write-Host "MemPro environment ready: $Venv"
Write-Host "Next: set MEMPRO_DATA_ROOT to the extracted 01_database_FORMAL directory, then run .\\run.ps1 doctor"
