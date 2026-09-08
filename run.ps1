[CmdletBinding()]
param(
    [Parameter(Mandatory = $true, Position = 0)]
    [ValidateSet('doctor', 'verify', 'verify-full', 'summary', 'snapshot-audit', 'figures', 'web')]
    [string]$Task,
    [string]$DataRoot,
    [string]$FigureOutputRoot,
    [string]$FigureAssetRoot,
    [switch]$Install,
    [switch]$SkipBuild
)

$ErrorActionPreference = 'Stop'
$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$VenvPython = Join-Path $RepoRoot '.venv\\Scripts\\python.exe'

if ($DataRoot) { $env:MEMPRO_DATA_ROOT = (Resolve-Path -LiteralPath $DataRoot).Path }
if ($Install -or -not (Test-Path $VenvPython)) { & (Join-Path $RepoRoot 'setup.ps1') }

switch ($Task) {
    'doctor'      { & $VenvPython -m mempro_repro doctor }
    'verify'      { & $VenvPython -m mempro_repro verify }
    'verify-full' { & $VenvPython -m mempro_repro verify --full }
    'summary'     { & $VenvPython -m mempro_repro summary }
    'snapshot-audit' { & $VenvPython -m mempro_repro snapshot-audit }
    'figures' {
        if (-not $env:MEMPRO_DATA_ROOT) { throw 'Set MEMPRO_DATA_ROOT or pass -DataRoot before regenerating figures.' }
        $Runner = Join-Path $RepoRoot 'figure_and_analysis_code\\24_MemPro_V7.2_NAR_analysis_and_figures_20260818\\01_scripts\\run_formal_figures.py'
        if (-not $FigureOutputRoot) { $FigureOutputRoot = Join-Path $RepoRoot 'outputs\\nar_figures_formal' }
        $Arguments = @($Runner, '--data-root', $env:MEMPRO_DATA_ROOT, '--output-root', $FigureOutputRoot)
        if ($FigureAssetRoot) { $Arguments += @('--asset-root', $FigureAssetRoot) }
        & $VenvPython @Arguments
    }
    'web' {
        $Website = Join-Path $RepoRoot 'website'
        if (-not (Test-Path (Join-Path $Website 'package.json'))) { throw 'Website source is missing.' }
        Push-Location $Website
        try {
            if (-not $SkipBuild) { npm ci; npm run build }
            npm run dev
        } finally { Pop-Location }
    }
}
