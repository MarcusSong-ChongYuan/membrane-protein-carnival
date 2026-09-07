param(
  [string]$Python = "python",
  [switch]$Install,
  [switch]$Full,
  [switch]$SkipPpt
)
$ErrorActionPreference = "Stop"
$Bundle = Split-Path -Parent $MyInvocation.MyCommand.Path
$Scripts = Join-Path $Bundle "08_reproducibility\scripts"
$env:MEMPRO_BUNDLE_ROOT = $Bundle
$env:MEMPRO_REVISION_ROOT = $Bundle

if ($Install) {
  & $Python -m venv (Join-Path $Bundle ".repro_env")
  $Python = Join-Path $Bundle ".repro_env\Scripts\python.exe"
  $Requirements = if ($Full) { "requirements-full-pagtn.txt" } else { "requirements-quick.txt" }
  & $Python -m pip install --upgrade pip
  & $Python -m pip install -r (Join-Path $Bundle $Requirements)
}

$Arguments = @((Join-Path $Scripts "agent_reproduce.py"))
if ($Full) { $Arguments += "--full" }
& $Python @Arguments

if (-not $SkipPpt) {
  $Node = Get-Command node -ErrorAction SilentlyContinue
  if ($Node -and $env:CODEX_NODE_MODULES) {
    & $Node.Source (Join-Path $Scripts "build_ppt.mjs")
  } else {
    Write-Warning "PPT rebuild skipped. Set CODEX_NODE_MODULES to node_modules containing @oai/artifact-tool. A verified PPT is already included."
  }
}
