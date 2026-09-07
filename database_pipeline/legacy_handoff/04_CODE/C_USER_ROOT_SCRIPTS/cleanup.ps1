$ErrorActionPreference = "SilentlyContinue"
$totalFreed = 0

# 1. Clean Windows Temp
$winTemp = "C:\Windows\Temp"
Write-Output "=== Cleaning Windows Temp ==="
$before = (Get-ChildItem $winTemp -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum).Sum
Get-ChildItem $winTemp -Recurse -ErrorAction Stop | Remove-Item -Force -Recurse -ErrorAction SilentlyContinue
$after = (Get-ChildItem $winTemp -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum).Sum
$freed = $before - $after
$totalFreed += $freed
Write-Output "Freed: $([math]::Round($freed/1MB, 1)) MB"

# 2. Clean User Temp
$userTemp = "$env:LOCALAPPDATA\Temp"
Write-Output "`n=== Cleaning User Temp ==="
$before = (Get-ChildItem $userTemp -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum).Sum
Get-ChildItem $userTemp -Recurse -ErrorAction Stop | Remove-Item -Force -Recurse -ErrorAction SilentlyContinue
$after = (Get-ChildItem $userTemp -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum).Sum
$freed = $before - $after
$totalFreed += $freed
Write-Output "Freed: $([math]::Round($freed/1MB, 1)) MB"

# 3. Clean .cache
$cacheDir = "C:\Users\Administrator\.cache"
if (Test-Path $cacheDir) {
    Write-Output "`n=== Cleaning .cache ==="
    $before = (Get-ChildItem $cacheDir -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum).Sum
    Get-ChildItem $cacheDir -Recurse -ErrorAction Stop | Remove-Item -Force -Recurse -ErrorAction SilentlyContinue
    $after = (Get-ChildItem $cacheDir -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum).Sum
    $freed = $before - $after
    $totalFreed += $freed
    Write-Output "Freed: $([math]::Round($freed/1MB, 1)) MB"
}

# 4. Remove Doubao (豆包)
$doubaoDir = "C:\Users\Administrator\AppData\Local\Doubao"
if (Test-Path $doubaoDir) {
    Write-Output "`n=== Removing Doubao ==="
    $doubaoSize = (Get-ChildItem $doubaoDir -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum).Sum
    Write-Output "Doubao size: $([math]::Round($doubaoSize/1GB, 2)) GB"
    Remove-Item $doubaoDir -Force -Recurse -ErrorAction SilentlyContinue
    if (-not (Test-Path $doubaoDir)) {
        $totalFreed += $doubaoSize
        Write-Output "Doubao removed successfully!"
    } else {
        Write-Output "Doubao removal partially done."
    }
}

# Summary
Write-Output "`n========================================"
Write-Output "TOTAL FREED: $([math]::Round($totalFreed/1GB, 2)) GB"

$free = [math]::Round((Get-PSDrive C).Free/1GB, 2)
Write-Output "C: Drive free space now: ${free}GB"
