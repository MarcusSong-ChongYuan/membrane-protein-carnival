$ErrorActionPreference = "SilentlyContinue"
$total = [math]::Round((Get-PSDrive C).Used/1GB, 2)
$free = [math]::Round((Get-PSDrive C).Free/1GB, 2)
Write-Output "C: Drive  Total=200GB  Used=${total}GB  Free=${free}GB"
Write-Output "========================================="

# Fast check: known large folders/files
$checks = @(
    @{Name="Hibernation File"; Path="C:\hiberfil.sys"},
    @{Name="Page File"; Path="C:\pagefile.sys"},
    @{Name="Swap File"; Path="C:\swapfile.sys"}
)
foreach ($c in $checks) {
    if (Test-Path $c.Path) {
        $f = Get-Item $c.Path -Force
        Write-Output "$([math]::Round($f.Length/1GB,2)) GB -- $($c.Name) ($($c.Path))"
    }
}

# Top-level directories (level 1 only, fast)
Write-Output "`n--- Top-level C:\ folders ---"
Get-ChildItem C:\ -Directory | ForEach-Object {
    try {
        $sub = Get-ChildItem $_.FullName -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum
        if ($sub.Sum -gt 1073741824) {
            Write-Output "$([math]::Round($sub.Sum/1GB,1)) GB -- $($_.FullName)"
        }
    } catch {}
} | Sort-Object

Write-Output "`n--- User profile top folders ---"
$userDir = "C:\Users\Administrator"
Get-ChildItem $userDir -Directory -Force | ForEach-Object {
    try {
        $sub = Get-ChildItem $_.FullName -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum
        if ($sub.Sum -gt 1073741824) {
            Write-Output "$([math]::Round($sub.Sum/1GB,1)) GB -- $($_.FullName)"
        }
    } catch {}
} | Sort-Object

Write-Output "`n--- Temp/cleanup locations ---"
$temps = @(
    "C:\Windows\Temp",
    "C:\Users\Administrator\AppData\Local\Temp",
    "C:\Windows\SoftwareDistribution\Download"
)
foreach ($t in $temps) {
    if (Test-Path $t) {
        try {
            $sub = Get-ChildItem $t -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum
            Write-Output "$([math]::Round($sub.Sum/1GB,1)) GB -- $t"
        } catch {
            Write-Output "Access Denied -- $t"
        }
    }
}

Write-Output "`n--- AppData subfolders ---"
$appdata = "C:\Users\Administrator\AppData"
Get-ChildItem $appdata -Directory | ForEach-Object {
    $level1 = Get-ChildItem $_.FullName -Directory
    foreach ($sub in $level1) {
        try {
            $files = Get-ChildItem $sub.FullName -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum
            if ($files.Sum -gt 536870912) {
                Write-Output "$([math]::Round($files.Sum/1GB,1)) GB -- $($sub.FullName)"
            }
        } catch {}
    }
} | Sort-Object -Descending

Write-Output "`nDone."
