# Check top-level folders on C drive
Get-ChildItem -Path C:\ -Directory -ErrorAction SilentlyContinue | ForEach-Object {
    $size = (Get-ChildItem -Path $_.FullName -Recurse -ErrorAction SilentlyContinue | Measure-Object -Property Length -Sum).Sum
    $sizeGB = [math]::Round($size/1GB, 2)
    Write-Output "$sizeGB GB  --  $($_.FullName)"
} | Sort-Object -Descending

Write-Output "`n=== User profile folders ==="
Get-ChildItem -Path "C:\Users\Administrator" -Directory -ErrorAction SilentlyContinue | ForEach-Object {
    $size = (Get-ChildItem -Path $_.FullName -Recurse -ErrorAction SilentlyContinue | Measure-Object -Property Length -Sum).Sum
    $sizeGB = [math]::Round($size/1GB, 2)
    Write-Output "$sizeGB GB  --  $($_.FullName)"
} | Sort-Object -Descending

Write-Output "`n=== Common cleanup locations ==="
$commonPaths = @(
    "C:\Windows\Temp",
    "C:\Users\Administrator\AppData\Local\Temp",
    "C:\Users\Administrator\Downloads",
    "C:\Windows\SoftwareDistribution\Download",
    "C:\Windows\Prefetch"
)
foreach ($p in $commonPaths) {
    if (Test-Path $p) {
        $size = (Get-ChildItem -Path $p -Recurse -ErrorAction SilentlyContinue | Measure-Object -Property Length -Sum).Sum
        $sizeGB = [math]::Round($size/1GB, 2)
        Write-Output "$sizeGB GB  --  $p"
    } else {
        Write-Output "NOT FOUND  --  $p"
    }
}

Write-Output "`n=== Recycle Bin ==="
$rb = "C:\`$Recycle.Bin"
if (Test-Path $rb) {
    $size = (Get-ChildItem -Path $rb -Recurse -Force -ErrorAction SilentlyContinue | Measure-Object -Property Length -Sum).Sum
    $sizeGB = [math]::Round($size/1GB, 2)
    Write-Output "$sizeGB GB  --  Recycle Bin"
}
