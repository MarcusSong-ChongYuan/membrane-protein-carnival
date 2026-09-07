Write-Output "Doubao exists: $(Test-Path C:\Users\Administrator\AppData\Local\Doubao)"
Write-Output ".cache exists: $(Test-Path C:\Users\Administrator\.cache)"

$t1 = (Get-ChildItem C:\Windows\Temp -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum).Sum
$t2 = (Get-ChildItem $env:LOCALAPPDATA\Temp -Recurse -File -ErrorAction Stop | Measure-Object Length -Sum).Sum
Write-Output "Windows Temp remaining: $([math]::Round($t1/1MB,1)) MB"
Write-Output "User Temp remaining: $([math]::Round($t2/1MB,1)) MB"

$free = [math]::Round((Get-PSDrive C).Free/1GB, 2)
Write-Output "C: free = ${free}GB"
