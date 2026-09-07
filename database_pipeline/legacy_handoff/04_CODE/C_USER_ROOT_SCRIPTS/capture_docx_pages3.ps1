param(
  [string]$DocxPath = 'C:\Users\Administrator\synbio_defense_build\docx-render\speech.docx',
  [string]$OutputDir = 'C:\Users\Administrator\synbio_defense_build\docx-pages-full'
)
Add-Type -AssemblyName System.Drawing
Add-Type @'
using System;
using System.Runtime.InteropServices;
public static class WindowCapture3 {
  [StructLayout(LayoutKind.Sequential)] public struct RECT { public int Left; public int Top; public int Right; public int Bottom; }
  [DllImport("user32.dll")] public static extern bool GetWindowRect(IntPtr hWnd, out RECT rect);
  [DllImport("user32.dll")] public static extern bool PrintWindow(IntPtr hWnd, IntPtr hdcBlt, uint nFlags);
}
'@
New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
$word = New-Object -ComObject Word.Application
$word.Visible = $true
$word.DisplayAlerts = 0
$word.AutomationSecurity = 3
$word.Options.SaveNormalPrompt = $false
foreach ($addin in @($word.COMAddIns)) { try { $addin.Connect = $false } catch {} }
try {
  $doc = $word.Documents.OpenNoRepairDialog($DocxPath, $false, $true, $false)
  $window = $word.ActiveWindow
  $window.WindowState = 1
  $window.View.Type = 3
  $window.View.Zoom.PageFit = 1
  try { $window.DocumentMap = $false } catch {}
  $pages = $doc.ComputeStatistics(2)
  Start-Sleep -Milliseconds 1200
  $process = Get-Process WINWORD | Where-Object { $_.MainWindowHandle -ne 0 } | Sort-Object StartTime -Descending | Select-Object -First 1
  if (-not $process) { throw 'Could not resolve the temporary Word window handle.' }
  $hwnd = [IntPtr]$process.MainWindowHandle
  for ($i = 1; $i -le $pages; $i++) {
    $range = $doc.GoTo(1, 1, $i)
    $window.ScrollIntoView($range, $true)
    Start-Sleep -Milliseconds 700
    $rect = New-Object WindowCapture3+RECT
    [WindowCapture3]::GetWindowRect($hwnd, [ref]$rect) | Out-Null
    $width = $rect.Right - $rect.Left
    $height = $rect.Bottom - $rect.Top
    $bitmap = New-Object System.Drawing.Bitmap($width, $height)
    $graphics = [System.Drawing.Graphics]::FromImage($bitmap)
    $hdc = $graphics.GetHdc()
    [WindowCapture3]::PrintWindow($hwnd, $hdc, 2) | Out-Null
    $graphics.ReleaseHdc($hdc)
    $graphics.Dispose()
    $file = Join-Path $OutputDir ('page-{0:D2}.png' -f $i)
    $bitmap.Save($file, [System.Drawing.Imaging.ImageFormat]::Png)
    $bitmap.Dispose()
    Write-Output $file
  }
  $doc.Close($false)
}
finally { $word.Quit() }
