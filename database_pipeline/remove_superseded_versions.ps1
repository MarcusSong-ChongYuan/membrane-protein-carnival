$ErrorActionPreference = 'Stop'

$formal = 'D:\finale\FORMAL'
$qa = Get-Content -LiteralPath (Join-Path $formal '07_QA\FORMAL_VALIDATION_REPORT.json') -Raw | ConvertFrom-Json
$negativeQa = Get-Content -LiteralPath (Join-Path $formal '07_QA\NEGATIVE_EVIDENCE_ARCHIVE_QA.json') -Raw | ConvertFrom-Json
if ($qa.status -ne 'PASS') { throw 'FORMAL QA is not PASS; deletion blocked.' }
if ($negativeQa.unique_content_files -lt 1) { throw 'Negative-evidence archive is empty; deletion blocked.' }

$targets = @(
    'D:\7.22',
    'D:\finale\00_交接说明',
    'D:\finale\01_正式数据_V6.2',
    'D:\finale\02_展示图表_V6.2',
    'D:\finale\02_V6.3_candidate_20260803',
    'D:\finale\03_Docking名单_V6.2',
    'D:\finale\03_V6.3.1_readiness_20260803',
    'D:\finale\04_复现代码',
    'D:\finale\04_V6.3.1_isoform_accession_audit_20260803',
    'D:\finale\05_质量控制',
    'D:\finale\05_V6.3.1_双页模块PPT重制_20260803',
    'D:\finale\06_清单与哈希',
    'D:\finale\07_M1-M8图表修订_20260806',
    'D:\finale\08_V6.4_database_freeze_and_figures_20260812',
    'D:\finale\09_V6.4.1_final_publication_triage_20260812',
    'D:\finale\09_V7_data_freeze_working_20260813',
    'D:\finale\10_MemPro_V7_data_freeze_20260813',
    'D:\finale\10_MemPro_V7_public_candidate_20260812',
    'D:\finale\11_MemPro_V7_automatic_data_release_20260813',
    'D:\finale\12_V7_final_completion_20260813',
    'D:\finale\13_MemPro_V7_final_data_release_20260813',
    'D:\finale\14_MemPro_V7.0.1_rescue_candidate_20260814',
    'D:\finale\15_MemPro_V7.0.1_rescue_final_20260814',
    'D:\finale\16_MemPro_V7.0.2_final_20260814',
    'D:\finale\17_MemPro_V7.1_candidate_20260814',
    'D:\finale\18_MemPro_V7.1_final_20260814',
    'D:\finale\20_MemPro_V7.1.1_final_20260814',
    'D:\finale\21_MemPro_V7.1.1_figures_20260814',
    'D:\finale\21_MemPro_V7.2_candidate_20260817',
    'D:\finale\22_MemPro_V7.2_final_20260817',
    'D:\finale\23_MemPro_V7.2_publication_figures_candidate_20260817',
    'D:\finale\MemPro_M1_M8_AGENT_REPRO_20260806',
    'D:\finale\MemPro_V6.2_数据库建立流程与科学意义',
    'D:\finale\_cleanup_pending_20260818',
    'C:\Users\Administrator\Desktop\FINAL',
    'C:\Users\Administrator\Desktop\MemPro_FULL_HANDOFF_20260817'
)

$allowedExact = @(
    [IO.Path]::GetFullPath('D:\7.22'),
    [IO.Path]::GetFullPath('C:\Users\Administrator\Desktop\FINAL'),
    [IO.Path]::GetFullPath('C:\Users\Administrator\Desktop\MemPro_FULL_HANDOFF_20260817')
)
$finaleRoot = [IO.Path]::GetFullPath('D:\finale') + [IO.Path]::DirectorySeparatorChar
$formalResolved = [IO.Path]::GetFullPath($formal)

$plan = [System.Collections.Generic.List[object]]::new()
foreach ($target in $targets) {
    $resolved = [IO.Path]::GetFullPath($target)
    $allowed = ($resolved.StartsWith($finaleRoot, [StringComparison]::OrdinalIgnoreCase) -and
                -not $resolved.StartsWith($formalResolved, [StringComparison]::OrdinalIgnoreCase)) -or
               ($allowedExact -contains $resolved)
    if (-not $allowed) { throw "Unsafe deletion target rejected: $resolved" }
    $exists = Test-Path -LiteralPath $resolved
    $bytes = if ($exists) {
        (Get-ChildItem -LiteralPath $resolved -Recurse -File -ErrorAction SilentlyContinue |
            Measure-Object Length -Sum).Sum
    } else { 0 }
    $plan.Add([pscustomobject]@{target=$resolved;exists=[int]$exists;bytes=$bytes;allowed=[int]$allowed})
}
$planPath = Join-Path $formal '05_metadata\SUPERSEDED_VERSION_DELETION_PLAN.tsv'
$plan | Export-Csv -LiteralPath $planPath -Delimiter "`t" -NoTypeInformation -Encoding utf8

foreach ($row in $plan | Where-Object exists -eq 1) {
    Remove-Item -LiteralPath $row.target -Recurse -Force
}

$result = [ordered]@{
    status = 'PASS'
    formal_qa = $qa.status
    negative_unique_content_files_preserved = $negativeQa.unique_content_files
    targets_planned = $plan.Count
    targets_removed = @($plan | Where-Object exists -eq 1).Count
    bytes_removed = ($plan | Measure-Object bytes -Sum).Sum
    deletion_plan = $planPath
}
$result | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $formal '07_QA\SUPERSEDED_VERSION_DELETION_REPORT.json') -Encoding utf8
$result | ConvertTo-Json
