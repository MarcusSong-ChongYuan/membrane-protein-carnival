param(
    [Parameter(Mandatory=$true)][string]$Python,
    [switch]$SkipScaffoldRecalculation
)

$ErrorActionPreference = "Stop"
$Release = Split-Path -Parent $PSScriptRoot
$Baseline = Join-Path $Release "01_baseline_v6_2"
$Sources = Join-Path $Release "11_source_snapshots"
$Scripts = Join-Path $Release "10_scripts"
$Work = Join-Path $Release "reproduced_v6_3"
$DiseaseOut = Join-Path $Work "disease"
$ProteinOut = Join-Path $Work "protein_annotation"
$Qa = Join-Path $Work "qa"
New-Item -ItemType Directory -Force -Path $DiseaseOut,$ProteinOut,$Qa | Out-Null

# Disease identity, hierarchy, therapeutic area and anatomy.
& $Python (Join-Path $Scripts "disease\build_disease_module.py") `
  --disease-relations (Join-Path $Baseline "protein_gene_disease_relations_v6_2.tsv") `
  --protein-master (Join-Path $Baseline "human_membrane_protein_master_v6_2.tsv") `
  --mondo-json (Join-Path $Sources "disease\mondo\mondo.json") `
  --mondo-release-json (Join-Path $Sources "disease\mondo\github_latest_release.json") `
  --doid-json (Join-Path $Sources "disease\doid\doid.json") `
  --uberon-obo (Join-Path $Sources "disease\uberon\uberon-basic.obo") `
  --opentargets-parquet (Join-Path $Sources "disease\opentargets\opentargets_26.06_disease.parquet") `
  --output-dir $DiseaseOut --qa-dir $Qa
& $Python (Join-Path $Scripts "disease\finalize_disease_module_v2.py") --output-dir $DiseaseOut --qa-dir $Qa

# Full 10,997-protein annotation and five-axis bridge.
& $Python (Join-Path $Scripts "protein\run_protein_cross_classification_v63.py") `
  --protein-master (Join-Path $Baseline "human_membrane_protein_master_v6_2.tsv") `
  --uniprot (Join-Path $Sources "protein\uniprot\human_reviewed_annotations.tsv.gz") `
  --uniprot-headers (Join-Path $Sources "protein\uniprot\human_reviewed_annotations.headers.txt") `
  --go-obo (Join-Path $Sources "protein\go\go-basic.obo") `
  --interpro-list (Join-Path $Sources "protein\uniprot\interpro_entry.list") `
  --pfam-clans (Join-Path $Sources "protein\uniprot\Pfam-A.clans.tsv.gz") `
  --reactome-pathways (Join-Path $Sources "protein\reactome\ReactomePathways.txt") `
  --reactome-relations (Join-Path $Sources "protein\reactome\ReactomePathwaysRelation.txt") `
  --reactome-version (Join-Path $Sources "protein\reactome\database_version.txt") `
  --output-dir $ProteinOut --qa-dir $Qa

# Evidence lineage and membership-preserving docking rerank.
& $Python (Join-Path $Scripts "protein\build_evidence_lineage_v63.py") `
  --evidence (Join-Path $Baseline "binding_evidence_master_v6_2.tsv.gz") `
  --sites (Join-Path $Baseline "binding_site_instances_v6_2.tsv.gz") `
  --pair-summary (Join-Path $Baseline "protein_compound_summary_v6_2.tsv.gz") `
  --docking-full (Join-Path $Sources "docking\docking_full_nonconflict_v2_v62.tsv.gz") `
  --docking-conflict (Join-Path $Sources "docking\docking_conflict_review_v2_v62.tsv.gz") `
  --protein-classification (Join-Path $ProteinOut "protein_cross_classification_summary_v0_1.tsv") `
  --output-dir $ProteinOut --qa-dir $Qa --work-db (Join-Path $Work "evidence_lineage.duckdb")

if (-not $SkipScaffoldRecalculation) {
  & $Python (Join-Path $Scripts "compound\run_m5e_scaffold_diversity_v2.py") `
    --compound-master (Join-Path $Baseline "small_molecule_master_v1_3.tsv") `
    --output-dir (Join-Path $Work "compound") --figure-dir (Join-Path $Work "figures") --qa-dir $Qa
}

$required = @(
  "DISEASE_MODULE_V0_2_VALIDATION.json",
  "PROTEIN_CROSS_CLASSIFICATION_V0_1_VALIDATION.json",
  "EVIDENCE_LINEAGE_V0_1_VALIDATION.json"
)
foreach ($file in $required) {
  $report = Get-Content -Raw -LiteralPath (Join-Path $Qa $file) | ConvertFrom-Json
  if ($report.status -ne "PASS") { throw "Validation failed: $file" }
}
Write-Output "MemPro V6.3 incremental reproduction: PASS"
