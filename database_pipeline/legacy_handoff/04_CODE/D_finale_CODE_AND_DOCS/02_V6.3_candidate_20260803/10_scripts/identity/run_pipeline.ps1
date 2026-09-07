param(
    [Parameter(Mandatory=$true)][string]$Root,
    [Parameter(Mandatory=$true)][string]$V62,
    [Parameter(Mandatory=$true)][string]$ComplexRoot,
    [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
$raw = Join-Path $Root "raw\uniprot"
$out = Join-Path $Root "output_v0_1"
$qa = Join-Path $Root "qa"
$scripts = Join-Path $Root "scripts"
New-Item -ItemType Directory -Force -Path $raw, $out, $qa | Out-Null

$bulkFasta = Join-Path $raw "UP000005640_human_with_isoforms.fasta.gz"
$headers = Join-Path $raw "response_headers.txt"
if (-not (Test-Path -LiteralPath $bulkFasta)) {
    curl.exe --fail --location --retry 5 --retry-all-errors `
      --dump-header $headers --output $bulkFasta `
      "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&includeIsoform=true&query=%28proteome%3AUP000005640%29"
}

& $Python (Join-Path $scripts "fetch_explicit_isoforms.py") `
  --complex-components (Join-Path $ComplexRoot "output_v0_2\complex_target_components_v0_2.tsv") `
  --output-fasta (Join-Path $raw "complex_explicit_isoforms.fasta") `
  --output-manifest (Join-Path $raw "complex_explicit_isoforms_manifest.json")

& $Python (Join-Path $scripts "run_build_identity_fixed.py") `
  --protein-master (Join-Path $V62 "human_membrane_protein_master_v6_2.tsv") `
  --binding-evidence (Join-Path $V62 "binding_evidence_master_v6_2.tsv.gz") `
  --binding-sites (Join-Path $V62 "binding_site_instances_v6_2.tsv.gz") `
  --disease-relations (Join-Path $V62 "protein_gene_disease_relations_v6_2.tsv") `
  --complex-components (Join-Path $ComplexRoot "output_v0_2\complex_target_components_v0_2.tsv") `
  --uniprot-fasta $bulkFasta --uniprot-headers $headers --output-dir $out

& $Python (Join-Path $scripts "run_extend_complex_identity_v2_fixed.py") `
  --protein-master (Join-Path $V62 "human_membrane_protein_master_v6_2.tsv") `
  --gene-entities (Join-Path $out "gene_entity_v0_1.tsv") `
  --canonical-entities (Join-Path $out "canonical_protein_entity_v0_1.tsv.gz") `
  --core-isoforms (Join-Path $out "protein_isoform_v0_1.tsv.gz") `
  --binding-evidence (Join-Path $V62 "binding_evidence_master_v6_2.tsv.gz") `
  --binding-resolution (Join-Path $out "binding_evidence_target_resolution_v0_1.tsv.gz") `
  --complex-components (Join-Path $ComplexRoot "output_v0_2\complex_target_components_v0_2.tsv") `
  --complex-master (Join-Path $ComplexRoot "output_v0_2\complex_target_master_v0_2.tsv") `
  --uniprot-fasta $bulkFasta `
  --explicit-isoform-fasta (Join-Path $raw "complex_explicit_isoforms.fasta") `
  --explicit-isoform-manifest (Join-Path $raw "complex_explicit_isoforms_manifest.json") `
  --output-dir $out --qa-dir $qa

& $Python (Join-Path $scripts "fetch_missing_external_canonical.py") `
  --external-entities (Join-Path $out "external_complex_protein_entity_v0_1.tsv.gz") `
  --output-fasta (Join-Path $raw "external_canonical_direct.fasta") `
  --output-manifest (Join-Path $raw "external_canonical_direct_manifest.json")

& $Python (Join-Path $scripts "finalize_identity_layer_v3.py") `
  --root $Root `
  --direct-fasta (Join-Path $raw "external_canonical_direct.fasta") `
  --direct-manifest (Join-Path $raw "external_canonical_direct_manifest.json")

$final = Get-Content -Raw -LiteralPath (Join-Path $qa "IDENTITY_LAYER_V0_3_VALIDATION.json") | ConvertFrom-Json
if ($final.status -ne "PASS") { throw "Final identity-layer validation failed." }
Write-Output "V6.3 candidate identity layer: PASS"
