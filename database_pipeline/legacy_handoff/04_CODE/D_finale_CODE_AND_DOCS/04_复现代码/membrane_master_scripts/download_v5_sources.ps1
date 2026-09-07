param(
    [switch]$Force
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$raw = Join-Path $root "raw"
New-Item -ItemType Directory -Force -Path $raw | Out-Null

$sources = @(
    @{
        id = "hpa_v25_1"
        url = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
        file = "hpa_proteinatlas_v25_1.tsv.zip"
        evidence = "experimental_subcellular_localization_and_hpa_predicted_membrane_class"
    },
    @{
        id = "opm_current"
        url = "https://opm-back.cc.lehigh.edu/opm-backend/primary_structures?fileFormat=csv"
        file = "opm_primary_structures_2026-07-24.csv"
        evidence = "experimentally_determined_membrane_structure_orientation"
    },
    @{
        id = "membranome_current"
        url = "https://opm-back.cc.lehigh.edu/membranome-backend/proteins?fileFormat=csv"
        file = "membranome_proteins_2026-07-24.csv"
        evidence = "curated_single_pass_membrane_protein_classification"
    },
    @{
        id = "unitmp_htp_d2_2"
        url = "https://htp.unitmp.org/data/HTP/data/d.2.2/sets/htp_all.xml"
        file = "unitmp_htp_all_d2_2.xml"
        evidence = "integrated_human_transmembrane_topology"
    },
    @{
        id = "unitmp_pdbtm_current"
        url = "https://pdbtm.unitmp.org/data/PDBTM/data/pdbtm_all.xml"
        file = "unitmp_pdbtm_all_2026-07-24.xml"
        evidence = "experimentally_determined_membrane_structure_orientation"
    },
    @{
        id = "tcdb_human"
        url = "https://www.tcdb.org/public/human.csv"
        file = "tcdb_human_2026-07-24.csv"
        evidence = "curated_transport_classification"
    },
    @{
        id = "gpcrdb_receptorlist_current"
        url = "https://gpcrdb.org/services/receptorlist/"
        file = "gpcrdb_receptorlist_2026-07-24.json"
        evidence = "curated_gpcr_class_family_and_subfamily"
    },
    @{
        id = "gpcrdb_proteinfamily_current"
        url = "https://gpcrdb.org/services/proteinfamily/"
        file = "gpcrdb_proteinfamily_2026-07-24.json"
        evidence = "curated_gpcr_family_hierarchy"
    },
    @{
        id = "opm_openapi"
        url = "https://opm.phar.umich.edu/swagger-opm-api.json"
        file = "opm_openapi_2026-07-24.json"
        evidence = "api_schema"
    },
    @{
        id = "membranome_openapi"
        url = "https://biomembhub.org/membranome/swagger-membranome-api.json"
        file = "membranome_openapi_2026-07-24.json"
        evidence = "api_schema"
    }
)

$manifest = @()
foreach ($source in $sources) {
    $destination = Join-Path $raw $source.file
    if ($Force -or -not (Test-Path -LiteralPath $destination)) {
        & curl.exe -L --fail --retry 3 --retry-delay 2 --max-time 900 `
            -o $destination $source.url
        if ($LASTEXITCODE -ne 0) {
            throw "Download failed: $($source.url)"
        }
    }

    $item = Get-Item -LiteralPath $destination
    $hash = (Get-FileHash -LiteralPath $destination -Algorithm SHA256).Hash.ToLowerInvariant()
    $manifest += [ordered]@{
        source_id = $source.id
        source_url = $source.url
        local_file = $source.file
        evidence_role = $source.evidence
        retrieved_at = (Get-Date).ToString("yyyy-MM-ddTHH:mm:ssK")
        bytes = $item.Length
        sha256 = $hash
    }
}

$manifestPath = Join-Path $raw "source_download_manifest.json"
$manifest | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $manifestPath -Encoding UTF8
Write-Output "Downloaded/verified $($manifest.Count) official source artifacts."
Write-Output $manifestPath
