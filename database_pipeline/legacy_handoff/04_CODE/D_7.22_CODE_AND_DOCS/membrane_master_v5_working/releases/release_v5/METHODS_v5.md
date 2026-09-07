# Methods

## Record unit and species

The public record unit is a canonical UniProtKB accession from the Homo sapiens reference proteome UP000005640. Alternative/unreviewed accessions are not silently merged: every TrEMBL candidate has a crosswalk row that records the identifiers used and whether the match is unique, ambiguous or absent.

## Evidence hierarchy

Evidence is stored source by source. The decision priority is:

1. Experimentally determined membrane structure/orientation (PDBTM or accession-mapped OPM integral segment).
2. Curated UniProt transmembrane/intramembrane annotation.
3. Specialist curated database classification (notably Membranome).
4. Integrated or multiple-source topology evidence (UniTmp HTP).
5. Prediction-only membrane evidence.
6. Subcellular localization evidence.
7. Name/family inference.

Location evidence is never treated as proof of membrane insertion. TCDB and pharmacology target membership are used for classification, not as sole membrane-membership evidence.

## Release tiers

- Tier A: integral membrane proteins supported by the v4 UniProt layer or strong specialist/structural evidence.
- Tier B: integral monotopic or lipid-anchored proteins.
- Tier C: peripheral membrane-associated proteins.
- Tier R: records with prediction-only, location-only, non-integral OPM association, or unresolved conflicts.
- Tier E: excluded/non-membrane decisions are retained in audit crosswalks rather than the published membrane union.

## Single-pass topology

For the 680 v4 `single_pass_unresolved_type` rows, UniTmp HTP terminal orientation was applied only when HTP reported one transmembrane segment. N-in/C-out supports type II. N-out/C-in plus a UniProt signal peptide supports type I; N-out/C-in with an N-terminal signal anchor supports type III. N-out/C-in without enough evidence to distinguish type I from III is labelled `single_pass_type_i_or_iii`, not forced.

## Functional classification

GPCRdb has highest precedence for GPCR hierarchy. GtoPdb supplies expert receptor, ion-channel, transporter and other target families. TCDB supplies TC identifiers for transporters. Existing v4 enzyme/receptor/channel assignments are retained when no specialist source overrides them. Conservative UniProt name/family rules provide a broad browsing category; `other_family_defined` means a family exists but a mature functional primary class is still unresolved.

## Reproducibility

Raw source URLs, retrieval dates, byte sizes and SHA-256 hashes are recorded in `source_download_manifest.json` and `SOURCE_REGISTRY_v5.tsv`. The build and validation scripts are included in the release.
