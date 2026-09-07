# Sequence policy v5.3

## Scope

The canonical sequence for each `target_uniprot_id` was retrieved from the
official UniProtKB REST API on 2026-07-24. Isoform-specific sequences are not
included in this release.

## Provenance

- UniProtKB release: 2026_02
- UniProtKB release date: 10-June-2026
- Retrieval endpoint: https://rest.uniprot.org/uniprotkb/search
- Returned fields: accession, entry name, reviewed status, sequence, length,
  and sequence version

## Integrity

`sequence_sha256_v53` is a SHA-256 digest calculated locally from the unwrapped
uppercase amino-acid sequence. `sequence_length_match_v53` compares the
retrieved sequence length with the length already stored in v5.2.

## Excel limitation

Excel permits at most 32,767 characters in one cell. The workbook therefore
splits sequences after 32,000 residues. Concatenating the two sequence-part
columns reconstructs the exact canonical sequence. TSV and FASTA outputs
always contain the unsplit full sequence.
