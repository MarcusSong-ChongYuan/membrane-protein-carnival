#!/usr/bin/env python3
"""Compatibility runner for the V0.2 extension builder."""

from __future__ import annotations

import extend_complex_identity_v2 as module


original_reader = module.reader


def compatible_reader(path):
    for row in original_reader(path):
        if "canonical_uniprot_accession" in row and "target_uniprot_id" not in row:
            row["target_uniprot_id"] = row["canonical_uniprot_accession"]
        if path.name == "protein_isoform_v0_1.tsv.gz":
            row.setdefault("entity_scope", "membrane_protein_master")
            row.setdefault("isoform_record_type", "bulk_reference_proteome_isoform")
        yield row


module.reader = compatible_reader
module.main()
