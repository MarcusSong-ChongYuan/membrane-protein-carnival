from __future__ import annotations

import HuMemLigDB_build_docx as base


_original_normalize = base.normalize_source


def normalize_source_v2(text: str) -> str:
    text = _original_normalize(text)
    text = text.replace(
        "The frozen V6.2 release contains 2,016,064 canonical compounds, 3,003,306 positive evidence records summarized into 1,502,456 protein-compound pairs, 95,598 binding-site instances, 10,264 protein-gene-disease relations and 7,784,223 mapped assay-inactive records.",
        "The frozen V6.2 release contains 2,016,064 canonical compounds, 3,003,306 positive evidence records summarized into 1,502,456 auditable protein-compound pairs, including a recommended subset of 1,305,791 pairs involving default E1-E2 proteins, together with 95,598 binding-site instances, 10,264 protein-gene-disease relations and 7,784,223 mapped assay-inactive records.",
    )
    text = text.replace(
        "The current audit universe contains 10,997 protein records, of which 7,904 satisfy the default E1-E2 high-confidence membrane criteria. The interaction layer contains more than 1.5 million unique protein-canonical compound pairs across 4,630 proteins, including 95,598 structure-resolved site records.",
        "The current audit universe contains 10,997 protein records, of which 7,904 satisfy the default E1-E2 high-confidence membrane criteria. The full interaction audit layer contains more than 1.5 million unique protein-canonical compound pairs across 4,630 proteins; applying the default protein criterion yields 1,305,791 pairs across 3,503 E1-E2 proteins. The release also includes 95,598 structure-resolved site records.",
    )
    text = text.replace(
        "Most pairs are supported by one contributing database label, while 295,207 pairs have two or more source labels; this field is retained as provenance support rather than interpreted as experimentally independent replication.",
        "Most pairs are supported by one contributing database label, while 295,207 pairs have two or more source labels; this field is retained as provenance support rather than interpreted as experimentally independent replication. Applying both the default protein layer and the default evidence policy yields 2,286,514 evidence records and 1,305,791 pairs across 3,503 E1-E2 proteins and 803,613 compounds; these values define the recommended core interaction subset.",
    )
    text = text.replace(
        "| Positive protein-compound pairs | 1,502,456 | Unique canonical protein-parent pairs |",
        "| Positive protein-compound pairs | 1,502,456 | Full auditable canonical protein-parent pair layer |\n| Core positive pairs | 1,305,791 | Default evidence on E1-E2 protein records |",
    )
    text = text.replace(
        "Third-party raw archives are not redistributed in the internal handoff package. Clean rebuilding requires reacquisition of the recorded source versions and configuration of machine-specific input paths.",
        "Third-party raw archives are not redistributed in the internal handoff package. Clean rebuilding requires reacquisition of the recorded source versions and configuration of machine-specific input paths. Several inherited sources, including BioLiP, sc-PDB, PDSP KiDatabase and DrugCentral, remain identifiable in the evidence table but require completion of their exact frozen-version metadata before submission.",
    )
    return text


base.normalize_source = normalize_source_v2


if __name__ == "__main__":
    base.main()
