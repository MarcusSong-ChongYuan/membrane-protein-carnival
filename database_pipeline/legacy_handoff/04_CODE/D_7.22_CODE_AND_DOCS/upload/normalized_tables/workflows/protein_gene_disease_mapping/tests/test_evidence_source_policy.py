from pathlib import Path

import yaml


WORKFLOW = Path(__file__).resolve().parents[1]
POLICY = yaml.safe_load((WORKFLOW / "config" / "evidence_source_policy.yaml").read_text(encoding="utf-8"))
REQUIRED = {
    "source_family", "source_level", "counts_toward_two_source_rule",
    "allowed_as_only_core_source", "rule_a_qualifies", "clinical_or_curated_genetic",
    "second_core_or_human_genetic", "high_authority_source", "source_note",
}


def test_policy_is_complete_and_conservative():
    defaults = POLICY["unknown_source_defaults"]
    assert REQUIRED <= defaults.keys()
    assert defaults["source_level"] == "unknown"
    assert defaults["counts_toward_two_source_rule"] is False
    assert defaults["rule_a_qualifies"] is False
    for datasource, entry in POLICY["source_families"].items():
        assert REQUIRED <= entry.keys(), datasource
        if entry["clinical_or_curated_genetic"]:
            assert entry["second_core_or_human_genetic"], datasource


def test_related_channels_collapse_to_one_family():
    sources = POLICY["source_families"]
    assert sources["crispr"]["source_family"] == sources["crispr_screen"]["source_family"]
    assert sources["crispr"]["source_family"] == "CRISPR"
    assert sources["eva"]["source_family"] == sources["eva_somatic"]["source_family"]
    assert sources["uniprot_literature"]["source_family"] == sources["uniprot_variants"]["source_family"]
