# Human Membrane-Protein Small-Molecule Database v1.0

This release is the canonical compound layer derived from the v4.1
small-molecule–membrane-protein evidence database.

## Release contents

- `small_molecule_master_v1_0.tsv`: 277,837 structure-normalized parent
  compounds.
- `compound_form_hierarchy_v1_0.tsv`: 279,870 exact forms, including 8,856
  salt or multicomponent forms.
- `compound_source_xref_v1_0.tsv`: all 302,884 source records and their
  mapping decisions.
- `compound_synonyms_v1_0.tsv`: preferred names and source-preserved
  synonyms.
- `protein_compound_summary_v4_2.tsv`: 416,766 default website
  protein–compound pairs.
- `compound_protein_summary_v1_0.tsv`: 243,899 compounds with default
  membrane-protein evidence.
- `binding_evidence_master_v4_2.tsv`: all 665,214 v4.1 evidence observations
  with canonical compound and form identifiers added.
- `human_membrane_audit_master_v5_5.tsv`: the v5.4 protein master with
  canonical compound counts appended.

## Compound scope

The website core contains 268,897 parent compounds. The extended chemical
layer contains 8,940 parent compounds, principally peptidic,
oligosaccharide, inorganic, or metal-containing entities. There are 518
unresolved source records and 3,399 excluded records. Excluded and unresolved
records remain in explicit audit tables and are not silently deleted.

## Identity policy

Source names are never used as the sole automatic merge key. Standard InChI
identity is used for structure grouping; stereo-specific full InChIKeys remain
separate. Salt and multicomponent forms are retained in the form table and
linked to a structure-derived parent. Connectivity-only families are provided
for search and review, not automatic merging.

Standard InChI can normalize tautomeric and protonation representations.
Groups with multiple standardized SMILES are retained with explicit QC flags
and source structures remain traceable in the cross-reference table.

