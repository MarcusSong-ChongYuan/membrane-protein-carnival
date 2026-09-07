# Methods: Small-Molecule Database v1.0

## Inputs

The input universe was `small_molecule_index_v4_1.tsv` plus five PubChem CIDs
that occurred in eight GtoPdb evidence rows but were missing from the v4.1
compound index. All 665,214 source-level binding observations were retained.

## Enrichment

- PubChem PUG REST supplied structures and identifiers for 6,157 missing
  PubChem records and the five repaired index-gap records.
- ChEMBL 37 official chemical representations supplied structures for
  ChEMBL source records. ChEMBL REST metadata supplied molecule type, parent
  hierarchy, preferred name, clinical phase, first approval, natural-product
  and probe flags.
- ChEBI July 2026 SDF and ontology files supplied stable ChEBI identifiers,
  nomenclature, direct chemical classes and `has role` relationships.
- DrugCentral, GtoPdb, BindingDB and legacy v3.1 identifiers and names were
  retained as source cross-references.

## Standardization

RDKit 2026.03.4 was used with the frozen
`MemProDB-compound-standardization-v1` policy:

1. parse the supplied SMILES;
2. perform RDKit cleanup and valence normalization;
3. preserve the standardized exact form;
4. derive a fragment parent for salts and multicomponent forms;
5. normalize parent charge where possible;
6. generate canonical isomeric SMILES, Standard InChI and InChIKey;
7. calculate formula, molecular weight, exact mass, XlogP, TPSA, hydrogen-bond
   counts, rotatable bonds, formal charge, heavy atoms and ring descriptors.

Original source structures are retained. Structure transformations and
failures are recorded per source record.

## Merge policy

- Same full Standard InChIKey: grouped, with multi-SMILES cases flagged.
- Same connectivity but different full InChIKey: kept separate.
- Salts and multicomponent forms: retained as distinct forms and linked to a
  parent compound.
- Same name but different structures: kept separate.
- Name-only records: not automatically merged.
- ChEMBL-parent disagreement: the reproducible structure-derived parent is
  retained and the discrepancy is audited.

## Binding remapping

`binding_evidence_master_v4_2.tsv` adds canonical parent and exact-form
identifiers. Default v4.2 website evidence requires:

1. the original evidence row to be default-included;
2. an unambiguous canonical compound mapping;
3. core small-molecule scope.

Extended, unresolved and excluded observations remain available for audit but
do not contribute to default website summaries.

