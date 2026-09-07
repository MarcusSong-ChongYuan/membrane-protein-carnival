# MemPro subunit and biological-assembly add-on

This add-on is built independently from the frozen V6.0 release. It must be
validated before its columns are merged into a future protein-master release.

## Evidence policy

- UniProtKB `Subunit structure` is retained verbatim. Oligomeric states are
  parsed conservatively and PubMed/ECO identifiers remain traceable.
- `ECO:0000269` in the UniProt comment is treated as direct experimental
  annotation.
- PDBe biological assemblies come from experimentally determined PDB entries.
  The assembly assignment is separately labelled `author_defined`,
  `software_defined`, or `unspecified`.
- A PDB biological assembly is **not** automatically labelled as directly
  experimentally confirmed physiological stoichiometry.
- Crystallographic chain count is never used on its own as the biological
  subunit count.
- Multiple observed states are retained. They are flagged rather than forced
  into one value.

## Main outputs

### `membrane_protein_subunit_summary_v1.tsv.gz`

One row per MemPro protein. Main fields:

- `subunit_state_summary`: monomer, homo/heterodimer, trimer, tetramer, etc.
- `subunit_count_values`, `subunit_count_min`, `subunit_count_max`
- `subunit_composition_summary`
- `biological_assembly_ids`: `PDB_ID:assembly_id`
- UniProt and PDBe evidence flags and identifiers
- `subunit_state_conflict_flag` and `subunit_state_conflict_reason`

### `membrane_protein_biological_assembly_v1.tsv.gz`

One row per MemPro protein–PDB–biological assembly. It retains assembly ID,
assignment class, oligomeric state, protein-subunit count, entity composition,
and source endpoint.

### `uniprot_subunit_annotations_v1.tsv.gz`

One row per MemPro protein containing the raw UniProt subunit comment, parsed
state/count values, PubMed IDs, ECO codes, and direct-experimental flag.
