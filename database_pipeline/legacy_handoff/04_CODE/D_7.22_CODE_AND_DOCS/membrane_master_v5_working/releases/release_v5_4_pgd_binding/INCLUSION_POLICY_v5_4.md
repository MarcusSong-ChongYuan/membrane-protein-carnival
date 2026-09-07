# v5.4 Protein–Disease and Small-Molecule Binding Inclusion Policy

Release date: 2026-07-24

## Protein universe

The protein universe is unchanged from v5.3: 10,997 audited canonical human
UniProt accessions. Disease and binding annotations never promote an E3
candidate or E0 excluded protein into the default membrane-protein website
release. Default website inclusion remains restricted to the v5.3 E1/E2 layer.

## Disease relations

The v3 Rule-B Open Targets relation table is retained at its original
protein–Ensembl-gene–disease granularity. Current UniProtKB disease comments are
added as a separate source layer and parsed only when UniProt explicitly labels
an `Involvement in disease` statement. OMIM identifiers and cited PubMed records
are retained when present.

Open Targets MONDO and UniProt/OMIM disease identifiers are not automatically
merged by name. Cross-ontology duplicates therefore remain separate relations
until an explicit MONDO–OMIM crosswalk is applied.

Disease evidence levels are separate from membrane evidence and binding
evidence:

- `medium`: supported v3 Rule-B association or curated UniProt disease statement.
- `high`: stronger v3 evidence or a UniProt causal-variant statement with a
  publication.
- `very_high`: v3 Rule-E evidence.

## Binding evidence

One row represents one source observation for one protein–compound context.

- `BE1`: compound-specific experimental structure or experimentally resolved
  binding-site evidence.
- `BE2`: direct quantitative binding evidence, including standardized `Ki` or
  `Kd`.
- `BE3`: functional/bioactivity evidence or a curated pharmacological target
  assertion that does not by itself prove a physical binding site.
- `BP`: predicted pocket. BP records are kept outside the experimental binding
  count and do not establish a protein–compound interaction.
- `BX`: excluded or unresolved evidence retained only for audit.

IC50 and EC50 are not promoted to BE2 solely because they are quantitative.
ChEMBL supplementation is restricted to single-protein binding assays with
standardized nM Ki/Kd values and `standard_flag=1`. BindingDB complex targets are
retained but excluded from default per-protein claims unless the target is a
single unambiguous protein.

Only defined, non-polymeric chemical entities contribute to the default
small-molecule count. Antibodies, proteins, nucleic acids, undefined polymers,
buffers, solvents, water and undefined mixtures are excluded. Peptide-like or
very large structures remain review records.

