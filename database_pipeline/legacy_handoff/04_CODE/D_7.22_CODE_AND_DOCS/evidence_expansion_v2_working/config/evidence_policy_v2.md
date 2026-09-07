# Evidence policy for the expansion workflow

## Positive evidence

- `BE1`: experimentally observed protein-small-molecule complex or explicit
  residue-level binding site.
- `BE2`: direct quantitative binding to a single mapped protein, including Kd
  and Ki.
- `BE3`: target-specific functional pharmacology or enzyme-ligand relation,
  including IC50, EC50, substrate, product, cofactor, activator and inhibitor.
- `MECH`: curated drug-target mechanism without a source-level quantitative
  binding experiment.

## Non-default evidence

- Target complexes without a unique human UniProt mapping.
- Cell-based or organismal assays without a unique target.
- Primary high-throughput screens without confirmatory or dose-response support.
- Cross-species evidence without a direct human experiment.
- Structure-only chemical mentions without protein-ligand contact evidence.
- Prediction-only relations.

## Negative evidence

Inactive and explicit non-binding observations are preserved in a separate
negative-evidence table. They never count as positive protein-compound
relations.

## Compound policy

Canonical parent and exact-form identifiers follow
`MemProDB-compound-standardization-v1`. Names alone never trigger an automatic
merge. Salts and multicomponent records remain linked exact forms.

