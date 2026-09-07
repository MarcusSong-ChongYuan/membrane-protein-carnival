# Methods

- V6.1 was treated as immutable.
- HPA RNA, IHC, MS and localization modalities retain their original units;
  missing values and measured zero are distinct.
- Automatic chemical identity merging requires a unique full InChIKey.
- Unassigned stereochemistry and uncurated salt/charge-parent relationships
  remain in review.
- Negative evidence is deduplicated by target UniProt, canonical compound,
  PubChem AID and SID, then compared with positive evidence at target-compound
  level.
- PDB asymmetric-unit chain counts are not used as physiological stoichiometry.
  Heteromer target copy number remains unknown without component-level mapping.
