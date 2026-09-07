# Revised figure captions — trial version 1

## Revised Figure M3D. HPA tissue-expression breadth with explicit mapping status

Expression breadth is shown only after separating HPA mapping status. The violin,
box and deterministic point sample use proteins with `hpa_mapping_status_v62=mapped`
and a numeric detected-tissue count; mapped zeros remain valid plotted observations.
The adjacent 100% bars separate mapped proteins with at least one detected tissue,
mapped proteins with zero tissues passing the HPA detection rule, and
unmapped/missing proteins. A mapped zero means that no tissue passed the chosen HPA
detection rule; it is not evidence that the protein is biologically absent.

## Revised Figure M4. Ontology-derived, multi-label disease landscape

Disease classification is taken from the Open Targets Platform 26.06
`therapeuticAreas` field rather than disease-name keywords. Open Targets entity IDs
are joined directly; UniProt/OMIM records are joined only through official Open
Targets cross-references. Ambiguous cross-references are retained as ambiguous and
are not used for automatic canonical disease merging. Therapeutic-area membership
is multi-label and non-exclusive. Raw membership counts can therefore exceed the
number of unique protein–disease pairs; diamond markers show fractional allocation,
where each pair contributes a total weight of one across its assigned areas.
Generic phenotype, biological-process, measurement, medical-procedure and non-human
branches are not treated as human organ systems. This figure represents clinical
therapeutic areas, not a DO/Uberon anatomical map.

## Revised Figure M5B. Biological-status intersections among annotated compounds

The UpSet plot includes core, QC-passed canonical compounds carrying at least one
of five status annotations: approved drug, clinical candidate, chemical probe,
endogenous ligand or natural product. Compounds with no status annotation are
excluded from the intersection denominator. Set and intersection axes use a log
scale. The approved-drug set size reflects the conservative annotation present in
the frozen HuMemLigDB release and is not a count of all approved drugs worldwide.
The coverage strip separates compounds with at least one positive status,
compounds whose five status fields were assessed but all negative, and compounds
for which the status fields are unavailable or incomplete.


## Revised Figure M5D. Physicochemical distributions of canonical compounds

Empirical cumulative distribution functions (ECDFs) summarize molecular weight,
XlogP, topological polar surface area (TPSA), and rotatable-bond count in the
deterministic V6.2 physicochemical cache. Each curve gives the fraction of observed
compounds at or below a given property value; dashed lines mark medians and the
annotations report the 10th, 50th and 90th percentiles. Display limits are applied
only for legibility and do not remove compounds from HuMemLigDB. These properties
describe size, lipophilicity, polarity and conformational flexibility. They are
not efficacy, approval, membrane permeability or universal docking criteria.
Lipinski Rule-of-Five counts are not used as the main chemical-space panel because
they were designed primarily for oral drug-likeness.
