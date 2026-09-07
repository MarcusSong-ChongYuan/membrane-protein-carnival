import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const outputDir = process.argv[2];
if (!outputDir) throw new Error("Usage: node build_review_workbook.mjs <output-directory>");
await fs.mkdir(outputDir, { recursive: true });

const wb = Workbook.create();
const COLORS = {
  navy: "#17324D",
  teal: "#1F7A8C",
  aqua: "#D9F0F2",
  pale: "#EEF5F7",
  green: "#DDF3E4",
  greenText: "#1E6B3A",
  amber: "#FFF2CC",
  red: "#FCE4E4",
  border: "#C8D5DB",
  white: "#FFFFFF",
  ink: "#23313A",
};

function title(sheet, text, subtitle, endCol = "H") {
  sheet.showGridLines = false;
  sheet.getRange(`A1:${endCol}2`).merge();
  sheet.getRange("A1").values = [[text]];
  sheet.getRange(`A1:${endCol}2`).format = {
    fill: COLORS.navy,
    font: { color: COLORS.white, bold: true, size: 18 },
    verticalAlignment: "center",
    wrapText: true,
  };
  sheet.getRange(`A3:${endCol}3`).merge();
  sheet.getRange("A3").values = [[subtitle]];
  sheet.getRange(`A3:${endCol}3`).format = {
    fill: COLORS.aqua,
    font: { color: COLORS.ink, italic: true, size: 10 },
    verticalAlignment: "center",
    wrapText: true,
  };
  sheet.getRange("A1:A3").format.rowHeight = 25;
}

function header(range) {
  range.format = {
    fill: COLORS.teal,
    font: { color: COLORS.white, bold: true },
    borders: { preset: "all", style: "thin", color: COLORS.border },
    verticalAlignment: "center",
    wrapText: true,
  };
}

function body(range) {
  range.format = {
    font: { color: COLORS.ink, size: 10 },
    borders: { preset: "all", style: "thin", color: COLORS.border },
    verticalAlignment: "top",
    wrapText: true,
  };
}

const summary = wb.worksheets.add("Summary");
title(summary, "MemPro V6.3 Candidate Identity Layer", "Gene → canonical protein → isoform → complex | UniProt 2026_02 | Final relational QA: PASS", "M");
summary.getRange("A5:C5").values = [["Entity / scope", "Count", "Interpretation"]];
header(summary.getRange("A5:C5"));
summary.getRange("A6:C11").values = [
  ["Gene entities", 10970, "HGNC/Ensembl/NCBI Gene identity layer"],
  ["Core canonical membrane proteins", 10997, "Frozen V6.2 count; one canonical UniProt accession per record"],
  ["Protein isoforms", 20211, "Exact UniProt isoform accessions with sequence and canonical FK"],
  ["Protein complexes", 19791, "Audit entities from the existing complex module"],
  ["External complex canonical proteins", 6126, "Context-only components; excluded from membrane-protein totals"],
  ["Website four-entity search rows", 61969, "Gene + core canonical + isoform + complex"],
];
body(summary.getRange("A6:C11"));
summary.getRange("B6:B11").format.numberFormat = "#,##0";
summary.getRange("A13:C13").values = [["QA item", "Result", "Meaning"]];
header(summary.getRange("A13:C13"));
summary.getRange("A14:C18").values = [
  ["Final validation", "PASS", "All primary-key, foreign-key, row-count and target-level checks passed"],
  ["Explicit complex isoforms", "168 / 168", "All explicit isoform IDs returned valid official UniProt FASTA"],
  ["Legacy binding isoform specificity", "Not recoverable", "V6.2 normalized records contain canonical accessions only; no isoform was inferred"],
  ["External accession review", 8, "Old/deleted external component IDs retained in a review table"],
  ["Core membrane-protein count changed?", "No", "Remains 10,997"],
];
body(summary.getRange("A14:C18"));
summary.getRange("B14").format = { fill: COLORS.green, font: { color: COLORS.greenText, bold: true } };
summary.getRange("E5:F10").values = [
  ["Entity", "Count"],
  ["Genes", 10970],
  ["Canonical", 10997],
  ["Isoforms", 20211],
  ["Complexes", 19791],
  ["External proteins", 6126],
];
header(summary.getRange("E5:F5"));
body(summary.getRange("E6:F10"));
const entityChart = summary.charts.add("bar", summary.getRange("E5:F10"));
entityChart.title = "Identity-layer entity counts";
entityChart.hasLegend = false;
entityChart.setPosition("H5", "M18");
summary.freezePanes.freezeRows(5);
summary.getRange("A:C").format.columnWidth = 24;
summary.getRange("C:C").format.columnWidth = 56;
summary.getRange("E:F").format.columnWidth = 18;

const model = wb.worksheets.add("Entity Model");
title(model, "Entity model and non-merging rules", "The same gene may have multiple canonical proteins; isoforms and complexes are independent entity layers.", "G");
model.getRange("A5:G5").values = [["Layer", "Primary key pattern", "Parent / FK", "Rows", "Automatic mapping", "Never do", "Website search"]];
header(model.getRange("A5:G5"));
model.getRange("A6:G10").values = [
  ["Gene", "HGNC:* / ENSG* / NCBIGene:*", "—", 10970, "Use stable gene IDs; keep multi-ID review flags", "Do not use gene symbol alone to merge distinct records", "Yes"],
  ["Canonical protein", "UNIPROT:<accession>", "gene_entity_id", 10997, "One frozen record per canonical UniProt accession", "Do not merge proteins because they share a gene", "Yes"],
  ["Protein isoform", "UNIPROT_ISOFORM:<accession-n>", "canonical_protein_entity_id", 20211, "Only exact, explicit UniProt isoform IDs", "Do not infer isoform from canonical accession", "Yes"],
  ["Protein complex", "HMCX-*", "component assertions", 19791, "Components may link canonical or isoform entities", "Do not force a complex into a single-protein row", "Yes"],
  ["External component protein", "UNIPROT:<accession>", "complex context", 6126, "Resolve sequence where possible; retain review state", "Do not count as a core membrane protein", "Via complex context"],
];
body(model.getRange("A6:G10"));
model.getRange("D6:D10").format.numberFormat = "#,##0";
model.getRange("A12:G12").merge();
model.getRange("A12").values = [["One gene → multiple canonical accessions is a valid one-to-many relation. The 10,997 canonical protein records are preserved even though they map to 10,970 gene entities."]];
model.getRange("A12:G12").format = { fill: COLORS.amber, font: { bold: true, color: COLORS.ink }, wrapText: true };
model.freezePanes.freezeRows(5);
model.getRange("A:G").format.columnWidth = 22;
model.getRange("E:F").format.columnWidth = 42;

const evidence = wb.worksheets.add("Evidence Mapping");
title(evidence, "Evidence target resolution", "Explicit isoform → isoform; canonical accession → canonical protein; gene-only → gene_product_unspecified.", "G");
evidence.getRange("A5:G5").values = [["Module", "Resolution level", "Rows", "% within module", "Target assertion", "Isoform inference", "Notes"]];
header(evidence.getRange("A5:G5"));
evidence.getRange("A6:G12").values = [
  ["Binding evidence", "canonical_protein", 2995192, null, "Canonical protein", "Forbidden", "Historical isoform specificity not recoverable"],
  ["Binding evidence", "complex_context_ambiguous", 8114, null, "Complex review", "Forbidden", "Preserved for complex-entity reassignment"],
  ["Binding sites", "canonical_protein", 95598, null, "Canonical protein", "Forbidden", "No explicit isoform ID in V6.2"],
  ["Disease", "canonical_protein", 4613, null, "UniProt accession", "Forbidden", "UniProtKB disease annotations"],
  ["Disease", "gene_product_unspecified", 5651, null, "Gene", "Forbidden", "Open Targets is gene-level"],
  ["Expression", "gene_product_unspecified", 10997, null, "Gene", "Forbidden", "HPA expression is never assigned to a specific isoform"],
  ["Binding evidence", "isoform", 0, null, "Isoform", "Exact only", "No explicit isoform IDs survived the frozen V6.2 normalization"],
];
evidence.getRange("D6").formulas = [["=C6/SUM($C$6:$C$7)"]];
evidence.getRange("D6:D7").fillDown();
evidence.getRange("D8").formulas = [["=1"]];
evidence.getRange("D9").formulas = [["=C9/SUM($C$9:$C$10)"]];
evidence.getRange("D9:D10").fillDown();
evidence.getRange("D11").formulas = [["=1"]];
evidence.getRange("D12").formulas = [["=0"]];
body(evidence.getRange("A6:G12"));
evidence.getRange("C6:C12").format.numberFormat = "#,##0";
evidence.getRange("D6:D12").format.numberFormat = "0.00%";
evidence.getRange("A14:G15").merge();
evidence.getRange("A14").values = [["Important: linking a legacy row to a canonical accession means the normalized target is the canonical protein record. It does not prove that the experiment used the canonical isoform sequence."]];
evidence.getRange("A14:G15").format = { fill: COLORS.amber, font: { bold: true, color: COLORS.ink }, wrapText: true, verticalAlignment: "center" };
evidence.freezePanes.freezeRows(5);
evidence.getRange("A:B").format.columnWidth = 24;
evidence.getRange("C:D").format.columnWidth = 16;
evidence.getRange("E:F").format.columnWidth = 20;
evidence.getRange("G:G").format.columnWidth = 48;

const complex = wb.worksheets.add("Complex Components");
title(complex, "Complex-component identity resolution", "All 97,993 component assertions are preserved; explicit isoforms were verified against UniProt.", "H");
complex.getRange("A5:D5").values = [["Component reference level", "Rows", "%", "Meaning"]];
header(complex.getRange("A5:D5"));
complex.getRange("A6:D10").values = [
  ["Core canonical protein", 59080, null, "Canonical protein in the 10,997 membrane master"],
  ["External canonical protein", 37807, null, "Non-core protein required to describe a complex"],
  ["Core isoform", 324, null, "Explicit isoform of a core membrane protein"],
  ["External isoform", 34, null, "Explicit isoform of an external complex component"],
  ["Non-protein / unresolved component type", 748, null, "RNA, small molecule, or source component without protein identity"],
];
complex.getRange("C6").formulas = [["=B6/SUM($B$6:$B$10)"]];
complex.getRange("C6:C10").fillDown();
body(complex.getRange("A6:D10"));
complex.getRange("B6:B10").format.numberFormat = "#,##0";
complex.getRange("C6:C10").format.numberFormat = "0.00%";
const complexChart = complex.charts.add("doughnut", complex.getRange("A5:B10"));
complexChart.title = "Component assertion composition";
complexChart.hasLegend = true;
complexChart.setPosition("F5", "M20");
complex.getRange("A13:C13").values = [["External sequence status", "Proteins", "Action"]];
header(complex.getRange("A13:C13"));
complex.getRange("A14:C16").values = [
  ["Resolved in reference proteome", 6102, "Use frozen bulk sequence"],
  ["Resolved by direct UniProt request", 16, "Use direct-accession sequence"],
  ["Old/deleted accession review", 8, "Retain source ID; do not silently remap"],
];
body(complex.getRange("A14:C16"));
complex.getRange("B14:B16").format.numberFormat = "#,##0";
complex.freezePanes.freezeRows(5);
complex.getRange("A:A").format.columnWidth = 34;
complex.getRange("B:C").format.columnWidth = 18;
complex.getRange("D:D").format.columnWidth = 52;

const files = wb.worksheets.add("File Guide");
title(files, "Authoritative files", "Use the newest version shown here. Earlier QA reports remain only as development audit history.", "E");
files.getRange("A5:E5").values = [["File", "Version", "Rows / scope", "Purpose", "Authoritative?"]];
header(files.getRange("A5:E5"));
files.getRange("A6:E19").values = [
  ["gene_entity_v0_1.tsv", "0.1", "10,970", "Gene entities", "Yes"],
  ["gene_identifier_v0_1.tsv", "0.1", "46,505", "HGNC/Ensembl/NCBI Gene aliases", "Yes"],
  ["canonical_protein_entity_v0_1.tsv.gz", "0.1", "10,997", "Core canonical membrane proteins", "Yes"],
  ["canonical_protein_gene_link_v0_1.tsv", "0.1", "10,997", "Protein-to-gene bridge", "Yes"],
  ["protein_isoform_v0_2.tsv.gz", "0.2", "20,211", "Isoform sequence and parent FK", "Yes"],
  ["binding_evidence_target_resolution_v0_1.tsv.gz", "0.1", "3,003,306", "Binding evidence target bridge", "Yes"],
  ["binding_site_target_resolution_v0_1.tsv.gz", "0.1", "95,598", "Binding-site target bridge", "Yes"],
  ["disease_relation_target_resolution_v0_1.tsv", "0.1", "10,264", "Disease target-level bridge", "Yes"],
  ["expression_target_resolution_v0_1.tsv", "0.1", "10,997", "HPA gene-level projection", "Yes"],
  ["complex_target_components_v0_4.tsv.gz", "0.4", "97,993", "Complex components with entity FKs", "Yes"],
  ["external_complex_protein_entity_v0_2.tsv.gz", "0.2", "6,126", "External component proteins", "Yes"],
  ["external_complex_protein_review_v0_2.tsv", "0.2", "8", "Old/deleted accession review", "Yes"],
  ["website_entity_search_index_v0_1.tsv.gz", "0.1", "61,969", "Four-entity search index", "Yes"],
  ["IDENTITY_LAYER_V0_3_VALIDATION.json", "0.3", "PASS", "Final validation report", "Yes"],
];
body(files.getRange("A6:E19"));
files.getRange("E6:E19").format = { fill: COLORS.green, font: { color: COLORS.greenText, bold: true }, borders: { preset: "all", style: "thin", color: COLORS.border } };
files.freezePanes.freezeRows(5);
files.getRange("A:A").format.columnWidth = 48;
files.getRange("B:C").format.columnWidth = 18;
files.getRange("D:D").format.columnWidth = 42;
files.getRange("E:E").format.columnWidth = 16;

const validation = wb.worksheets.add("Validation");
title(validation, "Final V0.3 validation checks", "Every check below passed. The old V0.1/V0.2 QA files are superseded.", "C");
validation.getRange("A5:C5").values = [["Check", "Result", "Interpretation"]];
header(validation.getRange("A5:C5"));
const checks = [
  ["Gene primary keys unique", "PASS", "No duplicate gene_entity_id"],
  ["Canonical primary keys unique", "PASS", "Exactly 10,997 canonical proteins"],
  ["Canonical → gene FK complete", "PASS", "Every core protein links to a gene entity"],
  ["Isoform primary keys unique", "PASS", "No duplicate isoform entity ID"],
  ["Isoform parent FK complete", "PASS", "Every isoform links to a core or external canonical protein"],
  ["Binding row count preserved", "PASS", "3,003,306 / 3,003,306"],
  ["Binding isoform not inferred", "PASS", "No canonical-only row was promoted to isoform"],
  ["Binding-site row count preserved", "PASS", "95,598 / 95,598"],
  ["Disease row count preserved", "PASS", "10,264 / 10,264"],
  ["Open Targets kept gene-level", "PASS", "5,651 gene_product_unspecified"],
  ["HPA kept gene-level", "PASS", "10,997 gene_product_unspecified projections"],
  ["Complex component count preserved", "PASS", "97,993 / 97,993"],
  ["Explicit complex isoforms resolved", "PASS", "168 / 168 official FASTA responses"],
  ["Review records retained", "PASS", "8 old/deleted external accessions queued"],
  ["Four website entity types present", "PASS", "Gene, canonical, isoform and complex"],
];
validation.getRange(`A6:C${5 + checks.length}`).values = checks;
body(validation.getRange(`A6:C${5 + checks.length}`));
validation.getRange(`B6:B${5 + checks.length}`).format = { fill: COLORS.green, font: { color: COLORS.greenText, bold: true }, borders: { preset: "all", style: "thin", color: COLORS.border } };
validation.freezePanes.freezeRows(5);
validation.getRange("A:A").format.columnWidth = 42;
validation.getRange("B:B").format.columnWidth = 14;
validation.getRange("C:C").format.columnWidth = 56;

for (const sheet of [summary, model, evidence, complex, files, validation]) {
  const used = sheet.getUsedRange();
  if (used) used.format.autofitRows();
}

const xlsx = await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(outputDir, "MemPro_V6_3_identity_layer_review.xlsx"));
const preview = await wb.render({ sheetName: "Summary", autoCrop: "all", scale: 1, format: "png" });
await fs.writeFile(path.join(outputDir, "MemPro_V6_3_identity_layer_review_preview.png"), new Uint8Array(await preview.arrayBuffer()));
const inspection = await wb.inspect({ kind: "sheet,region", sheetId: "Summary", range: "A1:M18", maxChars: 5000 });
await fs.writeFile(path.join(outputDir, "workbook_inspection.json"), inspection.ndjson ?? String(inspection), "utf8");
console.log(JSON.stringify({ workbook: path.join(outputDir, "MemPro_V6_3_identity_layer_review.xlsx"), preview: path.join(outputDir, "MemPro_V6_3_identity_layer_review_preview.png") }, null, 2));
