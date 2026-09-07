import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = "C:/Users/Administrator/7.22/membrane_master_v5_working";
const outputDir = "C:/Users/Administrator/outputs/membrane-v5-20260724";
const previewDir = path.join(outputDir, "previews");
await fs.mkdir(previewDir, { recursive: true });

function parseTsv(text) {
  const lines = text.replace(/^\uFEFF/, "").replace(/\r/g, "").split("\n");
  if (lines.at(-1) === "") lines.pop();
  const headers = lines[0].split("\t");
  const rows = lines.slice(1).map((line) => {
    const fields = line.split("\t");
    const row = {};
    headers.forEach((header, index) => {
      row[header] = fields[index] ?? "";
    });
    return row;
  });
  return { headers, rows };
}

async function readTsv(relativePath) {
  return parseTsv(await fs.readFile(path.join(root, relativePath), "utf8"));
}

function matrixFrom(rows, columns) {
  return [
    columns.map((column) => column.label),
    ...rows.map((row) =>
      columns.map((column) => {
        const raw = row[column.key] ?? "";
        if (column.type === "number" && raw !== "") {
          const numeric = Number(raw);
          return Number.isFinite(numeric) ? numeric : raw;
        }
        return raw;
      }),
    ),
  ];
}

function colLetter(indexZeroBased) {
  let value = indexZeroBased + 1;
  let result = "";
  while (value > 0) {
    const remainder = (value - 1) % 26;
    result = String.fromCharCode(65 + remainder) + result;
    value = Math.floor((value - 1) / 26);
  }
  return result;
}

function styleDataSheet(sheet, rowCount, colCount, tableName, widths = {}) {
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(2);
  const used = sheet.getRangeByIndexes(0, 0, rowCount + 1, colCount);
  used.format.font = { name: "Aptos", size: 9, color: "#172033" };
  used.format.verticalAlignment = "top";
  const header = sheet.getRangeByIndexes(0, 0, 1, colCount);
  header.format = {
    fill: "#17324D",
    font: { name: "Aptos Display", size: 10, bold: true, color: "#FFFFFF" },
    wrapText: true,
    verticalAlignment: "center",
    borders: { bottom: { style: "medium", color: "#0F766E" } },
  };
  header.format.rowHeight = 32;
  const table = sheet.tables.add(
    `A1:${colLetter(colCount - 1)}${rowCount + 1}`,
    true,
    tableName,
  );
  table.style = "TableStyleMedium2";
  table.showBandedColumns = false;
  table.showFilterButton = true;
  for (let index = 0; index < colCount; index += 1) {
    const width = widths[index] ?? 16;
    sheet.getRangeByIndexes(0, index, rowCount + 1, 1).format.columnWidth = width;
  }
  used.format.wrapText = false;
  return table;
}

async function writeMatrixInChunks(sheet, matrix, chunkRows = 750) {
  const columnCount = matrix[0].length;
  for (let start = 0; start < matrix.length; start += chunkRows) {
    const chunk = matrix.slice(start, start + chunkRows);
    sheet.getRangeByIndexes(start, 0, chunk.length, columnCount).values = chunk;
  }
}

const master = await readTsv("release/human_membrane_protein_master_v5.tsv");
const review = await readTsv("review/human_membrane_review_queue_v5.tsv");
const trembl = await readTsv("crosswalks/unreviewed_canonical_resolution_v5.tsv");
const legacy = await readTsv("crosswalks/legacy_377_resolution_v5.tsv");
const singlePass = await readTsv("crosswalks/single_pass_resolution_v5.tsv");
const sources = await readTsv("release/SOURCE_REGISTRY_v5.tsv");
const stats = JSON.parse(
  await fs.readFile(path.join(root, "reports/V5_BUILD_STATS.json"), "utf8"),
);

const masterColumns = [
  { key: "release_tier_v5", label: "Tier" },
  { key: "target_uniprot_id", label: "UniProt accession" },
  { key: "approved_symbol", label: "Gene symbol" },
  { key: "protein_name", label: "Protein name" },
  { key: "reviewed", label: "Reviewed" },
  { key: "record_status_v5", label: "Record status" },
  { key: "membrane_decision_v5", label: "Membrane decision" },
  { key: "membrane_confidence_v5", label: "Confidence" },
  { key: "membrane_scope_v5", label: "Membrane scope" },
  { key: "membrane_topology_v5", label: "Topology" },
  { key: "transmembrane_count_v5", label: "TM count", type: "number" },
  { key: "single_pass_resolution_v5", label: "Single-pass resolution" },
  { key: "functional_primary_class_v5", label: "Primary class" },
  { key: "functional_subclass_v5", label: "Functional subclass" },
  { key: "classification_source_v5", label: "Classification source" },
  { key: "classification_status_v5", label: "Classification status" },
  { key: "independent_membrane_source_count_v5", label: "Independent source count", type: "number" },
  { key: "independent_membrane_sources_v5", label: "Independent sources" },
  { key: "hpa_predicted_membrane_v5", label: "HPA predicted membrane", type: "number" },
  { key: "hpa_plasma_membrane_location_v5", label: "HPA plasma membrane", type: "number" },
  { key: "htp_present_v5", label: "UniTmp HTP", type: "number" },
  { key: "membranome_present_v5", label: "Membranome", type: "number" },
  { key: "opm_present_v5", label: "OPM", type: "number" },
  { key: "pdbtm_present_v5", label: "PDBTM", type: "number" },
  { key: "gpcrdb_class_v5", label: "GPCRdb class" },
  { key: "gpcrdb_family_v5", label: "GPCRdb family" },
  { key: "gtopdb_type_v5", label: "GtoPdb type" },
  { key: "gtopdb_family_name_v5", label: "GtoPdb family" },
  { key: "tcdb_tcids_v5", label: "TCDB IDs" },
  { key: "cross_source_conflict_v5", label: "Cross-source conflict" },
  { key: "review_reason_v5", label: "Review reason" },
];

const reviewColumns = [
  { key: "release_tier_v5", label: "Tier" },
  { key: "target_uniprot_id", label: "UniProt accession" },
  { key: "approved_symbol", label: "Gene symbol" },
  { key: "protein_name", label: "Protein name" },
  { key: "record_status_v5", label: "Record status" },
  { key: "membrane_decision_v5", label: "Membrane decision" },
  { key: "membrane_confidence_v5", label: "Confidence" },
  { key: "independent_membrane_sources_v5", label: "Independent sources" },
  { key: "hpa_predicted_membrane_v5", label: "HPA predicted membrane", type: "number" },
  { key: "hpa_plasma_membrane_location_v5", label: "HPA plasma membrane", type: "number" },
  { key: "htp_evidence_v5", label: "UniTmp evidence" },
  { key: "membranome_family_v5", label: "Membranome family" },
  { key: "opm_pdb_ids_v5", label: "OPM PDB IDs" },
  { key: "pdbtm_pdb_ids_v5", label: "PDBTM PDB IDs" },
  { key: "functional_primary_class_v5", label: "Primary class" },
  { key: "cross_source_conflict_v5", label: "Cross-source conflict" },
  { key: "review_reason_v5", label: "Review reason" },
];

const tremblColumns = [
  { key: "candidate_uniprot_id", label: "TrEMBL accession" },
  { key: "candidate_entry_name", label: "Entry name" },
  { key: "candidate_gene_symbol", label: "Gene symbol" },
  { key: "v4_gene_symbol_resolution_status", label: "Original gene-symbol status" },
  { key: "canonical_reviewed_uniprot_id", label: "Canonical reviewed accession" },
  { key: "canonical_match_methods", label: "Match methods" },
  { key: "canonical_resolution_status", label: "Resolution status" },
  { key: "external_membrane_support", label: "External membrane support" },
  { key: "htp_evidence", label: "UniTmp HTP evidence" },
  { key: "htp_num_tm", label: "UniTmp TM count", type: "number" },
];

const legacyColumns = legacy.headers.map((key) => ({
  key,
  label: key
    .replace(/_v5$/i, " v5")
    .replace(/_/g, " ")
    .replace(/\b\w/g, (letter) => letter.toUpperCase()),
  type: key.includes("count") ? "number" : undefined,
}));

const singlePassColumns = [
  { key: "target_uniprot_id", label: "UniProt accession" },
  { key: "approved_symbol", label: "Gene symbol" },
  { key: "original_topology", label: "Original topology" },
  { key: "resolved_topology_v5", label: "Resolved topology v5" },
  { key: "resolution_basis", label: "Resolution basis" },
  { key: "htp_evidence", label: "UniTmp evidence" },
  { key: "htp_num_tm", label: "UniTmp TM count", type: "number" },
  { key: "htp_n_terminal_side", label: "N-terminal side" },
  { key: "htp_c_terminal_side", label: "C-terminal side" },
  { key: "uniprot_signal_peptide", label: "UniProt signal peptide" },
];

const sourceColumns = sources.headers.map((key) => ({
  key,
  label: key.replace(/_/g, " ").replace(/\b\w/g, (letter) => letter.toUpperCase()),
}));

const workbook = Workbook.create();
const summarySheet = workbook.worksheets.add("Summary");
const masterSheet = workbook.worksheets.add("Master Index");
const reviewSheet = workbook.worksheets.add("Review Queue");
const tremblSheet = workbook.worksheets.add("TrEMBL Resolution");
const legacySheet = workbook.worksheets.add("Legacy 377");
const singlePassSheet = workbook.worksheets.add("Single-pass 680");
const sourceSheet = workbook.worksheets.add("Source Registry");

await writeMatrixInChunks(masterSheet, matrixFrom(master.rows, masterColumns));
styleDataSheet(masterSheet, master.rows.length, masterColumns.length, "MasterIndexTable", {
  0: 8, 1: 16, 2: 14, 3: 42, 4: 10, 5: 18, 6: 24, 7: 14, 8: 21,
  9: 22, 10: 10, 11: 24, 12: 24, 13: 36, 14: 24, 15: 20, 16: 14, 17: 34,
  24: 22, 25: 36, 26: 20, 27: 36, 28: 18, 29: 18, 30: 50,
});
masterSheet.getRange(`A2:A${master.rows.length + 1}`).conditionalFormats.add(
  "cellIs",
  { operator: "equal", formula: '"A"', format: { fill: "#DCFCE7", font: { color: "#166534", bold: true } } },
);
masterSheet.getRange(`A2:A${master.rows.length + 1}`).conditionalFormats.add(
  "cellIs",
  { operator: "equal", formula: '"B"', format: { fill: "#DBEAFE", font: { color: "#1E40AF", bold: true } } },
);
masterSheet.getRange(`A2:A${master.rows.length + 1}`).conditionalFormats.add(
  "cellIs",
  { operator: "equal", formula: '"C"', format: { fill: "#FEF3C7", font: { color: "#92400E", bold: true } } },
);
masterSheet.getRange(`A2:A${master.rows.length + 1}`).conditionalFormats.add(
  "cellIs",
  { operator: "equal", formula: '"R"', format: { fill: "#FEE2E2", font: { color: "#991B1B", bold: true } } },
);

await writeMatrixInChunks(reviewSheet, matrixFrom(review.rows, reviewColumns));
styleDataSheet(reviewSheet, review.rows.length, reviewColumns.length, "ReviewQueueTable", {
  0: 8, 1: 16, 2: 14, 3: 42, 4: 18, 5: 24, 6: 14, 7: 34, 10: 22,
  11: 30, 12: 28, 13: 28, 14: 24, 15: 20, 16: 54,
});
reviewSheet.getRange(`A2:A${review.rows.length + 1}`).format = {
  fill: "#FEE2E2",
  font: { color: "#991B1B", bold: true },
};

await writeMatrixInChunks(tremblSheet, matrixFrom(trembl.rows, tremblColumns));
styleDataSheet(tremblSheet, trembl.rows.length, tremblColumns.length, "TremblResolutionTable", {
  0: 18, 1: 24, 2: 15, 3: 28, 4: 22, 5: 38, 6: 36, 7: 28, 8: 22, 9: 14,
});

await writeMatrixInChunks(legacySheet, matrixFrom(legacy.rows, legacyColumns));
styleDataSheet(legacySheet, legacy.rows.length, legacyColumns.length, "LegacyResolutionTable", {
  0: 18, 1: 15, 2: 42, 3: 28, 4: 24, 5: 14, 6: 34, 7: 22, 8: 30,
  9: 32, 10: 12, 11: 34, 12: 54,
});

await writeMatrixInChunks(singlePassSheet, matrixFrom(singlePass.rows, singlePassColumns));
styleDataSheet(singlePassSheet, singlePass.rows.length, singlePassColumns.length, "SinglePassResolutionTable", {
  0: 18, 1: 15, 2: 28, 3: 28, 4: 42, 5: 22, 6: 14, 7: 16, 8: 16, 9: 26,
});

await writeMatrixInChunks(sourceSheet, matrixFrom(sources.rows, sourceColumns));
styleDataSheet(sourceSheet, sources.rows.length, sourceColumns.length, "SourceRegistryTable", {
  0: 24, 1: 20, 2: 46, 3: 52, 4: 48, 5: 58,
});
sourceSheet.getRange(`C2:F${sources.rows.length + 1}`).format.wrapText = true;
sourceSheet.getRange(`A2:F${sources.rows.length + 1}`).format.rowHeight = 48;

summarySheet.showGridLines = false;
summarySheet.getRange("A1:L2").merge();
summarySheet.getRange("A1").values = [["Human Membrane Protein Master v5"]];
summarySheet.getRange("A1:L2").format = {
  fill: "#102A43",
  font: { name: "Aptos Display", size: 22, bold: true, color: "#FFFFFF" },
  verticalAlignment: "center",
  horizontalAlignment: "left",
  borders: { bottom: { style: "thick", color: "#14B8A6" } },
};
summarySheet.getRange("A3:L3").merge();
summarySheet.getRange("A3").values = [[
  `Release ${stats.release_date} · Reviewed human protein union across UniProt and independent membrane resources`,
]];
summarySheet.getRange("A3:L3").format = {
  fill: "#E6FFFB",
  font: { size: 10, color: "#134E4A", italic: true },
};

summarySheet.getRange("A5:B5").values = [["All v5 records", null]];
summarySheet.getRange("D5:E5").values = [["Core A+B", null]];
summarySheet.getRange("G5:H5").values = [["Review tier R", null]];
summarySheet.getRange("J5:K5").values = [["Classified records", null]];
summarySheet.getRange("B5").formulas = [[`=COUNTA('Master Index'!$B$2:$B$${master.rows.length + 1})`]];
summarySheet.getRange("E5").formulas = [[
  `=COUNTIF('Master Index'!$A$2:$A$${master.rows.length + 1},"A")+COUNTIF('Master Index'!$A$2:$A$${master.rows.length + 1},"B")`,
]];
summarySheet.getRange("H5").formulas = [[`=COUNTIF('Master Index'!$A$2:$A$${master.rows.length + 1},"R")`]];
const classColumn = colLetter(masterColumns.findIndex((column) => column.key === "functional_primary_class_v5"));
summarySheet.getRange("K5").formulas = [[
  `=COUNTA('Master Index'!$B$2:$B$${master.rows.length + 1})-COUNTIF('Master Index'!$${classColumn}$2:$${classColumn}$${master.rows.length + 1},"unclassified")`,
]];
for (const block of ["A5:B6", "D5:E6", "G5:H6", "J5:K6"]) {
  summarySheet.getRange(block).format = {
    fill: "#F8FAFC",
    borders: { preset: "outside", style: "thin", color: "#CBD5E1" },
  };
}
for (const cell of ["A5", "D5", "G5", "J5"]) {
  summarySheet.getRange(cell).format.font = { size: 10, bold: true, color: "#475569" };
}
for (const cell of ["B5", "E5", "H5", "K5"]) {
  summarySheet.getRange(cell).format = {
    fill: "#FFFFFF",
    font: { size: 17, bold: true, color: "#0F766E" },
    numberFormat: "#,##0",
    horizontalAlignment: "right",
  };
}

summarySheet.getRange("A8:C8").values = [["Release tier", "Records", "Meaning"]];
summarySheet.getRange("A9:A12").values = [["A"], ["B"], ["C"], ["R"]];
summarySheet.getRange("B9").formulas = [[`=COUNTIF('Master Index'!$A$2:$A$${master.rows.length + 1},A9)`]];
summarySheet.getRange("B9:B12").fillDown();
summarySheet.getRange("C9:C12").values = [
  ["Integral membrane core"],
  ["Monotopic / lipid-anchored core"],
  ["Peripheral membrane extended set"],
  ["Prediction, location-only, or unresolved review"],
];
summarySheet.getRange("A8:C12").format.borders = {
  preset: "inside",
  style: "thin",
  color: "#E2E8F0",
};
summarySheet.getRange("A8:C8").format = {
  fill: "#17324D",
  font: { bold: true, color: "#FFFFFF" },
};
summarySheet.getRange("A9:A12").format.font = { bold: true };
summarySheet.getRange("B9:B12").format.numberFormat = "#,##0";
summarySheet.getRange("A9:C9").format.fill = "#DCFCE7";
summarySheet.getRange("A10:C10").format.fill = "#DBEAFE";
summarySheet.getRange("A11:C11").format.fill = "#FEF3C7";
summarySheet.getRange("A12:C12").format.fill = "#FEE2E2";

summarySheet.getRange("A15:B15").values = [["Independent source", "Master records supported"]];
const supportRows = [
  ["HPA predicted membrane", "hpa_predicted_membrane_v5"],
  ["HPA plasma-membrane location", "hpa_plasma_membrane_location_v5"],
  ["UniTmp HTP", "htp_present_v5"],
  ["Membranome", "membranome_present_v5"],
  ["OPM", "opm_present_v5"],
  ["PDBTM", "pdbtm_present_v5"],
];
summarySheet.getRange("A16:A21").values = supportRows.map(([label]) => [label]);
summarySheet.getRange("B16:B21").formulas = supportRows.map(([, key]) => {
  const column = colLetter(masterColumns.findIndex((item) => item.key === key));
  return [`=COUNTIF('Master Index'!$${column}$2:$${column}$${master.rows.length + 1},1)`];
});
summarySheet.getRange("A15:B21").format.borders = {
  preset: "inside",
  style: "thin",
  color: "#E2E8F0",
};
summarySheet.getRange("A15:B15").format = {
  fill: "#17324D",
  font: { bold: true, color: "#FFFFFF" },
};
summarySheet.getRange("B16:B21").format.numberFormat = "#,##0";

summarySheet.getRange("A24:D24").values = [["Audit workstream", "Input", "Resolved / assigned", "Residual review"]];
summarySheet.getRange("A25:D28").values = [
  ["TrEMBL canonicalization", stats.unreviewed_candidates, stats.unreviewed_resolution_status.absorbed_to_reviewed_canonical, stats.unreviewed_without_single_canonical],
  ["Legacy v3.1 records", stats.legacy_rows, stats.legacy_identifier_status.current_reviewed_accession + stats.legacy_identifier_status.mapped_by_approved_symbol, stats.legacy_membrane_scope_status.external_support_review + stats.legacy_membrane_scope_status.not_supported_as_membrane_in_v5],
  ["Single-pass topology", stats.single_pass_original_unresolved, stats.single_pass_resolution.single_pass_type_i + stats.single_pass_resolution.single_pass_type_ii + stats.single_pass_resolution.single_pass_type_iii, stats.single_pass_resolution.single_pass_type_i_or_iii + stats.single_pass_resolution.single_pass_unresolved],
  ["v4 unclassified proteins", stats.v4_unclassified, stats.v4_unclassified - stats.v5_unclassified_retained_v4, stats.v5_unclassified_retained_v4],
];
summarySheet.getRange("A24:D28").format.borders = {
  preset: "inside",
  style: "thin",
  color: "#E2E8F0",
};
summarySheet.getRange("A24:D24").format = {
  fill: "#17324D",
  font: { bold: true, color: "#FFFFFF" },
};
summarySheet.getRange("B25:D28").format.numberFormat = "#,##0";

summarySheet.getRange("A31:L31").merge();
summarySheet.getRange("A31").values = [["Interpretation and publication-use notes"]];
summarySheet.getRange("A31:L31").format = {
  fill: "#0F766E",
  font: { bold: true, color: "#FFFFFF", size: 11 },
};
summarySheet.getRange("A32:L35").merge(true);
summarySheet.getRange("A32:L35").values = [
  ["• “Complete” means a versioned union of the named sources under the documented v5 inclusion policy; it is not a claim of closed biological completeness."],
  ["• Tier A+B is the recommended core membrane-protein universe. Tier C is the peripheral extension. Tier R is intentionally excluded from public core counts until reviewed."],
  ["• HPA subcellular location supports membrane association but does not by itself promote a record into the integral-membrane core."],
  ["• GPCRdb, GtoPdb and TCDB are used for classification alignment, not as independent proof of membrane membership."],
];
summarySheet.getRange("A32:L35").format = {
  fill: "#F8FAFC",
  font: { size: 10, color: "#334155" },
  wrapText: true,
  verticalAlignment: "center",
};
summarySheet.getRange("A32:L35").format.rowHeight = 30;

const tierChart = summarySheet.charts.add("bar", summarySheet.getRange("A8:B12"));
tierChart.title = "v5 records by release tier";
tierChart.titleTextStyle.fontSize = 12;
tierChart.hasLegend = false;
tierChart.xAxis = { axisType: "textAxis", textStyle: { fontSize: 9 } };
tierChart.yAxis = { numberFormatCode: "#,##0" };
tierChart.setPosition("E8", "L22");

summarySheet.freezePanes.freezeRows(3);
summarySheet.getRange("A1:L35").format.font = {
  name: "Aptos",
  color: "#172033",
};
summarySheet.getRange("A1:L2").format.font = {
  name: "Aptos Display",
  size: 22,
  bold: true,
  color: "#FFFFFF",
};
summarySheet.getRange("A1:A35").format.columnWidth = 29;
summarySheet.getRange("B1:B35").format.columnWidth = 15;
summarySheet.getRange("C1:C35").format.columnWidth = 39;
for (const column of ["D", "E", "F", "G", "H", "I", "J", "K", "L"]) {
  summarySheet.getRange(`${column}1:${column}35`).format.columnWidth = 14;
}

const inspection = await workbook.inspect({
  kind: "workbook,sheet,table,formula,drawing",
  maxChars: 8000,
  tableMaxRows: 4,
  tableMaxCols: 8,
  options: { maxResults: 120 },
});
await fs.writeFile(path.join(outputDir, "workbook_inspection.ndjson"), inspection.ndjson ?? String(inspection), "utf8");

const renderSpecs = [
  ["Summary", "A1:L35"],
  ["Master Index", "A1:Q30"],
  ["Review Queue", "A1:Q30"],
  ["TrEMBL Resolution", "A1:J30"],
  ["Legacy 377", "A1:M30"],
  ["Single-pass 680", "A1:J30"],
  ["Source Registry", `A1:F${Math.min(sources.rows.length + 1, 30)}`],
];
for (const [sheetName, range] of renderSpecs) {
  const preview = await workbook.render({
    sheetName,
    range,
    autoCrop: "all",
    scale: 1,
    format: "png",
  });
  const safeName = sheetName.toLowerCase().replace(/[^a-z0-9]+/g, "_");
  await fs.writeFile(
    path.join(previewDir, `${safeName}.png`),
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(path.join(outputDir, "human_membrane_protein_master_v5.xlsx"));

const qa = {
  output: path.join(outputDir, "human_membrane_protein_master_v5.xlsx"),
  sheets: renderSpecs.map(([name]) => name),
  rows: {
    master: master.rows.length,
    review: review.rows.length,
    trembl: trembl.rows.length,
    legacy: legacy.rows.length,
    singlePass: singlePass.rows.length,
    sources: sources.rows.length,
  },
  expected: {
    master: 10997,
    review: 1109,
    trembl: 9823,
    legacy: 377,
    singlePass: 680,
  },
};
await fs.writeFile(path.join(outputDir, "workbook_build_qa.json"), JSON.stringify(qa, null, 2), "utf8");
console.log(JSON.stringify(qa, null, 2));
