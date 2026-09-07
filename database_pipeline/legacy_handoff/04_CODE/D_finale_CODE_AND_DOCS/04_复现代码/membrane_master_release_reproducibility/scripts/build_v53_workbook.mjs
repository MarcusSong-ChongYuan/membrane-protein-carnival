import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const releaseDir = "C:/Users/Administrator/7.22/membrane_master_v5_working/releases/release_v5_3_sequences";
const outputDir = "C:/Users/Administrator/outputs/membrane_v53_sequences_20260724";
const outputPath = `${outputDir}/human_membrane_protein_master_v5_3_with_sequences.xlsx`;

function parseTsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') {
        field += '"';
        i++;
      } else if (ch === '"') {
        quoted = false;
      } else {
        field += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === "\t") {
      row.push(field);
      field = "";
    } else if (ch === "\n") {
      row.push(field.replace(/\r$/, ""));
      rows.push(row);
      row = [];
      field = "";
    } else {
      field += ch;
    }
  }
  if (field || row.length) {
    row.push(field.replace(/\r$/, ""));
    rows.push(row);
  }
  return rows;
}

async function loadTsv(name) {
  const text = (await fs.readFile(`${releaseDir}/${name}`, "utf8")).replace(/^\uFEFF/, "");
  const matrix = parseTsv(text);
  const headers = matrix[0];
  return matrix.slice(1).filter((row) => row.length > 1).map((row) => {
    const record = {};
    headers.forEach((header, index) => {
      record[header] = row[index] ?? "";
    });
    return record;
  });
}

function colName(index) {
  let n = index + 1;
  let name = "";
  while (n > 0) {
    n--;
    name = String.fromCharCode(65 + (n % 26)) + name;
    n = Math.floor(n / 26);
  }
  return name;
}

function typed(value, field) {
  if (
    /(?:count|length|_flag|_present|website_default|length_match)$/i.test(field) &&
    value !== "" &&
    /^-?\d+$/.test(value)
  ) {
    return Number(value);
  }
  return value;
}

const palette = {
  navy: "#17324D",
  teal: "#1D6F78",
  green: "#2E7D32",
  paleBlue: "#E8F0F7",
  paleGreen: "#E2F0D9",
  paleGold: "#FFF2CC",
  white: "#FFFFFF",
  text: "#1F2937",
  grid: "#D7DEE5",
};

const fields = [
  "membrane_protein_id",
  "target_uniprot_id",
  "approved_symbol",
  "protein_name",
  "membrane_class_v52",
  "evidence_level_v52",
  "evidence_label_zh_v52",
  "sequence_length_v53",
  "sequence_version_v53",
  "sequence_status_v53",
  "sequence_length_match_v53",
  "canonical_sequence_part_1",
  "canonical_sequence_part_2",
  "sequence_sha256_v53",
  "sequence_source_url_v53",
  "uniprot_release_v53",
  "sequence_retrieval_date_v53",
  "direct_experimental_flag_v52",
  "release_scope_v52",
  "website_default_v52",
  "evidence_rule_v52",
  "evidence_basis_v52",
  "membrane_scope_v5",
  "membrane_topology_v5",
  "transmembrane_count_v5",
  "transmembrane_features",
  "intramembrane_features",
  "subcellular_location",
  "independent_membrane_source_count_v5",
  "independent_membrane_sources_v5",
  "hpa_reliability_if_v5",
  "hpa_subcellular_location_v5",
  "opm_present_v5",
  "membranome_present_v5",
  "decision_code_v51",
  "decision_basis_v51",
];

const labels = {
  membrane_protein_id: "HMP ID",
  target_uniprot_id: "UniProt accession",
  approved_symbol: "Approved symbol",
  protein_name: "Protein name",
  membrane_class_v52: "Membrane class",
  evidence_level_v52: "Evidence level",
  evidence_label_zh_v52: "证据等级",
  sequence_length_v53: "Sequence length",
  sequence_version_v53: "Sequence version",
  sequence_status_v53: "Sequence status",
  sequence_length_match_v53: "Length match",
  canonical_sequence_part_1: "Canonical sequence part 1",
  canonical_sequence_part_2: "Canonical sequence part 2",
  sequence_sha256_v53: "Sequence SHA-256",
  sequence_source_url_v53: "Sequence source URL",
  uniprot_release_v53: "UniProt release",
  sequence_retrieval_date_v53: "Sequence retrieval date",
  direct_experimental_flag_v52: "Direct experimental",
  release_scope_v52: "Release scope",
  website_default_v52: "Website default",
  evidence_rule_v52: "Evidence rule",
  evidence_basis_v52: "Evidence basis",
  membrane_scope_v5: "v5 membrane scope",
  membrane_topology_v5: "v5 topology",
  transmembrane_count_v5: "v5 TM count",
  transmembrane_features: "UniProt TM features",
  intramembrane_features: "UniProt intramembrane features",
  subcellular_location: "UniProt subcellular location",
  independent_membrane_source_count_v5: "Independent source count",
  independent_membrane_sources_v5: "Independent membrane sources",
  hpa_reliability_if_v5: "HPA reliability",
  hpa_subcellular_location_v5: "HPA location",
  opm_present_v5: "OPM present",
  membranome_present_v5: "Membranome present",
  decision_code_v51: "Previous decision code",
  decision_basis_v51: "Previous decision basis",
};

const widths = {
  membrane_protein_id: 18,
  target_uniprot_id: 16,
  approved_symbol: 18,
  protein_name: 42,
  membrane_class_v52: 14,
  evidence_level_v52: 14,
  evidence_label_zh_v52: 14,
  sequence_length_v53: 15,
  sequence_version_v53: 15,
  sequence_status_v53: 20,
  sequence_length_match_v53: 14,
  canonical_sequence_part_1: 48,
  canonical_sequence_part_2: 48,
  sequence_sha256_v53: 40,
  sequence_source_url_v53: 44,
  uniprot_release_v53: 15,
  sequence_retrieval_date_v53: 19,
  direct_experimental_flag_v52: 17,
  release_scope_v52: 15,
  website_default_v52: 16,
  evidence_rule_v52: 38,
  evidence_basis_v52: 58,
  membrane_scope_v5: 30,
  membrane_topology_v5: 26,
  transmembrane_count_v5: 13,
  transmembrane_features: 42,
  intramembrane_features: 42,
  subcellular_location: 62,
  independent_membrane_source_count_v5: 16,
  independent_membrane_sources_v5: 36,
  hpa_reliability_if_v5: 15,
  hpa_subcellular_location_v5: 34,
  opm_present_v5: 13,
  membranome_present_v5: 17,
  decision_code_v51: 38,
  decision_basis_v51: 55,
};

function valueFor(record, field) {
  if (field === "canonical_sequence_part_1") {
    return (record.canonical_sequence ?? "").slice(0, 32000);
  }
  if (field === "canonical_sequence_part_2") {
    return (record.canonical_sequence ?? "").slice(32000);
  }
  return typed(record[field] ?? "", field);
}

function addRecords(sheet, records, tableName) {
  const headers = fields.map((field) => labels[field]);
  const matrix = [
    headers,
    ...records.map((record) => fields.map((field) => valueFor(record, field))),
  ];
  const lastCol = colName(fields.length - 1);
  sheet.getRange(`A1:${lastCol}${matrix.length}`).values = matrix;
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(4);
  sheet.getRange(`A1:${lastCol}1`).format = {
    fill: palette.navy,
    font: { bold: true, color: palette.white, size: 10 },
    wrapText: true,
    verticalAlignment: "center",
    borders: { preset: "inside", style: "thin", color: "#35516B" },
  };
  sheet.getRange(`A1:${lastCol}1`).format.rowHeight = 36;
  sheet.getRange(`A2:${lastCol}${Math.max(matrix.length, 2)}`).format = {
    font: { color: palette.text, size: 9 },
    verticalAlignment: "top",
    borders: { insideHorizontal: { style: "thin", color: palette.grid } },
  };
  fields.forEach((field, index) => {
    const col = colName(index);
    sheet.getRange(`${col}:${col}`).format.columnWidth = widths[field] ?? 15;
  });
  const table = sheet.tables.add(`A1:${lastCol}${matrix.length}`, true, tableName);
  table.style = "TableStyleMedium2";
  table.showFilterButton = true;
}

function dataColumn(field) {
  return colName(fields.indexOf(field));
}

await fs.mkdir(outputDir, { recursive: true });
const e1 = await loadTsv("human_membrane_experimental_core_E1_v5_3.tsv");
const e2 = await loadTsv("human_membrane_strong_supported_E2_v5_3.tsv");
const defaults = await loadTsv("human_membrane_default_E1_E2_v5_3.tsv");
const e3 = await loadTsv("human_membrane_prediction_candidates_E3_v5_3.tsv");
const e0 = await loadTsv("human_membrane_excluded_E0_v5_3.tsv");

const workbook = Workbook.create();
const summary = workbook.worksheets.add("Summary");
const sheetE1 = workbook.worksheets.add("E1 Experimental");
const sheetE2 = workbook.worksheets.add("E2 Supported");
const sheetDefault = workbook.worksheets.add("E1+E2 Default");
const sheetE3 = workbook.worksheets.add("E3 Candidates");
const sheetE0 = workbook.worksheets.add("E0 Excluded");
const rules = workbook.worksheets.add("Sequence & Evidence Rules");

addRecords(sheetE1, e1, "E1ExperimentalSequenceTable");
addRecords(sheetE2, e2, "E2SupportedSequenceTable");
addRecords(sheetDefault, defaults, "DefaultSequenceTable");
addRecords(sheetE3, e3, "E3CandidateSequenceTable");
addRecords(sheetE0, e0, "E0ExcludedSequenceTable");

summary.showGridLines = false;
summary.getRange("A1:J2").merge();
summary.getRange("A1").values = [["Human Membrane Protein Master v5.3 — canonical sequences"]];
summary.getRange("A1:J2").format = {
  fill: palette.navy,
  font: { bold: true, color: palette.white, size: 20 },
  verticalAlignment: "center",
};
summary.getRange("A3:J3").merge();
summary.getRange("A3").values = [[
  "All 10,997 audited accessions include the UniProtKB 2026_02 canonical amino-acid sequence. The v5.2 ABC × E0–E3 classification is unchanged.",
]];
summary.getRange("A3:J3").format = {
  fill: palette.paleBlue,
  font: { italic: true, color: palette.text, size: 10 },
};
summary.getRange("A5:B11").values = [
  ["Evidence layer", "Count"],
  ["E1 — 实验证实", null],
  ["E2 — 强证据支持", null],
  ["Default website E1+E2", null],
  ["E3 — 预测候选", null],
  ["E0 — 排除", null],
  ["All audited records", null],
];
summary.getRange("B6").formulas = [[`=COUNTA('E1 Experimental'!$A$2:$A$${e1.length + 1})`]];
summary.getRange("B7").formulas = [[`=COUNTA('E2 Supported'!$A$2:$A$${e2.length + 1})`]];
summary.getRange("B8").formulas = [[`=COUNTA('E1+E2 Default'!$A$2:$A$${defaults.length + 1})`]];
summary.getRange("B9").formulas = [[`=COUNTA('E3 Candidates'!$A$2:$A$${e3.length + 1})`]];
summary.getRange("B10").formulas = [[`=COUNTA('E0 Excluded'!$A$2:$A$${e0.length + 1})`]];
summary.getRange("B11").formulas = [["=SUM(B8:B10)"]];
summary.getRange("A5:B5").format = {
  fill: palette.teal,
  font: { bold: true, color: palette.white },
};
summary.getRange("A6:B11").format.borders = {
  preset: "inside",
  style: "thin",
  color: palette.grid,
};
summary.getRange("B6:B11").format = {
  numberFormat: "#,##0",
  font: { bold: true, size: 12, color: palette.navy },
  horizontalAlignment: "right",
};

summary.getRange("A13:B17").values = [
  ["Default E1+E2 class", "Count"],
  ["A — integral", null],
  ["B — lipid-anchored/inserted", null],
  ["C — peripheral/associated", null],
  ["Total", null],
];
const classCol = dataColumn("membrane_class_v52");
summary.getRange("B14").formulas = [[`=COUNTIF('E1+E2 Default'!$${classCol}$2:$${classCol}$${defaults.length + 1},"A")`]];
summary.getRange("B15").formulas = [[`=COUNTIF('E1+E2 Default'!$${classCol}$2:$${classCol}$${defaults.length + 1},"B")`]];
summary.getRange("B16").formulas = [[`=COUNTIF('E1+E2 Default'!$${classCol}$2:$${classCol}$${defaults.length + 1},"C")`]];
summary.getRange("B17").formulas = [["=SUM(B14:B16)"]];
summary.getRange("A13:B13").format = {
  fill: palette.teal,
  font: { bold: true, color: palette.white },
};
summary.getRange("A14:B17").format.borders = {
  preset: "inside",
  style: "thin",
  color: palette.grid,
};
summary.getRange("B14:B17").format.numberFormat = "#,##0";

summary.getRange("D5:E10").values = [
  ["Sequence QA", "Result"],
  ["Sequences retrieved", null],
  ["Missing sequences", null],
  ["Length matches", null],
  ["Length mismatches", null],
  ["Excel split records", null],
];
const statusCol = dataColumn("sequence_status_v53");
const matchCol = dataColumn("sequence_length_match_v53");
const part2Col = dataColumn("canonical_sequence_part_2");
const sequenceSheets = [
  ["E1+E2 Default", defaults.length],
  ["E3 Candidates", e3.length],
  ["E0 Excluded", e0.length],
];
const countAcross = (column, criterion) =>
  sequenceSheets.map(([name, count]) => `COUNTIF('${name}'!$${column}$2:$${column}$${count + 1},"${criterion}")`).join("+");
const nonblankAcross = (column) =>
  sequenceSheets.map(([name, count]) => `COUNTA('${name}'!$${column}$2:$${column}$${count + 1})`).join("+");
summary.getRange("E6").formulas = [[`=${nonblankAcross(statusCol)}`]];
summary.getRange("E7").formulas = [["=B11-E6"]];
summary.getRange("E8").formulas = [[`=${countAcross(matchCol, "1")}`]];
summary.getRange("E9").formulas = [[`=${countAcross(matchCol, "0")}`]];
summary.getRange("E10").formulas = [[`=${countAcross(dataColumn("target_uniprot_id"), "Q8WZ42")}`]];
summary.getRange("D5:E5").format = {
  fill: palette.green,
  font: { bold: true, color: palette.white },
};
summary.getRange("D6:E10").format.borders = {
  preset: "inside",
  style: "thin",
  color: palette.grid,
};
summary.getRange("E6:E10").format = {
  numberFormat: "#,##0",
  font: { bold: true, size: 12, color: palette.navy },
  horizontalAlignment: "right",
};

summary.getRange("G5:H9").values = [
  ["Evidence level", "Count"],
  ["E1", null],
  ["E2", null],
  ["E3", null],
  ["E0", null],
];
summary.getRange("H6:H9").formulas = [["=B6"], ["=B7"], ["=B9"], ["=B10"]];
const chart = summary.charts.add("bar", summary.getRange("G5:H9"));
chart.title = "Evidence-layer distribution";
chart.hasLegend = false;
chart.yAxis = { numberFormatCode: "#,##0" };
chart.setPosition("G5", "J17");

summary.getRange("A19:J23").merge();
summary.getRange("A19").values = [[
  "Sequence storage: concatenate “Canonical sequence part 1” and “part 2” to reconstruct the exact canonical sequence. Part 2 is blank for every entry except Q8WZ42 (Titin), whose 34,350-aa sequence exceeds Excel's one-cell limit. TSV and FASTA retain every sequence unsplit.",
]];
summary.getRange("A19:J23").format = {
  fill: palette.paleGold,
  font: { size: 11, color: palette.text },
  wrapText: true,
  verticalAlignment: "center",
  borders: { preset: "outside", style: "thin", color: "#D6B656" },
};
summary.getRange("A:A").format.columnWidth = 35;
summary.getRange("B:B").format.columnWidth = 16;
summary.getRange("C:C").format.columnWidth = 4;
summary.getRange("D:D").format.columnWidth = 24;
summary.getRange("E:E").format.columnWidth = 16;
for (const column of ["F", "G", "H", "I", "J"]) {
  summary.getRange(`${column}:${column}`).format.columnWidth = 13;
}

rules.showGridLines = false;
rules.getRange("A1:E1").merge();
rules.getRange("A1").values = [["Sequence provenance and evidence model"]];
rules.getRange("A1:E1").format = {
  fill: palette.navy,
  font: { bold: true, color: palette.white, size: 16 },
};
rules.getRange("A3:E9").values = [
  ["Topic", "Value", "Scope", "Release behavior", "Source"],
  ["Sequence", "Canonical", "All 10,997 audited records", "Included in TSV, FASTA and Excel", "UniProtKB REST API"],
  ["UniProt release", "2026_02", "Retrieved 2026-07-24", "Sequence version and SHA-256 retained", "https://rest.uniprot.org/"],
  ["E1", "实验证实", "Direct experimental membrane evidence", "Strict experimental core", "Shown by default"],
  ["E2", "强证据支持", "Curated or multiple reliable sources", "Supported layer", "Shown by default"],
  ["E3", "预测候选", "Prediction, single source or unresolved ambiguity", "Candidate layer", "Hidden unless enabled"],
  ["E0", "排除", "Not retained as a membrane protein", "Internal audit blacklist", "Not shown"],
];
rules.getRange("A3:E3").format = {
  fill: palette.teal,
  font: { bold: true, color: palette.white },
};
rules.getRange("A4:E9").format = {
  wrapText: true,
  verticalAlignment: "top",
  borders: { preset: "inside", style: "thin", color: palette.grid },
  font: { size: 10, color: palette.text },
};
rules.getRange("A:A").format.columnWidth = 18;
rules.getRange("B:B").format.columnWidth = 20;
rules.getRange("C:C").format.columnWidth = 48;
rules.getRange("D:D").format.columnWidth = 34;
rules.getRange("E:E").format.columnWidth = 34;

const summaryCheck = await workbook.inspect({
  kind: "table",
  range: "Summary!A1:J23",
  include: "values,formulas",
  tableMaxRows: 23,
  tableMaxCols: 10,
  maxChars: 7000,
});
await fs.writeFile(`${outputDir}/summary_inspect.ndjson`, summaryCheck.ndjson, "utf8");
const errorScan = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 300 },
  summary: "v5.3 final formula error scan",
  maxChars: 3000,
});
await fs.writeFile(`${outputDir}/formula_errors.ndjson`, errorScan.ndjson, "utf8");

for (const [sheetName, range] of [
  ["Summary", "A1:J23"],
  ["E1 Experimental", "A1:P16"],
  ["E2 Supported", "A1:P16"],
  ["E1+E2 Default", "A1:P16"],
  ["E3 Candidates", "A1:P16"],
  ["E0 Excluded", "A1:P16"],
  ["Sequence & Evidence Rules", "A1:E9"],
]) {
  const preview = await workbook.render({ sheetName, range, scale: 1, format: "png" });
  await fs.writeFile(
    `${outputDir}/preview_${sheetName.replace(/[^A-Za-z0-9+]/g, "_")}.png`,
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);
console.log(JSON.stringify({
  outputPath,
  E1: e1.length,
  E2: e2.length,
  default: defaults.length,
  E3: e3.length,
  E0: e0.length,
}, null, 2));
