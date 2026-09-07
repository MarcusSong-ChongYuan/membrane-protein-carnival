import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const files = [
  "C:/Users/Administrator/Desktop/semifinal/normalized_tables_v3_1.xlsx",
  "C:/Users/Administrator/Desktop/semifinal/protein_gene_disease_rule_b_summary_v1.xlsx",
];

for (const file of files) {
  const wb = await SpreadsheetFile.importXlsx(await FileBlob.load(file));
  console.log(`\n### FILE ${file}`);
  console.log((await wb.inspect({
    kind: "sheet",
    include: "id,name",
    maxChars: 10000,
  })).ndjson);
  for (const sheet of wb.worksheets.items) {
    const used = sheet.getUsedRange(true);
    console.log(`\n## SHEET ${sheet.name} RANGE ${used?.address ?? "EMPTY"}`);
    if (used) {
      const rows = Math.min(8, used.rowCount);
      const cols = Math.min(50, used.columnCount);
      const sample = sheet.getRangeByIndexes(0, 0, rows, cols);
      console.log(JSON.stringify(sample.values));
    }
  }
}
