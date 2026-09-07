import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const workbookPath = process.argv[2];
const reportPath = process.argv[3];
if (!workbookPath || !reportPath) throw new Error("Usage: node verify_review_workbook.mjs <xlsx> <report-json>");

const blob = await FileBlob.load(workbookPath);
const wb = await SpreadsheetFile.importXlsx(blob);
const sheets = await wb.inspect({ kind: "sheet", include: "id,name", maxChars: 5000 });
const summary = await wb.inspect({ kind: "region", sheetId: "Summary", range: "A1:M18", maxChars: 10000 });
const formulas = await wb.inspect({ kind: "formula", maxChars: 10000, options: { maxResults: 200 } });
const joined = [sheets.ndjson ?? String(sheets), summary.ndjson ?? String(summary), formulas.ndjson ?? String(formulas)].join("\n");
const errorMatches = joined.match(/#REF!|#DIV\/0!|#VALUE!|#NAME\?|#N\/A/g) ?? [];
const report = {
  workbook: workbookPath,
  verified_at: new Date().toISOString(),
  expected_sheets: ["Summary", "Entity Model", "Evidence Mapping", "Complex Components", "File Guide", "Validation"],
  formula_error_count: errorMatches.length,
  status: errorMatches.length === 0 ? "PASS" : "FAIL",
  sheet_inspection: sheets.ndjson ?? String(sheets),
  summary_inspection: summary.ndjson ?? String(summary),
  formula_inspection: formulas.ndjson ?? String(formulas),
};
await fs.writeFile(reportPath, JSON.stringify(report, null, 2), "utf8");
console.log(JSON.stringify({ status: report.status, formula_error_count: report.formula_error_count }, null, 2));
