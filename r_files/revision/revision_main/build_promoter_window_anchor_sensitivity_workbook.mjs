import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const sourceDir = process.env.SENSITIVITY_SOURCE_DIR;
const outputPath = process.env.SENSITIVITY_WORKBOOK_PATH;
const previewDir = process.env.SENSITIVITY_PREVIEW_DIR;

if (!sourceDir || !outputPath || !previewDir) {
  throw new Error(
    "SENSITIVITY_SOURCE_DIR, SENSITIVITY_WORKBOOK_PATH, and " +
      "SENSITIVITY_PREVIEW_DIR are required.",
  );
}

const cases = [
  { width: 2001, mode: "interval", sheet: "W2001_Interval" },
  { width: 2001, mode: "midpoint", sheet: "W2001_Midpoint" },
  { width: 3001, mode: "interval", sheet: "W3001_Interval" },
  { width: 3001, mode: "midpoint", sheet: "W3001_Midpoint" },
  { width: 4001, mode: "interval", sheet: "W4001_Interval" },
  { width: 4001, mode: "midpoint", sheet: "W4001_Midpoint" },
];

function parseTsv(text, sourcePath) {
  const lines = text.trim().split(/\r?\n/);
  if (lines.length < 2) {
    throw new Error(`No gene rows found in ${sourcePath}`);
  }

  const header = lines[0].split("\t");
  const expected = ["ensembl_gene_id", "gene_symbol", "n"];
  if (header.length !== expected.length || header.some((v, i) => v !== expected[i])) {
    throw new Error(
      `Unexpected columns in ${sourcePath}: ${header.join(", ")}`,
    );
  }

  const rows = lines.slice(1).map((line, index) => {
    const fields = line.split("\t");
    const n = Number(fields[2]);
    if (fields.length !== 3 || !fields[0] || !Number.isInteger(n) || n < 1) {
      throw new Error(`Invalid row ${index + 2} in ${sourcePath}`);
    }
    return [fields[0], fields[1] || fields[0], n];
  });

  const ids = rows.map((row) => row[0]);
  if (new Set(ids).size !== ids.length) {
    throw new Error(`Duplicate Ensembl gene IDs found in ${sourcePath}`);
  }
  return [expected, ...rows];
}

const workbook = Workbook.create();
const rowCounts = {};

for (const sensitivityCase of cases) {
  const stem =
    `promoter_window_${sensitivityCase.width}bp_` +
    `anchor_${sensitivityCase.mode}_atac_overlap_50bp_` +
    "gene_counts_for_workbook.tsv";
  const sourcePath = path.join(sourceDir, stem);
  const matrix = parseTsv(await fs.readFile(sourcePath, "utf8"), sourcePath);
  const sheet = workbook.worksheets.add(sensitivityCase.sheet);
  const range = sheet.getRangeByIndexes(0, 0, matrix.length, 3);

  range.values = matrix;
  range.format.font = { name: "Aptos", size: 10, color: "#17212B" };
  range.format.rowHeight = 18;

  const header = sheet.getRange("A1:C1");
  header.format.fill = "#1F4E78";
  header.format.font = { name: "Aptos", size: 10, bold: true, color: "#FFFFFF" };
  header.format.horizontalAlignment = "center";
  header.format.verticalAlignment = "center";
  header.format.rowHeight = 24;
  header.format.borders = {
    bottom: { style: "medium", color: "#163A5C" },
  };

  sheet.getRange(`A2:B${matrix.length}`).format.horizontalAlignment = "left";
  sheet.getRange(`C2:C${matrix.length}`).format.horizontalAlignment = "right";
  sheet.getRange(`C2:C${matrix.length}`).format.numberFormat = "#,##0";
  sheet.getRange(`A1:A${matrix.length}`).format.columnWidth = 23;
  sheet.getRange(`B1:B${matrix.length}`).format.columnWidth = 20;
  sheet.getRange(`C1:C${matrix.length}`).format.columnWidth = 10;
  sheet.freezePanes.freezeRows(1);
  sheet.showGridLines = false;

  rowCounts[sensitivityCase.sheet] = matrix.length - 1;
}

const inspection = await workbook.inspect({
  kind: "sheet,table",
  include: "id,name,values",
  maxChars: 12000,
  tableMaxRows: 5,
  tableMaxCols: 3,
});
console.log(inspection.ndjson);

const errorScan = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(errorScan.ndjson);

await fs.mkdir(previewDir, { recursive: true });
for (const sensitivityCase of cases) {
  const preview = await workbook.render({
    sheetName: sensitivityCase.sheet,
    range: "A1:C20",
    scale: 1.5,
    format: "png",
  });
  const bytes = new Uint8Array(await preview.arrayBuffer());
  await fs.writeFile(
    path.join(previewDir, `${sensitivityCase.sheet}.png`),
    bytes,
  );
}

await fs.mkdir(path.dirname(outputPath), { recursive: true });
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(JSON.stringify({ outputPath, rowCounts }));
