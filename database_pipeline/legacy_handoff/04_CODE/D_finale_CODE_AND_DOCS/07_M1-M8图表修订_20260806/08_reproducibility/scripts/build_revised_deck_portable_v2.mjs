import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const source = path.join(scriptDir, "build_revised_deck.mjs");
const runtime = path.join(scriptDir, ".build_revised_deck_runtime.mjs");
const revisionRoot = process.env.MEMPRO_REVISION_ROOT ?? path.resolve(scriptDir, "..", "..");
const template = process.env.MEMPRO_PPT_TEMPLATE ?? path.resolve(scriptDir, "..", "reference_template.pptx");
const nodeModules = process.env.CODEX_NODE_MODULES ?? "C:/Users/Administrator/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules";
const artifactTool = pathToFileURL(path.resolve(nodeModules, "@oai", "artifact-tool", "dist", "artifact_tool.mjs")).href;

let code = await fs.readFile(source, "utf8");
code = code.replace('from "@oai/artifact-tool"', `from ${JSON.stringify(artifactTool)}`);
code = code.replace(
  'const sourcePptx = "C:/tmp/mempro_ppt_revision/source.pptx";',
  `const sourcePptx = ${JSON.stringify(template.replaceAll("\\", "/"))};`,
);
code = code.replace(
  'const revisionRoot = "D:/finale/07_M1-M8图表修订_20260806";',
  `const revisionRoot = ${JSON.stringify(revisionRoot.replaceAll("\\", "/"))};`,
);
code = code.replace(
  "if (!shape?.textFrame) throw new Error(`Text shape not found: ${id}`);\n  shape.textFrame.setText(value);",
  "if (!shape) throw new Error(`Text shape not found: ${id}`);\n  shape.text = value;",
);
code = code.replace(
  "await fs.writeFile(outputPptx, new Uint8Array(await pptx.arrayBuffer()));",
  "await pptx.save(outputPptx);",
);
await fs.writeFile(runtime, code, "utf8");
await import(`${pathToFileURL(runtime).href}?v=${Date.now()}`);
