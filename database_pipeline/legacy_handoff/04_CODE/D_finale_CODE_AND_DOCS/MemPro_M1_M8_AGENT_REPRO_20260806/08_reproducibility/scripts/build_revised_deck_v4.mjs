import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";

const source = "C:/tmp/mempro_ppt_revision/build_revised_deck.mjs";
const runtime = "C:/tmp/mempro_ppt_revision/build_revised_deck_runtime_v4.mjs";
let code = await fs.readFile(source, "utf8");
code = code.replace(
  "if (!shape?.textFrame) throw new Error(`Text shape not found: ${id}`);\n  shape.textFrame.setText(value);",
  "if (!shape) throw new Error(`Text shape not found: ${id}`);\n  shape.text = value;",
);
code = code.replace(
  "await fs.writeFile(outputPptx, new Uint8Array(await pptx.arrayBuffer()));",
  "await pptx.save(outputPptx);",
);
await fs.writeFile(runtime, code, "utf8");
await import(`${pathToFileURL(path.resolve(runtime)).href}?v=${Date.now()}`);
