import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const input = "MemPro_V6.3.1_数据库建立流程与科学意义_双页模块重制版.pptx";
const output = "MemPro_V6.3.1_数据库建立流程与科学意义_双页模块重制终版.pptx";
const deck = await PresentationFile.importPptx(await FileBlob.load(input));
const q = await deck.inspect({ kind: "slide,textbox,shape,notes", maxChars: 500000 });
const rows = q.ndjson.trim().split(/\n/).map((s) => JSON.parse(s));
const hit = rows.filter((r) => r.slide === 18 && r.name === "slide-title");
if (hit.length !== 1) throw new Error(`slide-title matches: ${hit.length}`);
const title = deck.resolve(hit[0].id);
title.text = "MemPro 已形成可追溯数据库；下一阶段聚焦复现与 HPC";
title.text.fontSize = 36;
title.text.bold = true;
title.text.color = "#12324A";
title.text.typeface = "Microsoft YaHei";

const slide = deck.slides.items[17];
const png = await deck.export({ slide, format: "png", scale: 1.5 });
await fs.writeFile("dual-module-final/slides/slide-18.png", new Uint8Array(await png.arrayBuffer()));
const layout = await slide.export({ format: "layout" });
await fs.writeFile("dual-module-final/layouts/slide-18.layout.json", await layout.text(), "utf8");
const inspect = await deck.inspect({ kind: "slide,textbox,shape,image,notes", maxChars: 1000000 });
await fs.writeFile("dual-module-final/final-inspect-fixed.ndjson", inspect.ndjson, "utf8");
const pptx = await PresentationFile.exportPptx(deck);
await pptx.save(output);
console.log(JSON.stringify({ output, fixedSlide: 18 }));
