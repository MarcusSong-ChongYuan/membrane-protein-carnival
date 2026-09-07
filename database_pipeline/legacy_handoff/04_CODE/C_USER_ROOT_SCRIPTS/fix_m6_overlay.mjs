import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const input = "MemPro_V6.3.1_M1-M8_full_figure.pptx";
const output = "MemPro_V6.3.1_M1-M8_full_figure_fixed.pptx";
const deck = await PresentationFile.importPptx(await FileBlob.load(input));

const m6 = deck.slides.items[13].images.add({
  blob: await fs.readFile("assets/M6.png"),
  fit: "contain",
  alt: "M6 binding evidence and binding-site landscape",
  name: "main-figure-M6-corrected",
});
m6.position = { left: 167, top: 122, width: 946, height: 564 };

await fs.mkdir("final-fixed/slides", { recursive: true });
await fs.mkdir("final-fixed/layouts", { recursive: true });
for (let i = 0; i < deck.slides.items.length; i += 1) {
  const slide = deck.slides.items[i];
  const n = String(i + 1).padStart(2, "0");
  const png = await deck.export({ slide, format: "png", scale: 1.5 });
  await fs.writeFile(`final-fixed/slides/slide-${n}.png`, new Uint8Array(await png.arrayBuffer()));
  const layout = await slide.export({ format: "layout" });
  await fs.writeFile(`final-fixed/layouts/slide-${n}.layout.json`, await layout.text(), "utf8");
}
const montage = await deck.export({ format: "webp", montage: true, scale: 0.55 });
await fs.writeFile("final-fixed/montage.webp", new Uint8Array(await montage.arrayBuffer()));
const inspect = await deck.inspect({ kind: "slide,textbox,shape,image,notes", maxChars: 1000000 });
await fs.writeFile("final-fixed/final-inspect.ndjson", inspect.ndjson, "utf8");
const pptx = await PresentationFile.exportPptx(deck);
await pptx.save(output);
console.log(JSON.stringify({ output, slides: deck.slides.items.length }));
