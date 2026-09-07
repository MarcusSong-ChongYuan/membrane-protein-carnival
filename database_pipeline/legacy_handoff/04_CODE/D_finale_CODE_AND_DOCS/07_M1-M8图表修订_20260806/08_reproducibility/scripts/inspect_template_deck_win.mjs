import fs from "node:fs/promises";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const workspace = path.resolve("C:/tmp/mempro_ppt_revision");
const pptx = path.join(workspace, "source.pptx");
const out = path.join(workspace, "template-inspect");

function tar(args, encoding) {
  const r = spawnSync("tar.exe", args, { encoding, maxBuffer: 100 * 1024 * 1024 });
  if (r.status !== 0) throw new Error(String(r.stderr || r.stdout || "tar failed"));
  return r.stdout;
}

async function writeBlob(file, blob) {
  await fs.writeFile(file, new Uint8Array(await blob.arrayBuffer()));
}

await fs.rm(out, { recursive: true, force: true });
await fs.mkdir(path.join(out, "source-slides"), { recursive: true });
await fs.mkdir(path.join(out, "layouts"), { recursive: true });
await fs.mkdir(path.join(out, "assets", "ppt", "media"), { recursive: true });

const presentation = await PresentationFile.importPptx(await FileBlob.load(pptx));
const slides = presentation.slides.items;
const artifacts = [];
for (let i = 0; i < slides.length; i += 1) {
  const n = String(i + 1).padStart(2, "0");
  const png = path.join(out, "source-slides", `source-slide-${n}.png`);
  const layout = path.join(out, "layouts", `source-slide-${n}.layout.json`);
  await writeBlob(png, await presentation.export({ slide: slides[i], format: "png", scale: 1 }));
  await writeBlob(layout, await presentation.export({ slide: slides[i], format: "layout" }));
  artifacts.push({ slide: i + 1, previewPath: png, layoutPath: layout });
}

const names = String(tar(["-tf", pptx], "utf8")).split(/\r?\n/).filter(Boolean);
const media = names.filter((x) => x.startsWith("ppt/media/") && !x.endsWith("/"));
const extracted = [];
for (const entry of media) {
  const bytes = tar(["-xOf", pptx, entry], undefined);
  const target = path.join(out, "assets", "ppt", "media", path.basename(entry));
  await fs.writeFile(target, bytes);
  extracted.push({ entry, path: target, bytes: bytes.length });
}

const fontSet = new Set();
for (const entry of names.filter((x) => /^ppt\/(slides|slideMasters|slideLayouts|theme)\/.*\.xml$/.test(x))) {
  const xml = String(tar(["-xOf", pptx, entry], "utf8"));
  for (const m of xml.matchAll(/\btypeface="([^"]+)"/g)) fontSet.add(m[1]);
}

const inspect = await presentation.inspect({ kind: "slide,textbox,shape,image,table,chart,notes,layout", maxChars: 500000 });
await fs.writeFile(path.join(out, "template-inspect.ndjson"), inspect.ndjson || "", "utf8");
const manifest = {
  sourcePptx: pptx,
  generatedAt: new Date().toISOString(),
  slideCount: slides.length,
  slideArtifacts: artifacts,
  extractedMedia: extracted,
  fonts: [...fontSet].sort(),
  packageParts: { mediaCount: media.length, slideXmlCount: names.filter((x) => /^ppt\/slides\/slide\d+\.xml$/.test(x)).length },
  inspectTruncated: Boolean(inspect.truncated),
};
await fs.writeFile(path.join(out, "template-manifest.json"), JSON.stringify(manifest, null, 2), "utf8");
console.log(JSON.stringify(manifest));
