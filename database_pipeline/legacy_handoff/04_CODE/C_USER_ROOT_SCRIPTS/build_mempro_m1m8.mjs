import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const INPUT = "starter_m1m8.pptx";
const OUTPUT = "MemPro_V6.3.1_M1-M8_full_figure.pptx";
const FIG = {
  4: "assets/M1.png",
  6: "assets/M2.png",
  8: "assets/M3.png",
  10: "assets/M4.png",
  12: "assets/M5.png",
  14: "assets/M6.png",
  18: "assets/M8.png",
};
const ALT = {
  4: "M1 database architecture, integration workflow, release layers, and quality control",
  6: "M2 membrane proteome landscape and classifications",
  8: "M3 expression and subcellular localization atlas",
  10: "M4 protein-disease association landscape",
  12: "M5 compound identity and chemical space",
  14: "M6 binding evidence and binding-site coverage",
  16: "M7 negative evidence, conflicts, review queues, and release quality",
  18: "M8 docking prioritization and readiness",
};

const deck = await PresentationFile.importPptx(await FileBlob.load(INPUT));
const inspected = await deck.inspect({
  kind: "slide,textbox,shape,image,notes",
  include: "id,slide,name,title,text,textPreview,bbox,bboxUnit,alt",
  maxChars: 1000000,
});
const rows = inspected.ndjson.trim().split(/\n/).filter(Boolean).map((s) => JSON.parse(s));

function rowFor(slide, name, kind = null) {
  const hits = rows.filter((r) => r.slide === slide && r.name === name && (!kind || r.kind === kind));
  if (hits.length !== 1) throw new Error(`Expected one ${kind ?? "object"} named ${name} on slide ${slide}, found ${hits.length}`);
  return hits[0];
}

function objectFor(slide, name, kind = null) {
  return deck.resolve(rowFor(slide, name, kind).id);
}

function setText(slide, name, value) {
  objectFor(slide, name).text = value;
}

function setNotes(slide, value) {
  const hits = rows.filter((r) => r.slide === slide && r.kind === "notes");
  if (hits.length !== 1) throw new Error(`Expected one notes object on slide ${slide}, found ${hits.length}`);
  deck.resolve(hits[0].id).setText(value);
}

function mainImage(slide) {
  const images = rows.filter((r) => r.slide === slide && r.kind === "image");
  if (!images.length) throw new Error(`No image on slide ${slide}`);
  images.sort((a, b) => (b.bbox?.[2] ?? 0) - (a.bbox?.[2] ?? 0));
  return deck.resolve(images[0].id);
}

for (const [slideText, path] of Object.entries(FIG)) {
  const slide = Number(slideText);
  const image = mainImage(slide);
  image.replace({
    blob: await fs.readFile(path),
    contentType: "image/png",
    alt: ALT[slide],
    fit: "contain",
  });
  image.frame = { left: 167, top: 122, width: 946, height: 564 };
  image.lockAspectRatio = false;
}

// Slides 14 and 16 were duplicated from the same source figure slide, so their
// inherited image parts are shared. Add M7 as a new top image to keep M6 intact.
const m7Image = deck.slides.items[15].images.add({
  blob: await fs.readFile("assets/M7.png"),
  fit: "contain",
  alt: ALT[16],
  name: "main-figure-M7",
});
m7Image.position = { left: 167, top: 122, width: 946, height: 564 };

// M7 prerequisite page: retain the exact six-box module geometry from slide 13.
setText(15, "section-kicker", "14 路 璐ㄩ噺鎺у埗");
setText(15, "slide-title", "姝ｈ礋璇佹嵁銆佸啿绐併€佹湭鏄犲皠璁板綍涓庡彂甯冭川閲忓垎鍒鐞?);
setText(15, "footer-page", "15");
const m7Boxes = [
  ["lim-head-0-0", "闃虫€ц瘉鎹?, "lim-body-0-0", "active銆佸畾閲忕粨鍚堟垨缁撴瀯浣嶇偣鏀寔鍏崇郴锛涙寜 BE1鈥揃E3/MECH 鍒嗗眰锛屽苟淇濈暀鍘熷疄楠屾潯浠躲€?],
  ["lim-head-0-1", "闃存€ц瘉鎹?, "lim-body-0-1", "inactive/negative 鐙珛淇濆瓨锛涘畠琛ㄧず鐗瑰畾瀹為獙鏉′欢涓嬫湭妫€鍑烘椿鎬э紝涓嶇瓑鍚屼簬姘镐箙鈥滀笉缁撳悎鈥濄€?],
  ["lim-head-0-2", "姝ｈ礋鍐茬獊", "lim-body-0-2", "鍚屼竴 canonical 铔嬬櫧鈥斿皬鍒嗗瓙鍑虹幇鐩稿弽缁撴灉鏃朵笉浜掔浉瑕嗙洊锛涙寜 assay銆佸墏閲忋€佹瀯寤轰綋鍜屾潯浠惰В閲娿€?],
  ["lim-head-1-0", "韬唤鏄犲皠鐘舵€?, "lim-body-1-0", "缁撴瀯鍞竴涓斿閿畬鏁磋€呰繘鍏ユ寮忓眰锛涙棤娉曞彲闈犺В鏋愮殑璁板綍杩涘叆 review queue锛屼笉寮鸿鍚堝苟涔熶笉闈欓粯鍒犻櫎銆?],
  ["lim-head-1-1", "鍘婚噸涓庤氨绯?, "lim-body-1-1", "鍚屼竴 PDB銆丄ID/SID/CID銆佽鏂囨垨瀹為獙琛ㄦ牸寤虹珛璋辩郴閿紱琚涓暟鎹簱杞浇涓嶉噸澶嶈浣滅嫭绔嬪疄楠屻€?],
  ["lim-head-1-2", "鍙戝竷 QA", "lim-body-1-2", "涓婚敭銆佸閿€佹潵婧愭暟閲忋€乵anifest銆佹枃浠跺搱甯屽拰鎶芥牱澶嶆牳鍏ㄩ儴閫氳繃鍚庯紝璁板綍鎵嶈繘鍏ュ喕缁?release銆?],
];
for (const [hn, ht, bn, bt] of m7Boxes) {
  setText(15, hn, ht);
  setText(15, bn, bt);
}
setNotes(15, `M7瑙勫垯椤靛厛鍖哄垎鍥涗欢浜嬶細闃虫€у拰闃存€ф槸瀹為獙缁撹锛涘啿绐佹槸鍚屼竴韬唤涓嬬殑鏉′欢宸紓锛涙湭鏄犲皠鏄韩浠戒笉纭畾锛決A鍐冲畾璁板綍鑳藉惁杩涘叆姝ｅ紡鍙戝竷灞傘€傛暟鎹簱淇濈暀杩欎簺宸紓锛屼笉鐢ㄤ竴涓粨璁鸿鐩栧彟涓€涓€俓n\n[Sources]\n- Internal: D:\\finale\\02_灞曠ず鍥捐〃_V6.2\\png\\M7_negative_conflicts_quality.png\n- Internal: MemPro negative-evidence, conflict, review-queue, manifest and QA tables`);

// M7 main figure page.
setText(16, "section-kicker", "15 路 璐ㄩ噺鎺у埗锝滀富鍥?);
setText(16, "slide-title", "璐熻瘉鎹槧灏勩€佹璐熷啿绐併€佸鏍搁槦鍒椾笌鍙戝竷璐ㄩ噺");
setText(16, "footer-page", "16");

// Docking pair moves two pages later; preserve all user-edited body copy.
setText(17, "section-kicker", "16 路 Docking");
setText(17, "footer-page", "17");
setText(18, "section-kicker", "17 路 Docking锝滀富鍥?);
setText(18, "footer-page", "18");

const notesBySlide = {
  4: `M1浠庡師濮嬫暟鎹簱杩涘叆鏁村悎灞傚紑濮嬶紝渚濇灞曠ず韬唤鏄犲皠銆乧anonical瀹炰綋銆佽法鏉ユ簮鍘婚噸銆佽瘉鎹垎灞傘€佸啿绐?澶嶆牳灞備互鍙婂喕缁撳彂甯冦€傞噸鐐规槸锛歩ntegration涓嶆槸绠€鍗曟嫾琛紝鑰屾槸鎶婃瘡鏉″叧绯诲彉鎴愬彲杩芥函銆佸彲鏍搁獙鐨勫彂甯冭褰曘€俓n\n[Sources]\n- Asset: D:\\finale\\02_灞曠ず鍥捐〃_V6.2\\png\\M1_database_architecture_quality.png\n- Internal: MemPro V6.2 frozen visualization snapshot; V6.3.1 narrative rules retained in this deck`,
  6: `M2闆嗕腑灞曠ず鑶滆泲鐧藉簱鐨凙BC鑶滅粨鍚堢被鍨嬨€丒1鈥揈3璇佹嵁鎶婃彙搴︺€佽法鑶滄鏁般€佸姛鑳界被鍒拰鏉ユ簮瑕嗙洊銆傝鍥炬椂鍏堣В閲夾BC涓嶦绛夌骇鏄袱鏉′笉鍚岃酱锛屽啀姣旇緝鍚勭被铔嬬櫧鐨勬暟閲忓拰鏈垎绫荤┖闂淬€俓n\n[Sources]\n- Asset: D:\\finale\\02_灞曠ず鍥捐〃_V6.2\\png\\M2_membrane_proteome_landscape.png\n- Internal: MemPro V6.2 frozen visualization snapshot`,
  8: `M3灞曠ずHPA缁勭粐/缁嗚優琛ㄨ揪銆佽〃杈惧箍搴︺€佷簹缁嗚優瀹氫綅鍜岀己澶辩姸鎬併€?鍊煎繀椤诲尯鍒唌apped zero涓巙nmapped/missing锛岀儹鍥鹃鑹茶〃绀虹浉瀵硅〃杈炬ā寮忚€屼笉鏄法铔嬬櫧缁濆琛ㄨ揪閲忋€俓n\n[Sources]\n- Asset: D:\\finale\\02_灞曠ず鍥捐〃_V6.2\\png\\M3_expression_localization_atlas.png\n- Human Protein Atlas 25.1; MemPro expression/localization module`,
  10: `M4灞曠ず铔嬬櫧鈥旂柧鐥呭叧绯昏妯°€佺柧鐥呯被鍒€佽瘉鎹笭閬撳強楂樿繛鎺ヨ妭鐐广€傛棫鍥剧敤浜庝繚鐣欏畬鏁寸殑瑙ｉ噴鎬х粺璁¤瑙掞紱姝ｅ紡绋胯В閲婃椂搴斾互MONDO exact-only缁熶竴銆丱pen Targets娌荤枟棰嗗煙鍜孌O/Uberon澶氭爣绛炬槧灏勪负鍑嗐€俓n\n[Sources]\n- Asset: D:\\finale\\02_灞曠ず鍥捐〃_V6.2\\png\\M4_disease_association_landscape.png\n- Open Targets 26.06; UniProtKB; MONDO/DO/Uberon mappings`,
  12: `M5灞曠ず灏忓垎瀛愯韩浠藉眰銆佹潵婧愪笌鐘舵€佷氦闆嗐€佺墿鍖栨€ц川鍙婂寲瀛︾┖闂淬€傚悕绉板彧鐢ㄤ簬妫€绱紱canonical parent涓巈xact form鐢辩粨鏋勬爣鍑嗗寲鍐冲畾銆傜墿鍖栨€ц川鏄暟鎹簱鎻忚堪鍜孌ocking鍑嗗鍙傝€冿紝涓嶆槸鑽晥鎴栧彛鏈嶆垚鑽€х殑纭€х粨璁恒€俓n\n[Sources]\n- Asset: D:\\finale\\02_灞曠ず鍥捐〃_V6.2\\png\\M5_chemical_identity_space.png\n- Internal: MemPro compound parent/form and QC tables`,
  14: `M6灞曠ず缁撳悎璇佹嵁鏉ユ簮銆丅E绛夌骇銆佺粨鏋勪綅鐐瑰彲鐢ㄦ€у強铔嬬櫧绫诲埆鈥斿寲瀛︾被鍒叧绯汇€備綅鐐硅鐩栧彧閽堝鍏锋湁PDB銆侀摼銆侀厤浣撳拰娈嬪熀鍧愭爣婧簮鐨刾air锛屼笉鑳借В閲婁负鍏ㄩ儴铔嬬櫧鈥斿皬鍒嗗瓙鍏崇郴閮芥湁鍙敤浜嶥ocking鐨勫彛琚嬨€俓n\n[Sources]\n- Asset: D:\\finale\\02_灞曠ず鍥捐〃_V6.2\\png\\M6_binding_evidence_sites.png\n- Internal: MemPro binding evidence and structural-site tables`,
  16: `M7鎶婅礋璇佹嵁鏄犲皠銆佹湭瑙ｆ瀽璁板綍銆佹璐熷啿绐佸拰鍙戝竷QA鏀惧湪鍚屼竴寮犲浘涓€傚啿绐佸苟涓嶈嚜鍔ㄨ〃绀烘煇涓€鏁版嵁搴撻敊璇紱涓嶅悓娴撳害銆佹瀯寤轰綋銆佺粓鐐广€佺墿绉嶆垨瀹為獙浣撶郴鍙骇鐢熺湡瀹炵殑鐘舵€佷緷璧栫粨鏋溿€俓n\n[Sources]\n- Asset: D:\\finale\\02_灞曠ず鍥捐〃_V6.2\\png\\M7_negative_conflicts_quality.png\n- Internal: MemPro negative evidence, conflicts, review queues and QA reports`,
  18: `M8灞曠ずDocking鍊欓€夌殑绛涢€夋紡鏂椼€丷0/D1/D2/D3浼樺厛绾с€佺粨鏋勪笌浣嶇偣鍑嗗搴︿互鍙婅绠楄妯°€傚€欓€夊悕鍗曞畬鎴愪笉绛変簬Docking宸叉墽琛岋紱浠嶉渶鍙椾綋銆侀厤浣撱€乥ox銆乺edocking鍜岀湡瀹濰PC璋冨害銆俓n\n[Sources]\n- Asset: D:\\finale\\02_灞曠ず鍥捐〃_V6.2\\png\\M8_docking_prioritization.png\n- Internal: MemPro docking pilot manifest and readiness rules`,
};
for (const [slide, note] of Object.entries(notesBySlide)) setNotes(Number(slide), note);

await fs.mkdir("final/slides", { recursive: true });
await fs.mkdir("final/layouts", { recursive: true });
for (let i = 0; i < deck.slides.items.length; i += 1) {
  const slide = deck.slides.items[i];
  const n = String(i + 1).padStart(2, "0");
  const png = await deck.export({ slide, format: "png", scale: 1.5 });
  await fs.writeFile(`final/slides/slide-${n}.png`, new Uint8Array(await png.arrayBuffer()));
  const layout = await slide.export({ format: "layout" });
  await fs.writeFile(`final/layouts/slide-${n}.layout.json`, await layout.text(), "utf8");
}
const montage = await deck.export({ format: "webp", montage: true, scale: 0.55 });
await fs.writeFile("final/montage.webp", new Uint8Array(await montage.arrayBuffer()));
const finalInspect = await deck.inspect({ kind: "slide,textbox,shape,image,notes", maxChars: 1000000 });
await fs.writeFile("final/final-inspect.ndjson", finalInspect.ndjson, "utf8");
const pptx = await PresentationFile.exportPptx(deck);
await pptx.save(OUTPUT);
console.log(JSON.stringify({ output: OUTPUT, slides: deck.slides.items.length }));

