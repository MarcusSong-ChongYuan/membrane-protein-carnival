import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const INPUT="template-starter.pptx";
const OUTPUT="MemPro_V6.3.1_数据库建立流程与科学意义_数据重绘一致版.pptx";
const FIG={
  4:"figures/M1_database_pipeline_redrawn.png",
  6:"figures/M2_protein_classification_redrawn.png",
  8:"figures/M3_expression_localization_redrawn.png",
  10:"figures/M4_disease_ontology_redrawn.png",
  12:"figures/M5_compound_upset_scaffold_redrawn.png",
  14:"figures/M6_binding_evidence_redrawn.png",
  16:"figures/M7_negative_conflict_qa_redrawn.png",
  18:"figures/M8_docking_pilot_redrawn.png",
};

const deck=await PresentationFile.importPptx(await FileBlob.load(INPUT));
const inspect=await deck.inspect({kind:"slide,textbox,shape,image,notes",include:"id,slide,name,text,textPreview,bbox",maxChars:1000000});
const rows=inspect.ndjson.trim().split(/\n/).filter(Boolean).map(JSON.parse);
function obj(slide,name){const h=rows.filter(r=>r.slide===slide&&r.name===name);if(h.length!==1)throw new Error(`target ${slide}/${name}: ${h.length}`);return deck.resolve(h[0].id);}
function setText(slide,name,value){obj(slide,name).text=value;}
function setNotes(slide,value){const h=rows.filter(r=>r.slide===slide&&r.kind==="notes");if(h.length!==1)throw new Error(`notes ${slide}: ${h.length}`);deck.resolve(h[0].id).setText(value);}
function addShape(slide,{x,y,w,h,fill="#00000000",line="#00000000",lineWidth=0,name}){return slide.shapes.add({geometry:"rect",name,position:{left:x,top:y,width:w,height:h},fill:{type:"solid",color:fill},line:{style:"solid",fill:line,width:lineWidth}});}
function addText(slide,value,{x,y,w,h,size=18,color="#12324A",bold=false,name}){const sh=addShape(slide,{x,y,w,h,name});sh.text=value;sh.text.fontSize=size;sh.text.color=color;sh.text.bold=bold;sh.text.typeface="Microsoft YaHei";sh.text.verticalAlignment="middle";sh.text.insets={left:4,right:4,top:2,bottom:2};return sh;}

for(const [sn,path] of Object.entries(FIG)){
  const slideNo=Number(sn), slide=deck.slides.items[slideNo-1];
  const im=slide.images.add({blob:await fs.readFile(path),fit:"contain",alt:`MemPro M${slideNo===4?1:slideNo/2} data-driven scientific figure`,name:`redrawn-M${slideNo===4?1:slideNo/2}`});
  im.position={left:36,top:128,width:1208,height:524};
}

// Exact title-to-figure alignment.
setText(10,"slide-title","疾病身份、本体层级、治疗领域、解剖映射与证据渠道");
setText(12,"slide-title","真实 UpSet 交集、Lipinski 描述、Bemis–Murcko 骨架与结构类别");
setText(13,"slide-title","BE等级、实验谱系、来源贡献和结构位点分别管理");
setText(13,"lim-head-1-2","位点可用性");
setText(13,"lim-body-1-2","只有具备PDB、链、配体、残基坐标或明确口袋来源的记录，才计入结构位点覆盖与Docking准备度。");
setText(14,"slide-title","来源、BE1–BE3、谱系覆盖与结合位点可用性");

// M7 prerequisite page.
setText(15,"section-kicker","14 · 质量控制");
setText(15,"slide-title","负证据、身份复核、正负冲突与发布质量分别管理");
setText(15,"footer-page","15");
const q=[
  ["lim-head-0-0","负证据","lim-body-0-0","inactive/negative独立保存；它表示特定实验条件下未检出活性，不等同于永久“不结合”。"],
  ["lim-head-0-1","身份解析","lim-body-0-1","full InChIKey或结构唯一时才自动映射；盐型、质子化和未定立体化学进入复核层。"],
  ["lim-head-0-2","正负冲突","lim-body-0-2","同一canonical蛋白—小分子出现相反结果时不覆盖；保留assay、浓度、构建体和实验条件。"],
  ["lim-head-1-0","复核队列","lim-body-1-0","无法可靠解析的记录不删除、不强并，按原因进入review queue并保留原始CID/AID/SID。"],
  ["lim-head-1-1","去重规则","lim-body-1-1","同一target＋canonical compound＋AID/SID只保留统一关系；转载记录不重复计为独立实验。"],
  ["lim-head-1-2","发布质量门","lim-body-1-2","主外键、来源对账、manifest、哈希、冲突标记和抽样复核全部通过后，才进入冻结release。"],
];
for(const [h,ht,b,bt] of q){setText(15,h,ht);setText(15,b,bt);}
setText(16,"section-kicker","15 · 质量控制｜主图");
setText(16,"slide-title","负证据解析漏斗、复核原因、正负冲突与发布质量门");
setText(16,"footer-page","16");

// Docking pair shifted after M7 and completed with the missing sixth rule.
setText(17,"section-kicker","16 · Docking");
setText(17,"slide-title","Docking名单分层与生产前门控");
setText(17,"footer-page","17");
const s17=deck.slides.items[16];
addShape(s17,{x:664,y:476,w:10,h:110,fill:"#775DA6",line:"#775DA6",name:"lim-accent-1-2-added"});
addText(s17,"生产前四道门",{x:692,y:476,w:500,h:38,size:21,color:"#12324A",bold:true,name:"lim-head-1-2-added"});
addText(s17,"受体链/组装体、配体微状态、搜索box与R0 redocking均通过后，才提交D1–D3生产任务。",{x:692,y:524,w:500,h:78,size:15.5,color:"#556575",name:"lim-body-1-2-added"});
setText(18,"section-kicker","17 · Docking｜主图");
setText(18,"slide-title","6,019个pilot pair的分层、结构等级与生产前门控");
setText(18,"footer-page","18");

const notes={
  4:`M1与规则页完全一致：从冻结来源开始，经身份标准化、canonical实体、证据/冲突处理，进入四张主表、复核队列、负证据表和可复现发布层。\n\n[Sources]\n- Internal: MemPro V6.3/V6.3.1 integration policy and frozen source registry`,
  6:`M2把ABC膜结合方式与E1–E3证据等级分开，同时展示五轴交叉注释覆盖和膜角色分布。“家族明确，角色未决”单独报告，不再作为模糊other。\n\n[Sources]\n- Internal: human_membrane_protein_master_v6_3_candidate.tsv.gz; n=10,997`,
  8:`M3明确区分mapped detected、mapped zero和unmapped/missing；组织广度采用箱线摘要；热图为行z-score，仅表示类内相对组织偏好。\n\n[Sources]\n- Human Protein Atlas 25.1\n- Internal: expression_location_summary_v2.tsv.gz and protein_subcellular_localization_v2.tsv.gz`,
  10:`M4只将direct MONDO或唯一official exact/equivalent映射进入canonical层。治疗领域来自Open Targets 26.06多标签；解剖映射来自MONDO/DO到Uberon，并区分asserted与inferred。\n\n[Sources]\n- Open Targets Platform 26.06\n- MONDO, Disease Ontology and Uberon frozen mappings\n- UniProtKB 2026_02`,
  12:`M5A是真正的UpSet：上方柱表示交集数量，下方点阵表示交集包含哪些状态。M5C/D使用全量Bemis–Murcko骨架，不再使用Morgan降维。Lipinski仅作描述，不作为收录或Docking硬门槛。\n\n[Sources]\n- Internal: small_molecule_master_v1_3.tsv and compound scaffold module\n- RDKit 2026.03.5; n=899,591 core-QC compounds`,
  14:`M6只展示正结合证据：来源×BE等级、谱系覆盖和结合位点可用性。数据库贡献数与独立实验数分开；位点覆盖仅限具有结构/链/配体/残基或明确口袋来源的记录。\n\n[Sources]\n- Internal: binding_evidence_master_v6_2.tsv.gz; V6.3 lineage audit; binding_site_instances_v6_2.tsv.gz`,
  15:`本页将负证据与正证据模块分开。未映射记录不删除，正负冲突不覆盖；结构唯一与发布QA共同决定记录进入正式层还是复核层。\n\n[Sources]\n- Internal: negative evidence identity policy and V6.3.1 review queues`,
  16:`M7展示10,246,948条原始负证据的解析结果：7,630,659条进入映射正式层，2,615,046条留在复核层，42,231条标记正负冲突。\n\n[Sources]\n- Internal: NEGATIVE_EVIDENCE_RESOLUTION_V62_VALIDATION.json\n- Internal: negative_unmapped_review_summary_v631.tsv`,
  18:`M8严格使用6,019条pilot清单，而非192,592条完整候选表。R0用于redocking，D1–D3用于生产优先级；当前执行状态均为not_run。\n\n[Sources]\n- Internal: docking_pilot_manifest_v631.tsv; exactly 6,019 pairs and 1,773 targets`,
};
for(const [n,v] of Object.entries(notes))setNotes(Number(n),v);

await fs.mkdir("final/slides",{recursive:true});await fs.mkdir("final/layouts",{recursive:true});
for(let i=0;i<deck.slides.items.length;i++){
  const slide=deck.slides.items[i],n=String(i+1).padStart(2,"0");
  const png=await deck.export({slide,format:"png",scale:1.5});await fs.writeFile(`final/slides/slide-${n}.png`,new Uint8Array(await png.arrayBuffer()));
  const lay=await slide.export({format:"layout"});await fs.writeFile(`final/layouts/slide-${n}.layout.json`,await lay.text(),"utf8");
}
const montage=await deck.export({format:"webp",montage:true,scale:.55});await fs.writeFile("final/montage.webp",new Uint8Array(await montage.arrayBuffer()));
const finalInspect=await deck.inspect({kind:"slide,textbox,shape,image,notes",maxChars:1000000});await fs.writeFile("final/final-inspect.ndjson",finalInspect.ndjson,"utf8");
const pptx=await PresentationFile.exportPptx(deck);await pptx.save(OUTPUT);
console.log(JSON.stringify({output:OUTPUT,slides:deck.slides.items.length}));
