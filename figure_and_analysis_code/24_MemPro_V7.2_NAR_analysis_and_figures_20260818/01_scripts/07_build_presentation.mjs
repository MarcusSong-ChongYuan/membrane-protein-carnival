import fs from "node:fs/promises";
import path from "node:path";
import { Presentation, PresentationFile } from "@oai/artifact-tool";

const ROOT = String.raw`D:\finale\24_MemPro_V7.2_NAR_analysis_and_figures_20260818`;
const OUT = path.join(ROOT, "07_presentation");
const MAIN = path.join(ROOT, "04_figures", "main", "png");
const GRAPHICS = path.join(ROOT, "05_graphics");
const RENDER = path.join(OUT, "rendered");
const W = 1280, H = 720;
const C = { dark: "#26323a", mid: "#69757d", blue: "#7b95c6", cyan: "#49c2d9", ltblue: "#a1d8e8", green: "#67a583", pale: "#d0e2c0", yellow: "#fded95", peach: "#ffc1a6", orange: "#f47254", red: "#c85e62", light: "#f4f7f8", white: "#ffffff" };

async function bytes(p) { const b = await fs.readFile(p); return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength); }
async function writeBlob(p, blob) { await fs.writeFile(p, new Uint8Array(await blob.arrayBuffer())); }

function box(slide, left, top, width, height, fill, radius="rounded-xl") {
  return slide.shapes.add({ geometry:"roundRect", position:{left,top,width,height}, fill, line:{style:"solid",fill:"none",width:0}, borderRadius:radius });
}
function text(slide, value, left, top, width, height, style={}) {
  const s=slide.shapes.add({geometry:"textbox",position:{left,top,width,height},fill:"none",line:{style:"solid",fill:"none",width:0}});
  s.text=value; s.text.style={fontFamily:"Arial",fontSize:style.fontSize??20,bold:style.bold??false,color:style.color??C.dark,alignment:style.alignment??"left",verticalAlignment:style.verticalAlignment??"middle"}; return s;
}
function title(slide, value, eyebrow) {
  text(slide,eyebrow.toUpperCase(),64,36,460,24,{fontSize:13,bold:true,color:C.blue});
  text(slide,value,64,62,1152,58,{fontSize:36,bold:true,color:C.dark});
  slide.shapes.add({geometry:"rect",position:{left:64,top:126,width:82,height:5},fill:C.orange,line:{style:"solid",fill:"none",width:0}});
}
function notes(slide, body, sources=[]) {
  slide.speakerNotes.textFrame.setText(`${body}\n\n[Sources]\n${sources.map(s=>`- ${s}`).join("\n")}`);
  slide.speakerNotes.setVisible(true);
}
async function addImage(slide,p,pos,alt) { slide.images.add({blob:await bytes(p),contentType:"image/png",alt,fit:"contain",position:pos}); }
function bulletBlock(slide, items, left, top, width, lineH=58) {
  items.forEach((it,i)=>{ slide.shapes.add({geometry:"ellipse",position:{left,top:top+i*lineH+10,width:12,height:12},fill:i===0?C.orange:C.blue,line:{style:"solid",fill:"none",width:0}}); text(slide,it,left+24,top+i*lineH,width-24,lineH-4,{fontSize:19,color:C.dark}); });
}
function addFooter(slide, n) { text(slide,`MemPro V7.2 · NAR analysis narrative`,64,684,700,20,{fontSize:11,color:C.mid}); text(slide,String(n),1180,684,36,20,{fontSize:11,color:C.mid,alignment:"right"}); }

async function main(){
  await fs.mkdir(RENDER,{recursive:true});
  const deck=Presentation.create({slideSize:{width:W,height:H}});
  let n=0;
  const add=()=>{const s=deck.slides.add();s.background.fill=C.white;n++;addFooter(s,n);return s;};

  let s=add();
  text(s,"MemPro V7.2",70,104,650,72,{fontSize:54,bold:true,color:C.blue});
  text(s,"An evidence-aware membrane-protein and small-molecule knowledge resource",70,184,780,92,{fontSize:30,bold:true,color:C.dark});
  text(s,"NAR database-manuscript analysis and figure narrative",72,300,640,40,{fontSize:20,color:C.mid});
  await addImage(s,path.join(GRAPHICS,"Section_Graphic_Integrated.png"),{left:740,top:105,width:455,height:455},"Integrated MemPro section graphic");
  text(s,"Frozen source: V7.2 · derivative analyses only",72,582,650,30,{fontSize:16,color:C.red});
  notes(s,"开场：这套图不再按‘数据库有多少行’罗列，而是按可追溯整合、功能组织、化学偏好、证据可靠性、疾病背景和发现假设推进。",["MemPro V7.2 frozen release"]);

  s=add(); title(s,"MemPro links traceable evidence to drug-discovery questions","Graphical abstract");
  await addImage(s,path.join(GRAPHICS,"Graphical_Abstract_MemPro_V72_internal_draft.png"),{left:70,top:150,width:1140,height:500},"MemPro graphical abstract internal draft");
  notes(s,"从左到右讲：公共数据层经过身份标准化、去重和证据分级，形成多层核心，最后产生偏好分析、证据解释和化学研究不足靶点三类输出。该图是内部设计稿，投稿前再核查版权与AI图像政策。",["MemPro V7.2","Figure contracts"]);

  s=add(); title(s,"Six figures build one cumulative argument","Manuscript logic");
  const steps=["1 · Traceable integration","2 · Five-axis protein organization","3 · Chemical and target preference","4 · Evidence independence","5 · Disease and anatomy","6 · Cross-module hypotheses"];
  steps.forEach((v,i)=>{const x=72+(i%3)*392,y=175+Math.floor(i/3)*190;box(s,x,y,340,125,[C.ltblue,C.pale,C.yellow,C.peach,C.orange,C.red][i]);text(s,v,x+22,y+18,296,42,{fontSize:21,bold:true});text(s,["Who entered and why?","How are proteins organized?","Which chemotypes prefer which targets?","How much support is genuinely distinct?","Where do disease links sit in the body?","Which targets deserve follow-up?"][i],x+22,y+64,296,40,{fontSize:16,color:C.mid});});
  notes(s,"这一页不是目录，而是论证链。每张图回答上一张自然提出的问题，最后落到可检验的药物发现假设。",["FIGURE_CONTRACTS.tsv"]);

  s=add(); title(s,"Identity, granularity and lineage are fixed before counting","Figure 1 · Technical background");
  await addImage(s,path.join(GRAPHICS,"Section_Graphic_Evidence.png"),{left:72,top:165,width:470,height:410},"Evidence section graphic");
  bulletBlock(s,["Canonical protein, isoform, processed form and complex are separate entities.","Canonical compound and parent/form identity prevent salt and stereochemistry conflation.","MONDO exact mapping preserves source disease IDs and ontology hierarchy.","Database count, evidence record and putative experiment lineage are different units."],590,170,590,82);
  notes(s,"先解释四个前置概念，再进入整图。重点是主键、外键和来源谱系；记录多不等于实验多。",["MemPro V7.2 manifests","F1C_source_evidence_tier_flow.tsv","F1D_module_protein_coverage.tsv"]);

  s=add(); title(s,"MemPro is a release contract, not a flat merge","Figure 1 · Main result");
  await addImage(s,path.join(MAIN,"Figure1_resource_architecture.png"),{left:120,top:140,width:1040,height:520},"Figure 1 resource architecture");
  notes(s,"按A到D讲。特别强调646,670是注册小分子，240,055才有正式pair；942,455是证据记录，不是独立实验。",["Figure1_resource_architecture.pdf","F1C_source_evidence_tier_flow.tsv","F1D_module_protein_coverage.tsv"]);

  s=add(); title(s,"A/B/C describes membrane attachment; five axes describe function","Figure 2 · Technical background");
  await addImage(s,path.join(GRAPHICS,"Section_Graphic_Protein.png"),{left:72,top:160,width:500,height:430},"Protein classification section graphic");
  bulletBlock(s,["A: transmembrane or strong integral-membrane evidence.","B: directly embedded without spanning the bilayer, including lipid anchors.","C: stable peripheral membrane association without membrane insertion.","Five independent axes retain structural family, molecular function, biological process, membrane role and specialist classification."],610,165,555,83);
  notes(s,"A/B/C是膜连接机制，不是功能等级；E1–E3是膜身份的证据强度。五轴解决‘家族明确但角色未决’以及一个蛋白多功能的问题。",["protein_analysis.tsv.gz","F2E_classification_completeness.tsv"]);

  s=add(); title(s,"Membrane role strongly structures molecular function","Figure 2 · Main result");
  await addImage(s,path.join(MAIN,"Figure2_membrane_functional_architecture.png"),{left:165,top:132,width:950,height:535},"Figure 2 membrane functional architecture");
  notes(s,"英雄面板是B：标准化残差定位具体富集组合，整体Cramér's V=0.619。D展示Tau；E明确未分类与无专科分类，不隐藏知识空缺。",["Figure2_membrane_functional_architecture.pdf","F2B_role_function_cell_statistics.tsv","F2D_tissue_specificity_tau.tsv"]);

  s=add(); title(s,"Chemical diversity needs structure-aware and target-aware summaries","Figure 3 · Technical background");
  await addImage(s,path.join(GRAPHICS,"Section_Graphic_Compound.png"),{left:72,top:160,width:500,height:430},"Compound section graphic");
  bulletBlock(s,["Bemis–Murcko scaffolds summarize reusable chemical cores.","Morgan/Tanimoto measures local substructure similarity; UMAP is only a projection.","Target breadth and membrane-role entropy separate selectivity from polypharmacology.","Physicochemical distributions are descriptive and are not a Lipinski pass/fail filter."],610,165,555,83);
  notes(s,"解释为什么不用一个拥挤UMAP作为主结论。骨架、物化性质、靶点广度和target-set Jaccard回答不同问题。",["compound_features.tsv.gz","F3A_scaffold_rank_abundance.tsv","F3E_chemical_vs_target_similarity.tsv"]);

  s=add(); title(s,"Chemical space is broad, but role-level preference is modest","Figure 3 · Main result");
  await addImage(s,path.join(MAIN,"Figure3_chemical_target_landscape.png"),{left:160,top:132,width:960,height:535},"Figure 3 chemical target landscape");
  notes(s,"A是骨架长尾；C把靶点数与跨角色熵分开；D的总体V=0.046，必须讲成弱偏好。E检验相似化合物是否共享靶点。",["Figure3_chemical_target_landscape.pdf","F3D_structure_role_cell_statistics.tsv","F3E_chemical_vs_target_similarity.tsv"]);

  s=add(); title(s,"Contributing sources and independent experiments are not synonyms","Figure 4 · Technical background");
  await addImage(s,path.join(GRAPHICS,"Section_Graphic_Evidence.png"),{left:72,top:165,width:475,height:405},"Evidence provenance section graphic");
  bulletBlock(s,["UpSet reports pair intersections across more than three databases.","Jaccard compares pair-set overlap, not biological agreement.","Putative experiment keys use publication, assay, target, compound and measurement context.","Site tiers and residue mapping quantify structural usability separately from interaction validity."],590,170,590,82);
  notes(s,"强调字段命名：distinct database count是贡献数据库数，putative experiment count是谱系推断，未做原文人工审核前不能称独立实验真值。",["pair_analysis.tsv.gz","F4A_source_upset.tsv","F4C_database_vs_putative_experiment.tsv"]);

  s=add(); title(s,"Coverage and evidence independence form different dimensions","Figure 4 · Main result");
  await addImage(s,path.join(MAIN,"Figure4_evidence_architecture.png"),{left:160,top:132,width:960,height:535},"Figure 4 evidence architecture");
  notes(s,"A/B说明来源重叠有限且结构不同；C直接证明数据库数和推定实验数并非一回事；F诚实展示绝大多数位点膜侧未知。",["Figure4_evidence_architecture.pdf","F4B_source_jaccard.tsv","F4F_membrane_side.tsv"]);

  s=add(); title(s,"Disease identity and anatomy require ontology boundaries","Figure 5 · Technical background");
  await addImage(s,path.join(GRAPHICS,"Section_Graphic_Disease.png"),{left:72,top:155,width:500,height:440},"Disease ontology section graphic");
  bulletBlock(s,["Only official exact/equivalent cross-references merge disease identity.","A disease can map to multiple anatomical systems; fractional counts avoid double counting.","Therapeutic area and anatomy are different ontology dimensions.","Unmapped formal diseases remain visible in every denominator audit."],610,165,555,83);
  notes(s,"解释为什么不再用关键词强迫疾病进入一个器官。疾病大类和亚型保留层级，多标签映射更接近本体真实结构。",["disease_relation_analysis.tsv.gz","F5_disease_mapping_coverage.tsv"]);

  s=add(); title(s,"Ontology mapping makes disease context interpretable","Figure 5 · Main result");
  await addImage(s,path.join(MAIN,"Figure5_disease_anatomy_context.png"),{left:160,top:132,width:960,height:535},"Figure 5 disease anatomy context");
  notes(s,"A是fractional人体图；B是治疗领域×解剖富集；D显示very high并非空值；E只展示阈值后的社区，避免毛线团。",["Figure5_disease_anatomy_context.pdf","F5B_therapeutic_anatomy_cell_statistics.tsv","F5D_disease_evidence_levels.tsv"]);

  s=add(); title(s,"Cross-module discovery must keep every dimension inspectable","Figure 6 · Technical background");
  await addImage(s,path.join(GRAPHICS,"Section_Graphic_Integrated.png"),{left:72,top:155,width:500,height:440},"Integrated discovery section graphic");
  bulletBlock(s,["Disease support, chemical coverage, anatomy and structure are displayed separately.","Degree-preserving permutations control network-degree bias in anatomy concordance.","Candidates require nonzero chemistry and a coordinate-ready structural route.","The output is a follow-up list—not disease causality or validated binding."],610,165,555,83);
  notes(s,"这一页先解释为什么不用单一总分。透明筛选能让审稿人看到每个候选为何进入，以及缺少什么。",["F6A_integrated_target_landscape.tsv","F6B_concordance_permutation_null.tsv"]);

  s=add(); title(s,"Integrated evidence reveals testable underexplored targets","Figure 6 · Main result");
  await addImage(s,path.join(MAIN,"Figure6_integrated_drug_discovery_insights.png"),{left:160,top:132,width:960,height:535},"Figure 6 integrated drug discovery insights");
  notes(s,"观察疾病—正常组织RNA一致性0.487，高于保持度数的零模型0.327，经验P=0.005。D列出透明候选，但仍需文献、实验或docking验证。",["Figure6_integrated_drug_discovery_insights.pdf","F6D_hypothesis_generation_candidates.tsv"]);

  s=add(); title(s,"The analysis package is ready; the website use case remains the next manuscript figure","Close");
  bulletBlock(s,["Six main figures, one supplementary chemical atlas and seven concept graphics are reproducible from V7.2.","Every panel has source data, denominator rules, statistical methods and provenance.","Figure 7 should be added only after search, download and API functions are accessible.","Manuscript claims must preserve the boundaries: concordance ≠ causality; docking ≠ validation."],115,175,1030,88);
  box(s,115,540,1030,70,C.pale);text(s,"Next deliverable: website/use-case evidence and manuscript Results text aligned to Figures 1–6",145,550,970,48,{fontSize:22,bold:true,alignment:"center"});
  notes(s,"收束：纯数据统计和图形叙事已经完成；网站可用后再增加Figure 7。不要以Thank you结束，而是给出下一步稿件动作。",["FINAL_QA_REPORT.json","MAIN_FIGURE_CAPTIONS_EN.md"]);

  for (const [i,slide] of deck.slides.items.entries()) {
    const stem=`slide-${String(i+1).padStart(2,"0")}`;
    await writeBlob(path.join(RENDER,`${stem}.png`),await deck.export({slide,format:"png",scale:1}));
    const layout=await slide.export({format:"layout"}); await fs.writeFile(path.join(RENDER,`${stem}.layout.json`),await layout.text());
  }
  await writeBlob(path.join(OUT,"MemPro_V72_NAR_story_montage.webp"),await deck.export({format:"webp",montage:true,scale:1}));
  const pptx=await PresentationFile.exportPptx(deck); await pptx.save(path.join(OUT,"MemPro_V72_NAR_story_and_figures.pptx"));
  console.log(JSON.stringify({slides:deck.slides.items.length,pptx:path.join(OUT,"MemPro_V72_NAR_story_and_figures.pptx")}));
}
main().catch(e=>{console.error(e);process.exitCode=1});
