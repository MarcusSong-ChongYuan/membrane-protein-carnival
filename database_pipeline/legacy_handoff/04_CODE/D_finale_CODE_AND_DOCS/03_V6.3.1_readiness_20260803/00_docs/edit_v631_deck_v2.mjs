import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const p=await PresentationFile.importPptx(await FileBlob.load("starter_v631.pptx"));
const q=await p.inspect({kind:"slide,textbox,shape,image,notes",maxChars:500000});
const rows=q.ndjson.trim().split(/\n/).map(s=>JSON.parse(s));
const find=(slide,name,kind)=>{
  const a=rows.filter(r=>r.slide===slide&&(!name||r.name===name)&&(!kind||r.kind===kind));
  if(a.length!==1) throw new Error(`target mismatch slide=${slide} name=${name} kind=${kind} n=${a.length}`);
  return p.resolve(a[0].id);
};
const tx=(s,n,t)=>{find(s,n).text=t};
const notes=(s,t)=>find(s,null,"notes").setText(t);
async function image(s,path,alt){
  const im=find(s,null,"image"), frame=im.frame, geometry=im.geometry, radius=im.borderRadius, rotation=im.rotation;
  im.replace({blob:await fs.readFile(path),contentType:"image/png",alt,fit:"contain"});
  im.frame=frame; im.crop={left:0,top:0,right:0,bottom:0}; im.fit="contain"; im.geometry=geometry; im.borderRadius=radius; im.rotation=rotation;
}

tx(1,"deck-title","MemPro V6.3.1");
tx(1,"deck-meta","V6.3候选冻结：2026-08-03  ·  V6.3.1就绪审计：2026-08-03");
tx(1,"deck-tagline","从可追溯数据库到可复现分析与HPC验证");
notes(1,"本稿以冻结的V6.3候选为数据基线；V6.3.1只增加工程就绪审计、复核队列与HPC准备。\n\n[Sources]\n- Internal: D:\\finale\\02_V6.3_candidate_20260803\\FREEZE_V6_3_CANDIDATE.json\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\05_qa\\V631_READINESS_VALIDATION.json");
for(let s=2;s<=18;s++) find(s,"footer-left").text.replace("MemPro V6.2","MemPro V6.3.1");

tx(5,"slide-title","V6.3整合14类主要来源，并区分可公开与受限衍生数据");
notes(5,"来源覆盖蛋白、结构、互作与疾病四层；每条记录保留来源版本和ID。\n\n[Sources]\n- UniProt: https://rest.uniprot.org/\n- HPA: https://www.proteinatlas.org/about/download\n- PDBe: https://www.ebi.ac.uk/pdbe/\n- PubChem: https://pubchem.ncbi.nlm.nih.gov/docs/bioassay\n- Open Targets: https://platform.opentargets.org/\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\01_engineering\\source_license_publication_gate.tsv");

tx(7,"section-kicker","05 · V6.3 / V6.3.1成果");
tx(7,"slide-title","V6.3候选已冻结；V6.3.1把剩余问题转为可执行队列");
tx(7,"metric-number-490-404","9,096");
tx(7,"metric-detail-490-404","其中 8,822 条进入默认疾病关系层");
tx(7,"pass-small","V6.3 candidate QA\n自动检查 = PASS");
tx(7,"overview-bottom-text","冻结数据不被改写；V6.3.1新增队列、环境锁和HPC门控，不等于人工复核已经完成。");
notes(7,"核心数量来自V6.3候选冻结包；人工文献、许可、第二设备复现与真实HPC仍开放。\n\n[Sources]\n- Internal: D:\\finale\\02_V6.3_candidate_20260803\\FREEZE_V6_3_CANDIDATE.json\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\05_qa\\V631_READINESS_VALIDATION.json");

tx(8,"slide-title","五轴交叉分类同时保留家族、功能、过程、膜角色与专科层级");
await image(8,"D:/finale/02_V6.3_candidate_20260803/08_figures/protein_annotation_v63_working/png/M2_five_axis_protein_classification_v63.png","V6.3五轴蛋白分类");
notes(8,"五轴分类避免用家族替代功能；自动注释覆盖全部蛋白，高价值未决项进入复核。\n\n[Sources]\n- Asset: D:\\finale\\02_V6.3_candidate_20260803\\08_figures\\protein_annotation_v63_working\\png\\M2_five_axis_protein_classification_v63.png\n- GO: https://current.geneontology.org/ontology/go-basic.obo\n- Reactome: https://reactome.org/download-data\n- InterPro: https://www.ebi.ac.uk/interpro/");

tx(9,"slide-title","表达与定位区分HPA已映射、测得为零与未映射/缺失");
notes(9,"0值不能一律解释为不表达；正式统计必须区分mapped zero和unmapped/missing。现有496条表达/定位复核项。\n\n[Sources]\n- HPA: https://www.proteinatlas.org/about/download\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\02_review_queues\\expression_location_review_summary_v631.tsv");

tx(10,"slide-title","疾病关系采用MONDO canonical、来源ID与多标签系统映射");
await image(10,"D:/finale/02_V6.3_candidate_20260803/08_figures/disease_v63_working/png/M4_disease_ontology_context_v63.png","V6.3疾病本体与多标签系统映射");
notes(10,"仅官方exact/equivalent映射自动合并；亚型保持层级，原始来源ID永久保留。仍有201个来源疾病ID待权威复核。\n\n[Sources]\n- Asset: D:\\finale\\02_V6.3_candidate_20260803\\08_figures\\disease_v63_working\\png\\M4_disease_ontology_context_v63.png\n- MONDO: https://mondo.monarchinitiative.org/pages/download/\n- Uberon: https://obofoundry.org/ontology/uberon.html");

tx(11,"slide-title","骨架多样性比低维化学空间更直接衡量结构覆盖");
await image(11,"D:/finale/02_V6.3_candidate_20260803/08_figures/disease_v63_working/png/M5E_scaffold_diversity.png","V6.3小分子骨架多样性");
notes(11,"本页用Bemis–Murcko骨架、环系和单例/重复骨架描述结构覆盖，替代拥挤的低维指纹嵌入。\n\n[Sources]\n- Asset: D:\\finale\\02_V6.3_candidate_20260803\\08_figures\\disease_v63_working\\png\\M5E_scaffold_diversity.png\n- Internal analysis: 899,591 core-QC compounds");

tx(12,"slide-title","数据库数、实验数、结构数与证据模态分别计数");
await image(12,"D:/finale/02_V6.3_candidate_20260803/08_figures/protein_annotation_v63_working/png/M7_evidence_lineage_v63.png","V6.3证据来源谱系");
notes(12,"数据库来源数不再被称为独立实验数；实验、结构、PubChem assay与证据模态分别建立谱系键。\n\n[Sources]\n- Asset: D:\\finale\\02_V6.3_candidate_20260803\\08_figures\\protein_annotation_v63_working\\png\\M7_evidence_lineage_v63.png\n- Internal: D:\\finale\\02_V6.3_candidate_20260803\\06_evidence_lineage\\binding_evidence_lineage_v0_1.parquet");

notes(13,"阳性和阴性可能因浓度、体系、构象与终点不同而同时成立。2,615,046条未映射负证据已形成原因汇总和分层样本。\n\n[Sources]\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\02_review_queues\\negative_unmapped_review_summary_v631.tsv\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\02_review_queues\\negative_unmapped_review_stratified_sample_v631.tsv");

tx(14,"slide-title","候选冻结与就绪审计均通过；人工和外部门槛保持开放");
tx(14,"qa-lead","PASS表示工程规则、主外键和审计产物可重放；不表示文献直证、授权或Docking已经完成。");
tx(14,"col-body-12","复合物直证、历史isoform、许可、第二设备复现和HPC执行仍有明确门槛。");
tx(14,"qa-gate-text","V6.3 candidate PASS  +  V6.3.1 automation PASS  →  Manual gates open");
notes(14,"V6.3.1状态为AUTOMATION_COMPLETE_MANUAL_GATES_OPEN，防止把待审候选误标为已确认。\n\n[Sources]\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\05_qa\\V631_READINESS_VALIDATION.json");

tx(16,"slide-title","当前不足已经量化，并拆成最小人工与外部执行任务");
tx(16,"lim-head-0-0","功能角色仍有未决项"); tx(16,"lim-body-0-0","E1/E2优先集仍有1,719个角色未决；已输出500个高价值P1项。");
tx(16,"lim-head-0-1","疾病映射存在边界项"); tx(16,"lim-body-0-1","201个来源疾病ID涉及1:N、obsolete、非疾病实体或source-only状态。");
tx(16,"lim-head-0-2","isoform历史证据不可强推"); tx(16,"lim-body-0-2","8个旧/删除accession及未指定isoform证据需权威记录或原文。");
tx(16,"lim-head-1-0","复合物直证仍需原文"); tx(16,"lim-body-1-0","24,954个候选中186个P1；完整复合物直证未核实前提升数为0。");
tx(16,"lim-head-1-1","负证据身份仍有复核层"); tx(16,"lim-body-1-1","2,615,046条未映射负证据已分层；仅结构身份唯一者可正式纳入。");
tx(16,"lim-head-1-2","公开与计算仍有外部门"); tx(16,"lim-body-1-2","许可、第二设备复现、受体/配体/box准备和真实HPC尚未完成。");
notes(16,"自动化只负责排序、抽样、对账与执行清单；科学判断和法律判断仍需相应证据。\n\n[Sources]\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\05_qa\\V631_READINESS_VALIDATION.json");

tx(17,"section-kicker","15 · 已完成与下一步");
tx(17,"slide-title","自动任务已落盘；后续按证据责任而不是文件数量推进");
tx(17,"road-head-0","已完成：工程与分层"); tx(17,"road-body-0","• 环境依赖锁\n• 路径与来源审计\n• 六类复核队列\n• 稳定ID与哈希清单");
tx(17,"road-head-1","人工：科学与许可"); tx(17,"road-body-1","• 复合物直证原文\n• 高价值功能未决项\n• 疾病/isoform边界\n• 再发布许可确认");
tx(17,"road-head-2","执行：复现与验证"); tx(17,"road-body-2","• 第二台设备一键复现\n• 受体/配体/box准备\n• redocking验收\n• HPC pilot后放行standard");
tx(17,"road-bottom-text","V6.3候选保持只读；V6.3.1是就绪增量，不把开放门槛写成完成。");
notes(17,"交接时先执行最小P1队列，并将决定写回review status、reviewer和evidence reference。\n\n[Sources]\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\README.md\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\01_engineering\\requirements-scientific-lock.txt");

tx(18,"slide-title","6,019对Pilot已整理为HPC数组，但尚未开始实际Docking");
tx(18,"dock-head","HPC三道门"); tx(18,"metric-detail-918-220","稳定ID、数组索引和批次映射已经生成");
tx(18,"metric-number-918-375","1,773"); tx(18,"metric-label-918-375","Unique targets"); tx(18,"metric-detail-918-375","逐靶点完成受体、box与参考配体验证");
tx(18,"dock-caution","Production gate\n受体 + 配体 + box + redocking RMSD ≤ 2 Å");
await image(18,"D:/finale/02_V6.3_candidate_20260803/08_figures/protein_annotation_v63_working/png/M8_docking_lineage_rerank_v63.png","V6.3 Docking谱系重排序");
notes(18,"Pilot成员保持6,019对和1,773个靶点不变。逐靶点完成结构、微状态、box和参考配体redocking后才能进入production。\n\n[Sources]\n- Asset: D:\\finale\\02_V6.3_candidate_20260803\\08_figures\\protein_annotation_v63_working\\png\\M8_docking_lineage_rerank_v63.png\n- Internal: D:\\finale\\03_V6.3.1_readiness_20260803\\04_docking_hpc\\docking_pilot_manifest_v631.tsv");

await fs.mkdir("final-previews-v2",{recursive:true}); await fs.mkdir("final-layouts-v2",{recursive:true});
for(let i=0;i<p.slides.items.length;i++){
  const slide=p.slides.items[i];
  const png=await p.export({slide,format:"png",scale:1});
  await fs.writeFile(`final-previews-v2/slide-${String(i+1).padStart(2,"0")}.png`,new Uint8Array(await png.arrayBuffer()));
  const layout=await slide.export({format:"layout"});
  await fs.writeFile(`final-layouts-v2/slide-${String(i+1).padStart(2,"0")}.layout.json`,await layout.text());
}
const montage=await p.export({format:"webp",montage:true,scale:0.5});
await fs.writeFile("final-montage-v2.webp",new Uint8Array(await montage.arrayBuffer()));
const inspect=await p.inspect({kind:"slide,textbox,image,notes",maxChars:500000});
await fs.writeFile("final-inspect-v2.json",JSON.stringify(inspect,null,2));
const pptx=await PresentationFile.exportPptx(p); await pptx.save("MemPro_V6.3.1_数据库建立流程与科学意义_更新版.pptx");
console.log(JSON.stringify({slides:p.slides.items.length,status:"exported"}));
