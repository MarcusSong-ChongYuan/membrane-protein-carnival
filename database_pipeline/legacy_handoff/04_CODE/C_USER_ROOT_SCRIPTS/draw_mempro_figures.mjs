import fs from "node:fs/promises";
import sharp from "file:///C:/Users/Administrator/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp/lib/index.js";

const W=2400,H=1040;
const C={navy:"#12324A",teal:"#159A92",coral:"#E66B5B",purple:"#775DA6",blue:"#3D7EA6",amber:"#E2A33B",green:"#3B8C6E",gray:"#526575",light:"#E8EFF2",pale:"#F7FAFA",white:"#FFFFFF",red:"#C94747"};
const data=JSON.parse(await fs.readFile("data/plot_data.json","utf8"));
const out="figures"; await fs.mkdir(out,{recursive:true});

const esc=s=>String(s??"").replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");
const fmt=n=>Number(n).toLocaleString("en-US");
const pct=(n,d)=>`${(100*n/d).toFixed(1)}%`;
function text(x,y,s,size=30,fill=C.navy,weight=400,anchor="start",extra="") {return `<text x="${x}" y="${y}" font-family="Microsoft YaHei,Arial,sans-serif" font-size="${size}" fill="${fill}" font-weight="${weight}" text-anchor="${anchor}" ${extra}>${esc(s)}</text>`;}
function rect(x,y,w,h,fill="none",stroke="none",sw=1,rx=0){return `<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${rx}" fill="${fill}" stroke="${stroke}" stroke-width="${sw}"/>`;}
function line(x1,y1,x2,y2,stroke=C.light,sw=2,dash=""){return `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${stroke}" stroke-width="${sw}" ${dash?`stroke-dasharray="${dash}"`:""}/>`;}
function circle(cx,cy,r,fill,stroke="none",sw=1){return `<circle cx="${cx}" cy="${cy}" r="${r}" fill="${fill}" stroke="${stroke}" stroke-width="${sw}"/>`;}
function panelHead(x,y,label,title,sub=""){return text(x,y,label,38,C.teal,700)+text(x+48,y,title,38,C.navy,700)+(sub?text(x+48,y+38,sub,23,C.gray):"");}
function base(content,foot="MemPro V6.3.1 · data-driven redraw"){return `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">${rect(0,0,W,H,C.white)}${content}${line(40,1004,2360,1004,"#D7E2E5",2)}${text(40,1030,foot,20,C.gray)}</svg>`;}
function arrow(x1,y,x2,color=C.teal){return line(x1,y,x2-18,y,color,6)+`<polygon points="${x2-20},${y-11} ${x2},${y} ${x2-20},${y+11}" fill="${color}"/>`;}
function hbars(x,y,w,h,items,{color=C.teal,labelW=260,log=false,max=null,valueFmt=fmt,barH=34,gap=18}={}){
  let s=""; const vals=items.map(d=>+d[1]); const mx=max??Math.max(...vals); const plotW=w-labelW-100;
  items.forEach((d,i)=>{const yy=y+i*(barH+gap); const v=+d[1]; const frac=log?(Math.log10(v+1)/Math.log10(mx+1)):v/mx; s+=text(x+labelW-14,yy+barH*0.72,d[0],24,C.gray,400,"end"); s+=rect(x+labelW,yy,plotW,barH,"#EDF3F4","none",0,5); s+=rect(x+labelW,yy,Math.max(2,plotW*frac),barH,color,"none",0,5); s+=text(x+labelW+plotW*frac+10,yy+barH*0.72,valueFmt(v),23,C.gray,500);});
  return s;
}
function stackedH(x,y,w,h,rows,series,colors){
  let s=""; const max= Math.max(...rows.map(r=>series.reduce((a,k)=>a+(r[k]||0),0))); const labelW=100, plotW=w-labelW-80;
  rows.forEach((r,i)=>{const yy=y+i*94; const total=series.reduce((a,k)=>a+(r[k]||0),0); s+=text(x+labelW-14,yy+35,r.label,31,C.navy,700,"end"); let xx=x+labelW; for(let j=0;j<series.length;j++){const v=r[series[j]]||0; const ww=plotW*v/max; if(ww>0){s+=rect(xx,yy,ww,48,colors[j],C.white,2,3); if(ww>70)s+=text(xx+ww/2,yy+32,fmt(v),21,C.white,600,"middle");}xx+=ww;} s+=text(x+labelW+plotW+12,yy+32,fmt(total),23,C.gray,600);}); return s;
}
function heat(v){const z=Math.max(-2.5,Math.min(2.5,+v)); if(z<0){const t=(z+2.5)/2.5; return `rgb(${Math.round(61+(255-61)*t)},${Math.round(126+(255-126)*t)},${Math.round(166+(255-166)*t)})`;} const t=z/2.5; return `rgb(${Math.round(255+(230-255)*t)},${Math.round(255+(107-255)*t)},${Math.round(255+(91-255)*t)})`;}
async function save(name,svg){await fs.writeFile(`${out}/${name}.svg`,svg,"utf8"); await sharp(Buffer.from(svg)).png().toFile(`${out}/${name}.png`);}

// M1: exact architecture used by the prerequisite slide.
{
  let s=panelHead(50,62,"M1","数据库建立与发布链","每条记录保留来源、身份、证据和发布决策");
  const stages=[
    {x:55,w:400,t:"① 冻结原始来源",c:"#E8F1F7",items:["蛋白：UniProt / HPA / HTP","结构：PDBe / OPM / PDBTM","互作：ChEMBL / PubChem / BRENDA","疾病：Open Targets / UniProtKB"]},
    {x:505,w:385,t:"② 身份标准化",c:"#DDF2F0",items:["UniProt / gene / isoform","CID + full InChIKey","MONDO exact/equivalent","保留全部 source ID"]},
    {x:940,w:385,t:"③ Canonical实体",c:"#E5F3EC",items:["protein / isoform / complex","compound parent / exact form","canonical disease / hierarchy","名称相同不自动合并"]},
    {x:1375,w:385,t:"④ 证据与冲突",c:"#FFF3D8",items:["E1–E3 与 BE1–BE3","实验 / 结构 / 文献谱系键","阳性、阴性分别保存","冲突按实验条件解释"]},
    {x:1810,w:535,t:"⑤ 冻结发布",c:"#FBE8E4",items:["四张主表 + 扩展实体表","review queue / negative table","主键、外键、来源对账","manifest / hash / API / HPC"]},
  ];
  stages.forEach((g,i)=>{s+=rect(g.x,150,g.w,710,g.c,"#C9D9DE",3,18);s+=rect(g.x,150,g.w,72,i===4?C.coral:C.navy,"none",0,18);s+=text(g.x+g.w/2,198,g.t,29,C.white,700,"middle");g.items.forEach((v,j)=>{s+=rect(g.x+28,260+j*132,g.w-56,92,C.white,"#D4E1E4",2,12);s+=text(g.x+g.w/2,315+j*132,v,24,C.navy,600,"middle");}); if(i<stages.length-1)s+=arrow(g.x+g.w+8,505,stages[i+1].x-8);});
  s+=text(1200,928,"可追溯链：来源版本 → 原始记录 → 身份映射 → canonical实体 → 证据等级 → 发布/复核决策",29,C.teal,700,"middle");
  await save("M1_database_pipeline_redrawn",base(s,"Sources: frozen source registry; MemPro V6.3/V6.3.1 integration policy"));
}

// M2: ABC x evidence and five-axis coverage.
{
  const m=data.M2; let s=panelHead(45,60,"A","ABC × E1–E3/E0","膜结合方式与证据把握度是两条独立轴");
  const rows=["A","B","C","unknown"].map(k=>({label:k==="unknown"?"E0排除":k,...(m.class_evidence[k]||{})}));
  s+=stackedH(45,150,1040,420,rows,["E1","E2","E3","E0"],[C.teal,C.blue,C.purple,C.gray]);
  [[C.teal,"E1"],[C.blue,"E2"],[C.purple,"E3"],[C.gray,"E0"]].forEach((d,i)=>{s+=rect(160+i*180,565,26,26,d[0]);s+=text(195+i*180,587,d[1],24,C.gray);});
  s+=text(55,670,`总计 ${fmt(m.total)} 个蛋白`,34,C.navy,700); s+=text(55,716,"默认确认层：E1 + E2；E3保留为候选，E0为审计排除",25,C.gray);
  s+=panelHead(1190,60,"B","五轴交叉注释覆盖","家族、功能、过程、膜角色、专科层级并行保留");
  s+=hbars(1190,145,1160,400,Object.entries(m.axis_coverage),{color:C.teal,labelW:290,max:100,valueFmt:v=>`${v.toFixed(1)}%`,barH:34,gap:18});
  s+=panelHead(1190,620,"C","膜角色分布","“家族明确，角色未决”单列，不再伪装成other");
  const roleNames={membrane_associated_enzyme:"膜相关酶",family_defined_membrane_role_unresolved:"家族明确，角色未决",receptor:"受体",transporter:"转运体",membrane_scaffold_or_linker:"膜支架/连接",signaling_regulator:"信号调节",ion_channel:"离子通道",adhesion_recognition:"黏附/识别",immune_or_cell_recognition:"免疫/细胞识别"};
  s+=hbars(1190,700,1160,260,m.top_roles.slice(0,7).map(d=>[roleNames[d[0]]||d[0],d[1]]),{color:C.blue,labelW:300,barH:26,gap:10});
  await save("M2_protein_classification_redrawn",base(s,"Data: human_membrane_protein_master_v6_3_candidate.tsv.gz; n=10,997"));
}

// M3: mapping status, breadth, localization and z-score heatmap.
{
  const m=data.M3,total=Object.values(m.mapping).reduce((a,b)=>a+b,0); let s=panelHead(45,58,"A","HPA映射状态","把mapped zero与unmapped/missing分开");
  const mapItems=[["mapped detected",m.mapping.mapped_detected,C.teal],["mapped zero",m.mapping.mapped_zero,C.amber],["unmapped / missing",m.mapping.unmapped_missing,C.gray]];
  let xx=60; const sw=1010; mapItems.forEach(d=>{const ww=sw*d[1]/total;s+=rect(xx,145,ww,70,d[2],C.white,2); if(ww>120)s+=text(xx+ww/2,188,`${fmt(d[1])} · ${pct(d[1],total)}`,24,C.white,700,"middle");xx+=ww;}); mapItems.forEach((d,i)=>{s+=rect(80+i*320,250,24,24,d[2]);s+=text(114+i*320,271,d[0],23,C.gray);});
  s+=panelHead(45,370,"B","组织检出广度","箱体=四分位区间；竖线=中位数；仅mapped记录");
  const classes=["A","B","C","unknown"], labels={A:"A 跨膜",B:"B 嵌膜",C:"C 外周",unknown:"E0/未定"};
  classes.forEach((k,i)=>{const q=m.breadth[k], y=465+i*92, x0=230, pw=800, scale=v=>x0+pw*v/51;s+=text(200,y+22,labels[k],25,C.gray,600,"end");s+=line(x0,y+12,x0+pw,y+12,"#DDE7E9",12);s+=line(scale(q.min),y+12,scale(q.max),y+12,C.blue,4);s+=rect(scale(q.q1),y-5,Math.max(4,scale(q.q3)-scale(q.q1)),34,"#D9EEF0",C.teal,3,5);s+=line(scale(q.median),y-7,scale(q.median),y+31,C.coral,5);s+=text(1040,y+22,`median ${q.median.toFixed(0)} / 51`,22,C.gray);});
  s+=panelHead(1190,58,"C","亚细胞定位（前8）","一个蛋白可有多个定位术语");
  s+=hbars(1190,145,1160,360,m.top_localizations.slice(0,8),{color:C.blue,labelW:310,barH:28,gap:13});
  s+=panelHead(1190,560,"D","组织偏好热图","行z-score：同一功能类内部的相对组织偏好");
  const rows=m.zscore.filter(r=>["enzyme","receptor","transporter","ion_channel","adhesion_recognition","other_family_defined"].includes(r.functional_class)); const tissues=["liver","kidney","skeletal muscle","thyroid gland","skin","bone marrow","heart muscle","blood vessel"];
  const rn={enzyme:"酶",receptor:"受体",transporter:"转运体",ion_channel:"离子通道",adhesion_recognition:"黏附/识别",other_family_defined:"家族明确"}; const tn={liver:"肝",kidney:"肾","skeletal muscle":"骨骼肌","thyroid gland":"甲状腺",skin:"皮肤","bone marrow":"骨髓","heart muscle":"心肌","blood vessel":"血管"};
  const hx=1445,hy=650,cw=102,ch=48; tissues.forEach((t,j)=>s+=text(hx+j*cw+cw/2,632,tn[t],20,C.gray,500,"middle")); rows.forEach((r,i)=>{s+=text(hx-18,hy+i*ch+31,rn[r.functional_class],21,C.gray,500,"end");tissues.forEach((t,j)=>{const v=+r[t];s+=rect(hx+j*cw,hy+i*ch,cw-4,ch-4,heat(v),C.white,1,3);});});
  s+=text(2240,962,"蓝：相对低  ·  白：平均  ·  珊瑚：相对高",21,C.gray,400,"end");
  await save("M3_expression_localization_redrawn",base(s,"Data: HPA 25.1 expression/location module; n=10,997 proteins"));
}

// M4: exact-only disease identity and ontology-derived multi-label classifications.
{
  const m=data.M4; let s=panelHead(45,58,"A","疾病身份统一","仅direct MONDO或唯一official exact/equivalent自动进入canonical层");
  const canonNames={unique_official_exact_mapping:"官方唯一exact",direct_mondo_id:"直接MONDO ID",non_disease_mondo_entity_review:"非疾病实体复核",unmapped_source_only:"source-only",unique_mapping_to_obsolete_mondo_review:"obsolete复核"};
  s+=hbars(45,150,1100,320,Object.entries(m.canonicalization).sort((a,b)=>b[1]-a[1]).map(d=>[canonNames[d[0]]||d[0],d[1]]),{color:C.teal,labelW:330,barH:36,gap:18});
  s+=rect(70,520,310,88,"#DDF2F0","none",0,10)+text(225,555,fmt(m.source_diseases),34,C.teal,700,"middle")+text(225,588,"源疾病实体",22,C.gray,400,"middle");
  s+=rect(410,520,310,88,"#E8F1F7","none",0,10)+text(565,555,fmt(m.relations),34,C.blue,700,"middle")+text(565,588,"canonical蛋白–疾病关系",22,C.gray,400,"middle");
  s+=rect(750,520,330,88,"#FBE8E4","none",0,10)+text(915,555,"201",34,C.coral,700,"middle")+text(915,588,"需复核映射",22,C.gray,400,"middle");
  s+=panelHead(1195,58,"B","Open Targets治疗领域","官方多标签；同一疾病可进入多个领域");
  s+=hbars(1195,145,1150,410,m.therapeutic_areas.slice(0,8),{color:C.blue,labelW:430,barH:28,gap:15});
  s+=panelHead(45,690,"C","DO/MONDO → Uberon解剖映射","区分直接asserted与本体祖先推断");
  const ai=m.anatomy_assertion,at=Object.values(ai).reduce((a,b)=>a+b,0); const aitems=[["MONDO祖先推断",ai.inferred_from_mondo_ancestor,C.teal],["DOID祖先推断",ai.inferred_from_doid_ancestor,C.blue],["直接asserted",ai.asserted,C.coral]]; let ax=60; aitems.forEach(d=>{const ww=1030*d[1]/at;s+=rect(ax,780,ww,62,d[2],C.white,2);if(ww>120)s+=text(ax+ww/2,820,pct(d[1],at),23,C.white,700,"middle");ax+=ww;}); aitems.forEach((d,i)=>{s+=rect(65+i*340,875,22,22,d[2]);s+=text(98+i*340,895,`${d[0]} ${fmt(d[1])}`,21,C.gray);});
  s+=panelHead(1195,650,"D","疾病证据渠道","来自Open Targets聚合证据与UniProtKB人工疾病注释");
  s+=hbars(1195,740,1150,230,m.evidence_channels.slice(0,6),{color:C.purple,labelW:300,barH:24,gap:10});
  await save("M4_disease_ontology_redrawn",base(s,"Data: MONDO/DO/Uberon frozen mappings; Open Targets 26.06; UniProtKB 2026_02"));
}

// M5: true UpSet, Lipinski description, scaffold diversity and structural classes.
{
  const m=data.M5; let s=panelHead(40,58,"A","生物学状态交集（UpSet）","柱=交集数量；点阵=该交集包含的状态");
  const sets=["Approved","Clinical","Probe","Endogenous","Natural"], combos=m.status_intersections.slice(0,9); const ux=250,uy=190,barTop=105,barH=250,colW=90,max=Math.max(...combos.map(d=>d.count));
  combos.forEach((d,j)=>{const x=ux+j*colW,h=barH*Math.log10(d.count+1)/Math.log10(max+1);s+=rect(x,barTop+barH-h,54,h,C.blue,"none",0,4);s+=text(x+27,barTop+barH-h-10,fmt(d.count),20,C.gray,500,"middle");});
  sets.forEach((set,i)=>{const y=uy+250+i*54;s+=text(215,y+8,set,23,C.gray,500,"end"); combos.forEach((d,j)=>{const on=d.combination.split("+").includes(set),x=ux+j*colW+27;s+=circle(x,y,9,on?C.navy:"#DCE5E8");});});
  combos.forEach((d,j)=>{const ys=sets.map((set,i)=>d.combination.split("+").includes(set)?uy+250+i*54:null).filter(v=>v!==null);if(ys.length>1)s+=line(ux+j*colW+27,Math.min(...ys),ux+j*colW+27,Math.max(...ys),C.navy,4);});
  s+=panelHead(1190,58,"B","Lipinski性质满足度","描述性统计，不作为收录或Docking硬门槛");
  const lip=[ ["HBD ≤ 5",98.152],["HBA ≤ 10",97.667],["≤1项违反",85.293],["XlogP ≤ 5",77.985],["MW ≤ 500 Da",70.945] ];
  s+=hbars(1190,150,1160,330,lip,{color:C.teal,labelW:250,max:100,valueFmt:v=>`${v.toFixed(1)}%`,barH:34,gap:18});
  s+=panelHead(40,745,"C","Bemis–Murcko骨架多样性","替代无区分度的Morgan降维散点");
  const stats=[["核心QC化合物",m.core_qc,C.navy],["独特骨架",m.scaffold.unique_murcko_scaffolds,C.teal],["singleton骨架",m.scaffold.singleton_scaffolds,C.purple],["骨架多样性比",m.scaffold.scaffold_diversity_ratio,C.coral]];
  stats.forEach((d,i)=>{const x=55+i*270;s+=text(x,830,d[0],21,C.gray);s+=text(x,884,i===3?d[1].toFixed(3):fmt(d[1]),36,d[2],700);}); s+=text(55,940,`Top 10骨架仅覆盖 ${(100*m.scaffold.top10_scaffold_compound_share).toFixed(1)}% 化合物；结构空间高度分散`,23,C.gray);
  s+=panelHead(1190,585,"D","结构类别与骨架比","条长=化合物数；末端数字=独特骨架/化合物");
  const clsNames={polycyclic_organic:"多环有机",polycyclic_organic_compound:"多环有机化合物",cyclic_organic_compound:"环状有机化合物",cyclic_organic:"环状有机",macrocyclic_compound:"大环",carbohydrate_like:"糖类样",lipid_like:"脂质样"};
  s+=hbars(1190,680,1160,280,m.scaffold_by_class.slice(0,7).map(d=>[`${clsNames[d.class]||d.class} · ${d.ratio.toFixed(3)}`,d.compounds]),{color:C.blue,labelW:380,log:true,barH:25,gap:11});
  await save("M5_compound_upset_scaffold_redrawn",base(s,"Data: 899,591 core-QC compounds; RDKit 2026.03.5; full-data UpSet and Bemis–Murcko analysis"));
}

// M6: source-by-tier, lineage and binding-site readiness.
{
  const m=data.M6; let s=panelHead(45,58,"A","来源 × BE等级","每个来源按BE1/BE2/BE3分列；横轴为log10记录数");
  const sources=Object.entries(m.source_tier).slice(0,8); const sx=360,sy=145,pw=760,colors={BE1:C.coral,BE2:C.teal,BE3:C.blue};
  sources.forEach(([src,counts],i)=>{const y=sy+i*72;s+=text(330,y+32,src,22,C.gray,500,"end");["BE1","BE2","BE3"].forEach((tier,j)=>{const v=counts[tier]||0, ww=pw*Math.log10(v+1)/6.4;s+=rect(sx,y+j*17,ww,13,colors[tier],"none",0,2);});s+=text(1140,y+31,fmt(Object.values(counts).reduce((a,b)=>a+b,0)),20,C.gray,500);});
  [[C.coral,"BE1结构/位点"],[C.teal,"BE2定量直接结合"],[C.blue,"BE3功能/酶配体"]].forEach((d,i)=>{s+=rect(400+i*245,735,22,22,d[0]);s+=text(432+i*245,754,d[1],20,C.gray);});
  s+=panelHead(1260,58,"B","来源谱系覆盖","数据库数不等于独立实验数");
  const lin=[["PubChem assay键",m.lineage.pubchem_assay_lineage],["多证据模态pair",m.lineage.multi_modality_pairs],["文献实验代理键",m.lineage.literature_experiment_proxy],["独立结构键",m.lineage.structure_lineage]];
  s+=hbars(1260,150,1080,340,lin,{color:C.purple,labelW:280,barH:36,gap:22});
  s+=panelHead(1260,590,"C","结合位点可用性","仅统计有PDB/链/配体/残基或明确位点来源的记录");
  s+=text(1290,680,fmt(m.site_total),46,C.navy,700)+text(1510,680,"位点实例",25,C.gray);
  s+=text(1290,735,fmt(m.site_specificity.compound_specific),38,C.teal,700)+text(1510,735,"compound-specific",23,C.gray);
  const stNames={experimental_structure_residue_contact:"结构残基接触",experimental_structure_match:"结构匹配",experimental_complex_pocket:"复合物口袋",bindingdb_ligand_target_complex:"BindingDB复合物",uniprot_curated_binding_site:"UniProt人工位点"};
  s+=hbars(1260,780,1080,170,m.site_type.slice(0,5).map(d=>[stNames[d[0]]||d[0],d[1]]),{color:C.teal,labelW:300,barH:22,gap:8});
  await save("M6_binding_evidence_redrawn",base(s,"Data: binding_evidence_master_v6_2 + V6.3 lineage audit + 95,598 binding-site instances"));
}

// M7: negative evidence, conflicts, review queue and QA.
{
  const m=data.M7,n=m.negative_table; let s=panelHead(45,58,"A","负证据身份解析漏斗","负证据不会被阳性覆盖；未映射记录进入复核层");
  const funnel=[["原始负证据",n.source_rows,C.navy],["成功映射正式层",n.mapped_release_rows,C.teal],["未映射复核层",n.unmapped_review_rows,C.amber],["正负冲突标记",n.positive_negative_conflict_rows,C.coral]];
  const max=n.source_rows; funnel.forEach((d,i)=>{const y=150+i*115,ww=960*d[1]/max;s+=text(240,y+38,d[0],26,C.gray,600,"end");s+=rect(270,y,960,62,"#EDF3F4","none",0,6);s+=rect(270,y,Math.max(5,ww),62,d[2],"none",0,6);s+=text(1245,y+40,`${fmt(d[1])} · ${pct(d[1],max)}`,23,C.gray,600);});
  s+=text(275,655,`精确去重移除 ${fmt(n.duplicates_removed)} 条；正负冲突率=${pct(n.positive_negative_conflict_rows,n.mapped_release_rows)}（相对正式负证据）`,23,C.gray);
  s+=panelHead(1320,58,"B","未映射复核原因","不确定立体化学和盐/母体关系是主要来源");
  const rrNames={unassigned_stereochemistry:"未指定立体化学",salt_parent_not_already_curated_as_exact_form:"盐型/母体未唯一",record_has_no_parseable_CID:"无可解析CID",charge_or_protonation_parent_uncertain:"电荷/质子化不确定","canonicalization / kekulization":"结构规范化失败","full InChIKey mismatch":"full InChIKey不一致"};
  s+=hbars(1320,150,1030,340,m.review_reasons.map(d=>[rrNames[d[0]]||d[0],d[1]]),{color:C.amber,labelW:330,log:true,barH:32,gap:18});
  s+=panelHead(1320,590,"C","发布质量门","7项规则全部通过才冻结release");
  const qLabels={full_inchikey_required_for_exact_form_merge:"exact form必须full InChIKey",unassigned_stereochemistry_not_auto_merged:"未定立体化学不自动合并",uncurated_salt_or_charge_parent_not_auto_merged:"盐型/电荷不强并母体",dedup_key_target_canonical_aid_sid:"按target+canonical+AID+SID去重",positive_negative_conflicts_flagged:"正负冲突显式标记",unmapped_records_retained:"未映射记录保留",v61_immutable:"V6.1冻结层只读"};
  Object.entries(m.qa).forEach(([k,v],i)=>{const x=1340+(i%2)*490,y=680+Math.floor(i/2)*74;s+=circle(x,y,16,v?C.teal:C.red);s+=text(x,y+7,v?"✓":"×",22,C.white,700,"middle");s+=text(x+32,y+8,qLabels[k]||k,21,C.gray,500);});
  await save("M7_negative_conflict_qa_redrawn",base(s,"Data: NEGATIVE_EVIDENCE_RESOLUTION_V62_VALIDATION.json; V6.3.1 review queue summary"));
}

// M8: the validated 6,019-pair pilot, not the 192,592-pair full candidate table.
{
  const m=data.M8; let s=panelHead(45,58,"A","Pilot分层","6,019个pair / 1,773个靶点");
  const tierOrder=["R0","D1","D2","D3"],tc={R0:C.purple,D1:C.teal,D2:C.blue,D3:C.amber};
  tierOrder.forEach((k,i)=>{const x=85+i*265,y=155,h=390*m.tiers[k]/Math.max(...Object.values(m.tiers));s+=rect(x,540-h,150,h,tc[k],"none",0,8);s+=text(x+75,575,k,30,C.navy,700,"middle");s+=text(x+75,530-h,fmt(m.tiers[k]),25,C.gray,600,"middle");});
  s+=text(65,650,"R0：回对接验证；D1→D3：由强证据/高准备度到探索层",24,C.gray);
  s+=panelHead(1190,58,"B","受体结构等级","S1为pair-specific共晶位点");
  const sn={S1_pair_specific_experimental_site:"S1 pair-specific位点",S2_membrane_experimental_OPM_PDBTM:"S2 膜实验结构",S3_other_experimental_PDB:"S3 其他实验PDB"};
  s+=hbars(1190,150,1150,250,m.structure_grade.map(d=>[sn[d[0]]||d[0],d[1]]),{color:C.teal,labelW:360,barH:46,gap:28});
  s+=panelHead(45,760,"C","生产前四道门","名单完成≠Docking已执行");
  const gates=[
    ["受体",m.receptor_prep[1]?.[1]||1473,"S1需复核链/组装体",C.coral],
    ["配体",m.pairs,"全部需3D微状态准备",C.blue],
    ["Box",m.box_status[1]?.[1]||1473,"R0由共晶配体定义",C.purple],
    ["执行",m.execution.not_run,"当前全部not_run",C.gray],
  ];
  gates.forEach((d,i)=>{const x=70+i*555;s+=rect(x,830,500,105,"#F4F8F8","#D6E2E5",2,12);s+=rect(x,830,12,105,d[3],"none",0,6);s+=text(x+38,870,d[0],28,C.navy,700);s+=text(x+150,870,fmt(d[1]),27,d[3],700);s+=text(x+38,912,d[2],22,C.gray);});
  s+=text(1195,600,"HPC执行顺序",27,C.navy,700); const flow=["R0 redocking","RMSD ≤ 2 Å","D1生产","D2扩展","D3探索"]; flow.forEach((v,i)=>{const x=1195+i*225;s+=rect(x,650,185,62,i===0?"#EDE6F5":"#E2F1F0",i===0?C.purple:C.teal,2,10);s+=text(x+92,689,v,21,C.navy,600,"middle");if(i<flow.length-1)s+=arrow(x+190,681,x+220,C.teal);});
  await save("M8_docking_pilot_redrawn",base(s,"Data: docking_pilot_manifest_v631.tsv; exactly 6,019 pairs; execution status=not_run"));
}

console.log("generated",(await fs.readdir(out)).filter(f=>f.endsWith(".png")));
