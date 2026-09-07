# 02｜端到端工作流程

## A. 数据来源

蛋白：UniProtKB、HPA、HTP、Membranome、OPM、PDBTM。  
结构：PDBe、PDBbind、OPM、PDBTM、BioLiP、sc-PDB。  
互作：ChEMBL、BindingDB、PubChem BioAssay、BRENDA、IUPHAR/GtoPdb、PDSP KiDatabase、DrugCentral。  
疾病：Open Targets、UniProtKB疾病注释、MONDO、DO、Uberon。  
DrugBank因许可不可得，不纳入。

## B. 蛋白层

1. canonical UniProt accession作为主键。
2. 保存canonical序列和SHA-256。
3. 基因层使用HGNC、Ensembl Gene和NCBI Gene。
4. A：跨膜/膜内嵌；B：脂质锚或其他直接嵌膜；C：稳定外周膜结合。
5. E1/E2默认确认核心；E3候选。
6. 同一蛋白同时保留secondary membrane modes。
7. isoform、processed form和state独立建表。
8. 形式未指明的互作指向canonical，并标记ISOFORM_UNSPECIFIED。
9. 免疫球蛋白/TCR V/D/J片段不作为独立膜蛋白。
10. 分泌蛋白和腔内蛋白若无膜连接机制则排除。

## C. 五轴交叉分类

structural_family、molecular_function、biological_process、membrane_role、specialist_classification。  
来源包括GO、InterPro/Pfam、Reactome、GPCRdb、IUPHAR、TCDB、EC及MemPro规则。家族与功能不能互相替代。

## D. 小分子层

1. compound_internal_id为内部主键。
2. CID是PubChem编号；InChIKey用于跨来源结构身份。
3. canonical母体与form分层。
4. form包括盐、多组分、立体、质子化和电荷形式。
5. 名称不触发自动合并；结构唯一才合并。
6. standard SMILES、InChI、InChIKey、分子式和物化性质保留。
7. approved、clinical、endogenous、natural product、probe为独立生物状态标签。

## E. 互作证据

1. evidence记录是单条来源/实验记录；pair是去重后的protein—canonical compound关系。
2. BE1结构直接证据；BE2直接定量结合；BE3功能/药理/活性。
3. 保存原始值、单位、standard nM、activity type、assay、文献、PDB和来源ID。
4. v3.1_mixed_sources的74,932条全部恢复到具体来源。
5. PubMed+assay+target+compound+measurement构建实验谱系键。
6. PDB+chain+ligand构建结构谱系键。
7. AID+SID+CID构建PubChem谱系键。
8. 数据库数量、独立实验、独立结构和证据模态分开。
9. 阳性、阴性和冲突不互相覆盖。
10. 非冲突负证据仅归档；冲突负证据用于上下文提示。

## F. 结合位点

S1：exact cocrystal且坐标完整。  
S2：坐标完整但查询配体身份未完全解决。  
S3：来源残基或部分坐标。  
无residue仍是有效互作。来源残基与三维映射坐标分开。链不唯一时列候选链。膜侧可为intramembrane、interface、cytoplasmic、non-cytoplasmic、mixed或unknown。

## G. 疾病

仅使用官方exact/equivalent映射进行身份合并；broad/narrow/related只建关系。名称相同不自动合并。保留source disease ID。疾病允许多个治疗领域和多个器官系统。

## H. 表达与定位

RNA、IHC、MS分层；正常与疾病组织分层；mapped+measured才进入比例分母；missing/unmapped不等同于not detected。RNA不得称作蛋白表达。

## I. 复合物

独立protein_complex实体；组件可指canonical或isoform；保存assembly、copy count、stoichiometry、状态和required/optional字段。全部统一发布，但证据范围必须区分complex-specific与component-context-only。

## J. 冻结

从candidate生成正式release，执行主键、外键、来源、数量和规则QA；生成manifest和SHA-256；确认blocking_errors为空；设置只读。不得原地修改正式冻结目录。
