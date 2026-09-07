# MemPro V7.2 网页数据库条目、分类与层次规划

## 一、首页统计必须采用双口径

- 可检索 canonical 膜相关蛋白：**7,800**。
- 默认确认核心 E1/E2：**7,374**；E3候选：**426**。
- 有正式小分子互作的蛋白：**2,829**。
- canonical小分子总目录：**646,700**；进入正式pair的小分子：**240,055**。
- 正式蛋白—小分子pair：**529,168**。
- 正式正证据：**942,455**。
- canonical疾病：**3,739**；蛋白—疾病关系：**6,753**，涉及蛋白**2,347**。
- 蛋白复合物：**4,584**；组分断言：**14,825**；正式复合物结合证据：**24,954**。
- 结合位点：**55,610**。
- HPA表达测量：**3,046,789**；亚细胞定位：**19,091**。

首页必须把“总实体目录”和“有互作、疾病或位点证据的子集”分开，不能让用户误以为7,800个蛋白全部都有小分子互作。

## 二、数据库层次

### L0 版本与治理层
来源注册表、来源版本、release rules、manifest、SHA-256和validation。对应网页About、Methods、Downloads和Release history。

### L1 核心实体层
1. Gene：HGNC/Ensembl/NCBI Gene与canonical蛋白汇总。
2. Canonical protein：UniProt accession、名称、序列、A/B/C和E等级。
3. Protein isoform：isoform ID、序列哈希及形式特异分类。
4. Processed form：signal peptide、propeptide、chain等加工形式及坐标。
5. Protein complex：复合物、组分、assembly、stoichiometry及状态。
6. Canonical compound：canonical母体、结构、物化性质和身份键。
7. Compound form：盐型、多组分、立体异构体、质子化/电荷形式。
8. Canonical disease：MONDO锚点及保留的来源疾病实体。

### L2 分类与注释层
- 膜连接分类：A/B/C。
- 膜身份证据：E1/E2/E3。
- 五轴交叉分类：structural_family、molecular_function、biological_process、membrane_role、specialist_classification。
- 状态依赖膜结合。
- 组织RNA、cell-type RNA、IHC、质谱及亚细胞定位。
- 疾病治疗领域、父子层级及解剖/器官多标签。

### L3 关系层
Protein—Compound pair；Protein—Disease；Complex—Component；Complex—Compound；Disease hierarchy/xref；Disease—Anatomy/Therapeutic area。

### L4 证据与谱系层
每条正向证据、BE1/BE2/BE3、定量值、单位、assay、PDB、文献、来源记录ID、experiment/structure lineage key、位点S1/S2/S3，以及上下文负证据和冲突。

### L5 应用派生层
统计图、API、批量下载和Docking readiness。Docking readiness是结构计算条件，不等同于互作证据等级。

## 三、蛋白条目和分类

### ABC分类
- A：integral/intramembrane，**5,680**。
- B：脂质锚或其他直接嵌膜，**613**。
- C：稳定外周膜结合，**1,507**。

### 主膜连接方式（全部）
1. transmembrane_or_intramembrane
2. structure_supported_integral_membrane
3. beta_barrel_integral_membrane
4. explicit_uniprot_covalent_lipid_anchor
5. covalent_lipid_anchor
6. uniprot_peripheral_membrane_annotation
7. validated_direct_or_state_dependent_peripheral
8. structure_supported_peripheral
9. validated_complex_mediated_membrane_association
10. immunoglobulin_heavy_chain_membrane_isoform

### E等级
E1 5,808；E2 1,566；E3 426。蛋白浏览默认E1/E2，E3通过开关显示。

### 五轴交叉分类
- structural_family
- molecular_function
- biological_process
- membrane_role
- specialist_classification

全部实际标签、来源、PRIMARY/ADDITIONAL和是否推断见 `MemPro_V7.2_FIVE_AXIS_LABELS.tsv`。

## 四、小分子条目和分类

### canonical层
名称、standard SMILES/InChI/InChIKey、结构键、分子式、MW、XlogP、TPSA、氢键供受体、可旋转键、环、结构类别、来源ID、生物状态、靶点数及最佳BE。

### form层（全部）
parent_form；salt_or_multicomponent_form；stereochemically_defined_parent；charge_or_standardized_form；charged_or_protonation_form。

### 推荐网页过滤
身份置信度；结构类别；approved drug；clinical candidate；endogenous ligand；natural product；chemical probe；开发阶段；靶点数；最佳BE。

## 五、互作和证据条目

- pair层：529,168个唯一protein—canonical compound pair。
- evidence层：942,455条，允许一个pair拥有多条实验或数据库记录。
- BE1：结构直接证据。
- BE2：直接定量结合证据。
- BE3：功能、药理或生物活性证据。
- 主要activity：Ki、Kd、IC50、EC50、AC50、Km及其他来源原始类型。
- 数据库数量与独立实验、独立结构、独立模态分开显示。
- 无reported residue仍是有效互作。

## 六、结合位点

- S1 exact cocrystal + complete coordinates：442。
- S2 coordinate complete, ligand unresolved：48,655。
- S3 source-only or partial：6,513。
- residue状态：complete、partial、source-only、no residue。
- 链状态：unique、multiple equivalent、not reported。
- 膜侧：intramembrane、membrane interface/mixed、cytoplasmic、non-cytoplasmic、both/mixed、extramembrane unresolved、unknown。
- 页面必须分别显示source residue和mapped 3D coordinate。

## 七、疾病

- canonical疾病实体：3,739。
- 关系：6,753。
- 证据：very_high、high、medium。
- 22个Open Targets治疗领域；18个解剖系统标签；允许一个疾病属于多个系统。
- source ID、canonical ID、mapping predicate和本体层级同时保留。

## 八、表达与定位

- RNA与protein分层。
- measurement type：nCPM、nTPM、TPM、pTPM、IHC_level、MS_intensity。
- detection：detected、not_detected、missing_or_not_measured。
- mapped和unmapped单列。
- RNA、IHC、MS不得相加；missing不得解释为not detected。

## 九、复合物

- 复合物：4,584。
- heteromeric：4,275；homomer/single-component：309。
- 全部统一发布，不分formal/beta。
- 必须显示component set、membrane components、assembly、copy count、stoichiometry、状态和证据范围。
- component-context-only不能冒充complex-specific binding。

## 十、推荐网页栏目

1. Home
2. Browse proteins
3. Protein detail：Overview｜Forms｜Classification｜Interactions｜Sites/Structures｜Expression/Localization｜Diseases｜Complexes｜Provenance
4. Browse compounds
5. Compound detail
6. Browse complexes
7. Disease explorer
8. Evidence explorer
9. Site/Structure explorer
10. Expression atlas
11. Download/API
12. About/Methods/QC

## 十一、第一版MVP

优先完成Protein、Compound、Pair/Evidence、Disease、Site五类页面，以及全局搜索和下载。第二轮增加Expression atlas和Complex explorer。底层冻结表已经支持这些模块，不需要改写V7.2。

全局搜索应支持UniProt、gene symbol、HGNC、Ensembl、NCBI Gene、compound name、CID、ChEMBL ID、InChIKey、disease name、MONDO/OMIM/Orphanet/EFO、PDB和complex ID。
