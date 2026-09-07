# MemPro V7.2 投稿图手工文字替换表

本文件只列出 SVG 中因 `textwrap.shorten()` 而显示为省略号的文字。完整名称来自对应 Figure Source Data 和绘图脚本；未修改任何图片。

## Figure 1｜Data provenance and release framework

### Panel b：Source → evidence modality → BE tier

按图中从上到下顺序：

| 当前显示 | 完整显示文字 | Source Data 原值 |
|---|---|---|
| IUPHAR/BPS Guide to… | IUPHAR/BPS Guide to PHARMACOLOGY | IUPHAR/BPS Guide to PHARMACOLOGY |
| Target-specific… | Target-specific activity | pubchem_target_specific_quantitative_activity |
| … | Functional/pharmacology | functional_or_bioactivity_measurement / binding_assay_activity |

### Panel c：Protein-level module coverage

| 当前显示 | 完整显示文字 | Source Data 原值 |
|---|---|---|
| expression mapped… | expression mapped and measured | expression_mapped_measured |

## Figure 2｜Five-axis organization of membrane proteins

### Panel b：Role × molecular-function enrichment

横轴从左到右：

| 当前显示 | 完整显示文字 | Source Data 原值 |
|---|---|---|
| binding… | binding activity | binding_activity |
| catalytic… | catalytic activity | catalytic_activity |
| ion channel… | ion channel activity | ion_channel_activity |
| molecular… | molecular function regulator | molecular_function_regulator |
| other… | other molecular function annotated | other_molecular_function_annotated |
| signaling… | signaling receptor activity | signaling_receptor_activity |
| structural… | structural molecule activity | structural_molecule_activity |
| transmembrane… | transmembrane transporter activity | transmembrane_transporter_activity |

纵轴中被省略的项目，按图中从上到下顺序：

| 当前显示 | 完整显示文字 | Source Data 原值 |
|---|---|---|
| membrane… | membrane-associated enzyme | membrane_associated_enzyme |
| signaling… | signaling regulator | signaling_regulator |
| family defined… | family-defined membrane role unresolved | family_defined_membrane_role_unresolved |
| membrane… | membrane scaffold or linker | membrane_scaffold_or_linker |
| junction or… | junction or adhesion | junction_or_adhesion |
| adhesion… | adhesion recognition | adhesion_recognition |
| membrane… | membrane trafficking | membrane_trafficking |
| immune or cell… | immune or cell recognition | immune_or_cell_recognition |
| membrane… | membrane organizer | membrane_organizer |

### Panel c：Primary-category flow across the classification system

每一列均按图中从上到下顺序。

#### Membrane role

| 当前显示 | 完整显示文字 |
|---|---|
| membrane… | membrane-associated enzyme |
| family defined… | family-defined membrane role unresolved |
| signaling… | signaling regulator |
| membrane scaffold… | membrane scaffold or linker |
| adhesion… | adhesion recognition |
| immune or cell… | immune or cell recognition |
| membrane… | membrane trafficking |

#### Molecular function

| 当前显示 | 完整显示文字 |
|---|---|
| signaling… | signaling receptor activity |
| transmembrane… | transmembrane transporter activity |
| molecular… | molecular function regulator |
| ion channel… | ion channel activity |
| other molecular… | other molecular function annotated |

#### Biological process

| 当前显示 | 完整显示文字 |
|---|---|
| transport and… | transport and localization |
| signal… | signal transduction |
| developmental… | developmental process |
| response to… | response to stimulus |
| protein metabolic… | protein metabolic process |
| lipid metabolic… | lipid metabolic process |
| immune system… | immune system process |
| other biological… | other biological process annotated |

### Panel e：Axis completeness

按图中从上到下顺序：

| 当前显示 | 完整显示文字 | Source Data 原值 |
|---|---|---|
| specialist… | specialist classification | specialist_classification |
| molecular… | molecular function | molecular_function |
| biological… | biological process | biological_process |
| structural… | structural family | structural_family |

## Figure 3｜Chemical diversity and membrane-target preference

### Panel d：Structure-class preference

横轴中被省略的项目，按图中从左到右顺序：

| 当前显示 | 完整显示文字 | Source Data 原值 |
|---|---|---|
| membrane… | membrane-associated enzyme | membrane_associated_enzyme |
| signaling… | signaling regulator | signaling_regulator |
| membrane… | membrane scaffold or linker | membrane_scaffold_or_linker |
| family defined… | family-defined membrane role unresolved | family_defined_membrane_role_unresolved |
| adhesion… | adhesion recognition | adhesion_recognition |
| immune or cell… | immune or cell recognition | immune_or_cell_recognition |

纵轴：

| 当前显示 | 完整显示文字 | Source Data 原值 |
|---|---|---|
| … | peptidic/peptidomimetic | peptidic/peptidomimetic |

## Figure 4｜Evidence provenance and site annotation completeness

### Panels a–c：Evidence-source labels

以下三个省略形式均应替换为同一个完整名称：

| 当前显示 | 完整显示文字 |
|---|---|
| IUPHAR/BPS Guide to… | IUPHAR/BPS Guide to PHARMACOLOGY |
| IUPHAR/BPS… | IUPHAR/BPS Guide to PHARMACOLOGY |
| IUPHAR/BPS Guide… | IUPHAR/BPS Guide to PHARMACOLOGY |

### Panel e：Residue-coordinate mapping by site tier

| 当前显示 | 建议完整显示文字 | Source Data 原值 |
|---|---|---|
| S1 EXACT COCRYSTAL… | S1 exact co-crystal coordinate complete | S1_EXACT_COCRYSTAL_COORDINATE_COMPLETE |
| S2 COORDINATE… | S2 coordinate complete; ligand unresolved | S2_COORDINATE_COMPLETE_LIGAND_UNRESOLVED |

### Panel f：Membrane-side completeness

| 当前显示 | 完整显示文字 | Source Data 原值 |
|---|---|---|
| extramembrane side… · 4,390 | extramembrane side unresolved · 4,390 | extramembrane_side_unresolved |
| membrane interface or… · 1,328 | membrane interface or mixed · 1,328 | membrane_interface_or_mixed |

## Figure 5｜Disease, therapeutic-area and anatomy context

### Panel b：Therapeutic area × anatomy enrichment

横轴从左到右（只列被省略者）：

| 当前显示 | 完整显示文字 |
|---|---|
| … | musculoskeletal system |
| entire sense… | entire sense organ system |
| cardiovascular… | cardiovascular system |
| circulatory… | circulatory system |
| digestive… | digestive system |
| alimentary… | alimentary part of gastrointestinal system |
| hematopoietic… | hematopoietic system |
| integumental… | integumental system |
| reproductive… | reproductive system |

纵轴从上到下（只列被省略者）：

| 当前显示 | 完整显示文字 |
|---|---|
| genetic, familial or… | genetic, familial or congenital disease |
| musculoskeletal or… | musculoskeletal or connective tissue disease |
| disorder of visual… | disorder of visual system |
| integumentary system… | integumentary system disorder |
| endocrine system… | endocrine system disorder |

### Panel e：Selected shared-target disease communities

| 当前显示 | 完整显示文字 | Canonical disease ID |
|---|---|---|
| catecholaminergic… | catecholaminergic polymorphic ventricular tachycardia 4 | MONDO:0013966 |

## 无需替换的图

- Figure 6：未检测到省略号文字。
- Supplementary Figure S1：未检测到省略号文字。
- Graphical Abstract：未检测到省略号文字。
- 六张章节概念 Graphic：未检测到省略号文字。

## 手工修改提示

- 请编辑 `04_figures/main/svg/` 中的 SVG，不要编辑 PNG。
- 若完整名称太长，优先手动换成两行或三行，不要再使用省略号。
- `membrane…` 在不同位置对应不同术语，必须依据本文件中的面板和排列顺序替换。
- 修改后另存新文件，保留原始投稿图和对应 Source Data 不变。
