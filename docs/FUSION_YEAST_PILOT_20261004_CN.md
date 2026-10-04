# 融合试点：酵母 GEM 注册 + kcat 约束 FBA（GECKO 式）首次闭环

日期：2026-10-04。性质：工程试点（engineering pilot），不是验证实验，不产生发表级结论。对应讨论见 `docs/10H2DA_YEASTMETWIN_COMPARISON_20261004_CN.md` 第 5 节的融合差距矩阵；本试点落实其中的"双向耦合"两项。

## 1. 完成了什么

1. **第三个 FBA 上下文模型注册**：用户 "Metabolic model prediction" 项目的酵母基因组尺度模型（Yeast-MetaTwin，Yeast-GEM 风格 r_ 标识符）注册进 `src/biocandidate/fba_context.py`，模型文件以 SHA256 钉死（`configs/fba_context_model_yeast_metatwin.json`，只读来源，重哈希副本在 gitignored `data/`）。注册表从 E. coli 硬编码泛化为每模型声明式 spec（交换反应、特征反应、预设、默认培养基）。
2. **GECKO 式 kcat 约束 FBA 试点**（`scripts/fusion_yeast_kcat_fba_pilot.py`）：把 BioCandidateRanker 预测的 kcat 以酶容量上界（ub = kcat × E × 3600）注入酵母 GEM 的 10H2DA 候选路线。
3. **上下文附加演示**（`scripts/fusion_yeast_context_attach_demo.py`）：酵母 FBA 上下文向量经 workbench 预测路径附加到 BCR 推理。

## 2. 注册细节

- 模型身份：SHA256 `7b2576232bd093001c2f173e07fb9d29c0ae565b5fc3d4be0d88a4b7d444deec`，5,435,712 字节，加载时强制校验。
- 8 维特征契约（与所有已冻结检查点的 `fba_context_dim=8` 一致）：生物量目标 `r_2111`、葡萄糖交换 `r_1714`、氧交换 `r_1992`、糖酵解 PFK `r_0886`、PPP 入口 G6PDH `r_0466`、TCA 柠檬酸合酶 `r_0300`、ATP 维持 `r_4046`、POX1 癸酰-CoA 氧化 `r_0120`（10H2DA 前体链的关口）。酵母无 ED 途径与膜转氢酶，E. coli 特征槽中的 EDD/THD2pp 由氧交换与 POX1 替代（映射在代码注释与 config 中声明）。
- 默认培养基：葡萄糖 10 mmol/gDW/h、氧 1000（模型文件自身上限）；`run_fba_context` 的默认值语义从模块级 E. coli 常量改为按注册模型解析，显式传参行为不变。
- 基线解：生长 0.9241 h⁻¹（status optimal），condition_id 绑定模型哈希+场景。
- 测试：`tests/test_fba_context.py` 注册表契约扩为三模型并新增酵母基线/确定性测试；`test_fba_context` 20 项、`test_fba_surrogate`/`test_cli_predict_fba`/`test_workbench` 36 项全部通过。

## 3. 试点一：kcat 约束 FBA

外部 10H2DA 候选层（4 条 cand_* 反应 + demand）从目标工作流配置（SHA256 已记录）精确重建。

**锚点复现**：无 kcat 约束时最大 10H2DA demand 通量 = **0.274011**，与 Yeast-MetaTwin 9 月 2 日报告一致到小数点后 6 位。复现链（模型副本→候选反应→生物量 10% 下限→demand 最大化）完整有效。

**kcat 来源**：上一轮 18 候选预测（`artifacts/10h2da-yeast-prediction/task1_10h2da_pathway-predictions.json`，SHA256 已记录）：CYP52M1 路线 A 2.01 s⁻¹、路线 B 3.72 s⁻¹；TES1 6.51/7.66 s⁻¹；POX1 9.37 s⁻¹（均为未校准预测值）。

**酶量敏感性扫描**（酶量 E 单位 mmol/gDW，路线专属扫描时封堵另一条路线）：

| E (mmol/gDW) | 路线 A（游离酸） | 路线 B（CoA 结合态） | 双路线+POX1 |
| --- | --- | --- | --- |
| 1e-5 | 0.0723 | 0.1338 | 0.2062 |
| 3e-5 | 0.2170 | 0.274011 | 0.274011 |
| ≥1e-4 | 0.274011 | 0.274011 | 0.274011 |

约束值与理论上限逐点吻合（如路线 A 在 E=1e-5 时 2.01×1e-5×3600=0.0723）。

**维持锚点通量所需最低酶量**（kcat × E × 3600 ≥ 0.274011）：

| 步骤 | 酶 | kcat (s⁻¹) | 最低酶量 (mg/gDW) |
| --- | --- | --- | --- |
| 路线 A ω-羟化 | CYP52M1 | 2.01 | **2.24** |
| 路线 B ω-羟化 | CYP52M1 | 3.72 | **1.21** |
| 路线 A 硫酯水解 | TES1/PTE1 | 6.51 | 0.45 |
| 路线 B 硫酯水解 | TES1/PTE1 | 7.66 | 0.38 |
| 前体氧化 | POX1 | 9.37 | 0.67 |

**试点读数**（假设成立时的可检验推论，不是结论）：

1. 0.274011 的模型通量在现实酶量（个位数 mg/gDW，对照酵母总蛋白约 500 mg/gDW）下**不受 kcat 限制**；约束只在酶量降到 ~0.4–2.2 mg/gDW 以下才生效。FBA 的"化学计量上限"与动力学可行域之间存在可量化的余量。
2. **动力学瓶颈是 CYP52M1 ω-羟化步**，且同等酶投入下路线 B 的极限通量约为路线 A 的 1.85 倍（0.134 vs 0.072 @1e-5），或同等通量下路线 B 只需一半的 CYP 酶量（1.21 vs 2.24 mg/gDW）。注意与 MetaTwin 综合评分的微弱 A 优（24.884 vs 24.681）方向相反——两法分歧本身就是决定性实验（分别以两路线构建株比较 CYP 表达量-产量关系）的选题依据。
3. 上述 A/B 比较建立在未校准 kcat 上（sd≈1.3 log10），差异在噪声内，只能作为假设进入实验设计，不能作为设计定论。

## 4. 试点二：上下文附加演示

酵母模型 baseline 场景（生长 0.9241，condition_id `sha256:68a56620…`）的 8 维向量附加到 11 个任务 1 候选的 BCR 推理：管线打通，log10(kcat) 位移 −0.04～−0.14。**按既有声明，FBA 输入通路在所有已发布检查点上未训练，附加后的数值不具备科学效力**；原始预测与附加预测并列存档于 `artifacts/fusion-yeast-pilot/context_attached_prediction_demo.json`。

## 5. 声明边界

- kcat 为未校准模型预测（validation-only 校准为负结果，conformal 区间为空）；酶量为假设值，未测量。
- FBA 通量为模拟值，不是实验通量；锚点复现只说明同一模型+同一协议得到同一数字，不构成对 0.274011 的实验支持。
- 外部项目（Metabolic model prediction）全程只读；目标工作流配置与预测文件的 SHA256 已记录在输出中。
- 本试点不改变任何已冻结结果、协议或发布清单（`fba_context.py` 不在 v6 冻结清单内）。

## 6. 下一步（按价值排序）

1. **分歧触发实验**：HFD1 排序分歧（MetaTwin 证据第 1 vs BCR kcat 第 4）与路线 A/B 酶量-产量关系，两项并入 MetaTwin 侧构建体设计矩阵。
2. **湿实验回流**：任何一块发酵/酶活数据都同时解冻 BCR 的 `governed Km training` / `experimental-flux training` 工作流（当前 `blocked_pending_external_data`）。
3. **校准工件**：为 ESM-2 基线重配 validation-only 校准，使酶量敏感性扫描的置信区间可报告。

## 7. 可追溯产物

- `configs/fba_context_model_yeast_metatwin.json`（模型身份）
- `scripts/fusion_yeast_kcat_fba_pilot.py`、`scripts/fusion_yeast_context_attach_demo.py`
- `artifacts/fusion-yeast-pilot/pilot_results.json`（锚点、扫描、酶成本；含全部 SHA256）
- `artifacts/fusion-yeast-pilot/context_attached_prediction_demo.json`（原始/附加预测并列）
