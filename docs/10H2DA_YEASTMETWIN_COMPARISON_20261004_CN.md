# 10H2DA/10HDA 任务：BioCandidateRanker 与 Yeast-MetaTwin 预测对比

日期：2026-10-04。性质：跨项目对比练习（research exercise），不是基准评估，不产生任何发表级结论，也不是独立验证。

## 1. 目的与范围

用本项目（BioCandidateRanker，`artifacts/full-homology-baseline/best.pt`，ESM-2 多模态 log10(kcat) 模型，历史同源冷分离测试 RMSE 1.3916 / Pearson 0.381）对 "D:\biological\Metabolic model prediction" 中 2026-09-02 的两份 Yeast-MetaTwin 预测做同任务重演，然后逐项对比：

1. **任务 1**：10H2DA（(E)-10-hydroxydec-2-enoic acid，蜂王酸）酿酒酵母代谢路径预测；
2. **任务 2**：10HDA（10-hydroxydecanoic acid）末端醇氧化候选酶排序。

Yeast-MetaTwin 侧基准来自 `10H2DA_酿酒酵母_Yeast-MetaTwin代谢路径预测报告.docx`、`10HDA_末端醇氧化候选酶分析.md` 与 `10HDA_末端醇氧化候选酶排序.csv`。

## 2. 执行方式

- 序列：15 个候选酶的 reviewed UniProt FASTA（2026-10-04 通过 REST 下载），存于 `artifacts/10h2da-yeast-prediction/sequences/`。
- 底物：10 种（含 3 个 CoA 硫酯），全部通过 RDKit 解析校验；SMILES 与原子计数见 `build_candidates.py`。
- 输入构建：`artifacts/10h2da-yeast-prediction/build_candidates.py`（任务 1 共 11 个酶-底物对，任务 2 共 7 个）。
- 推理：torch 2.12 下 CLI `predict` 因冻结 checkpoint 先于 torch 2.6 而失败（README 已知问题），改走 workbench 安全加载路径（`run_predictions.py` 直接复用 `scripts/serve_workbench.py` 的 `safe_globals` 允许清单与 `run_prediction` 流程），CUDA 执行。
- 原始输出与输入身份哈希：`task1_10h2da_pathway-predictions.json`、`task2_10hda_terminal_oxidation-predictions.json`。
- 无校准工件（README：当前无匹配校准文件），conformal 区间为空；下述不确定性均为原始（未校准）模型标准差。

## 3. 任务 2 对比：10HDA 末端醇氧化候选酶排序

预测值换算为线性 kcat（s⁻¹）便于阅读；sd 为 log10 尺度。

### 第一步：10-hydroxydecanoate → 10-oxodecanoate（EC 1.1.1.1）

| 候选酶 | MetaTwin 排名/评分 | MetaTwin 依据 | BCR 预测 log10(kcat)（线性 kcat） |
| --- | --- | --- | --- |
| SFA1 (P32771) | 1 / 100 | 直接底物证据：纯化酶对 10-hydroxydecanoate 的 Km=7.2 mM | **1.375（23.7）** |
| ADH2 (P00331) | 2 / 45 | 同家族泛活性 | 1.228（16.9） |
| ADH1 (P00330) | 3 / 25 | 泛活性，生理方向常为还原 | 1.264（18.4） |

**第一步一致**：两法均把 SFA1 排第一。但 ADH1/ADH2 的次序两法相反（差距 ~0.04 log，远小于 sd≈1.28，无统计意义）。

### 第二步：10-oxodecanoate → sebacate（EC 1.2.1.3/1.2.1.4）

| 候选酶 | MetaTwin 排名/评分 | MetaTwin 依据 | BCR 预测 log10(kcat)（线性 kcat） |
| --- | --- | --- | --- |
| HFD1 (Q04458) | 1 / 80 | 长链脂肪醛脱氢酶类别证据；无精确 10-oxodecanoate 测定 | **0.700（5.0）——四个中最低** |
| ALD2 (P47771) | 2 / 45（与 ALD3 并列） | 泛醛氧化 | 1.366（23.2） |
| ALD3 (P54114) | 2 / 45（与 ALD2 并列） | 泛醛氧化 | 1.352（22.5） |
| ALD6 (P54115) | 3 / 30 | 以乙醛为主 | 1.260（18.2） |

**第二步实质性分歧**：MetaTwin 凭策展类别证据把 HFD1 排第一；BCR 把 HFD1 排最后（与 ALD2 差 ~0.67 log ≈ 4.7 倍）。该差距仍小于未校准 sd（~1.28），排序本身不可判优。可能原因：(a) EC 1.2.1.3 上下文哈希桶稀释了 HFD1 的底物特异性信号；(b) kcat 先验对长链脂肪醛脱氢酶家族校准不足；(c) 线性 kcat 本就不是 MetaTwin 评分所依据的量（其证据核心是底物类别匹配）。**此分歧上，MetaTwin 的证据驱动排序更可靠。**

## 4. 任务 1 对比：10H2DA 酿酒酵母路径

### 4.1 BCR 对路径末端酶与消耗侧酶的动力学预测

| 酶-底物对 | 角色 | log10(kcat)（线性 kcat） |
| --- | --- | --- |
| CYP52M1 × (E)-dec-2-enoic acid | 路线 A ω-羟化 | 0.303（2.0） |
| CYP52M1 × (E)-dec-2-enoyl-CoA | 路线 B ω-羟化 | 0.570（3.7） |
| CYP52M1 × decanoic acid | 家族活性锚点 | -0.001（1.0） |
| TES1 × (E)-dec-2-enoyl-CoA | 硫酯水解（两条路线共用） | 0.814（6.5） |
| TES1 × (E)-10-hydroxydec-2-enoyl-CoA | 路线 B 产物释放 | 0.884（7.7） |
| CYP505 × (E)-dec-2-enoic acid | 区域选择性对照（亚末端） | 0.454（2.8） |
| POX1 × decanoyl-CoA | 前体氧化/竞争 | 0.972（9.4） |
| FOX2 × (E)-10-hydroxydec-2-enoyl-CoA | 消耗入口（绕过 POX1 假说） | **1.457（28.6）** |
| POT1 × 3-oxodecanoyl-CoA | 链缩短硫解酶 | 1.176（15.0） |
| FAA2 × decanoic acid | 前体活化（双重角色） | 0.332（2.1） |
| FAA2 × 10H2DA | 消耗再活化（双重角色） | 0.433（2.7） |

（NCP1/CPR1 是 P450 电子供给蛋白，底物不是目标分子，不属于本排序器可评分的酶-底物对，未纳入。）

### 4.2 与 MetaTwin 报告结论的对照

- **路线 A vs B 不可分，两法一致**。MetaTwin：FBA 通量同为 0.274011，综合分 24.884 vs 24.681，建议并行验证。BCR：路线 B ω-羟化 kcat 略高于路线 A（0.570 vs 0.303），差距远小于 sd≈1.33-1.35，同样不可分。
- **消耗风险获得量化支持**。MetaTwin 定性警告：第一轮 β-氧化可能经 FOX2/POT1 绕过 POX1，单敲 POX1 不充分；FAA2 有促进生产与促进消耗的双重角色。BCR 的预测把这一警告定量化：消耗侧 FOX2（28.6）、POT1（15.0）、POX1（9.4）的预测 kcat 全面高于生产侧 CYP52M1（1.0-3.7）与 TES1（6.5-7.7）；FAA2 对 10H2DA 的再活化 kcat（2.7）与对 decanoate 的生产活化（2.1）相当——"双重角色"得到数值印证。
- **Tier-1 组合不冲突**。BCR 未发现任何比 CYP52M1/TES1 更优的替代末端酶（CYP505 对照的 0.454 介于两者之间，与"仅作区域对照"的定位一致）。

## 5. 能力差距矩阵

| 维度 | Yeast-MetaTwin | BioCandidateRanker |
| --- | --- | --- |
| 路径枚举/发现（逆合成、FBA 可行性） | 有：两条候选路线 + 质量电荷平衡 + 通量上限 0.274011 | 无：只能对给定酶-底物对打分，不做路径推理 |
| 可能方式途径的覆盖面 | 宽：生产路线、降解/消耗假说、区室衔接、副产物二酸支路 | 窄：仅在给定候选内排序；途径层面依赖外部给定 |
| 消耗/降解路径 | 有定性假说分级（β-氧化最可信，POX1 可绕过） | 无路径推理，但能给消耗侧酶定量 kcat 先验 |
| 酶级动力学定量 | 间接：UniKP 类分数仅用于排序，明确声明不等价于 kcat | 直接：log10(kcat)+不确定性；但无 Km 头（未训练），给不了 kcat/Km |
| 调控手段（工程杠杆） | 具体且分级：ACC1/FAS1/FAS2、FAA2 双重角色警示、ZWF1/GND1 NADPH、POX1/FOX2/POT1 敲除条件、外排工程、区室衔接方案 | 无工程设计能力；只能按 kcat 先验提示"哪些酶消耗快需要关注" |
| 不确定性 | 证据等级 A/B/C 分级 | 模型 sd，但未校准（README：validation-only 标量校准为负结果），区间尚不可用于决策 |
| 证据治理 | UniProt/Rhea/PubMed 检索 + 精确底物证据标注 | 输入哈希固定、无校准工件则不出区间、非验证警告 |
| 关键弱点 | 无 A 级精确底物证据；FBA 无酶容量/毒性/热力学 | 上下文哈希桶对未见 EC（1.14.14.80 等）信号弱；单 checkpoint；排序在噪声内 |

## 6. 局限性

- 单一冻结 checkpoint、单次推理；未做多种子或集成。
- EC 1.14.14.80、1.3.3.6 等上下文与 CoA 硫酯底物在 UniKP 开发语料中稀少，哈希桶化上下文对这些候选的信号贡献有限，预测主要受序列与分子图驱动。
- 不确定性未校准：所有排序差异均不可单独作为决策依据。
- kcat 只是催化能力的一个维度；MetaTwin 第一步 SFA1 的核心证据是 Km=7.2 mM（亲和性），本模型当前无法复核该维度（Km 头未训练）。
- 本对比不是外部验证：两者都无 10H2DA/10HDA 精确底物的 A 级实测数据可比。

## 7. 结论

1. **可能方式途径的覆盖面：Yeast-MetaTwin 明显更广**。路径枚举、两条候选路线、降解假说、区室衔接、FBA 通量全部只有 MetaTwin 侧具备；BioCandidateRanker 在路径层面没有独立预测能力，其贡献限定在"给定候选的酶级动力学排序"。
2. **调控手段的可靠性：MetaTwin 更可用，BCR 提供补充定量**。MetaTwin 的工程杠杆带证据分级和条件警示（如 FAA2 双重角色、POX1 单敲不充分），直接可执行；BCR 的预测与这些警示方向一致并把消耗风险量化（消耗侧酶预测 kcat 高于生产侧），但其未校准的不确定性和 kcat 单维度（无 Km）使它不能独立支撑调控决策。
3. **分歧点明确一处**：第二步醛脱氢酶排序（HFD1）。MetaTwin 的策展证据排序比 BCR 的 kcat 先验更可靠；BCR 把 HFD1 排最后提示其对该家族的 kcat 先验可能偏保守，值得在获得实测数据后复核。
4. 两法对"路线 A/B 并行验证、Tier-1 用 CYP52M1+TES1"和"消耗风险真实存在"这两个核心结论一致。

## 8. 可追溯产物

- `artifacts/10h2da-yeast-prediction/sequences/*.fasta`（15 个 UniProt reviewed 序列）
- `artifacts/10h2da-yeast-prediction/build_candidates.py`、`run_predictions.py`
- `artifacts/10h2da-yeast-prediction/task1_10h2da_pathway.json` / `task2_10hda_terminal_oxidation.json`（输入，含 RDKit 校验 SMILES）
- `artifacts/10h2da-yeast-prediction/task1_10h2da_pathway-predictions.json` / `task2_10hda_terminal_oxidation-predictions.json`（原始输出，含 checkpoint 身份与输入 sha256）
