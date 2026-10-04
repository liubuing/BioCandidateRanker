# 融合设计矩阵增补件：把两个分歧变成实验（2026-10-04）

配套文件：`artifacts/fusion-yeast-pilot/design_matrix_addendum.csv`（机器可读，字段对齐 Yeast-MetaTwin 部署的 `06_evaluation/10h2da_construct_design_matrix.csv` 的 24 列 schema 并追加融合专属列）。上游文档：`docs/FUSION_YEAST_PILOT_20261004_CN.md`、`docs/10H2DA_YEASTMETWIN_COMPARISON_20261004_CN.md`。

外部项目（Metabolic model prediction / Integrated_Yeast_MetaTwin_Deployment）保持只读：本增补件不写入其目录，供其采用。原矩阵 80 行覆盖候选构建体的构建优先级（Tier1 = 双路线 CYP52M1+TES1）；本增补件不重复构建行，只增加**裁决性实验族**（F1–F3）与**标签回流规范**（F4）。

## 增补件的三个来源

| 来源 | 结论 | 转化为 |
| --- | --- | --- |
| 同任务对比的排序分歧（第二步醛脱氢酶） | MetaTwin 证据排序 HFD1 第 1 vs BCR kcat 先验把 HFD1 排第 4（5.0 vs ALD2 23.2 s⁻¹） | F1a/F1b：酶学测定 + 敲除面板 |
| kcat 约束 FBA 试点的核心读数 | 同等 CYP 酶量下路线 B 极限通量 ≈1.85×路线 A；维持 0.274 通量路线 B 只需 0.54× 的 CYP 酶量（1.21 vs 2.24 mg/gDW） | F2a：路线 × CYP 表达滴定双因素实验 |
| 对比与试点的共同缺口 | CYP52M1 对两个精确 C10 底物无任何 A 级证据；消耗侧只有定性假说 | F2b：精确底物体外酶学；F3a/F3b：示踪与双角色定量 |

## F1：第二步醛脱氢酶裁决（HFD1 vs ALD2/ALD3）

- **F1a（酶学，P0）**：纯化或 competent lysate 的 HFD1/ALD2/ALD3/ALD6 对 10-oxodecanoate 直接测 kcat、Km（NADH 340 nm 偶联，需验证偶联酶不限速）。这是纯裁决实验：任一方向的结果都有明确解释——HFD1 更快则 MetaTwin 证据排序正确、BCR 对 EC 1.2.1.3 长链底物的 kcat 先验偏保守（需复核）；反之则 kcat 先验正确、类别证据在该底物上不外推。
- **F1b（体内，P1）**：WT / hfd1Δ / ald2Δald3Δ / 三敲 + 回补，喂 10-hydroxydecanoate 测 10-oxodecanoate 积累与 sebacate 生成速率。三敲表型才等于单 hfd1Δ 时说明 ALD 有体内冗余（MetaTwin 已把 ALD2/3 标为 redundancy candidate）。
- 回流：F1a 产出 8 个直接标签（4 对 × kcat+Km，另有 4 个推导 kcat/Km），直接服务 BCR 的 `governed Km training` 与独立动力学基准扩充（均处 `blocked_pending_external_data`）。

## F2：路线 A/B 的酶量-产量关系（融合试点预测的正面检验）

- **F2a（体内，P0）**：两条 Tier1 路线背景 × CYP52M1 表达滴定（启动子强度/拷贝数梯度），每株定量 CYP52M1 丰度（表位标签或靶向蛋白质组学），LC-MS 测 10H2DA 滴度与路线中间体。判读看"滴度-CYP 丰度"初始斜率之比：B/A > 1.5 支持试点预测（超出未校准噪声的可辨识阈值）；< 1.2 不可判；< 1 支持原矩阵综合评分的微弱 A 优。
- **F2b（体外，P0）**：CYP52M1 + NCP1/CPR1（或同源 CPR）重构体系对 trans-2-decenoate 与 trans-dec-2-enoyl-CoA 测 kcat/Km，LC-MS 确认 ω 位产物（若主要产物是 ω-1/ω-2，则无论速率大小，路线比较作废——这正是原矩阵把 CYP505 只作区域对照的原因）。BCR 的预测值 2.01/3.72 s⁻¹ 首次获得被实测检验的机会。
- 回流：F2b 最多 4 个直接标签；F2a 的表观活性为情境参考，不入直接标签。

## F3：消耗侧定量化

- **F3a（示踪，P1）**：标记 10H2DA 喂 WT / pox1Δ / fox2Δ / pot1Δ / faa2Δ，定量 8-hydroxyoctanoate、乙酰-CoA 衍生物与 CO₂ 掺入。裁决"第一轮绕过 POX1"假说：消耗在 fox2Δ/pot1Δ 中下降即证实入口；全敲仍有的残余指向外排/酯化汇。链缩短产率即生产侧风险的定量。
- **F3b（双角色，P2）**：faa2Δ / 过表达 / WT 三比，分别测癸酸活化供给与 10H2DA 再活化消耗——两者之比决定 FAA2 该过表达、不动还是敲除（两法都只给出了"双重角色"的定性警告或接近的 kcat 对，2.1 vs 2.7 s⁻¹）。
- 回流：F3a 若以单位生物量速率报告，构成 `experimental-flux training` 工作流要求的"直接实验通量行"。

## F4：标签回流规范（数据飞轮的接口）

所有合格测量按 `docs/RESEARCH_PROTOCOL.md` 与 strict v1 normalized-kinetics 适配器入库：`evidence: "direct"`、正线性值 + 单位（kcat s⁻¹、Km M）、PubChem CID 底物、UniProt 序列、provenance 对象（实验记录/文稿哈希）。验收遵循 `docs/EXTERNAL_RESPONSE_ACCEPTANCE.md`。F1a+F2b 预期 ≥12 个直接动力学标签——这一个实验季就能把 BCR 三个被阻塞工作流推进到可评分。

## 判读纪律

- F2a/F2b 的 A/B 比较建立在未校准 kcat（sd≈1.3 log10）之上，方向一致才升级为设计定论；方向冲突时以实测为准并复核 BCR 上下文编码。
- 所有读数以 LC-MS 等靶向定量为准，不以生长或模型通量代替（沿用原矩阵的读出纪律）。
- 本增补件是设计清单，不是湿实验 protocol；剂量、培养与检测条件由执行方按其规范确定。
