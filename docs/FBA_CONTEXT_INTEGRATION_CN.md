# 基因组尺度代谢通量融合（FBA 上下文）

日期：2026-09-25。本模块将两套外部技术栈融合进 BioCandidateRanker 的既有接口，作为功能拓展。

## 融合了什么

| 来源 | 技术栈 | 本项目落点 |
|---|---|---|
| 辅因子改造论文（Wang et al.，D-泛酸生产，浙江工业大学） | **E. coli iML1515** 基因组尺度模型 + FBA/FVA + 基因敲除/过表达边界模拟 + 葡萄糖 15 / 氧气 20 mmol/gDW/h 培养基 | `--preset baseline / ptsG_knockout / thd2pp_overexpression / nadK_overexpression / gltA_overexpression / paper_stack` |
| 用户精氨酸/固碳项目（`D:\biological\大肠\05_FBA计算验证\融合FBA计算.py`） | **E. coli iJO1366** + 8 条合成反应（CCR_EMA、NOG_F6P、PTXD 等）+ NOG/亚磷酸/CO₂ 场景上限 | `--preset fusion_M2_NOG / fusion_M3_PtxD / fusion_M5_NOG_PtxD`（仅 iJO1366） |

两个模型都输出**同一契约**：8 维 `fba_context` 向量（`ModelConfig.fba_context_dim = 8`，与所有已冻结检查点一致），特征反应为各模型生物量目标 + 葡萄糖交换 + EMP（PFK）/ PPP（G6PDH2r）/ ED（EDD）/ TCA（CS）/ 转氢酶（THD2pp）/ ATP 维持（ATPM），变换为 signed log10(1+|v|)。

## 文件与身份

- `src/biocandidate/fba_context.py`：多模型注册表、预设、FBA/FVA 运行、特征提取
- `data/iML1515.xml`、`data/iJO1366.json`：模型资产（gitignored `data/`），身份记录于
  `configs/fba_context_model.json` 与 `configs/fba_context_model_ijo1366.json`，加载时强制校验 SHA256
- `scripts/serve_workbench.py`：`GET /api/fba/status`、`POST /api/fba/run`；`POST /api/predict` 接受可选 `fba` 字段
- 工作台新增"代谢通量"页签：选模型、选场景、跑 FBA、查看通量与 FVA 区间，一键"附加到预测"

## 两个如实声明

1. **附加后的预测不具备科学效力**。所有随附检查点均在无通量标签的数据上训练，模型的 FBA 输入通路处于未训练状态；附加上下文改变输出数值只说明"管线打通"，不构成任何验证。预测响应中的 `fba.note` 字段与页面提示均已写明这一点。
2. **FBA 输出是模拟值，不是实验通量**。论文自身的通量也是 iML1515 预测（作者明确指出未来需 13C 标记验证）；大肠项目的融合场景同为模型计算。二者均不满足 `experimental_flux_training` 工作流"直接实验通量行"的解冻门槛。

## 实现中发现并修复的两个模型陷阱

1. **共享模型对象的结构性污染**：cobra 的上下文管理器不能完全回滚"新增反应/代谢物"，跨调用复用同一模型并重复添加融合层会破坏质量平衡，得到看似 optimal 实则错误的解（表现为目标值 0.0094）。修复：每次运行在 `model.copy()` 上进行，且融合层添加改为幂等。
2. **iJO1366 的双生物量反应**：模型文件默认目标是 `BIOMASS_Ec_iJO1366_core_53p95M`，而精氨酸项目使用 `BIOMASS_Ec_iJO1366_WT_53p95M`。若不显式设定目标，对 WT 施加生长下限而最大化 core 会产生方向完全相反的"最优解"。修复：`run_fba_context` 强制 `scoped.objective = 注册表声明的目标`，目标函数成为场景条件的一部分（进入 condition_id 哈希）。

## 使用

```bash
# 命令行
python -m biocandidate.fba_context --model-id iML1515 --preset paper_stack --fva
python -m biocandidate.fba_context --model-id iJO1366 --preset fusion_M5_NOG_PtxD

# 工作台
python scripts/serve_workbench.py    # http://127.0.0.1:8808/ → 代谢通量页签
```

`run_fba_context` 返回 `condition_id`（对模型哈希、预设、培养基、生长下限的规范化哈希），使任意上下文向量可追溯到其产生场景；该 id 随 `FBAFeatureMetadata` 一同进入预测请求，满足既有的身份校验链。


## 代理蒸馏：让原模型本身获得 GEM 能力（2026-09-25 增补）

按用户澄清的需求——"融合"指**原模型（ranker 架构）本身具备两个 GEM 的能力**——新增 `src/biocandidate/fba_surrogate.py`：

- **场景网格**：培养基（葡萄糖 ×6 档 × 氧气 ×9 档，覆盖论文 15/20 工作点与大肠项目 10/20 场景）与遗传场景（敲除/过表达/融合途径）全交叉；
- **批量求解**：单副本逐场景求解（边界改动逐场景回滚），iML1515 240 场景约 12 秒；
- **蒸馏训练**：以 8 维场景向量为 `fba_context` 输入、FBA 生物量通量为 `log10_flux` 标签，从零训练 ranker 架构检查点（蛋白/分子模态关闭，与既有 simulated-flux smoke 同款做法）；
- **效果**：iML1515 代理 test RMSE(log10 生长) 0.0116（均值基线 0.2799，约 24 倍）；iJO1366 代理 0.0149（基线 0.3361，约 23 倍）。线性尺度约 2–3% 生长误差；
- **推理接口**：`POST /api/fba/surrogate_predict` 以场景参数直接返回代理生长速率，工作台"代谢通量"页签内置"代理 vs FBA 真值"对比面板。

推理侧的一个关键细节：训练载体记录带有固定的 organism/EC/酶类型/反应上下文字符串，`use_context=True` 时这些字段参与前向哈希——推理若用默认空字段会查询模型从未见过的上下文，产生无意义输出。服务端代理端点已在内部构造规范载体，调用方无需关心。

### 声明边界（增补）

代理检查点是对**模拟输出**的蒸馏：在采样分布内复现 FBA 结果，分布外不保证；不是实验通量模型，不替代 GEM，不构成对真实细胞的任何验证。检查点与数据集清单中的 `claim_boundary` 字段均已声明。
