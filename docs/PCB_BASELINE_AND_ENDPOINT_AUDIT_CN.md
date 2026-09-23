# PLOS Computational Biology 第一工作包执行结果

日期：2026-09-17。已执行源代码核查、数据输入重叠审计、本机资源检查和 CataPro 预测头运行检查。未训练完整新基线，未产生新的真实数据预测；独立外部池仍未评分。

## 关键发现与执行决定

1. **CataPro 公开权重不能直接作为当前内部测试的独立对照。** 官方 kcat 表的 27,658 条记录与本项目 1,646 条测试记录中，1,232 条序列相同，1,127 条序列与规范化 SMILES 对相同（约 68.5%）。这是官方十折数据并集的输入重叠，不是逐个 checkpoint 的训练成员审计，也未比较标签。官方推理平均十个 fold 的模型；应优先采用本项目固定划分重新训练。
2. **现有 context 不是 pH／温度。** 原始语料只有 ECNumber、Organism、Sequence、Smiles、Substrate、Type、Unit、Value 八列。collator 的四个 context 输入为 organism、ec、enzyme_type、reaction；UniKP 适配器不填 reaction，后者为空。原报告中的实验条件假设需要其他有条件元数据的数据集才能检验。
3. **当前是序列冷划分，不是底物冷划分。** 重新核验三个分区哈希均与冻结清单一致，跨分区完全相同序列与序列—SMILES 对为零；训练与测试共享 337 个规范化 SMILES，测试共 512 个不同 SMILES。这里计数的是底物种类，不是样本比例；不能据此宣称新底物泛化。
4. **聚合时上下文一致性需单列诊断。** 原始精确 Sequence／Smiles 分组内，有 192 组 Organism 不一致、204 组 ECNumber 不一致、6 组 Type 不一致。这些是包含后续被拒绝行的原始字符串分组统计，并非清洗后数量。现有聚合取 log10 标签中位数并保留首行元数据，需进一步在已接受训练行上审计，不直接改变已报告的冻结结果。

## 基线复现可行性

| 项目 | CataPro | CatPred |
|---|---|---|
| 固定官方代码 | cc89b2c81768665cf6fd76dfda607ce88691f601 | b314a28a84d388755237de90d9f336f951b291b3 |
| 代码许可 | MIT；已保存 LICENSE | MIT；已保存 LICENSE 与 LICENSE.txt |
| 输入与表示 | 蛋白序列、SMILES；ProtT5 1024＋MolT5 768＋MACCS 167 | 序列、SMILES；蛋白记录及相关依赖，具体训练配置待逐项锁定 |
| 目标差异 | kcat 为 log10(s^-1)；Km 为 log10(mM)，效率为 log10(s^-1 mM^-1) | kcat／Km／Ki；各配置的单位与输出转换须在适配时核验 |
| 上游运行方式 | 原脚本按给定 fold 做十折训练；推理平均十折模型 | 有论文复现脚本与预训练数据包；推荐 Linux GPU |
| 本机已完成 | 原版 kcat 预测头在 CUDA 上完成合成输入前向、反向和一步优化 | 依赖与脚本静态核查，尚未运行模型 |
| 主要剩余条件 | MolT5 权重与完整特征链；固定分区训练入口；编码器身份和预处理核对 | 隔离 Linux 环境、缺失依赖、数据包清单／体积、训练配置和重叠审计 |
| 决定 | 优先实现同划分重训，公开十折权重不进入独立主比较 | 并行准备环境与配置，不能用未审计网页预测代替基线 |

MIT 为代码许可，不自动解决上游数据库数据或独立基础模型权重的再分发条件。CataPro 公共 kcat 表中没有原始测量文献 DOI、pH、温度列，不能仅因表格较新就当作符合时间和来源门槛的独立外部数据。

## 已验证资源与运行边界

本机：Windows 11，RTX 3080，20,480 MiB 显存；当前 Python 3.14.5，PyTorch 2.12.0+cu130，CUDA 可用。CatPred environment.yml 指定 Python 3.9，并依赖当前环境缺少的 torch-scatter、progres、descriptastorus 等。需要独立环境，不能直接往现有科研环境混装依赖。

WSL 列举在当前沙箱返回 E_ACCESSDENIED；这不表示机器未安装 Linux。旧 MMseqs 日志指向 Ubuntu-24.04，环境是否可用于新基线须单独核查。

CataPro kcat 头有 506,113 个参数。8 条合成特征完成前向、反向和一步优化，耗时约 0.62 秒、峰值分配显存约 26 MiB。该测量不含编码器、真实特征、完整训练和 GPU 驱动保留显存，不能用来估计全流程速度或显存需求。

本地已有 ProtT5 文件，但尚未验证与 CataPro 所需模型完全一致。CataPro 对超过 1,000 氨基酸的序列保留前后各 500 位，现有 UniKP 缓存不能未核对预处理就复用。其 CLI 的 batch_size 参数也未传入示例特征生成链，不能把改 CLI 数字视为改变编码器批量。示例 CSV 以 index_col=0 读取，输入适配须保留索引列，避免误吞 Enzyme_id。

## 数据终点台账

机器可读表：artifacts/pcb-readiness/endpoint-register.csv。

| 数据 | 可支持的用途 | 尚未解决的问题 |
|---|---|---|
| 现有 16,443 个聚合样本 | 内部 log10 kcat 回归、探索性元数据诊断 | 逐条文献来源、再分发权限、聚合元数据一致性 |
| EnzEngDB 6,423 条／51 campaigns | 已观察的 campaign 内排序诊断 | fitness 单位跨任务不可直接合并；新规则需另找独立评价 |
| IMDH 512 基因型／1,024 辅因子观测 | 已观察的辅因子内排序 | 绝对单位常数未解决，且与训练序列同源 |
| 时间外部池 192 条／25 家族／54 底物 | 受治理候选池，暂不评价 | 仍缺至少 108 条、5 个家族等门槛 |
| 新前瞻性排序任务 | 尚未获得可评价独立名单 | 完整 roster、独立保管、冻结和标签交付 |

## 下一步已明确的实施范围

优先实现 CataPro 的固定 train／validation 接口与特征缓存。沿用原模型、损失与优化器，验证集用于早停；明确将“官方十折”改为“本项目冻结划分重训”的协议差异。测试标签不用于调参；编码器身份、长序列处理和随机种子写入补充协议后，再进行正式比较。

CatPred 先补齐包级配置和所需文件列表，探测独立 Linux 环境，再做官方示例运行。当前不能承诺全量训练耗时或宣称完成复现。

先对已接受训练行量化聚合元数据冲突，决定是否需要开发集上的敏感性分析。所有与已有测试结果相关的新解释均按事后探索处理。冻结外部池门槛不变。

## 复现与证据文件

运行 `python scripts/audit_pcb_local_evidence.py` 可重做本地哈希、元数据、输入重叠和资源审计。外部语料路径在脚本中显式指定且只读。

`scripts/audit_pcb_baseline_sources.py` 保存上游源文件、版本、哈希及获取错误；`scripts/audit_catapro_overlap.py` 检验固定数据 blob 后审计输入重叠；`scripts/smoke_catapro_head.py` 只运行合成特征检查。

输出目录：artifacts/pcb-readiness。主要文件：source-receipt.json、local-evidence-audit.json、catapro-data-audit.json、catapro-head-smoke.json、endpoint-register.csv。源代码与代码许可保存在对应方法的 source 子目录。

官方来源：[CataPro](https://github.com/zchwang/CataPro)、[CatPred](https://github.com/maranasgroup/CatPred)。审计使用的具体版本见上表及 source-receipt.json；当前主分支不能替代固定版本复现。
