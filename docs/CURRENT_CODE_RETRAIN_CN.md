# 当前实现三种子重训

日期：2026-09-28。

最新状态：三种子训练与内部诊断均已完成。当前重训平均RMSE为1.41730，未建立强基线优势。详情见 `docs/CURRENT_RETRAIN_RESULTS_CN.md`；下文启动状态为历史执行记录。

## 本次执行

已启动独立后台队列 `scripts/run_current_code_retrain.py`，按7、42、123顺序从头训练完整模型。使用当前源码的冻结副本 `artifacts/current-code-retrain/source/biocandidate`，不加载历史模型权重。ESM2预训练权重除外：其本来就是冻结的输入编码器。

协议：`configs/current_code_retrain_protocol.json`。设置继承既有 ESM2 实验记录：ESM2-t6-8M冻结，d_model 256，8头，protein/molecule/fusion层数3/4/2，chunk128，context4096，task_query融合、heteroscedastic损失，AdamW学习率0.0003，batch4、累积8步，最多10轮、warmup2、cosine调度、验证损失早停patience4。以最低验证损失选模型。三种子不按测试成绩筛选。

此次是在已观察内部测试成绩之后冻结的重训，不构成独立确认性实验。原始数据与同源划分沿用历史哈希；本次队列只运行训练和验证，不输出测试成绩。三种子全部完成后，再冻结检查点并作一次内部诊断测试。

## 状态与产物

- `artifacts/current-code-retrain/queue-stdout.log`、`queue-stderr.log`：队列日志。
- `artifacts/current-code-retrain/seed{7,42,123}/train.log`：每轮训练/验证记录。
- 同目录 `execution.json`、`best.pt`、`latest.pt`：状态、最佳模型与最后轮次模型。
- `completed-training-manifest.json`：仅全部训练成功后生成。
- `environment-freeze.txt`：Python依赖记录；运行环境为Windows Python3.14、Torch2.12.0+cu130，与CatPred隔离环境不同。

启动已确认：seed7子进程存在、验证均值基线记录已输出、GPU正在使用。尚未宣称任何一组已完成。

## 中断边界

原训练CLI的resume分支未完整恢复scheduler、历史最佳验证值与早停计数。本队列因此拒绝对未完成目录自动续跑；中断时保留全部产物，需单独审计恢复策略。不得用该resume选项宣称精确恢复。队列持有Windows文件锁，防止重复启动。

## 外部验证准备

现有全局家族审计的FASTA和映射文件哈希已复核。合格池为192条记录、25个全局家族、54种底物。最小缺口是108条记录和5个新家族；不能仅从5个全新家族各取20条，因为这最多增加100条，仍需至少8条来自其他未达到上限的家族。更稳妥的搜集目标为至少6个新增家族和108条最终合格记录；候选数量需高于此数以容纳同源排除和来源核查。

补充顺序：

1. 先按论文和数据仓库核查2023年以后原始实验、测量表许可、绝对kcat与单位；不采集模型预测或相对活性充数。
2. 每行补齐DOI、表格/单元格、全序列或明确变体、底物结构、温度/pH（有则记录），并与既有来源注册表去重。
3. 对全部开发序列执行同源排除，再与已有候选一起做全局家族聚类；每家族最多20条，取舍不得依据标签或模型得分。
4. 达到300条、30家族、50底物，并冻结来源、成员、检查点和指标后才能评分。

本轮完成的是现有证据和缺口复核，未新增文献测量记录，未对外部候选评分。训练完成与外部验证完成是两个独立里程碑。
