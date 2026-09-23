# CatPred 固定划分复现：环境与开发试运行

更新：2026-09-21。

## 最新：逐轮恢复队列

attempt2 已保留 1 个完成成员，第二成员中断。新增逐轮恢复机制，经真实 CPU 连续训练/中断恢复对照验证，最终权重、Torch/采样器随机状态和最佳验证分数一致。独立后台恢复队列已提交；后续以 `full-recoverable` 成员收据和 `recoverable-queue-stdout.log` 为准。剩余成员采用预先固定的独立种子，属于已披露的随机流修订。详见 `docs/CATPRED_RECOVERY_CN.md`。尚未完成全部成员或真实测试，不沿用早期完成时间估算。

## 当前执行状态（2026-09-20）

v3 首次作业日志最后停在 2026-09-19 22:41，首成员 epoch 10（第 11 轮）中途。核查时 WSL 中没有匹配的队列、runner 或训练进程；无异常堆栈，原因未确定。原 `running` 收据已审计为 `interrupted`，原收据副本及日志保留。完成成员数为 0，种子 42/123 未启动。

上游最佳检查点不含完整训练优化器状态，故未从该权重进行冒充精确续训的暖启动。按相同初始种子、配置和数据在 `full-v3-attempt2` 从头执行，重试协议为 `configs/catpred_mode_b_protocol_v3_attempt2.json`。这属于执行重试，无统计配置变更，不将中断结果计入正式重复次数。

新队列通过隐藏的独立 Windows WSL 进程启动，启动命令退出后已另行确认队列与 runner 存活。此方式减少对当前工具终端的依赖，但不保证机器关机或 WSL 被关闭后仍运行。日志：`artifacts/catpred-mode-b/full-v3-attempt2-stdout.log` 与 `full-v3-attempt2-stderr.log`；各组训练日志和状态在 `full-v3-attempt2/seed*/`。真实测试尚未执行。

## 全量训练进度

### v3 主存修正

**全量首轮已实测通过。** seed 7 的第一个成员完成 epoch 0，进入 epoch 1；首轮训练加验证约 122 秒，验证 RMSE 1.320997、MAE 1.013127、R² 0.100254，检查点已生成。这里只是一个未完成成员的开发指标，不能与测试表比较，也不能据此宣布完整基线完成。

全量训练采样中，进程主存高水位 15,744,616 KiB（约 16.1 GB），交换占用 0；GPU 总已用显存采样为 15,320 MiB（含其他进程与运行时，不等同于本模型 tensor allocated）。因此 556 MB 的映射后 RSS 不应写作训练峰值。进度证据：`full-v3/seed7/first-epoch-progress.json`。按首轮 122 秒粗算，三种子共 900 轮约 30 小时，仅为早期估算，需随完整成员耗时更新。

v2 在加载 ESM 缓存时出现 `OSError: [Errno 12] Cannot allocate memory`，436 秒后退出，未生成检查点，后续种子按队列规则未启动。GPU 显存压力测试不能替代全量主存检查。

新增 `scripts/catpred_memory_entry.py`，仅将 CPU ESM 缓存读取改为 `torch.load(weights_only=True, mmap=True)`；保留路径边界校验，其他检查点加载保持原路径。模型、样本、特征值、优化器及种子设置不变，仍使用已披露的非标准残基适配。

完整加载验证通过：6,992 条，逻辑特征 15,027,773,440 字节，全部映射后进程 RSS 556,220,416 字节，映射用时 145.5 秒。9、383、2,047 残基的三个代表性张量与普通读取逐值完全一致。该 RSS 是映射后的测量，不是完整训练的峰值保证。证据：`artifacts/catpred-mode-b/mmap-verification.json`。

已冻结 v3 协议并提交顺序队列：`scripts/run_catpred_full_queue_v3.sh`。当前日志：`artifacts/catpred-mode-b/full-v3-queue-relaunch.log`；训练与状态目录：`artifacts/catpred-mode-b/full-v3/seed7`。v1、v2 的协议、失败日志和文件均保留。当前仍须确认完整训练轮次与最终检查点，不能将加载验证通过计作模型实验完成。

### 2026-09-19 输入兼容修正与重启

首个原版全量作业在 epoch 0 遇到 `KeyError: 'X'`，退出码 1，耗时约 607 秒。上游序列注意力分支只为 20 种标准残基建立字典；开发集实际包含 X/U/B/O。旧目录 `full/seed7` 完整保留，其中初始检查点不属于训练完成的结果。

训练集 40 行、验证集 3 行含非标准残基。适配入口 `scripts/catpred_noncanonical_entry.py` 将该分支中的所有非标准残基映射到既有 index 20（零 embedding），保留序列长度、原始输入及 ESM2 特征；标准残基映射保持一致。不删除样本，不增加可训练参数，不修改上游源码。此策略与 padding 共用 embedding，属于明确的方法适配；正式报告必须标为“适配版 CatPred Mode B”，不可声称原版实现完全不变。预测时也必须使用相同入口适配。

验证完成：两项单元测试通过；包含全部 43 行非标准开发记录的 smoke（训练 56 行、验证 19 行、3 epochs）成功完成，日志时间 21 秒，指标为有限值。该验证仍使用验证集作为 upstream test 别名，其成绩不进入正式测试表。

已冻结 `configs/catpred_mode_b_protocol_v2.json`，保留 v1 协议及失败日志。训练设置未改变，seed 7 从头重跑；新结果目录为 `full-v2/seed7`。队列 `scripts/run_catpred_full_queue.sh` 将依次执行 7、42、123，任意失败即停止；不会自动执行真实测试。当前队列日志为 `artifacts/catpred-mode-b/full-v2-queue.log`，各组状态收据为 `full-v2/seed*/execution.json`。尚未获得完整全量成绩。

完整训练/验证特征缓存已生成并通过逐张量形状与有限值检查：6,992 条不同序列，1,084.8 秒（约 18 分钟），峰值 GPU allocated 3,596,895,744 字节。每个缓存文件均已记录 SHA256。

最长序列开发压力测试也已完成：32 条最长训练记录、16 条最长验证记录，batch 16、3 epochs；最长原序列 3,712 aa，上游 ESM2 分支按其既定长度上限截断，未改算法。日志时间 65 秒。其验证别名指标不属于正式测试结果。

全量配置已冻结到 `configs/catpred_mode_b_protocol.json`：种子 7/42/123，每种子 10 个成员、30 epochs、batch 16、MVE；每成员按最低验证 RMSE 选检查点。runner 核对源码、数据、依赖清单、特征 manifest 及全部缓存哈希后才启动。

原 v1 seed 7 已失败，见上方修正记录。v2 队列已提交，日志为 `artifacts/catpred-mode-b/full-v2/seed7/training.log`，状态收据为同目录 `execution.json`（输入校验完成后创建）。种子 42/123 需等待 seed 7 成功完成才会启动。不能把启动、首个检查点或微型开发指标当成完整基线结果。

正式入口：WSL 执行 `scripts/run_catpred_full_queue.sh`，或隔离 Python 运行 `scripts/run_catpred_full_v2.py --seed 7`。已有目录会被拒绝覆盖；中断运行须先审计再制定续跑方案。训练期间 upstream test 参数仍仅指向验证集。真正内部测试须等全部成员冻结后另行执行。

## 已完成

- 官方源代码固定到 `b314a28a84d388755237de90d9f336f951b291b3`，91 个选定源文件按 Git blob SHA1 核验。下载清单：`artifacts/catpred-mode-b/source-receipt.json`。
- Ubuntu-24.04-D WSL，项目内隔离 Python 3.12 环境，PyTorch 2.4.1+cu124，实测可用 RTX 3080。完整依赖保存在 `environment-freeze.txt`。
- 采用官方 kcat 管线：分子 DMPNN、序列注意力、650M ESM2、MVE 损失；未修改上游模型代码。650M 编码器与本项目 8M 编码器容量不同，应在正式比较中披露。
- 开发冒烟测试已成功结束：固定训练集前 32 行、验证集前 8 行，共 24 条不同序列；3 epochs、1 个模型、batch 16、seed/pytorch_seed 均为 7。完成特征缓存、训练、验证选模及检查点保存。
- 上游日志记录运行时间 1 分 20 秒；不包含首次环境安装、导入准备等完整墙钟成本，不能据此线性估计全量训练时间。

## 结果解释

上游训练入口要求评估分区，因此冒烟测试将同一份 8 行验证数据同时传入 `separate_val_path` 与 `separate_test_path`。日志中的 `test rmse=0.707919` 实际是这个微型开发子集指标，**不是真实测试集成绩，不能进入论文结果表**。真实内部测试集和外部候选池均未用于本次 CatPred 运行。

文件中的 pdbpath 是稳定的序列标识符，不代表已生成结构；本次 kcat 管线只使用序列记录，不使用结构坐标。

## 复现入口

1. Windows：`python scripts/prepare_catpred_smoke.py`。
2. WSL：`bash /mnt/d/biological/BioCandidateRanker/scripts/run_catpred_smoke.sh`。

已有输出应保留，后续正式实验使用新输出目录。试运行入口会拒绝覆盖已有模型。上游源码及依赖清单须与本次收据核对后再正式使用。

## 正式实验前仍需完成

1. 完成已经冻结的三种子训练。上游复现 shell 的循环变量未实际传给 seed 参数；当前 runner 显式传入 seed 与 pytorch_seed。
2. 长序列显存测试已通过；仍需记录全量单轮时间、主存占用和训练完成状态。
3. 只依训练/验证集选检查点，冻结完整集成后单次读取内部测试集。明确上游训练尾部开发集评估的用途，避免与正式测试混淆。
4. BioCandidateRanker 历史指标已通过历史 forward 与 batch 64 恢复。正式比较须分别报告历史版本与修复后重新训练版本，不能把旧权重的当前代码回放等同于重训结果。

当前状态是“开发验证与完整特征缓存完成，全量首种子已提交”，尚未完成 CatPred 全量公平对照。
