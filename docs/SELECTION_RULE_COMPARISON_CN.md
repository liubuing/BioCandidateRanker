# 选模规则对照执行记录

日期：2026-09-28。

## 设计

仅训练/验证数据，种子7/42/123。沿用此前模型、冻结ESM2权重、优化器、学习率、batch4、累积8步、warmup2和10轮cosine调度。为避免某一种选模指标主导停止时间，将patience从4改为0，三组均完成预定10轮。此停止规则变更在运行前冻结。

每组共享同一训练轨迹：`best.pt`保存最低验证高斯损失的模型，`best_rmse.pt`保存最低验证RMSE的模型；严格小于才更新，相等时保留最早轮次。增加的保存操作不改变反向传播、梯度或优化器更新。

源码副本位于 `artifacts/selection-rule-comparison/source/biocandidate`。只在副本中增加RMSE检查点保存逻辑，既有重训源代码及模型不改动。协议为 `configs/selection_rule_comparison_protocol.json`。

## 状态与运行

已提交三种子后台队列，顺序执行。状态见 `artifacts/selection-rule-comparison/seed*/execution.json`，训练过程见各目录 `train.log`；只有全部成功后才生成 `completed-training-manifest.json`。没有自动测试集推理步骤。

执行脚本：`scripts/run_selection_rule_comparison.py`。存在未完成目录时拒绝静默重跑或不完整恢复；进程持锁防止重复启动。

全部完成后运行 `scripts/summarize_selection_rules.py`：核验30轮日志、两类检查点哈希、选模轮次是否等于对应指标的最小值，并生成 `selection-summary.json`。

## 判断边界

比较验证RMSE、MAE、Pearson、高斯损失、覆盖率、选模轮次和种子波动。不因为某个种子结果好就单独推广。按验证RMSE选模后，该验证RMSE更低具有选择上的必然性，不能当成独立泛化证据。两种规则的最终泛化比较需要训练分区内的家族嵌套划分或达标的独立外部数据；不得继续使用已观察内部测试集筛选规则。

预计量级：前一轮三种子共17个epoch约59分钟；本轮30个epoch在类似负载下约1.5–2小时，仅为排程估算。以实际进程、日志和完成清单为准。
