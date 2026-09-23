# CataPro 固定划分训练记录与运行说明

2026-09-17，面向 PLOS Computational Biology 的现代基线工作包。

**2026-09-18 更新：三种子完整训练及单次内部测试均已完成。** 新增单独评分入口 scripts/evaluate_catapro_once.py，受 configs/catapro_internal_test_freeze.json 约束；完成后再次调用只校验保存记录，不重复推理。结果见 [内部对照报告](CATAPRO_INTERNAL_COMPARISON_CN.md)。

## 已实现

- 入口：scripts/run_catapro_mode_b.py。支持 cache 和 train 两个动作，只载入训练与验证分区，没有测试集评分入口。
- 协议：configs/catapro_mode_b_protocol.json，在新训练前记录源代码版本、原始分区清单哈希、编码器身份、优化器及模型选择规则。
- 预测头：原版 CataPro kcat 模型，506,113 个参数，输入为 ProtT5 1024 维、MolT5 768 维和 MACCS 167 维。
- 缓存：按编码器文件哈希、协议和输入身份分目录保存。训练前检查分区文件哈希、样本顺序、特征形状与文件哈希；已完成 checkpoint 再使用前也检查哈希。
- 输出：artifacts/catapro-mode-b 下的协议关联特征清单、编码器一致性记录、训练日志和最佳 checkpoint。

## 特征来源与一致性

MolT5 来自官方 laituan245/molt5-base-smiles2caption，固定 revision 7b7d4b0ab8b66b351e669b1f66272418ba15c3d9。下载 990,402,637 字节权重，SHA256 为 111fc277b175d7b0586905424448e20f1dce1dc1b1e912c6e447921b7d8621b1。文件清单见 models/molt5-base-smiles2caption/download-receipt.json。

复用本地 UniKP 的固定 ProtT5 编码器缓存。权重哈希与旧运行记录一致，预处理均为 UZOB 替换 X、超过 1000 氨基酸保留首尾各 500、去除 EOS 后均值池化。按训练输入长度选择 9、383、3712 氨基酸三个样本重新提取，最大绝对误差约 2.90e-7，低于预先设置的 1e-4 容差。该检查结合代码核查支持复用，不代表对每条缓存逐一重新提取。

MolT5 取三个固定训练样本，与固定上游 utils.py 的 CPU 原版提取结果比较；GPU 缓存最大绝对误差约 3.58e-7，MACCS 指纹完全相同。正式训练必须检查这份一致性记录与当前特征身份匹配。

MolT5 保留官方逐条输入和不额外截断的处理。长 SMILES 可触发分词器默认长度提示；本轮完整缓存生成成功且所有特征均通过维度和有限值检查。没有为了消除提示改变分子表示规则。

## 固定划分与训练规则

训练 13,157 条、验证 1,640 条；绝对参数目标为 log10(kcat / s^-1)。原版十折程序改为本项目固定 train／validation，属于明确披露的协议适配，不能称为原论文十折结果复现。

保留原版预测头、RMSE loss（epsilon 1e-6）、Adam、学习率 0.001、weight decay 0.01、batch size 64、dropout 0、最多 150 epochs、patience 20、min delta 0.001。按原版平均验证 batch RMSE loss 早停和选模型，同时单独记录整个验证集的 RMSE 与 MAE。三随机种子为 7、42、123。

实现上的差异：使用 GPU randperm 生成训练顺序；若出现单样本末 batch 则合并到前一个 batch，以满足 BatchNorm；冻结编码器只提取特征，不参与优化；均值池化使用向量化计算并检查数值一致性。当前训练规模末 batch 不为 1。未设置完全确定性内核保证，不宣称跨硬件逐位重现。

与原 SOTA 协议的区别：CataPro 使用验证集早停；不修改旧 UniKP／DLKcat 的冻结运行规则。不同模型的选择方式须在最终比较表中披露。

## 已执行检查

32 条训练和 32 条验证样本的完整特征与两轮训练已跑通。这是接口检查，不用于选择超参数，也不构成模型性能结论。

全部 train／validation 特征缓存已生成，包含 2,510 个不同底物的分子表示。本次完整缓存步骤耗时约 111 秒，依赖已有且经核验的 ProtT5 缓存，不能视为从零提取两种编码器的总耗时。

5 项测试覆盖协议改动、缓存内容改动、样本身份错配、缺失测试文件时仍能只读 train／validation，以及长序列／非标准氨基酸预处理。

## 复现命令

在仓库根目录使用已验证的科研 Python 环境：

```powershell
python scripts/prepare_catapro_molt5.py
python scripts/run_catapro_mode_b.py cache --limit 32
python scripts/run_catapro_mode_b.py train --mode smoke
python scripts/verify_catapro_molecule_features.py
python scripts/run_catapro_mode_b.py cache
python scripts/run_catapro_mode_b.py train --mode full
python -X pycache_prefix=artifacts/catapro-mode-b/pycache -m pytest tests/test_catapro_mode_b.py -q
```

下载需要网络，其他模型加载为本地模式。ProtT5 路径显式保留现有外部只读部署位置。程序拒绝身份不一致的缓存或已完成运行；不要通过覆盖旧协议消除报错。重新生成 feature-manifest 会改变其运行记录哈希，已有训练需要使用原特征清单才能匹配。

## 结果边界

验证集分数参与早停，不能拿它替代内部测试或外部表现。训练入口不读取原 1,646 条测试集；该测试随后由单独冻结入口完成评价。时间外部池和前瞻性排序数据未被访问。

三种子 checkpoint 和单次测试入口现已固定。与其他模型比较时报告输入、模型容量、训练选择方式和不确定性；下一步补 CatPred 和家族层面的配对分析，独立外部结论仍须满足既有门槛。
