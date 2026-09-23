# CataPro 现代基线实测结果

2026-09-18。完整特征缓存、三种子固定划分训练及单次内部测试已完成。

| 方法 | RMSE | MAE | Pearson |
|---|---|---|---|
| BioCandidateRanker ESM-2 (ours, reference) | 1.3916 ± 0.0183 | 1.0708 ± 0.0117 | 0.3810 ± 0.0222 |
| UniKP Mode B (retrained) | 1.4093 ± 0.0094 | 1.0785 ± 0.0041 | 0.3777 ± 0.0222 |
| DLKcat Mode B (retrained) | 1.7138 ± 0.1415 | 1.3351 ± 0.1251 | 0.1463 ± 0.0764 |
| CataPro fixed-split retraining | 1.3876 ± 0.0195 | 1.0650 ± 0.0105 | 0.4320 ± 0.0004 |

**判断：CataPro 的平均 RMSE 和 MAE 数值略低、Pearson 数值更高；本项目尚未证明优于该现代基线。** RMSE 仅相差约 0.0040，不能据三种子标准差认定统计显著。所有结果为原有 1,646 条内部测试集，均不是独立外部验证。

训练集 13,157 条、验证集 1,640 条。CataPro 采用官方预测头与默认优化参数，并用固定验证集早停，替代官方十折训练；UniKP 与 DLKcat 按其既有冻结协议训练。模型容量、输入模态及模型选择方式不同，不能把这个比较解释为单独检验融合结构优劣。

CataPro 三种子训练分别运行 28、23、31 轮，最佳 epoch（从 0 开始）为 7、2、10。总计约 72 秒为缓存特征后的预测头训练时间，不含编码器下载或从零蛋白表示提取。

复用 ProtT5 缓存经过权重、预处理和短／中／长序列数值检查；MolT5 与原版 CPU 提取器及 MACCS 也通过检查。7 项针对缓存、数据分区和冻结文件的测试通过。保存预测的行顺序、标签、哈希和重算指标均核验一致。

下一步应补 CatPred 并进行按家族的配对差值分析；在没有新证据前，不因看到测试结果而调整 CataPro 或本项目的超参数。论文重点仍需落在可独立复验的生物学问题与泛化边界。

记录：artifacts/catapro-mode-b/internal-test/summary.json、comparison.json、seed*-predictions.csv；冻结协议：configs/catapro_internal_test_freeze.json；运行说明：docs/CATAPRO_MODE_B_RUNBOOK_CN.md。
