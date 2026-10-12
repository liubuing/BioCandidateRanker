# 植物 CCD 序列预筛与单位核查

日期：2026-10-03。来源：10.3390/biology14070752，原始全文 PMC12292525。

## 完成结果

已核实两条公开蛋白参考序列，并完成冻结阈值下的开发集同源预筛：

| 酶 | 参考登录号 | 长度 | 合格开发集命中 |
|---|---|---:|---:|
| OeCCD1 | XP_022895849.1 | 542 aa | 0 |
| InCCD1 | XP_019150854.1 | 567 aa | 0 |

版本与参数：MMseqs `5d152c612b6ad2a56f657b7a02c127eceaea2a75`，identity 30%、coverage 80%、cov-mode 0。开发 FASTA 哈希与冻结值一致。
查询为论文给出的 CCD 催化组分，不是完整 GST 融合蛋白。零命中仅表示通过本次组分预筛，不表示已完成实测构建审计或已增加两个全局家族。

## 为什么暂缓纳入

全文第 2.7 节定义 1 U 为每分钟产生 **1 nmol** β-ionone；第 3.3 节报告两种 GST 融合蛋白分子量约 87.43 和 89.55 kDa。
根据表 2 的 Vmax 换算如下：

| 酶 | Vmax（原文 U/mg） | 按原文 nmol 定义换算 kcat（s⁻¹） | 原文 kcat（s⁻¹） |
|---|---:|---:|---:|
| OeCCD1 | 2.30 | 0.00335148 | 3.35 |
| InCCD1 | 1.22 | 0.00182085 | 1.82 |

公式为 `kcat = Vmax × 分子量(Da) / 60000000`，其中 U 按文中 nmol/min 定义。
原表 kcat 与把 U 改作 μmol/min 后的换算值相符，差异约为 1000 倍，远大于四舍五入可解释的范围。
这可能是方法部分的单位笔误，也可能涉及其他未交代的归一化；当前证据不足以判断哪个数值应修正。

因此两条记录标为 `held_unit_definition_conflict`，保留原值，不生成更正标签，也不加入正式池。
动力学部分只笼统写最适条件，未明确列出浓度序列；完整融合链及其他测量条件的核查尚未完成。
在这一单位问题澄清前，暂停该来源的进一步逐项整理，转向其他来源。

## 可核对的问题草稿

> 第 2.7 节的酶活单位是否应为 μmol/min，而不是 nmol/min？表 2 的 Vmax 和 kcat 分别使用什么单位、酶分子量和浓度归一化方式？可否提供实际底物浓度序列以及完整 GST 融合链？

仅保存草稿，未联系作者。

## 产物与当前总量

目录：`artifacts/external/temporal-absolute-kinetics/ccd-prescreen-20261003/`。
包含两条参考 FASTA、`catalytic-components.fasta`、`development-hits.tsv`、`mmseqs.log`、`homology-audit.json`、`unit-consistency.json` 和带处置状态的 `source-records.csv`。
复现脚本：`scripts/prescreen_plant_ccd.py`。

本轮新增正式记录 **0**。现有池仍为 **192 条、25 家族、54 底物**，未运行模型评分。
原检索清单保留为初筛快照；本报告和本来源审计回执更新该来源的后续状态。
