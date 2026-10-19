# Purdue Cdc14（10.4231/429k-qe94）本机获取与暂存

## 1. 为什么现在做

外部请求 #1（Purdue Cdc14）此前记录为 `acquisition_blocked`，原因是本机访问
`purr.purdue.edu` 时 TLS 握手失败。2026-10-19 复核时该端点已可访问，因此直接从
公开仓库获取，不再等待邮件回复。数据许可为 **CC0 1.0**，测量文章为 JBC 2024
（DOI `10.1016/j.jbc.2024.107644`）。

## 2. 已获取的内容

- 页面：`https://purr.purdue.edu/publications/4552`
- 下载：`https://purr.purdue.edu/publications/4552/serve/1?render=archive`
- 外层 zip：52,298,036 字节，SHA256
  `4736b6b250f44e9fa5a06f96396112bee357aff03c43d4084a8d41c5222cfbca`
- 内层 bundle.zip：52,340,648 字节，SHA256
  `3cb195611eefe84e16b074fecff55057ae9acfe56fde7bed5e35845a2b50ac8c`
- 解包后 42 个数据文件，合计 126,206,610 字节，另有数据字典
  `Data Dictionary and Data Files Description.docx`（27,566 字节）。
- 本地路径：
  `artifacts/external/temporal-absolute-kinetics/purdue-429k-qe94/raw/purr-extract/10_4231_429K-QE94/bundle/`

## 3. 已暂存的候选行（未入池）

`Figure5/Fig5E_S3D_hCdc14_DiFMUP_Normalization.xlsx` 的 `Sheet1` 直接给出
3 次独立试验的 `kcat (s-1)`、`KM (µM)`、`kcat/KM (M-1 s-1)`，覆盖 6 个
hCdc14 构建体（hCdc14A WT/Q395A/3A；hCdc14B WT/QA/3A）。

已暂存为逐试验 CSV：

- `artifacts/external/temporal-absolute-kinetics/purdue-429k-qe94/candidate-kcat-km-20261019.csv`
- 18 行（6 构建体 × 3 试验），`admitted_records = 0`
- 底物：DiFMUP，**PubChem CID 2786976**（`C10H7F2O6P`）
- 测定条件（文章 Methods）：25 mM Hepes pH 7.5、150 mM NaCl、2 mM EDTA、
  2 mM TCEP、30 °C；连续荧光 DiFMUP 法
- 许可：CC0 1.0
- 采集与暂存收据：
  `artifacts/external/temporal-absolute-kinetics/purdue-429k-qe94/acquisition-receipt-20261019.json`

**注意**：此前本地 `raw/pubchem-difmup.json` 里的 CID 219176 是错误匹配
（对应一个含氟硅烷，不是 DiFMUP）。本次改用 PubChem CID **2786976**；
旧文件保留为历史，不应再作为底物映射依据。

## 4. 为什么还不能直接入池

这些行仍缺四项才能通过冻结协议：

1. **精确构建序列**：文章只写 hCdc14A 1-413 / hCdc14B 1-411、N 端 GST
   （pDEST15），没有给出 GST-linker-attB1 的确切连接序列，也没有说 GST 是否切除。
   协议要求最终测定链含标签、linker、切割 scar。
2. **hCdc14B 的 3A 位点**：补充表 S1 明确 QA 为 **Q390A**；但 3A 的三个位点在
   本地正文与补充表里只写“3A”，没有逐位点说明。
3. **逐 assay 酶浓度**：xlsx 给出的是已换算成 s-1 的速率与拟合 kcat；Methods 说
   蛋白浓度用 SDS-PAGE/ Bradford 测定，但未在表里逐行给出。
4. **同源与家族审计**：已对 hCdc14A 1-413 与 hCdc14B 1-411 跑 pinned MMseqs2
   （30% identity / 80% coverage / coverage mode 0），对 UniKP reference 的
   **命中数为 0**（7,822 条目标序列）。这支持“人 Cdc14 是一个新的同源冷家族”，
   但还没有做 20 条/家族 cap 与全局聚类。

因此本轮只把“运输失败”改成了“已获取、待策展”，没有把 18 行提升为接受记录，
也没有对任何模型跑预测。

## 5. 它对池子的实际意义

- 同源预检：hCdc14A 1-413 与 hCdc14B 1-411 对 UniKP reference（7,822 条）
  在 30% identity / 80% coverage / coverage mode 0 下 **0 命中**。这支持它作为
  一个新的同源冷家族，但最终家族划分仍要跑全局聚类。
- 这是 **1 个候选家族**（人 Cdc14A/B），不是 2 个；即使全部通过，也只能补
  1 个家族，仍缺至少 1 个家族。
- 6 条构建体级记录（或 18 条逐试验记录，按协议聚合规则取其一）相对于
  99 条缺口只是小头；但它把“Purdue 不可达”这个外部依赖从关键路径上拿掉了。
- 下一步是在本机完成 exact construct、hCdc14B 3A 位点、逐 assay 酶浓度与
  全局家族审计；若这几项闭合，再把行写入 strict registry，并重新跑 global audit。

## 6. 边界

- 全部为 CC0 公开数据的本机副本；没有发送邮件、没有使用预测值、没有推断
  缺失的构建序列或酶浓度。
- 未入池、未评分、未训练；`predictions_permitted` 仍为 `false`。
