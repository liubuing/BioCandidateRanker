# 2026-10-19 本机可补项收尾

## 范围

只做本机可复核、不新增外部数据、不放宽任何冻结门槛的工作。没有新增前瞻池记录、
没有新增家族、没有生成任何针对前瞻池的预测、没有发送任何外部请求。

## 已补

### 1. 版本与文档口径

- `README.md` 的校准工件说明已更正：`artifacts/uncertainty-calibration-v2/grouped-cv.json`
  与 `artifacts/full-homology-baseline/best.pt` 匹配（identity，1,640 条 validation），
  工作台可返回 conformal 区间；其他 checkpoint 仍无匹配工件。
- `README.md` 的 `training.py` 冻结说明从"v6 manifest"改为"v6–v9 manifests"。
- 前瞻池当前口径统一为：原始 **209 / 28 / 66 / 21 sources**；质量暂挂后
  **201 / 28 / 65**；门槛 300 / 30 / 50；原始缺口 91 条 / 2 家族，
  质量暂挂口径缺口 99 条 / 2 家族。
- `docs/CURRENT_RESULTS.md` 增加当前口径说明，并把发布列表补到 v8 / v9。
- `docs/RELEASE_AND_GOVERNANCE_TOOLING.md` 的缺口从 108 / 5 更正为 99 / 2
  （原始口径 91 / 2）。
- `configs/external_data_blockers.json` 的 `absolute_kinetics_expansion` 当前值更新为
  201 / 28 / 65，并显式记录 historical_pool 209 / 28 / 66、8 条 quality hold、
  carry-forward 缺口 99 / 2 与 raw 缺口 91 / 2。

### 2. 前瞻池面板

- `scripts/prospective_pool_dashboard.py` 默认 release manifest 从 v7 改为 v9。
- 重新生成 `artifacts/prospective-pool/dashboard.json` 与 `dashboard.md`：
  读取 2026-10-04 targeted audit，显示 209 / 28 / 66，`predictions_permitted = false`，
  仍差 91 条 / 2 家族。

### 3. 外部队列登记

此前 8 份草稿只有 4 份在队列里。现在 8 份全部登记（仍然 0 封发出）：

- `purdue-cdc14-files-v1`
- `sabio-rk-permission-sequences-v1`
- `yeast-gem-flux-metadata-v1`
- `prospective-campaign-collaborator-v1`
- `multicopper-oxygen-v1`（新增）
- `phage-fusion-construct-v1`（新增）
- `vapal2-active-core-v1`（新增）
- `absolute-kcat-panel-v1`（新增，收件人仍待定）

`tests/test_external_request_queue.py` 已同步为 8 条；`sent_at` 与
`message_reference` 仍全部为空，`messages_sent` 仍为 0。

### 4. 归档测试

`tests/test_project_state_archive.py` 原先把 `RESUME_HERE.md` 当作不可变归档文件，
但后续版本会正常更新这个"当前发布指针"，导致 1 项失败。现改为：

- 收据中的不可变状态文件仍逐字节校验；
- `RESUME_HERE.md` 按"当前发布指针"校验其确实指向某个 release manifest；
- 2026-07-27 的历史身份仍保留在收据与 git 历史中。

## 本机补不了的部分

前瞻池仍缺 **99 条记录 / 2 个新家族**（原始口径 91 / 2）。本机没有可合法新增的
合格直接动力学记录：第四、五轮候选来源没有一个同时通过 exact construct、底物、
饱和/端点、许可与 MMseqs 同源五道门；8 条 CtNDT 肌苷记录仍因来源一致性暂挂。
因此没有为了凑数而放宽门槛，也没有把工业过程数据混入动力学池。

要真正解阻，仍然只能走：

1. 定向作者/仓库请求（Purdue、SABIO-RK、yeast-GEM、三篇定向论文作者）；
2. 5 个无关家族的湿实验 panel（至少 108 条 attempted，20 条/家族上限）；
3. 若目标是工业阶段，另开过程数据线（titer、产率、时间序列、溶氧、酶丰度），
   不进入当前前瞻池。

## 验证

```powershell
python scripts\verify_software_release.py --manifest configs\software_implementation_release_v9.json
python -m pytest -q
python -m ruff check .
```

本机结果：

- release v9：`valid: true`；
- 全量 pytest：**425 passed, 0 failed**（8 subtests passed）；
- ruff：`All checks passed!`。

## 结论

本机能补的只是口径、面板、队列、校准说明和归档测试；**前瞻池的科学缺口没有被
"补"掉，也不应该用本机数据伪造掉**。当前仍是 `blocked_pending_external_data`。
