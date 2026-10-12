# BioCandidateRanker → Yeast-MetaTwin 平台插件接入报告

日期：2026-10-12
接入方向：BioCandidateRanker（BCR）作为 kcat 预测插件接入 Integrated_Yeast_MetaTwin_Deployment
状态：**已完成接入并验证；两侧各一个提交，均未推送**

---

## 一、结论摘要

BCR 已作为 `kcat_prediction` 插件注册进平台，与既有的 CLEAN、DLKcat、UniKP 并列。平台自身的插件测试 8/8 通过，端到端 smoke 返回 `ready`。插件清单如实标注为 `blocked`（原因见第七节）。

接入前的状态是**零接入**：全平台检索 `biocandidate` 零命中，`PLUGIN_SPECS` 只有三个插件。此前的 10-04 融合是单向的——BCR 侧注册了平台的酵母 GEM 并做了对比与试点，平台对此一无所知。

本插件带来的是「一个有完整审计链的对照预测器」，**不是更准的预测器**：BCR 自己的复评结论是对 CatPred、CataPro、UniKP 的配对区间全部跨 0，未建立优越性。

---

## 二、交付物

### BCR 侧（提交 `009db7c`）

| 文件 | 字节 | 作用 |
|---|---:|---|
| `integrations/yeast_metatwin_plugin/biocandidate.py` | 15,358 | 平台侧适配器 |
| `integrations/yeast_metatwin_plugin/benchmark.py` | 10,003 | 清单生成器（镜像 UniKP 同名脚本） |
| `integrations/yeast_metatwin_plugin/catalog_entry.py` | 2,722 | catalog 注册数据（纯数据） |
| `integrations/yeast_metatwin_plugin/install.py` | 9,019 | 安装/卸载器 |
| `integrations/yeast_metatwin_plugin/README.md` | — | 采用说明 |
| `tests/test_plugin_adapter.py` | 15,307 | 12 项契约测试 |

### 平台侧（提交 `bf859d1`）

| 文件 | 作用 |
|---|---|
| `04_prediction_plugins/plugin_runtime/biocandidate.py` | 适配器落点 |
| `04_prediction_plugins/BioCandidateRanker/benchmark.py` | 清单生成器 |
| `04_prediction_plugins/BioCandidateRanker/benchmark_manifest.json` | 基准清单（blocked） |
| `04_prediction_plugins/BioCandidateRanker/readiness_manifest.json` | 就绪清单（blocked） |
| `04_prediction_plugins/BioCandidateRanker/LICENSE` | MIT |
| `04_prediction_plugins/plugin_runtime/catalog.py` | 追加一段 marker 界定的注册块（+33 行） |
| `.gitignore` | 新增三项机器本地状态的忽略规则 |

---

## 三、插件契约符合情况

平台契约由 `plugin_runtime/schema.py` 定义。逐项对照：

| 契约字段 | BCR 提供 | 状态 |
|---|---|---|
| `capability` | `kcat_prediction` | 符合 |
| `sequence` + `substrate_smiles` 输入 | 接受；其余字段用项目默认值 | 符合 |
| `PluginResult.predictions` | `log10_kcat`（log10 尺度）+ `kcat`（线性 s⁻¹） | 符合 |
| `transforms` | 头部截断（>512 残基）与非标准残基替换均申报 | 符合 |
| `uncertainty` | 异方差高斯标准差；区间为未校准 1σ | 符合，`calibrated=false` |
| `applicability` | 无 OOD 模型 → `not_available`、`ood=null` | **能力缺失，如实申报** |
| `benchmark` | 无公开可重建基准 → `not_run` | **如实申报** |
| `PluginSpec.required_assets` | 检查点 + 适配器，均带 SHA256 | 符合 |
| `readiness_manifest` | 已生成，状态 `blocked` | 符合 |

---

## 四、安装与使用

安装（幂等，可重复执行）：

```bash
python integrations/yeast_metatwin_plugin/install.py ^
  --platform "D:/biological/Metabolic model prediction/Integrated_Yeast_MetaTwin_Deployment"
```

预览不落盘：加 `--dry-run`。卸载（`catalog.py` 字节级还原）：加 `--uninstall`。

从平台侧调用：

```python
from plugin_runtime.schema import PluginInput
from plugin_runtime import biocandidate

result = biocandidate.predict(
    PluginInput(
        request_id="req-1", capability="kcat_prediction",
        sequence="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANL",
        substrate_smiles="CC(=O)O",
        context={"organism": "Saccharomyces cerevisiae", "ec": "1.1.1.1"},
    ),
    plugin_root=PLATFORM / "04_prediction_plugins",
)
```

`plugin_root` 遵循平台约定，即 `04_prediction_plugins` 目录（与 UniKP 适配器一致），因此资产路径与 catalog 条目逐字对应。

### 安装与卸载的可回滚性（已实测）

- 重复安装：exit 0，文件内容不变（幂等）
- 卸载后 `catalog.py` **字节级还原**：sha256 一致，5,489 字节
- 对平台的改动仅限文件末尾一段 `# --- BEGIN/END BioCandidateRanker plugin ---` 界定的块，其余部分一字节未动

---

## 五、环境要求

| 依赖 | 说明 |
|---|---|
| torch | 平台主环境钉 `2.13.0`；适配器用 `safe_globals` 白名单解决旧检查点加载问题（见下） |
| fair-esm | 冻结的 ESM-2 骨干；首次运行从 ESM 模型库下载 `esm2_t6_8M_UR50D` 并缓存 |
| rdkit | 底物 SMILES 解析 |
| BCR 代码 | 由 `BIOCANDIDATE_RANKER_ROOT` 环境变量或 `ranker_root.txt` 指向 |

### torch 版本问题的处理

平台 `constraints.txt` 钉死 `torch==2.13.0`，而 BCR 检查点早于 torch 2.6：`biocandidate.training.load_checkpoint` 的 `weights_only=True` 会拒绝加载。`training.py` 属 BCR 的 `software-implementation-v6` 字节级冻结清单，不可就地修改，因此适配器在调用点包了 `torch.serialization.safe_globals` 白名单（只放行可信的 numpy/torch 类型），与 BCR 自己 `scripts/serve_workbench.py` 使用的是同一套。**torch < 2.6 时该包装是空操作。**

隔离环境示例（与平台「plugin environments are intentionally separate」的做法一致）：

```bash
set BIOCANDIDATE_RANKER_ROOT=D:/biological/BioCandidateRanker
set PYTHONPATH=%BIOCANDIDATE_RANKER_ROOT%/src;%PYTHONPATH%
```

---

## 六、验证结果

| 验证项 | 结果 |
|---|---|
| 平台插件注册表 | 列出 `['CLEAN', 'DLKcat', 'UniKP', 'BioCandidateRanker']` |
| 平台插件/运行时测试 | **8/8 通过** |
| 平台全量测试 | **58 通过 / 1 失败**（失败为预存在，见第八节） |
| 端到端 smoke（平台自己的 `FIXED_SMOKE_INPUTS`） | `status=ready`；`log10_kcat=0.4098`、`kcat=2.5694 s⁻¹`、`sd=1.5554` |
| 资产核验（按 `min_size_bytes`） | 检查点 95,441,381 字节、适配器 15,358 字节，均通过 |
| BCR 侧契约测试 | **12/12 通过** |
| BCR 全量测试 | **315 通过 / 1 失败**（失败为预存在的 v6 冻结漂移） |

smoke 结果的不确定性报为 `available` 且 `calibrated=false`（无匹配校准工件）；`applicability.ood=null`；`benchmark=not_run`——三项均与第七节的申报一致。

---

## 七、三个诚实申报

插件清单 `readiness_manifest.json` 的状态是 **blocked**，不是 ready。原因真实且不掩盖：

1. **benchmark gate 阻塞**。BCR 的内部同源冷测试是训练语料的已观察划分，不是独立基准；开发语料不可再分发，第三方无法从公开产物重建该评估；对 CatPred、CataPro、UniKP 的配对家族 bootstrap 区间全部跨 0。清单中 `independent=false`、`metrics_published=false`、`metrics={}`。
2. **无适用域（OOD）能力**。BCR 不提供 OOD 模型，`applicability.ood` 永远为 `null`，标注 `not_implemented`。平台契约允许这样申报，但 UniKP 插件是真做了 OOD 评估的——这个能力差距不应被掩盖。
3. **只有 log10_kcat 有标签**。主检查点声明了四个任务，但 `log10_km`、`log10_activity`、`log10_flux` 的验证样本数为 0，适配器不输出它们，清单的 `untrained_declared_tasks` 明确列出。

### 一个容易踩的坑：上下文

检查点 `use_context=True`，四个上下文字段（organism、EC、酶类型、反应）会被哈希进前向传播。平台契约只强制 `sequence` + `substrate_smiles`，调用方若不填上下文，适配器会使用 BCR CLI 的默认值（空串 / `unknown`），并在 `messages` 与 `applicability.reasons` 中**明确报告**哪些字段是默认值——因为空上下文相对训练语料是分布漂移。建议调用方从 `PluginInput.context` 传入全部四字段。

---

## 八、已知问题与遗留

| 项 | 状态 | 说明 |
|---|---|---|
| 平台预存在测试失败 | 未修 | `tests/test_generic_target_workflow.py::test_reaction_evidence_readiness_verifies_both_versioned_manifests`：两份冻结的反应证据清单记录的 `constraints.txt` 字节数与实际（306 字节）不符。**与插件无关**，插件未触碰该文件 |
| BCR 预存在测试失败 | 未修 | `tests/test_software_release.py`：`src/biocandidate/data/homology.py` 的字节与 v6 冻结清单不符（因未提交的健壮性修复）。正确做法是重冻 v7 |
| 92 MB 检查点 | 已忽略 | 加入平台 `.gitignore`，可由 `install.py` 从 BCR 检出重建。若希望像 UniKP 那样入库跟踪，去掉该规则即可 |
| 两侧提交未推送 | 待确认 | BCR `009db7c` → `origin/master`；平台 `bf859d1` → `origin/main` |
| 校准区间为空 | 已知 | 无与检查点匹配的校准工件，保形区间功能空转 |

---

## 九、验证过程中发现并修复的三个真实缺陷

三个缺陷各自都已补测试覆盖：

1. **非标准残基计数错误**。原实现用「去重残基种类数」当「替换位置数」，导致**每条全标准序列都被误报含非标准残基**。修正为按位置计数。
2. **Windows 行尾破坏**。`Path.write_text` 会把 `\n` 翻译成本地行尾（CRLF），第一版安装器把平台整个 LF 的 `catalog.py` 从 5,489 字节改写成 5,607 字节。改为字节级读写，卸载才能字节级还原。
3. **`ranker_root.txt` 路径约定错位**。适配器曾按 `plugin_runtime` 解释 `plugin_root`，而平台约定是 `04_prediction_plugins`（UniKP 适配器即如此）。检查点路径侥幸正确，指针查找错了一级；**BCR 侧测试因 `sys.path` 泄漏而侥幸通过**，是「平台视角」的验证才暴露。测试已改为从纯净 catalog 起算，并在断言指针解析前移除可导入的 `biocandidate`。

另有一次操作失误已修正：为查清注册流程，我用当前解释器跑了平台的 `audit_phase2_plugin_runtime_compatibility.py`，覆盖了三份已跟踪的审计报告。原报告在平台的专用 UniKP 环境（Python 3.10 + sklearn 1.2）下生成、模型 `load_ok`，我那份显示 `load_failed`。已 `git checkout` 还原，未提交任何降级内容。

---

## 十、资产身份

| 资产 | 字节 | SHA256（前 16 位） |
|---|---:|---|
| 检查点 `esm2-t6-opt-seed7.best.pt` | 95,441,381 | `5a85d9613cceac0d` |
| 适配器 `biocandidate.py` | 15,358 | `eda2930def6da4dc` |
| 清单生成器 `benchmark.py` | 10,003 | `c6186853f42441ec` |
| 安装器 `install.py` | 9,019 | `2223924a795e859a` |

插件版本标识：`software-implementation-v6+esm2-t6-seed7`
检查点来源：`artifacts/esm2-t6-opt-seed7/best.pt`（BCR 仓库）
上游仓库：`https://github.com/liubuing/BioCandidateRanker`
