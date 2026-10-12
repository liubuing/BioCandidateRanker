# BioCandidateRanker plugin for the Yeast-MetaTwin Platform

把 BioCandidateRanker（BCR）的 kcat 预测能力作为插件接入
`Integrated_Yeast_MetaTwin_Deployment` 平台。本目录是**可审查的贡献件**，安装前不改动目标项目任何文件。

## 它给平台带来什么

一个 `kcat_prediction` 能力的插件，与平台已有的 UniKP / DLKcat / CLEAN 并列：

- 输入：`sequence` + `substrate_smiles`（+ 可选四元上下文）
- 输出：`log10_kcat` 与线性 `kcat`（s⁻¹），附异方差高斯标准差
- 模型：BCR 的冻结 ESM-2 多模态检查点（`software-implementation-v6`）

**不带来精度优势。** BCR 自己的复评结论是对 CatPred、CataPro、UniKP 的配对区间全部跨 0，
未建立优越性。接入的价值是「多一个有完整审计链的对照预测器」，不是「更准的预测器」。

## 安装

```bash
# 1. 目标项目侧需要一个能跑 BCR 的环境（见下方"环境"一节）
# 2. 安装（幂等，可重复执行）
python integrations/yeast_metatwin_plugin/install.py \
  --platform "D:/biological/Metabolic model prediction/Integrated_Yeast_MetaTwin_Deployment"

# 预览将写什么（不落盘）
python integrations/yeast_metatwin_plugin/install.py --platform ... --dry-run

# 卸载（catalog.py 字节级还原，插件目录整体删除）
python integrations/yeast_metatwin_plugin/install.py --platform ... --uninstall
```

安装会写入平台侧这些路径：

| 路径 | 内容 |
|---|---|
| `04_prediction_plugins/plugin_runtime/biocandidate.py` | 适配器 |
| `04_prediction_plugins/BioCandidateRanker/benchmark.py` | 清单生成器 |
| `04_prediction_plugins/BioCandidateRanker/benchmark_manifest.json` | 生成 |
| `04_prediction_plugins/BioCandidateRanker/readiness_manifest.json` | 生成 |
| `04_prediction_plugins/BioCandidateRanker/ranker_root.txt` | 指向本仓库 |
| `04_prediction_plugins/BioCandidateRanker/assets/esm2-t6-opt-seed7.best.pt` | 检查点 |
| `04_prediction_plugins/BioCandidateRanker/LICENSE` | MIT |
| `04_prediction_plugins/plugin_runtime/catalog.py` | **仅追加一个 marker 界定的块** |

对 `catalog.py` 的改动是文件末尾一段由 `# --- BEGIN/END BioCandidateRanker plugin ---` 界定的块，
文件其余部分一个字节都不动——这是卸载能字节级还原的原因。

## 环境

平台主环境 `constraints.txt` 钉死 `torch==2.13.0`。BCR 的检查点早于 torch 2.6，
`biocandidate.training.load_checkpoint` 的 `weights_only=True` 会拒绝加载。

`training.py` 属 BCR 的 `software-implementation-v6` 字节级冻结清单，不可就地修改，
因此适配器在调用点包了 `torch.serialization.safe_globals` 白名单（只放行可信的
numpy/torch 类型，与 BCR 自己的 `scripts/serve_workbench.py` 同一套）。
**torch < 2.6 时该包装是空操作。**

另外需要：

- `fair-esm`（冻结的 ESM-2 骨干），首次运行会从 ESM 模型库下载 `esm2_t6_8M_UR50D` 权重并缓存
- `rdkit`（底物 SMILES 解析）
- 一个可导入的 `biocandidate` 包：由 `BIOCANDIDATE_RANKER_ROOT` 环境变量、
  或插件目录旁的 `ranker_root.txt`（安装器已写入）指向

```bash
# 隔离环境示例（推荐，与平台的"plugin environments are intentionally separate"一致）
export BIOCANDIDATE_RANKER_ROOT="D:/biological/BioCandidateRanker"
export PYTHONPATH="$BIOCANDIDATE_RANKER_ROOT/src:$PYTHONPATH"
```

## 诚实申报（这三条必须让平台看到）

插件清单里的状态是 **blocked**，不是 ready。原因是真实的：

1. **benchmark gate 阻塞**：BCR 的内部同源冷测试是训练语料的已观察划分，不是独立基准；
   开发语料不可再分发，第三方无法从公开产物重建该评估；对三个现代基线的配对区间全部跨 0。
   `benchmark_manifest.json` 的 `independent=false`、`metrics_published=false`、`metrics={}`。
2. **无适用域（OOD）能力**：BCR 不提供 OOD 模型，所以 `applicability.ood` 永远报 `null`，
   `assessment_status="not_available"`。平台契约允许这样报，但 UniKP 插件是真做了 OOD 评估的——
   这个差距不该被掩盖。
3. **只有 log10_kcat 有标签**：主检查点声明了 4 个任务，但 `log10_km` / `log10_activity` /
   `log10_flux` 的验证样本数为 0，适配器不会输出它们。

## 一个容易踩的坑：上下文

检查点 `use_context=True`，四个上下文字段（**organism、EC、酶类型、反应**）会被哈希进前向传播。
平台契约只强制 `sequence` + `substrate_smiles`，所以调用方若不填上下文，适配器会使用
BCR CLI 的默认值（空串 / `unknown`），并在 `messages` 和 `applicability.reasons` 里**明确报告**
哪些字段是默认值——因为空上下文相对训练语料是分布漂移。

建议调用方从 `PluginInput.context` 传入：

```python
PluginInput(
    request_id="req-1", capability="kcat_prediction",
    sequence="...", substrate_smiles="...",
    context={"organism": "Saccharomyces cerevisiae", "ec": "1.1.1.1",
             "enzyme_type": "wildtype", "reaction": "..."},
)
```

## 已验证的行为

- 用平台**真实的** `plugin_runtime/schema.py` 构造 `PluginResult` 不抛错（契约合规）
- 三种输入情形正确：全标准残基（无误报）、含 X/Z/U（报替换数）、>512 残基（报头部截断）
- 错误能力与非法残基 → `RunStatus.UNSUPPORTED`，不抛异常
- 安装幂等；卸载后 `catalog.py` 字节级还原（sha256 相同）
- 清单为 LF（避免 Windows 上 `write_text` 把整个文件转成 CRLF）

## 文件

| 文件 | 作用 |
|---|---|
| `biocandidate.py` | 平台侧适配器（落点 `plugin_runtime/`） |
| `benchmark.py` | 生成 benchmark / readiness 清单（镜像 UniKP 的同名脚本） |
| `catalog_entry.py` | catalog 注册数据（纯数据，便于审查） |
| `install.py` | 安装/卸载器 |
| `tests/test_plugin_adapter.py` | BCR 侧契约测试（用平台真实 schema） |
