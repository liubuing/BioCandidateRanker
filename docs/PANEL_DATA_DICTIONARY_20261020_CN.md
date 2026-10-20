# 前瞻 kcat Panel 数据字典（草案）

版本：0.1  
日期：2026-10-20  
配套协议：`docs/PANEL_PROTOCOL_20261020_CN.md`  
状态：**草案，未冻结**。字段冻结前不得开始正式标签采集。

## 1. 通用约定

| 项 | 约定 |
|---|---|
| 编码 | UTF-8，逗号分隔 CSV；FASTA 用 ASCII；JSON 用 UTF-8 |
| 表头 | 第一行必须是字段名，不得合并单元格、不得有多行表头 |
| 日期 | ISO 8601，`YYYY-MM-DD`；时间戳用 `YYYY-MM-DDTHH:MM:SS+08:00` |
| 单位 | 写在字段名或独立字段中，不在数值里混单位 |
| 缺失 | 用空值 + `status`/`reason_code` 说明；**禁止推断、插值、按均值填补** |
| 删失 | 用 `censored`，不要写成 0 或检测限 |
| 未检出 | 用 `not_detected`，保留检测限和原始读数 |
| 失败 | 用 `failed`，写入 `panel_attempts.csv`，不能删除 |
| 数值 | 不得使用 `"<0.01"`、`"ND"` 等字符串替代数值；检测限单独字段 |
| 重复 | 技术孔、同批重复、独立蛋白批次必须分开记录，不能混为一条 |
| 标识 | 所有 ID 全局唯一、稳定、不含空格；不得用行号当 ID |
| 哈希 | 每个交付文件在 `panel_provenance.json` 中记录 SHA256 和字节数 |
| 本地路径 | 对外交付版本不得包含绝对本地路径 |
| 版本 | 字段变更必须升版本并写变更记录；不得静默改列含义 |

## 2. 文件清单

| 文件 | 粒度 | 用途 |
|---|---|---|
| `panel_families.csv` | 家族 | 家族选择、同源筛查、备用顺序 |
| `panel_constructs.csv` | 构建体 | 最终测定链、标签、突变、载体 |
| `panel_assays.csv` | 构建体 × 底物 × 条件 | 测定条件、酶浓度、饱和/偶联证据 |
| `panel_records.csv` | 构建体 × 底物（聚合后） | 正式 kcat/Km 记录 |
| `panel_attempts.csv` | 每次尝试 | 失败、失活、删失、未检出、排除 |
| `panel_sequences.fasta` | 最终链 | 精确氨基酸序列 |
| `panel_provenance.json` | 数据集 | 来源、许可、文件哈希、角色 |
| `custodian_labels.csv` | 记录 | 独立保管人持有的标签，建模方在盲存前不可见 |

## 3. `panel_families.csv`

| 字段 | 类型 | 必填 | 说明 / 允许值 |
|---|---|---|---|
| `family_id` | string | 是 | 主键，如 `fam-001` |
| `family_name` | string | 是 | 人类可读名称 |
| `source_accession` | string | 是 | 代表序列的 UniProt/GenBank/PDB accession |
| `family_definition` | string | 是 | 序列簇、结构域或文献定义的家族边界 |
| `mmseqs_target_sha256` | string(64) | 是 | UniKP 开发集 FASTA 的 SHA256 |
| `mmseqs_min_identity` | number | 是 | 固定 0.30 |
| `mmseqs_coverage` | number | 是 | 固定 0.80 |
| `mmseqs_coverage_mode` | integer | 是 | 固定 0 |
| `mmseqs_hit_count` | integer | 是 | 该家族对开发集的命中数；纳入必须为 0 |
| `selection_rank` | integer | 是 | 1–12，label-independent 排序 |
| `backup_order` | integer | 否 | 备用家族顺序；主家族留空 |
| `feasibility_notes` | string | 是 | 表达、纯化、底物、读出可行性 |
| `status` | enum | 是 | `primary` / `backup` / `excluded` |

## 4. `panel_constructs.csv`

| 字段 | 类型 | 必填 | 说明 / 允许值 |
|---|---|---|---|
| `construct_id` | string | 是 | 主键，如 `fam-001-c01` |
| `family_id` | string | 是 | 外键 → `panel_families.family_id` |
| `parent_accession` | string | 是 | 亲本 UniProt/GenBank accession |
| `parent_sequence_id` | string | 是 | FASTA 中的亲本 ID |
| `source_organism` | string | 是 | 亲本来源物种（不是表达宿主） |
| `ec_number` | string | 是 | 如 `3.1.3.16`；多 EC 用分号分隔 |
| `reaction` | string | 是 | 反应式，含底物/产物 |
| `construct_name` | string | 是 | 实验内部名称 |
| `mutations` | string | 是 | `WT` 或 HGVS/位点列表，如 `Q395A;L399A` |
| `tag_n_terminal` | string | 是 | 如 `GST`、`6xHis`、`none` |
| `tag_c_terminal` | string | 是 | 如 `6xHis`、`none` |
| `linker_sequence` | string | 是 | 精确氨基酸序列；无则写 `none` |
| `cleavage_site` | string | 是 | 如 `thrombin LVPRGS`、`none` |
| `processed_termini` | string | 是 | 成熟链的 N/C 端边界；未知写 `unknown` 并说明 |
| `final_chain_id` | string | 是 | 外键 → `panel_sequences.fasta` header |
| `expression_vector` | string | 是 | 如 `pDEST15` |
| `expression_host` | string | 是 | 如 `E. coli BL21(DE3)` |
| `purification_method` | string | 是 | 如 `GST affinity`、`Ni-NTA` |
| `is_tag_cleaved` | enum | 是 | `yes` / `no` / `unknown` |
| `stoichiometry` | string | 是 | 单体/二聚体/复合体及计量 |
| `selection_basis` | string | 是 | 为什么选它；不得写“模型预测高” |
| `status` | enum | 是 | `planned` / `expressed` / `purified` / `failed` |

## 5. `panel_assays.csv`

| 字段 | 类型 | 必填 | 说明 / 允许值 |
|---|---|---|---|
| `assay_id` | string | 是 | 主键，如 `fam-001-c01-s01-a01` |
| `construct_id` | string | 是 | 外键 → `panel_constructs.construct_id` |
| `substrate_id` | string | 是 | 内部底物 ID |
| `substrate_name` | string | 是 | 化学名/常用名 |
| `substrate_pubchem_cid` | integer | 是 | PubChem CID；无 CID 时留空并给结构 |
| `substrate_isomeric_smiles` | string | 是 | 精确异构 SMILES |
| `assay_type` | string | 是 | 如 `continuous fluorescence`、`discontinuous endpoint` |
| `readout` | string | 是 | 如 `DiFMU fluorescence 460 nm`、`NADH A340` |
| `instrument` | string | 是 | 型号/平台 |
| `buffer` | string | 是 | 完整组成，含浓度 |
| `pH` | number | 是 | 实测 pH |
| `temperature_C` | number | 是 | 实测温度 |
| `reaction_volume_uL` | number | 是 | 反应体积 |
| `enzyme_concentration_nM` | number | 是 | 活性酶浓度；不能只给总蛋白 |
| `enzyme_concentration_basis` | enum | 是 | `A280` / `Bradford` / `active_site_titration` / `other` |
| `active_fraction` | number | 否 | 0–1；未知留空并说明 |
| `substrate_range_uM` | string | 是 | 逗号分隔的 8 个非零浓度 + 0 |
| `n_nonzero_points` | integer | 是 | 必须 ≥8 |
| `max_substrate_uM` | number | 是 | 最高非零浓度 |
| `km_fitted_uM` | number | 是 | 用于饱和判断的拟合 Km |
| `max_substrate_multiple_of_km` | number | 是 | 必须 ≥5 |
| `fixed_components` | string | 是 | 辅底物/辅因子/金属/氧气/偶联酶及浓度 |
| `fixed_component_saturation_evidence` | string | 是 | 饱和曲线/文献/全局拟合证据；无则写 `not_established` |
| `coupling_enzyme` | string | 否 | 偶联酶及浓度 |
| `coupling_non_rate_limiting_evidence` | string | 否 | 必需时填写；无则写 `not_applicable` 或 `not_established` |
| `plate_id` | string | 是 | 板号/随机化区块 |
| `assay_date` | date | 是 | 测定日期 |
| `operator` | string | 是 | 操作者 ID（非姓名也可） |
| `qc_status` | enum | 是 | `pass` / `fail` / `review` |
| `qc_flags` | string | 否 | 分号分隔，如 `cv_gt_15;turnover_gt_10pct` |

## 6. `panel_records.csv`

这是与现有 `candidate_records.csv` 对接的主表。一条记录 = 一个构建体 × 一个底物 × 一组匹配条件，跨 3 个独立蛋白批次按预注册规则聚合。

| 字段 | 类型 | 必填 | 说明 / 允许值 |
|---|---|---|---|
| `candidate_id` | string | 是 | 主键，如 `panel-fam001-c01-s01` |
| `family_id` | string | 是 | 外键 |
| `construct_id` | string | 是 | 外键 |
| `assay_id` | string | 是 | 外键 |
| `endpoint` | enum | 是 | `kcat_s-1` / `km_uM` / `kcat_per_km_M-1_s-1` |
| `kcat_s-1` | number | 是 | 正数；不可用 `<`/`ND` |
| `kcat_sd_s-1` | number | 是 | 独立批次 SD；单批次时写 `NA` 并说明 |
| `km_uM` | number | 是 | 正数 |
| `km_sd_uM` | number | 是 | 独立批次 SD |
| `kcat_per_km_M-1_s-1` | number | 是 | 由 kcat/Km 计算，单位 M⁻¹ s⁻¹ |
| `kcat_per_km_sd` | number | 否 | 传播或独立拟合 SD |
| `error_type` | string | 是 | 如 `SD of independent protein batches` |
| `biological_replicates` | integer | 是 | 独立蛋白批次数，目标 3 |
| `technical_replicates` | integer | 是 | 每浓度技术孔数，目标 2 |
| `fit_model` | string | 是 | `Michaelis-Menten` / `substrate-inhibition` / `global` |
| `fit_r2` | number | 否 | 拟合 R² |
| `fit_notes` | string | 是 | 加权、剔除、边界、异常说明 |
| `initial_rate_constraints` | string | 是 | 线性区间、底物消耗 ≤10%、时间点数 |
| `status_at_normalization` | enum | 是 | `accepted_prospective_panel` / `inactive` / `censored` / `failed` / `excluded` |
| `exclusion_reason` | string | 否 | 非 accepted 时必填 |
| `source_row` | string | 是 | 原始文件 + 行/单元格 |
| `provenance_id` | string | 是 | 外键 → `panel_provenance.json` 的文件 ID |

**一致性检查**：`kcat_per_km_M-1_s-1` 应等于 `kcat_s-1 / (km_uM × 1e-6)`，相对误差 ≤1%；超出时必须在 `fit_notes` 说明是独立拟合还是传播导致。

## 7. `panel_attempts.csv`

| 字段 | 类型 | 必填 | 说明 |
|---|---|---|---|
| `attempt_id` | string | 是 | 主键 |
| `family_id` | string | 是 | 外键 |
| `construct_id` | string | 是 | 外键 |
| `assay_id` | string | 否 | 未进入测定时留空 |
| `stage` | enum | 是 | `cloning` / `expression` / `purification` / `assay` / `fit` / `review` |
| `status` | enum | 是 | `failed` / `inactive` / `censored` / `not_detected` / `excluded` |
| `reason_code` | string | 是 | 受控词表，如 `no_expression`、`insoluble`、`no_saturation`、`coupling_rate_limiting` |
| `reason_detail` | string | 是 | 具体说明，不推断 |
| `attempt_date` | date | 是 | 日期 |
| `operator` | string | 是 | 操作者 ID |
| `source_file` | string | 是 | 原始文件/行 |

## 8. `panel_sequences.fasta`

- Header 格式：`>{final_chain_id}|{construct_id}|{family_id}|{parent_accession}|{mutations}`
- 序列必须是**最终测定链**，含标签、linker、切割 scar、突变、加工末端；
- 不允许空格、数字、`*`、`X`（除非有明确说明并单独标记）；
- 每个 header 必须在 `panel_constructs.final_chain_id` 中出现；
- 文件 SHA256 写入 `panel_provenance.json`。

## 9. `panel_provenance.json`

最低字段：

```json
{
  "schema_version": 1,
  "dataset_id": "panel-2026-01",
  "title": "Prospective kcat panel",
  "license": "CC0-1.0 or CC-BY-4.0",
  "license_url": "https://...",
  "institution": "...",
  "contributors": [{"name": "...", "role": "..."}],
  "protocol_version": "0.1",
  "protocol_sha256": "...",
  "retrieval_or_release_date": "2026-10-20",
  "files": [
    {"path": "panel_records.csv", "sha256": "...", "size_bytes": 0, "role": "records"}
  ],
  "data_use_terms": "...",
  "claim_boundary": "Prospective wet-lab panel; not the temporal literature pool."
}
```

## 10. `custodian_labels.csv`

由独立保管人持有，建模方在 `blind-predictions.csv` 存入前不可见。

| 字段 | 类型 | 必填 | 说明 |
|---|---|---|---|
| `candidate_id` | string | 是 | 外键 → `panel_records.candidate_id` |
| `disposition` | enum | 是 | `active` / `inactive` / `censored` |
| `label` | number | 条件 | active 时必填；censored 可空 |
| `label_unit` | string | 条件 | 如 `s-1` |
| `label_sd` | number | 否 | 独立批次 SD |
| `censoring_reason` | string | 条件 | censored 时必填 |
| `notes` | string | 否 | 其他说明 |

## 11. 与现有 `candidate_records.csv` 的字段映射

| 现有字段 | Panel 来源 | 备注 |
|---|---|---|
| `candidate_id` | `panel_records.candidate_id` | 直接 |
| `article_doi` | `panel_provenance.dataset_id` / 论文 DOI | 未发表时用数据集 ID |
| `source_table` | `panel_records.source_row` | 文件 + sheet/行 |
| `source_cell` | `panel_records.source_row` | 可合并 |
| `organism` | `panel_constructs.source_organism` | 亲本来源物种；表达宿主另记 |
| `sequence_id` | `panel_constructs.final_chain_id` | 直接 |
| `accession` | `panel_constructs.parent_accession` | 直接 |
| `construct` | `panel_constructs.construct_name` + 标签/突变 | 人类可读 |
| `fusion_topology` | `panel_constructs.tag_n_terminal`/`tag_c_terminal`/`linker_sequence` | 合并描述 |
| `catalytic_activity` | `panel_records.endpoint` | 统一为 kcat |
| `ec_number` | `panel_constructs.ec_number` | 直接 |
| `reaction` | `panel_constructs.reaction` | 直接 |
| `variable_substrate` | `panel_assays.substrate_name` | 直接 |
| `substrate_pubchem_cid` | `panel_assays.substrate_pubchem_cid` | 直接 |
| `substrate_isomeric_smiles` | `panel_assays.substrate_isomeric_smiles` | 直接 |
| `endpoint` | `panel_records.endpoint` | 直接 |
| `kcat_s-1` | `panel_records.kcat_s-1` | 直接 |
| `kcat_sd_s-1` | `panel_records.kcat_sd_s-1` | 直接 |
| `km_uM` | `panel_records.km_uM` | 直接 |
| `km_sd_uM` | `panel_records.km_sd_uM` | 直接 |
| `error_type` | `panel_records.error_type` | 直接 |
| `biological_replicates` | `panel_records.biological_replicates` | 直接 |
| `fit` | `panel_records.fit_model` + `fit_notes` | 合并 |
| `assay_type` | `panel_assays.assay_type` | 直接 |
| `assay_buffer` | `panel_assays.buffer` | 直接 |
| `assay_pH` | `panel_assays.pH` | 直接 |
| `assay_temperature_C` | `panel_assays.temperature_C` | 直接 |
| `reaction_volume_uL` | `panel_assays.reaction_volume_uL` | 直接 |
| `enzyme_concentration_nM` | `panel_assays.enzyme_concentration_nM` | 直接 |
| `substrate_range_uM` | `panel_assays.substrate_range_uM` | 直接 |
| `maximum_substrate_uM` | `panel_assays.max_substrate_uM` | 直接 |
| `maximum_substrate_multiple_of_km` | `panel_assays.max_substrate_multiple_of_km` | 直接 |
| `saturation_threshold_multiple_of_km` | 固定 5.0 | 协议固定 |
| `saturation_passes` | `panel_assays.max_substrate_multiple_of_km >= 5` | 布尔 |
| `initial_rate_constraints` | `panel_records.initial_rate_constraints` | 直接 |
| `readout` | `panel_assays.readout` | 直接 |
| `status_at_normalization` | `panel_records.status_at_normalization` | 仅当日后按文献池协议审计并并入时，才映射为 `accepted_homology_cold_pool` |

## 12. 校验规则

正式冻结前必须全部通过：

1. 主键唯一；所有外键存在；
2. `family_id` 在 `panel_families` 中唯一，`construct_id`/`assay_id`/`candidate_id` 同理；
3. `panel_records` 的每个 `candidate_id` 在 `custodian_labels` 中恰好一次；
4. `kcat`、`Km` 为正数；`kcat/Km` 一致性 ≤1%；
5. `max_substrate_multiple_of_km >= 5`；`n_nonzero_points >= 8`；
6. `enzyme_concentration_nM > 0`；不能是总蛋白当活性酶；
7. `status_at_normalization` 与 `custodian_labels.disposition` 一致（`accepted_prospective_panel` ↔ `active`，`inactive` ↔ `inactive`，`censored` ↔ `censored`）；
8. 每个 `construct_id` 在 FASTA 中有且仅有一条最终链；
9. `panel_attempts` 覆盖所有未进入 `panel_records` 的构建体/测定；
10. 每个文件有 SHA256；对外版本无绝对本地路径；
11. 技术重复不产生额外 `candidate_id`；
12. 不出现预测值、相对活性、fitness、%WT、RFU/min 作为标签。

## 13. 示例（合成，不是真实测量）

```csv
candidate_id,family_id,construct_id,assay_id,endpoint,kcat_s-1,kcat_sd_s-1,km_uM,km_sd_uM,kcat_per_km_M-1_s-1,kcat_per_km_sd,error_type,biological_replicates,technical_replicates,fit_model,fit_r2,fit_notes,initial_rate_constraints,status_at_normalization,exclusion_reason,source_row,provenance_id
panel-fam001-c01-s01,fam-001,fam-001-c01,fam-001-c01-s01-a01,kcat_s-1,2.50,0.18,24.0,2.1,104166.7,12000,SD of 3 independent protein batches,3,2,Michaelis-Menten,0.998,"weighted; no outlier removal",4 time points; substrate consumption <10%; 8 nonzero concentrations,accepted_prospective_panel,,records.csv:row 2,prov-001
```

示例数值仅用于格式演示，不得当作实验数据。

## 14. 版本与变更

| 版本 | 日期 | 变更 |
|---|---|---|
| 0.1 | 2026-10-20 | 初稿；字段、校验、映射建立 |
