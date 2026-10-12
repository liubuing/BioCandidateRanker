# 投稿元数据待填清单

> 日期：2026-09-25（当日执行批次见文末执行日志）。这份清单是"机器人扫不出来的字段"的唯一权威列表：面板的占位符扫描覆盖稿件正文，本清单覆盖元数据。每填一项在复选框打勾；全部完成后才具备提交条件。

## 1. 稿件正文中的占位符（双盲投稿的处理方式特殊）

双盲审稿要求稿件**不含**作者信息，所以正文占位符不是"填上真实姓名"，而是"改成匿名指向"：

- [x] ~~`manuscript/Bioinformatics_Main_Final_EN.md` 与 `_CN.md` 的 `[repository URL]`~~
      已改为投稿口径："public archive (DOI to be inserted upon acceptance)"，与 Cover
      Letter 口径一致（2026-09-25 执行）
- [x] 两份主稿的 `[corresponding author]` / `[通讯作者]` 行
      已按双盲要求删除；署名信息仅录入投稿系统（2026-09-25 执行）
- [x] 两份主稿的 `[Authors anonymized for review]` / `[作者匿名审稿]`
      保持原样（这就是双盲要求的状态），补充材料已核对一致
- [x] ~~`*Bioinformatics* online` 期刊名残留~~ 已改为期刊中性表述
      "Supporting information is provided with the submission"（2026-09-25 执行）
- [x] 占位符全量扫描：`[corresponding|repository URL|TBD|TODO|待填]` 零残留（2026-09-25）

## 2. 投稿系统元数据（正文不出现，但必须真实；仅能由作者本人填写）

- [ ] 通讯作者姓名、机构、通讯地址、邮箱
- [ ] 全体作者姓名与机构（与致谢一致；"作者和致谢须真实"来自 PLOS 计划）
- [ ] 利益冲突声明（当前口径：无）
- [ ] 资助信息（Funding statement；PLOS 要求即使"无资助"也要声明）
- [ ] Author Summary 最终版（150–200 词，从 `Author_Summary_CN_EN.md` 粘贴并按需微调）

## 3. 数据与代码可用性（PLOS 政策硬性字段）

- [ ] 归档 DOI：把路径清洗版归档 `artifacts/release-archive/software-implementation-v6-redacted.zip` 上传 Zenodo
      → 获得 DOI 后同步更新：稿件 Data Availability、Cover Letter、本清单
      （归档 zip 与校验脚本均已就绪，仅剩上传这一人工步骤）
- [x] 数据可用性声明终稿三件事（对照 `artifacts/license-audit/data-availability-statement.md` 起草版）：
      训练语料不可再分发的表述已写入 EN/CN 主稿摘要 Availability 段；
      EnzEngDB（CC BY 4.0）/IMDH（CC0）引用与署名方式已写入正文数据节；
      冻结划分重建指引（本地获取副本重建脚本）已声明（2026-09-25 执行）
- [x] 代码许可确认：根目录已补 `LICENSE`（MIT + 第三方条款说明，与
      `artifacts/license-audit/development-corpus-license-audit.md` 一致）（2026-09-25 执行）

## 4. 归档内容终检（上传前 5 分钟）

- [x] `python scripts/verify_release_archive.py artifacts/release-archive/software-implementation-v6-redacted.zip` 返回 valid: true（2026-09-25 运行）
- [x] 归档内不含任何绝对本地路径（验证器已内置该检查，valid 即通过）
- [ ] 归档 manifest 的 `built_from_commit` 与投稿时仓库 HEAD 一致（投稿当天核对）

## 5. 前瞻池路线确认（影响 Cover Letter 措辞）

- [x] 路线确认：**C（端点迁移贡献改写）**。理由：(a) 出版台账决策门明确
      "池子仍被阻断 → 准备更窄的基准/负结果文章，不放宽冻结门槛"，路线 B 直接
      违反 fail-closed 铁律；(b) 路线 A（继续采集）需要未发出的外部请求获得
      响应，时间上与 85 天投稿窗口冲突；(c) 现有 Cover Letter "What we did" (c)
      与稿件结论（端点不可迁移性 = 贡献 3、"scientifically informative negative
      result"）已经是按 C 写的，经核对**无需改动措辞**（2026-09-25 确认）

## 6. 导出与格式

- [x] 四份 DOCX（EN/CN 主稿 + EN/CN 补充材料）已从修正后的 .md 源重新导出
      （`scripts/export_manuscript_docx.py`，2026-09-25；旧的 7/30 版已覆盖）
- [ ] 投稿版 PDF：用 Word 从最终 DOCX 导出（本环境无 pandoc/LaTeX；PLOS 初投
      接受 DOCX，PDF 为可选项）
- [ ] 格式备注：主稿摘要当前为结构式（Motivation/Results/Availability），PLOS
      初投实行 format-free 政策，可按现状提交；若编辑部要求再改为非结构化摘要

## 执行日志（2026-09-25）

| 动作 | 文件 |
|---|---|
| 补代码许可 | `LICENSE`（新建，MIT + 第三方条款） |
| 清理 EN 主稿 | `manuscript/Bioinformatics_Main_Final_EN.md`（Availability 重写、Contact 行删除、期刊名中性化） |
| 清理 CN 主稿 | `manuscript/Bioinformatics_Main_Final_CN.md`（同上中文版） |
| 重导出 DOCX ×4 | `manuscript/*.docx`（由 `scripts/export_manuscript_docx.py` 生成，脚本入库） |
| 归档校验 | `artifacts/release-archive/software-implementation-v6-redacted.zip` → valid: true |
| 路线决策 | C，Cover Letter 无需改动 |
