# 投稿元数据待填清单

> 日期：2026-09-25。这份清单是"机器人扫不出来的字段"的唯一权威列表：面板的占位符扫描覆盖稿件正文，本清单覆盖元数据。每填一项在复选框打勾；全部完成后才具备提交条件。

## 1. 稿件正文中的占位符（双盲投稿的处理方式特殊）

双盲审稿要求稿件**不含**作者信息，所以正文占位符不是"填上真实姓名"，而是"改成匿名指向"：

- [ ] `manuscript/Bioinformatics_Main_Final_EN.md` 与 `_CN.md` 的 `[repository URL]`
      → 改为 "will be deposited upon acceptance; DOI to be inserted"（与 Cover
      Letter 口径一致），或投稿前已上传则填匿名化仓库链接
- [ ] 两份主稿的 `[corresponding author]`
      → 双盲版删除该行；署名信息仅录入投稿系统
- [ ] 两份主稿的 `[Authors anonymized for review]` / `[作者匿名审稿]`
      → 保持原样（这就是双盲要求的状态），但确认 supplementary 一致
- [ ] 两份补充材料的对应字段与主稿同步

## 2. 投稿系统元数据（正文不出现，但必须真实）

- [ ] 通讯作者姓名、机构、通讯地址、邮箱
- [ ] 全体作者姓名与机构（与致谢一致；"作者和致谢须真实"来自 PLOS 计划）
- [ ] 利益冲突声明（当前口径：无）
- [ ] 资助信息（Funding statement；PLOS 要求即使"无资助"也要声明）
- [ ] Author Summary 最终版（150–200 词，从 `Author_Summary_CN_EN.md` 粘贴并按需微调）

## 3. 数据与代码可用性（PLOS 政策硬性字段）

- [ ] 归档 DOI：把路径清洗版归档 `artifacts/release-archive/software-implementation-v6-redacted.zip` 上传 Zenodo
      → 获得DOI 后同步更新：稿件 Data Availability、Cover Letter、本清单
- [ ] 数据可用性声明终稿：从 `artifacts/license-audit/data-availability-statement.md`
      起草版核对三件事——训练语料不可再分发的表述、EnzEngDB/IMDH 的引用方式、
      冻结划分重建指引是否存在
- [ ] 代码许可确认：仓库现有许可是否覆盖投稿（若仓库无 LICENSE 文件，这是投稿前必须补的空缺）

## 4. 归档内容终检（上传前 5 分钟）

- [ ] `python scripts/verify_release_archive.py artifacts/release-archive/software-implementation-v6-redacted.zip` 返回 valid: true
- [ ] 归档内不含任何绝对本地路径（验证器已内置该检查，valid 即通过）
- [ ] 归档 manifest 的 `built_from_commit` 与投稿时仓库 HEAD 一致

## 5. 前瞻池路线确认（影响 Cover Letter 措辞）

- [ ] 确认路线：C（端点迁移贡献改写）/ A（继续采集）/ B（放宽门槛并披露）
      → 确认后 Cover Letter 第二段 "What we did" 的 (c) 与稿件结论需按路线微调
