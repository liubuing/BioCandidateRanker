# Cover Letter — PLOS Computational Biology（草稿，投稿前替换占位符）

> 收件：PLOS Computational Biology 编辑部（Senior Editors）
> 栏目：Research Article, Biological Macromolecules
> 状态：草稿 v1（2026-09-25）。占位符见文末清单，与 SUBMISSION_METADATA_CHECKLIST.md 一致。

---

Dear Editors,

We are pleased to submit our manuscript, "Pretrained Protein Representations
Improve Enzyme Kinetic Prediction under Homology-Cold Evaluation," for
consideration as a Research Article in PLOS Computational Biology.

**The biological question.** Computational predictions of enzyme turnover
number (kcat) increasingly drive which enzyme variants get measured in the
laboratory. A generation of multimodal architectures assumes that jointly
encoding protein sequence, substrate structure, and reaction context is what
makes these predictions work. Whether that assumption survives when a model
must face enzyme families unlike its training data is both a practical
question for enzyme engineering and a testable hypothesis about where
biomolecular signal actually lives.

**What we did.** Under a frozen, homology-cold evaluation (MMseqs2 at 30%
identity, no test cluster sharing detectable homology with training data), we
(a) disentangled the contribution of each input modality with per-seed
controls, (b) retrained the strongest published competitors on the identical
split, and (c) asked whether paper-accurate kcat predictions translate into
ranking real enzyme-engineering campaigns (EnzEngDB, 51 homology-cold
campaigns) and a classical mutation landscape (512-genotype IMDH).

**What we found.** (1) Multimodal gains depend almost entirely on the protein
representation: adding substrate and reaction context helps only once a
pretrained protein encoder is in place. (2) After retraining competitors on
the same split, our model does not establish superiority — the honest point
estimate is nearly identical to UniKP's. (3) Internal predictive accuracy
fails to transfer to experiment-ranking endpoints, with near-random Spearman
correlations on both external tasks.

**Why this belongs in PLOS Computational Biology.** The field publishes many
multimodal architectures and few controlled tests of what the multimodality
contributes. Our results replace an assumed mechanism with a measured one,
bound the transferability of kcat models to real engineering decisions, and
quantify how much of a published gain survives an honest retraining
comparison. We believe negative and boundary-mapping results of this kind are
exactly what rigorous methodology sections of computational biology should
contain, and they align with the journal's emphasis on methods that advance
biological understanding rather than leaderboard position.

**Reproducibility.** Every number in the paper derives from frozen protocols
and SHA256-identified inputs. A deterministic release archive (code, frozen
split manifests, evaluation protocols, and a verifier that recomputes every
hash) will be deposited with a DOI upon acceptance; per-seed results,
paired-family bootstrap intervals, and the audit scripts are included. The
development corpus itself cannot be redistributed pending upstream permission,
and the data availability statement states this explicitly together with
reconstruction instructions.

This manuscript is not under consideration elsewhere and all authors have
approved the submission. The authors declare no competing interests.

Thank you for your consideration.

Sincerely,

`[CORRESPONDING AUTHOR NAME]`
on behalf of all authors

`[INSTITUTION]` — `[EMAIL]`

---

## 提交前必须替换的占位符

| 占位符 | 位置 | 内容 |
|---|---|---|
| `[CORRESPONDING AUTHOR NAME]` | 署名 | 真实通讯作者姓名 |
| `[INSTITUTION]` | 署名 | 所属机构 |
| `[EMAIL]` | 署名 | 可公开联系邮箱 |
| 归档 DOI | "Reproducibility" 段 | Zenodo DOI（路径清洗版归档上传后获得） |

## 语言校准说明

- 全文避免"novel / superior / unlock"类词；主张限定为"replace an assumed
  mechanism with a measured one"——与稿件标题和摘要的克制口径一致
- "will be deposited with a DOI upon acceptance"：PLOS 接受投稿时归档即可，
  无需在投稿前公开 DOI；若选择投稿前上传，把句子改为 "is deposited at
  [DOI]" 并同步更新稿件 Data Availability
- 若编辑部要求建议审稿人，另行准备 3–5 名无利益冲突人选（方法学方向：
  同源划分/蛋白质表征/酶动力学数据集各一）
