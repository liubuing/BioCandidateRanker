# Author Summary（PLOS Computational Biology 必投稿件）

> 要求：150–200 词，避免行话，面向非本领域读者。以下英文为投稿版本；中文为对照参考，不随稿件提交。

## English（投稿用，186 词）

Enzymes accelerate the chemical reactions of life, and knowing how fast each
enzyme works — its turnover number, kcat — helps biologists design everything
from drugs to industrial fermentation. Because measuring kcat in the
laboratory is slow and expensive, computer models have been built to predict
it from an enzyme's amino-acid sequence and the small molecule it acts on.
Recent models combine several kinds of information at once, and it is widely
assumed that this "multimodal" combination is what makes them work.

We tested that assumption under strict conditions that prevent the model from
recognizing enzymes similar to its training data. Three findings emerged.
First, when models must face unfamiliar enzyme families, adding substrate and
reaction information helps only after a strong protein representation is in
place. Second, a simpler retrained competitor performs as well as our
multimodal model, so no superiority claim survives honest comparison. Third,
predictions that are accurate on paper do not translate into ranking real
enzyme-engineering experiments. We release frozen data splits, evaluation
protocols, and audit tools so that multimodal claims in enzyme kinetics can be
checked rather than assumed.

*Word count: 186.*

## 中文对照（内部参考）

酶加速生命体内的化学反应，知道每种酶的工作速度（转换数 kcat）能帮助生物学家设计从药物到工业发酵的各种方案。由于在实验室测量 kcat 既慢又贵，研究者建立了从氨基酸序列和小分子底物预测它的计算模型。近年模型同时组合多种信息，并普遍假设这种"多模态"组合是其有效的原因。

我们在防止模型认出与训练数据相似酶的严格条件下检验了这一假设，得到三个发现。第一，面对陌生酶家族时，只有先有强蛋白质表征，加入底物和反应信息才有帮助。第二，重新训练的更简单竞争者与我们的多模态模型表现相当，因此任何优越性声明都经不起诚实比较。第三，纸面准确的预测并不能转化为对真实酶工程实验的排序能力。我们发布冻结的数据划分、评估协议与审计工具，使酶动力学中的多模态声明可以被检验而非被默认。

## 写作依据（逐句对应稿件）

| 句子 | 出处 |
|---|---|
| 严格条件 | 主稿 2.3 节同源冷划分（30% 一致性 MMseqs2） |
| 发现一 | S4a 模态诊断：1.5025→1.4454→1.3916 |
| 发现二 | 表 1：UniKP 1.4093 vs 本模型 1.3916，差异小 |
| 发现三 | EnzEngDB Spearman 0.071、IMDH −0.166 |
| 发布承诺 | 冻结协议、归档、verify 工具链 |
