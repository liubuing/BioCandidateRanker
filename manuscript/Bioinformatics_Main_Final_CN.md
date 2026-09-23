# 预训练蛋白质表征改善同源冷评估下的酶动力学预测

[作者匿名审稿]

## 摘要

**动机：** 联合编码蛋白质序列、底物分子图和反应上下文的多模态深度学习架构已被提出用于酶动力学预测，但其在同源冷评估条件下相对于单模态基线的优势尚未得到充分表征。

**结果：** 本文在冻结的MMseqs2同源冷划分（30%一致性，16,838对）上评估BioCandidateRanker的log10(kcat)预测。从头训练蛋白质编码器时，Morgan指纹MLP几乎匹配多模态模型（RMSE 1.4756 vs. 1.4699）。使用冻结ESM-2后，三个种子的平均RMSE依次为：仅蛋白质1.5025、蛋白质＋底物1.4454、完整模型1.3916。同一划分上重新训练的UniKP达到RMSE 1.4093。这些模态比较是在已观察过的内部测试集上进行的事后诊断。模型向EnzEngDB适应度排序（Spearman 0.071）和IMDH景观（Spearman −0.166）的零样本迁移较差。

**数据与代码：** 项目仓库已备有代码和评估清单；投稿前仍须建立公开归档版本和数据可用性声明，并解决开发语料库的再分发许可。

**通讯作者：** [通讯作者]

**补充信息：** 补充数据可在*Bioinformatics*在线版获取。

---

## 1 引言

酶动力学参数，特别是催化转换数kcat，是代谢工程、酶发现和系统生物学的基础（Bar-Even et al., 2011; Kroll et al., 2023）。kcat的实验测定劳动密集，推动了计算预测方法的发展以优先筛选湿实验验证候选（Schomburg et al., 2004）。

近期方法将kcat预测构建为多模态学习问题，联合编码蛋白质序列、底物分子结构和反应上下文（Li et al., 2022; Kroll et al., 2023; Boorla and Maranas, 2025）。其隐含假设是多模态融合提供超越任何单一模态的互补信号——这一原则在其他领域已得到充分验证（Baltrušaitis et al., 2019），但在酶动力学中尚未在防止序列同源信息泄漏的评估条件下得到严格检验。

在进化尺度序列数据库上预训练的蛋白质语言模型已经变革了蛋白质表征学习（Rives et al., 2021; Lin et al., 2023; Elnaggar et al., 2022）。ESM-2（Lin et al., 2023）和ProtTrans（Elnaggar et al., 2022）产生捕获结构和功能信息的逐残基嵌入，无需显式结构测定（Jumper et al., 2021）。这些表征是否为动力学参数预测提供了超越简单序列特征的可操作信号，仍是一个开放问题（Hie et al., 2022）。

我们识别出一个关键混淆因素：当蛋白质编码器在有限数据（约13K训练对）上从头训练时，蛋白质模态贡献的信号可忽略不计，仅使用底物的Morgan指纹基线（Rogers and Hahn, 2010）即可达到可比性能。多模态架构看似有效，仅因为它没有显著低于基线——而非因为它提供了真正的融合收益。这与分子性质预测文献中的关切相呼应，其中简单基线经常匹配或超越复杂架构（Wu et al., 2018; Yang et al., 2019）。

我们观察到，预训练蛋白质语言模型表征改善了内部测试结果。将从头蛋白质编码器替换为冻结ESM-2嵌入（Lin et al., 2023）将多模态信号从可忽略转变为实质性，在冻结同源冷测试集上三个随机种子均产生一致改进。

现有外部评估显示，这种内部预测增益未能迁移到所测试的候选排序任务。在EnzEngDB工程适应度数据（Long et al., 2026）和Lunzer古老适应性IMDH景观（Lunzer et al., 2005）上的零样本评估均产生与随机排序不可区分或更差的结果。这种端点不可迁移性——绝对kcat预测与相对工程适应度之间——本身是一个具有领域启示的科学性负面结果（Hie et al., 2022; Starr and Thornton, 2016）。

本文贡献为：(1) 在冻结同源冷划分上评估预训练蛋白质表征和模态事后对照；(2) 在同一划分上比较重新训练的UniKP与DLKcat；(3) 表征kcat预测向工程适应度排序迁移不佳的现象。

---

## 2 方法

### 2.1 架构

BioCandidateRanker将三种模态编码到共享d维空间后进行任务查询融合，遵循多模态Transformer架构的一般范式（Vaswani et al., 2017; Baltrušaitis et al., 2019）：

**蛋白质序列。** 我们评估两种蛋白质编码器。*分块Transformer*将氨基酸序列分割为固定大小的块（默认128残基），在每个块内应用局部Transformer编码器（Vaswani et al., 2017），复杂度为O(L × chunk_size)，每块产生一个池化标记。该设计平衡了长序列的感受野与计算可行性，类似于长文档NLP中的分块注意力策略（Beltagy et al., 2020）。*ESM-2编码器*将20种标准氨基酸重映射到ESM-2词表，运行冻结预训练骨干网络（esm2_t6_8M_UR50D，6层，320维；Lin et al., 2023），通过学习线性-SiLU-LayerNorm模块将逐残基表征投影到d维，并应用相同的分块池化策略。超过512残基的序列在ESM-2骨干前截断。

**底物分子图。** 原子从六个分类特征（原子序数、度、形式电荷、芳香性、杂化、氢计数）嵌入，键从类型、共轭和环成员嵌入。具有GRU更新（Cho et al., 2014）和均值聚合的稀疏消息传递神经网络（Gilmer et al., 2017）为每个分子产生单个图级标记。这遵循MPNN（Gilmer et al., 2017）建立的分子图神经网络范式，并由SchNet（Schütt et al., 2018）和DimeNet（Gasteiger et al., 2020）扩展。

**反应上下文。** 生物体、EC号、酶类型和反应字符串被哈希到固定词表（4096桶）并以学习的字段位置嵌入，遵循高基数分类数据的特征哈希策略（Weinberger et al., 2009）。可选FBA上下文向量（通量平衡分析集成的预留接口；Orth et al., 2010）投影到共享维度但在所有报告实验中未激活。

**任务查询融合。** 三个模态标记序列与可学习任务特定查询向量拼接。交叉注意力层（Vaswani et al., 2017）后接Transformer编码器融合信息，逐任务MLP头产生均值和对数方差输出用于异方差高斯预测（Nix and Weigend, 1994; Kendall and Gal, 2017）。

### 2.2 数据

开发语料库来源于UniKP/DLKcat数据集（Li et al., 2022; Kroll et al., 2023）：17,010原始行，过滤非正值和多组分SMILES后接受16,838行。重复酶-底物对（395）在log10尺度上以中位数聚合。数据源仅有语料库级来源；记录级引用不可用。底层动力学数据主要来源于BRENDA（Schomburg et al., 2004），由UniKP项目汇编。

### 2.3 同源冷划分

蛋白质序列使用MMseqs2（Steinegger and Söding, 2017；版本5d152c612b6ad2a56f657b7a02c127eceaea2a75）以30%最小一致性、80%覆盖率、覆盖模式0进行聚类。该阈值遵循同源感知蛋白质基准的既定实践（Rost, 1999）。冻结划分将2,204个簇分配到训练/验证/测试分区：中位数聚合后13,157 / 1,640 / 1,646行。无测试簇与任何训练簇共享超过阈值的同源性。划分清单通过SHA256绑定到源文件身份，确保可重复性（Pineau et al., 2021）。

### 2.4 训练方案

所有模型使用AdamW（Loshchilov and Hutter, 2019），权重衰减0.01，梯度裁剪范数1.0，异方差高斯NLL损失（Nix and Weigend, 1994）支持掩码多任务。ESM-2实验使用批大小4，梯度累积8个微批（有效批32），余弦退火学习率调度（Loshchilov and Hutter, 2017）带2轮线性预热，验证损失早停耐心4。从头基线使用批大小32训练3轮遵循原始方案。三个种子（7, 42, 123）在执行前冻结的方案下运行，遵循多种子评估最佳实践（Pineau et al., 2021; Bouthillier et al., 2021）。

### 2.5 基线

(1) 训练集均值预测器。(2) 氨基酸组成MLP（20维，128→64，40轮）。(3) Morgan半径-2 2048位指纹MLP（Rogers and Hahn, 2010；相同架构）。(4) 两种特征集的岭回归变体（alpha从{0.1, 1, 10, 100}在验证集上选择）。(5) 具有相同编码器和训练预算的晚拼接多模态模型，代表常见的早融合vs晚融合比较（Baltrušaitis et al., 2019）。(6) DLKcat（Li et al., 2022）在GPL-3.0-only下于IMDH景观上评估。

### 2.6 外部评估

**EnzEngDB v1**（Long et al., 2026；Zenodo DOI 10.5281/zenodo.17310823，CC BY 4.0）：462,092实验CSV行，245,945严格接受跨160个活动。冻结同源冷选择移除所有与UniKP有MMseqs2命中的序列，保留6,423行跨51个活动。指标在每个活动内计算并宏平均，遵循酶工程数据的活动感知评估方案（Hie et al., 2022）。

**Lunzer IMDH景观**（Lunzer et al., 2005；Dryad DOI 10.5061/dryad.7nd70，CC0）：512个完整六突变基因型的*E. coli* 3-异丙基苹果酸脱氢酶，具有拟合的ln(Km)和ln(kcat/Km)用于NAD和NADP。所有512个变体具有UniKP同源物（67–69%一致性），使其成为突变敏感性测试而非同源冷验证。该景观已被广泛用作适应度预测方法的基准（Starr and Thornton, 2016）。

### 2.7 可重复性

所有数据源、划分、检查点和评估输出通过SHA256文件身份绑定。如果源文件哈希改变，训练将失败。软件实现以文件级清单冻结。MMseqs2版本、聚类参数和划分分配记录在机器可读清单中。该方法解决了机器学习（Pineau et al., 2021）和计算生物学（Hutson, 2018）中记录的可重复性问题。

---

## 3 结果

### 3.1 从头多模态悖论

在冻结同源冷划分下，具有从头分块Transformer蛋白质编码器的任务查询多模态模型达到三种子平均测试RMSE 1.4699 ± 0.0064，Pearson 0.2471 ± 0.0282。Morgan指纹MLP——一个完全忽略蛋白质序列的单模态模型——达到1.4756 ± 0.0024 RMSE和0.2566 ± 0.0188 Pearson。多模态模型的RMSE优势为0.006 log10(s⁻¹)，远在种子变异的一个标准差之内。

该结果表明，从头蛋白质编码器在约13K对、64隐藏维和3轮训练下，未能提取补充底物图的信号。多模态架构无害，但其融合机制相比仅底物预测未提供可测量收益。这一发现与分子性质预测中的观察一致，其中图神经网络在中等规模数据集上经常无法超越Morgan指纹基线（Wu et al., 2018; Yang et al., 2019）。

### 3.2 预训练表征改善内部预测

将蛋白质编码器替换为冻结ESM-2（esm2_t6_8M_UR50D）表征并将共享维度增至256，结果发生质变（表1）。三种子ESM-2多模态模型达到RMSE 1.3916 ± 0.0183，MAE 1.0708 ± 0.0117，Pearson 0.3810 ± 0.0222。

**表1.** 冻结同源冷划分上的三种子测试结果（1,646条记录）。

| 模型 | RMSE（均值 ± SD） | MAE（均值 ± SD） | Pearson（均值 ± SD） |
|---|---|---|---|
| 训练集均值 | 1.5101 | 1.1795 | 0.0000 |
| 氨基酸组成MLP | 1.5777 ± 0.0121 | 1.2185 ± 0.0104 | 0.0428 ± 0.0064 |
| Morgan MLP | 1.4756 ± 0.0024 | 1.1469 ± 0.0022 | 0.2566 ± 0.0188 |
| 分块Transformer多模态 | 1.4699 ± 0.0064 | 1.1268 ± 0.0063 | 0.2471 ± 0.0282 |
| 晚拼接多模态 | 1.4801 ± 0.0230 | 1.1396 ± 0.0114 | 0.2821 ± 0.0278 |
| ESM-2仅蛋白质 | 1.5025 ± 0.0064 | 1.1656 ± 0.0051 | 0.2174 ± 0.0410 |
| ESM-2蛋白质＋底物 | 1.4454 ± 0.0156 | 1.1200 ± 0.0149 | 0.2791 ± 0.0294 |
| **ESM-2多模态（本文）** | **1.3916 ± 0.0183** | **1.0708 ± 0.0117** | **0.3810 ± 0.0222** |
| UniKP Mode B（重新训练） | 1.4093 ± 0.0094 | 1.0785 ± 0.0041 | 0.3777 ± 0.0222 |
| DLKcat Mode B（重新训练） | 1.7138 ± 0.1415 | 1.3351 ± 0.1251 | 0.1463 ± 0.0764 |

与从头多模态模型相比，RMSE改进为0.078 RMSE（5.3%）和0.134 Pearson（54%）。相比已测试的Morgan指纹单模态基线的改进为0.084 RMSE（5.7%）和0.124 Pearson（48%）。三个种子均显示一致方向；无种子交叉（逐种子结果见补充表S2）。改进幅度与蛋白质功能预测中用预训练表征替换学习嵌入所报告的增益相当（Rives et al., 2021; Hie et al., 2022）。

UniKP与DLKcat在冻结训练集上按已发表超参数各重新训练三个种子；DLKcat的RMSE为1.7138 ± 0.1415、Pearson为0.1463 ± 0.0764。已发表检查点因测试序列与其训练集重叠而未纳入公平比较。ESM-2模态对照是在该内部测试结果已被观察后加入的，因此属于诊断，不能视作融合收益的独立确认。

### 3.3 架构消融

固定ESM-2骨干及训练设置后，在仅蛋白质模型中加入底物输入使平均RMSE降低0.0571；继续加入反应上下文再降低0.0538。三个种子的配对差异方向一致（补充表S3）。这是逐步加入输入的消融，不能单独识别生物机制，也不能证明在独立kcat数据上的泛化。

从头模型上的事后单种子消融（种子42，基线RMSE 1.4638；补充表S1）提示部分组件在该训练预算下可能有贡献：移除反应上下文增加RMSE +0.041，将任务特定查询替换为共享向量增加+0.027，将分块Transformer替换为全局均值池化增加+0.027。这些效应适度但方向一致，但不足以确认预训练模型中的融合增益。

异方差高斯NLL目标（Nix and Weigend, 1994; Kendall and Gal, 2017）相比固定方差MSE在三种子上改进点RMSE 0.032 ± 0.005，确认了学习不确定性对点预测的价值，即使校准仍不完美。

### 3.4 不确定性校准

仅验证集标量校准（平均尺度1.002 ± 0.051）未实质性改进测试NLL或将覆盖率推向标称值。一sigma覆盖率0.66–0.70对标称0.68；二sigma覆盖率0.89–0.95对标称0.95。不确定性输出应解释为相对置信度而非校准概率。此校准偏差与深度学习不确定性估计中的发现一致（Guo et al., 2017; Kuleshov et al., 2018），表明确保可靠区间估计可能需要保形预测方法（Vovk et al., 2005; Angelopoulos and Bates, 2023）。

### 3.5 外部迁移失败

**EnzEngDB零样本排序。** kcat检查点作为活动内适应度排序代理，在51个同源冷活动上达到宏平均Spearman 0.071和成对准确率0.524——仅略高于随机（0.000和0.500）。使用成对logistic损失（Burges et al., 2005）训练的专用活动感知排序器达到三种子平均Spearman −0.026 ± 0.080，置信区间覆盖零。无训练模型超越随机十分位富集。

**IMDH突变景观。** 在512个六突变变体（与训练数据67–69%一致性）上，BioCandidateRanker达到辅因子内Spearman −0.166，DLKcat（Li et al., 2022）达到−0.259。两个预测器均差于随机排序。该景观对进一步模型选择关闭。

这些失败不能归因于编码器质量：ESM-2模型更优的内部Pearson（0.38 vs. 0.25）未能挽救外部排序。瓶颈在于端点语义——kcat预测不意味着工程适应度排序，绝对动力学预测不意味着相对突变效应排序。这一观察与更广泛的发现一致：蛋白质语言模型表征虽然对某些任务强大，但并非普遍迁移到所有下游应用（Hie et al., 2022; Madani et al., 2023）。

---

## 4 讨论

### 4.1 表征质量作为约束瓶颈

核心观察是：预训练蛋白质表征改善了冻结内部划分上的结果。同条件ESM-2模态对照显示，加入底物、再加入反应上下文均降低了内部测试误差。这些事后结果尚不能证明前瞻性融合收益；从头模型与ESM-2模型的比较还同时改变了共享维度。

这对领域具有实践启示。使用从头编码器在中等数据集上关于多模态架构优越性的声明可能在测量噪声而非信号。Morgan MLP平齐结果应作为强制健全性检查：如果指纹基线（Rogers and Hahn, 2010）匹配你的多模态模型，你的蛋白质编码器没有贡献。我们呼应分子机器学习中严格基线的呼吁（Wu et al., 2018; Bouthillier et al., 2021）。

### 4.2 端点不可迁移性差距

外部评估失败约束了kcat预测模型的实际效用。EnzEngDB适应度值（Long et al., 2026）反映活动特定工程目标（表达、稳定性、工艺条件下活性），不可还原为kcat。IMDH景观（Lunzer et al., 2005）测量特定辅因子结合上下文中突变的相对适应度效应，其中排序取决于Km和kcat/Km权衡，仅kcat预测器无法捕获。

这一差距表明候选排序应用需要：(a) 多目标预测（kcat、Km、表达、稳定性）配合任务特定排序头，或(b) 从活动特定训练数据直接预测适应度，如蛋白质适应度景观方法所探索（Starr and Thornton, 2016; Hie et al., 2022; Starr and Thornton, 2016）。以当前公共数据的规模和来源质量，两者均无法实现独立验证。

### 4.3 局限性

ESM-2实验使用最小变体（t6，8M参数，6层），受GPU内存约束（8 GB）。另有一个t12模型完成训练但未做测试评估；更大变体也未评估。512残基截断影响少数长酶。开发语料库缺乏记录级引用，无法审计与现有预测器训练数据的重叠——聚合动力学数据库的已知局限（Schomburg et al., 2004; Kroll et al., 2023）。前瞻性独立基准（192/300条记录，25/30个家族）仍未完成；未对该池生成模型预测。

### 4.4 走向诚实基准

我们在观察结果前冻结了前瞻性外部评估方案；部分内部稳健性和消融分析属于事后诊断。我们记录负面结果而不进行事后挽救尝试，并在评估后关闭外部基准的进一步模型选择。前瞻性时间基准方案要求30个家族300条记录且与训练零重叠，方允许任何预测——当前公共数据无法满足的标准。我们将此差距作为发现而非需要默默绕过的局限报告，遵循计算生物学中负责任基准的建议（Pineau et al., 2021; Hutson, 2018）。

---

## 5 结论

冻结同源冷测试表明，预训练蛋白质表征改善本模型的预测。ESM-2事后模态对照显示，加入底物和反应上下文后内部误差依次降低；完整模型与重新训练的UniKP表现接近。内部预测增益尚未证明工程候选排序效用：在已测试的两个排序端点上迁移较差。提出更广泛的泛化主张前，仍需独立的同终点kcat验证。

---

## 致谢

作者感谢匿名审稿人的建设性反馈。计算资源由[机构匿名]提供。

## 资助

本工作未获得外部资助。

## 利益冲突

无声明。

---

## 参考文献

Angelopoulos,A.N. and Bates,S. (2023) Conformal prediction: a gentle introduction. *Found. Trends Mach. Learn.*, **16**, 494–591.

Baltrušaitis,T., Ahuja,C. and Morency,L.P. (2019) Multimodal machine learning: a survey and taxonomy. *IEEE Trans. Pattern Anal. Mach. Intell.*, **41**, 423–443.

Bar-Even,A., Noor,E., Savir,Y., Liebermeister,W., Davidi,D., Tawfik,D.S. and Milo,R. (2011) The moderately efficient enzyme is evolutionary and physicochemical trend shaping kinetic parameters. *Biochemistry*, **50**, 4402–4410.

Boorla,V.S. and Maranas,C.D. (2025) CatPred: a comprehensive framework for deep learning in vitro enzyme kinetic parameters. *Nat. Commun.* DOI: 10.1038/s41467-025-57215-9.

Beltagy,I., Peters,M.E. and Cohan,A. (2020) Longformer: the long-document transformer. *arXiv*, arXiv:2004.05150.

Burges,C., Shaked,T., Renshaw,E., Lazier,A., Deeds,M., Hamilton,N. and Hullender,G. (2005) Learning to rank using gradient descent. In *Proceedings of the 22nd International Conference on Machine Learning*, pp. 89–96.

Cho,K., van Merriënboer,B., Gulcehre,C., Bahdanau,D., Bougares,F., Schwenk,H. and Bengio,Y. (2014) Learning phrase representations using RNN encoder-decoder for statistical machine translation. In *Proceedings of EMNLP*, pp. 1724–1734.

Elnaggar,A., Heinzinger,M., Dallago,C., Rehawi,G., Wang,Y., Jones,L. et al. (2022) ProtTrans: toward understanding the language of life through self-supervised learning. *IEEE Trans. Pattern Anal. Mach. Intell.*, **44**, 7112–7127.

Gasteiger,J., Groß,J. and Günnemann,S. (2020) Directional message passing for molecular graphs. In *International Conference on Learning Representations*.

Gilmer,J., Schoenholz,S.S., Riley,P.F., Vinyals,O. and Dahl,G.E. (2017) Neural message passing for quantum chemistry. In *Proceedings of the 34th International Conference on Machine Learning*, pp. 1263–1272.

Guo,C., Pleiss,G., Sun,Y. and Weinberger,K.Q. (2017) On calibration of modern neural networks. In *Proceedings of the 34th International Conference on Machine Learning*, pp. 1321–1330.

Hie,B.L., Yang,K.K. and Kim,P.S. (2022) Evolutionary velocity with protein language models predicts evolutionary dynamics of diverse proteins. *Cell Syst.*, **13**, 274–285.e9. DOI: 10.1016/j.cels.2022.03.001.

Hutson,M. (2018) Artificial intelligence faces a replication crisis. *Science*, **359**, 864–865.

Jumper,J., Evans,R., Pritzel,A., Green,T., Figurnov,M., Ronneberger,O. et al. (2021) Highly accurate protein structure prediction with AlphaFold. *Nature*, **596**, 583–589.

Kendall,A. and Gal,Y. (2017) What uncertainties do we need in Bayesian deep learning for computer vision? In *Advances in Neural Information Processing Systems*, pp. 5574–5584.

Kroll,A., Ranjan,S., Engqvist,M.K.M. and Lercher,M.J. (2023) Turnover number predictions for kinetically uncharacterized enzymes using machine and deep learning. *Nat. Commun.*, **14**, 4138. DOI: 10.1038/s41467-023-39840-4.

Kuleshov,V., Fenner,N. and Ermon,S. (2018) Accurate uncertainties for deep learning using calibrated regression. In *Proceedings of the 35th International Conference on Machine Learning*, pp. 2796–2804.

Li,F., Yuan,L., Lu,H. et al. (2022) Deep learning-based kcat prediction enables improved enzyme-constrained model reconstruction. *Nat. Catal.*, **5**, 662–672. DOI: 10.1038/s41929-022-00798-z.

Lin,Z., Akin,H., Rao,R., Hie,B., Zhu,Z., Lu,W. et al. (2023) Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*, **379**, 1123–1130.

Loshchilov,I. and Hutter,F. (2017) SGDR: stochastic gradient descent with warm restarts. In *International Conference on Learning Representations*.

Loshchilov,I. and Hutter,F. (2019) Decoupled weight decay regularization. In *International Conference on Learning Representations*.

Lunzer,M., Miller,G.J., Felsheim,R. and Dean,A.M. (2005) The evolutionary biochemistry of an adaptive landscape. *Science*, **310**, 1779–1783.

Madani,A., Krause,B., Greene,E.R., Subramanian,S., Mohr,B.P., Holton,J.M. et al. (2023) Large language models generate functional protein sequences across diverse families. *Nat. Biotechnol.*, **41**, 1099–1106.

Bouthillier,X., Delaunay,P., Bronzi,M. et al. (2021) Accounting for variance in machine learning benchmarks. In *Proceedings of Machine Learning and Systems*, **3**.

Nix,D.A. and Weigend,A.S. (1994) Estimating the mean and variance of the target probability distribution. In *Proceedings of the IEEE International Conference on Neural Networks*, pp. 55–60.

Orth,J.D., Thiele,I. and Palsson,B.O. (2010) What is flux balance analysis? *Nat. Biotechnol.*, **28**, 245–248.

Pineau,J., Vincent-Lamarre,P., Sinha,K., Larivière,V., Beygelzimer,A., d'Alché-Buc,F. et al. (2021) Improving reproducibility in machine learning research. *J. Mach. Learn. Res.*, **22**, 1–20.

Rives,A., Meier,J., Sercu,T., Goyal,S., Lin,Z., Liu,J. et al. (2021) Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. *Proc. Natl Acad. Sci.*, **118**, e2016239118.

Rogers,D. and Hahn,M. (2010) Extended-connectivity fingerprints. *J. Chem. Inf. Model.*, **50**, 742–754.

Rost,B. (1999) Twilight zone of protein sequence alignments. *Protein Eng.*, **12**, 85–94.

Schomburg,I., Chang,A., Ebeling,C., Gremse,M., Heldt,C., Huhn,G. and Schomburg,D. (2004) BRENDA, the enzyme database: updates and major new developments. *Nucleic Acids Res.*, **32**, D431–D433.

Schütt,K.T., Unke,O.T. and Gastegger,M. (2018) Equivariant message passing for the prediction of tensorial properties and molecular spectra. In *Proceedings of the 38th International Conference on Machine Learning*, pp. 9377–9388.

Starr,T.N. and Thornton,J.W. (2016) Epistasis in protein evolution. *Protein Sci.*, **25**, 1204–1218.

Steinegger,M. and Söding,J. (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat. Biotechnol.*, **35**, 1026–1028.

Long,Y., Abbasinejad,F., Li,F.-Z. et al. (2026) Enzyme Engineering Database (EnzEngDB): a platform for sharing and interpreting sequence–function relationships across protein engineering campaigns. *Nucleic Acids Res.*, **54**, D564–D571. DOI: 10.1093/nar/gkaf1142. Dataset: DOI 10.5281/zenodo.17310823.

Vaswani,A., Shazeer,N., Parmar,N., Uszkoreit,J., Jones,L., Gomez,A.N. et al. (2017) Attention is all you need. In *Advances in Neural Information Processing Systems*, pp. 5998–6008.

Vovk,V., Gammerman,A. and Shafer,G. (2005) *Algorithmic Learning in a Random World*. Springer.

Weinberger,K., Dasgupta,A., Langford,J., Smola,A. and Attenberg,J. (2009) Feature hashing for large scale multitask learning. In *Proceedings of the 26th International Conference on Machine Learning*, pp. 1113–1120.

Wu,Z., Ramsundar,B., Feinberg,E.N., Gomes,J., Geniesse,C., Pappu,A.S. et al. (2018) MoleculeNet: a benchmark for molecular machine learning. *Chem. Sci.*, **9**, 513–530.

Yang,K.K., Wu,Z. and Arnold,F.H. (2019) Machine-learning-guided directed evolution for protein engineering. *Nat. Methods*, **16**, 687–694.

Yu,H. et al. (2023) UniKP: a unified framework for the prediction of enzyme kinetic parameters. *Nat. Commun.*, **14**, 8505. DOI: 10.1038/s41467-023-44113-1.
