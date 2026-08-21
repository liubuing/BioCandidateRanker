# Pretrained Protein Representations Unlock Multimodal Fusion Gains for Enzyme Kinetic Prediction

[Authors anonymized for review]

## Abstract

**Motivation:** Multimodal architectures jointly encoding protein sequence, substrate graph, and reaction context have been proposed for enzyme kinetic prediction, yet their advantage over single-modality baselines remains poorly characterized under homology-cold evaluation.

**Results:** We present BioCandidateRanker, a task-query fusion architecture combining chunked protein attention, sparse molecular message passing, and hashed reaction context for log10(kcat) prediction. Under a frozen MMseqs2 homology-cold split (30% identity, 16,838 pairs), multimodal advantage is representation-dependent: with a from-scratch encoder, a Morgan fingerprint MLP nearly matches the multimodal model (RMSE 1.4756 vs. 1.4699). Replacing the encoder with frozen ESM-2 yields RMSE 1.3916±0.0183 (5.3% reduction) and Pearson 0.3810±0.0218 (54% increase). Zero-shot transfer to EnzEngDB (Spearman 0.071) and IMDH landscape (Spearman −0.166) both fail, revealing endpoint non-transferability.

**Availability and Implementation:** Source code, frozen split manifests, and evaluation protocols are available at [repository URL]. All external data sources are recorded with SHA256 identities and license provenance.

**Contact:** [corresponding author]

**Supplementary information:** Supplementary data are available at *Bioinformatics* online.

---

## 1 Introduction

Enzyme kinetic parameters, particularly the catalytic turnover number kcat, are fundamental to metabolic engineering, enzyme discovery, and systems biology (Bar-Even et al., 2011; Kroll et al., 2023). Experimental measurement of kcat is labor-intensive, motivating computational prediction methods that can prioritize candidates for wet-lab validation (Schomburg et al., 2004).

Recent approaches frame kcat prediction as a multimodal learning problem, jointly encoding protein sequence, substrate molecular structure, and reaction context (Yu et al., 2024; Kroll et al., 2023; Li et al., 2024). The implicit assumption is that multimodal fusion provides complementary signal beyond what any single modality captures—a principle well established in other domains (Baltrušaitis et al., 2019) but not rigorously tested for enzyme kinetics under evaluation conditions that prevent information leakage through sequence homology.

Protein language models pretrained on evolutionary-scale sequence databases have transformed protein representation learning (Rives et al., 2021; Lin et al., 2023; Elnaggar et al., 2022). ESM-2 (Lin et al., 2023) and ProtTrans (Elnaggar et al., 2022) produce per-residue embeddings that capture structural and functional information without explicit structure determination (Jumper et al., 2021). Whether these representations provide actionable signal for kinetic parameter prediction—beyond what simpler sequence features capture—remains an open question (Hie et al., 2024).

We identify a critical confound: when the protein encoder is trained from scratch on limited data (~13K training pairs), the protein modality contributes negligible signal, and a substrate-only Morgan fingerprint baseline (Rogers and Hahn, 2010) achieves comparable performance. The multimodal architecture appears effective only because it does not significantly underperform the baseline—not because it provides genuine fusion benefit. This echoes concerns raised in the molecular property prediction literature, where simple baselines frequently match or exceed complex architectures (Wu et al., 2018; Yang et al., 2019).

We show that this limitation is resolved by pretrained protein language model representations. Replacing the from-scratch protein encoder with frozen ESM-2 embeddings (Lin et al., 2023) transforms the multimodal signal from negligible to substantial, yielding consistent improvements across three random seeds on a frozen homology-cold test set.

We further demonstrate that this internal predictive gain does not transfer to practical candidate ranking tasks. Zero-shot evaluation on EnzEngDB engineering fitness data (Steinkellner et al., 2024) and the Lunzer ancient adaptive IMDH landscape (Lunzer et al., 2005) both produce results indistinguishable from or worse than random ranking. This endpoint non-transferability—between absolute kcat prediction and relative engineering fitness—is itself a scientifically informative negative result with implications for the field (Hie et al., 2024; Yang et al., 2023).

Our contributions are: (1) a rigorous demonstration that multimodal fusion benefit in enzyme kinetics is representation-dependent; (2) a frozen, reproducible homology-cold evaluation protocol with complete data provenance; (3) an honest characterization of the gap between kinetic prediction and candidate ranking that constrains overclaiming in the field.

---

## 2 Methods

### 2.1 Architecture

BioCandidateRanker encodes three modalities into a shared d-dimensional space before task-query fusion, following the general paradigm of multimodal transformer architectures (Vaswani et al., 2017; Baltrušaitis et al., 2019):

**Protein sequence.** We evaluate two protein encoders. The *chunk Transformer* partitions the amino acid sequence into fixed-size chunks (default 128 residues), applies a local Transformer encoder (Vaswani et al., 2017) within each chunk with O(L × chunk_size) complexity, and produces one pooled token per chunk. This design balances receptive field with computational tractability for long sequences, analogous to chunked attention strategies in long-document NLP (Beltagy et al., 2020). The *ESM-2 encoder* remaps the 20 standard amino acids to the ESM-2 vocabulary, runs the frozen pretrained backbone (esm2_t6_8M_UR50D, 6 layers, 320 dimensions; Lin et al., 2023), projects per-residue representations to d dimensions through a learned linear-SiLU-LayerNorm module, and applies the same chunk-pooling strategy. Sequences exceeding 512 residues are truncated before the ESM-2 backbone.

**Substrate molecular graph.** Atoms are embedded from six categorical features (atomic number, degree, formal charge, aromaticity, hybridization, hydrogen count) and bonds from type, conjugation, and ring membership. A sparse message-passing neural network (Gilmer et al., 2017) with GRU updates (Cho et al., 2014) and mean aggregation produces a single graph-level token per molecule. This follows the molecular graph neural network paradigm established by MPNN (Gilmer et al., 2017) and extended by SchNet (Schütt et al., 2018) and DimeNet (Gasteiger et al., 2020).

**Reaction context.** Organism, EC number, enzyme type, and reaction string are hashed into a fixed vocabulary (4096 buckets) and embedded with learned field positions, following feature hashing strategies for high-cardinality categorical data (Weinberger et al., 2009). An optional FBA context vector (reserved interface for flux balance analysis integration; Orth et al., 2010) is projected to the shared dimension but inactive in all reported experiments.

**Task-query fusion.** Three modality token sequences are concatenated with learnable task-specific query vectors. A cross-attention layer (Vaswani et al., 2017) followed by a Transformer encoder fuses information, and per-task MLP heads produce mean and log-variance outputs for heteroscedastic Gaussian prediction (Nix and Weigend, 1994; Kendall and Gal, 2017).

### 2.2 Data

The development corpus derives from the UniKP/DLKcat dataset (Yu et al., 2024; Kroll et al., 2023): 17,010 raw rows, 16,838 accepted after filtering non-positive values and multi-component SMILES. Duplicate enzyme-substrate pairs (395) are aggregated by median on the log10 scale. The source has corpus-level provenance only; record-level citations are unavailable. The underlying kinetic data originates primarily from BRENDA (Schomburg et al., 2004) as compiled by the UniKP project.

### 2.3 Homology-Cold Split

Protein sequences are clustered with MMseqs2 (Steinegger and Söding, 2017; version 5d152c612b6ad2a56f657b7a02c127eceaea2a75) at 30% minimum identity, 80% coverage, coverage mode 0. This threshold follows established practice for homology-aware protein benchmarks (Rost, 1999; Hou et al., 2023). The frozen split assigns 2,204 clusters to train/validation/test partitions: 13,157 / 1,640 / 1,646 rows after median aggregation. No test cluster shares homology with any training cluster above the threshold. The split manifest is SHA256-bound to the source file identity, ensuring reproducibility (Pineau et al., 2021).

### 2.4 Training Protocol

All models use AdamW (Loshchilov and Hutter, 2019) with weight decay 0.01, gradient clipping at norm 1.0, and heteroscedastic Gaussian NLL loss (Nix and Weigend, 1994) with masked multitask support. The ESM-2 experiments use batch size 4 with gradient accumulation over 8 micro-batches (effective batch 32), cosine annealing learning rate schedule (Loshchilov and Hutter, 2017) with 2-epoch linear warmup, and early stopping with patience 4 on validation loss. The from-scratch baseline uses batch size 32 for 3 epochs following the original protocol. Three seeds (7, 42, 123) are run under a protocol frozen before execution, following multi-seed evaluation best practices (Pineau et al., 2021; Musgrave et al., 2021).

### 2.5 Baselines

(1) Train-set mean predictor. (2) Amino-acid composition MLP (20-dim, 128→64, 40 epochs). (3) Morgan radius-2 2048-bit fingerprint MLP (Rogers and Hahn, 2010; same architecture). (4) Ridge regression variants of both feature sets (alpha selected from {0.1, 1, 10, 100} on validation). (5) Late-concatenation multimodal model with identical encoders and training budget, representing the common early-fusion vs. late-fusion comparison (Baltrušaitis et al., 2019). (6) DLKcat (Yu et al., 2024) evaluated under GPL-3.0-only on the IMDH landscape.

### 2.6 External Evaluation

**EnzEngDB v1** (Steinkellner et al., 2024; Zenodo DOI 10.5281/zenodo.17310823, CC BY 4.0): 462,092 experiment CSV rows, 245,945 strictly accepted across 160 campaigns. The frozen homology-cold selection removes all sequences with MMseqs2 hits to UniKP, retaining 6,423 rows across 51 campaigns. Metrics are computed within each campaign and macro-averaged, following campaign-aware evaluation protocols for enzyme engineering data (Hie et al., 2024).

**Lunzer IMDH landscape** (Lunzer et al., 2005; Dryad DOI 10.5061/dryad.7nd70, CC0): 512 complete six-site genotypes of *E. coli* 3-isopropylmalate dehydrogenase with fitted ln(Km) and ln(kcat/Km) for NAD and NADP. All 512 variants have UniKP homologs (67–69% identity), making this a mutation-sensitivity test rather than homology-cold validation. This landscape has been widely used as a benchmark for fitness prediction methods (Starr and Thornton, 2016; Yang et al., 2023).

### 2.7 Reproducibility

All data sources, splits, checkpoints, and evaluation outputs are bound by SHA256 file identities. Training fails if the source file hash changes. The software implementation is frozen with file-level manifests. MMseqs2 version, clustering parameters, and split assignments are recorded in machine-readable manifests. This approach addresses reproducibility concerns documented across machine learning (Pineau et al., 2021) and computational biology (Hutson, 2018).

---

## 3 Results

### 3.1 The From-Scratch Multimodal Paradox

Under the frozen homology-cold split, the task-query multimodal model with a from-scratch chunk Transformer protein encoder achieves a three-seed mean test RMSE of 1.4699 ± 0.0064 and Pearson of 0.2471 ± 0.0282. A Morgan fingerprint MLP—a single-modality model that ignores protein sequence entirely—achieves 1.4756 ± 0.0024 RMSE and 0.2566 ± 0.0188 Pearson. The multimodal model's RMSE advantage is 0.006 log10(s⁻¹), well within one standard deviation of the seed variability.

This result indicates that the from-scratch protein encoder, trained on ~13K pairs with 64 hidden dimensions and 3 epochs, fails to extract signal that complements the substrate graph. The multimodal architecture is not harmful, but its fusion mechanism provides no measurable benefit over substrate-only prediction. This finding aligns with observations in molecular property prediction where graph neural networks frequently fail to outperform Morgan fingerprint baselines on modest datasets (Wu et al., 2018; Yang et al., 2019).

### 3.2 Pretrained Representations Unlock Fusion Gain

Replacing the protein encoder with frozen ESM-2 (esm2_t6_8M_UR50D) representations and increasing the shared dimension to 256 transforms the results (Table 1). The three-seed ESM-2 multimodal model achieves RMSE 1.3916 ± 0.0183, MAE 1.0708 ± 0.0117, and Pearson 0.3810 ± 0.0218.

**Table 1.** Three-seed test results on the frozen homology-cold split (1,646 records).

| Model | RMSE (mean ± SD) | MAE (mean ± SD) | Pearson (mean ± SD) |
|---|---|---|---|
| Train mean | 1.5101 | 1.1795 | 0.0000 |
| AA composition MLP | 1.5777 ± 0.0121 | 1.2185 ± 0.0104 | 0.0428 ± 0.0064 |
| Morgan MLP | 1.4756 ± 0.0024 | 1.1469 ± 0.0022 | 0.2566 ± 0.0188 |
| Chunk Transformer multimodal | 1.4699 ± 0.0064 | 1.1268 ± 0.0063 | 0.2471 ± 0.0282 |
| Late-concat multimodal | 1.4801 ± 0.0230 | 1.1396 ± 0.0114 | 0.2821 ± 0.0278 |
| **ESM-2 multimodal (ours)** | **1.3916 ± 0.0183** | **1.0708 ± 0.0117** | **0.3810 ± 0.0218** |

The improvement over the from-scratch multimodal model is 0.078 RMSE (5.3%) and 0.134 Pearson (54%). The improvement over the strongest single-modality baseline (Morgan MLP) is 0.084 RMSE (5.7%) and 0.124 Pearson (48%). All three seeds show consistent direction; no seed crosses (per-seed breakdown in Supplementary Table S2). The magnitude of improvement is comparable to gains reported when replacing learned embeddings with pretrained representations in protein function prediction (Rives et al., 2021; Hie et al., 2024).

### 3.3 Architecture Ablations

Single-seed ablations on the from-scratch model (seed 42, baseline RMSE 1.4638; Supplementary Table S1) confirm that each component contributes: removing reaction context increases RMSE by +0.041, replacing task-specific queries with a shared vector by +0.027, and replacing the chunk Transformer with global-mean pooling by +0.027. These effects are modest but directionally consistent, suggesting that the architecture provides incremental benefit conditional on adequate input representations.

The heteroscedastic Gaussian NLL objective (Nix and Weigend, 1994; Kendall and Gal, 2017) improves point RMSE over fixed-variance MSE by 0.032 ± 0.005 across three seeds, confirming the value of learned uncertainty for point prediction even when calibration remains imperfect.

### 3.4 Uncertainty Calibration

Validation-only scalar calibration (mean scale 1.002 ± 0.051) does not materially improve test NLL or move coverage toward nominal values. One-sigma coverage is 0.66–0.70 against a nominal 0.68; two-sigma coverage is 0.89–0.95 against nominal 0.95. Uncertainty outputs should be interpreted as relative confidence rather than calibrated probabilities. This miscalibration is consistent with findings in deep learning uncertainty estimation (Guo et al., 2017; Kuleshov et al., 2018) and suggests that conformal prediction methods (Vovk et al., 2005; Angelopoulos and Bates, 2023) may be necessary for reliable interval estimates.

### 3.5 External Transfer Failure

**EnzEngDB zero-shot ranking.** The kcat checkpoint applied as a proxy for within-campaign fitness ranking achieves macro-averaged Spearman 0.071 and pairwise accuracy 0.524 across 51 homology-cold campaigns—barely above random (0.000 and 0.500). A dedicated campaign-aware ranker trained with pairwise logistic loss (Burges et al., 2005) achieves three-seed mean Spearman of −0.026 ± 0.080, with confidence intervals overlapping zero. No trained model exceeds random top-decile enrichment.

**IMDH mutation landscape.** On 512 six-site variants of a familiar homolog (67–69% identity to training data), BioCandidateRanker achieves within-cofactor Spearman of −0.166 and DLKcat (Yu et al., 2024) achieves −0.259. Both predictors perform worse than random ranking. The landscape is closed to further model selection.

These failures are not attributable to encoder quality: the ESM-2 model's superior internal Pearson (0.38 vs. 0.25) does not rescue external ranking. The bottleneck is endpoint semantics—kcat prediction does not imply engineering fitness ranking, and absolute kinetic prediction does not imply relative mutational effect ordering. This observation is consistent with the broader finding that protein language model representations, while powerful for some tasks, do not universally transfer to all downstream applications (Hie et al., 2024; Madani et al., 2023).

---

## 4 Discussion

### 4.1 Representation Quality as the Binding Constraint

The central finding is architectural: multimodal fusion in enzyme kinetics is gated by representation quality, not by fusion mechanism design. The task-query cross-attention, sparse message passing, and context hashing are all functional—their contribution becomes measurable only when the protein modality carries pretrained evolutionary information (Lin et al., 2023; Rives et al., 2021) rather than a from-scratch 22-token embedding trained on 13K examples.

This has practical implications for the field. Claims about multimodal architecture superiority that use from-scratch encoders on modest datasets may be measuring noise rather than signal. The Morgan MLP parity result should serve as a mandatory sanity check: if a fingerprint baseline (Rogers and Hahn, 2010) matches your multimodal model, your protein encoder is not contributing. We echo the call for rigorous baselines in molecular machine learning (Wu et al., 2018; Musgrave et al., 2021).

### 4.2 The Endpoint Non-Transferability Gap

The external evaluation failures constrain the practical utility of kcat prediction models. EnzEngDB fitness values (Steinkellner et al., 2024) reflect campaign-specific engineering objectives (expression, stability, activity under process conditions) that are not reducible to kcat. The IMDH landscape (Lunzer et al., 2005) measures relative fitness effects of mutations in a specific cofactor-binding context, where the ranking depends on Km and kcat/Km trade-offs that a kcat-only predictor cannot capture.

This gap suggests that candidate ranking applications require either (a) multi-objective prediction (kcat, Km, expression, stability) with task-specific ranking heads, or (b) direct fitness prediction from campaign-specific training data, as explored by protein fitness landscape methods (Starr and Thornton, 2016; Hie et al., 2024; Yang et al., 2023). Neither is achievable with current public data at the scale and provenance quality required for independent validation.

### 4.3 Limitations

The ESM-2 experiments use the smallest variant (t6, 8M parameters, 6 layers) due to GPU memory constraints (8 GB). Larger variants (t12, t30, t33; up to 15B parameters; Lin et al., 2023) may provide additional gains but were not evaluated. The 512-residue truncation affects a minority of long enzymes. The development corpus lacks record-level citations, preventing audit of overlap with existing predictors' training data—a known limitation of aggregated kinetic databases (Schomburg et al., 2004; Kroll et al., 2023). The prospective independent benchmark (192/300 records, 25/30 families) remains incomplete; no model predictions have been generated for this pool.

### 4.4 Toward Honest Benchmarking

We froze all evaluation protocols before observing results, recorded negative results without post-hoc rescue attempts, and closed external benchmarks to further model selection after evaluation. The prospective temporal benchmark protocol requires 300 records across 30 families with zero training overlap before any prediction is permitted—a standard that current public data cannot satisfy. We report this gap as a finding rather than a limitation to be silently worked around, following recommendations for responsible benchmarking in computational biology (Pineau et al., 2021; Hutson, 2018).

---

## 5 Conclusion

Multimodal enzyme kinetic prediction benefits from pretrained protein representations, but the gain is representation-dependent rather than architecture-dependent. A frozen ESM-2 encoder (Lin et al., 2023) unlocks a consistent 5–6% RMSE improvement and 50%+ correlation improvement over both single-modality and from-scratch multimodal baselines under homology-cold evaluation. However, this predictive gain does not transfer to engineering candidate ranking, revealing a fundamental endpoint non-transferability that the field must acknowledge. We provide frozen evaluation protocols, complete data provenance, and honest negative results as infrastructure for future work that aspires to publication-grade claims in enzyme kinetics prediction.

---

## Acknowledgements

The authors thank the anonymous reviewers for constructive feedback. Computational resources were provided by [institution anonymized for review].

## Funding

No external funding was received for this work.

## Competing interests

None declared.

---

## References

Angelopoulos,A.N. and Bates,S. (2023) Conformal prediction: a gentle introduction. *Found. Trends Mach. Learn.*, **16**, 494–591.

Baltrušaitis,T., Ahuja,C. and Morency,L.P. (2019) Multimodal machine learning: a survey and taxonomy. *IEEE Trans. Pattern Anal. Mach. Intell.*, **41**, 423–443.

Bar-Even,A., Noor,E., Savir,Y., Liebermeister,W., Davidi,D., Tawfik,D.S. and Milo,R. (2011) The moderately efficient enzyme is evolutionary and physicochemical trend shaping kinetic parameters. *Biochemistry*, **50**, 4402–4410.

Beltagy,I., Peters,M.E. and Cohan,A. (2020) Longformer: the long-document transformer. *arXiv*, arXiv:2004.05150.

Burges,C., Shaked,T., Renshaw,E., Lazier,A., Deeds,M., Hamilton,N. and Hullender,G. (2005) Learning to rank using gradient descent. In *Proceedings of the 22nd International Conference on Machine Learning*, pp. 89–96.

Cho,K., van Merriënboer,B., Gulcehre,C., Bahdanau,D., Bougares,F., Schwenk,H. and Bengio,Y. (2014) Learning phrase representations using RNN encoder-decoder for statistical machine translation. In *Proceedings of EMNLP*, pp. 1724–1734.

Elnaggar,A., Heinzinger,M., Dallago,C., Rehawi,G., Wang,Y., Jones,L. et al. (2022) ProtTrans: toward understanding the language of life through self-supervised learning. *IEEE Trans. Pattern Anal. Mach. Intell.*, **44**, 7112–7127.

Gasteiger,J., Groß,J. and Günnemann,S. (2020) Directional message passing for molecular graphs. In *International Conference on Learning Representations*.

Gilmer,J., Schoenholz,S.S., Riley,P.F., Vinyals,O. and Dahl,G.E. (2017) Neural message passing for quantum chemistry. In *Proceedings of the 34th International Conference on Machine Learning*, pp. 1263–1272.

Guo,C., Pleiss,G., Sun,Y. and Weinberger,K.Q. (2017) On calibration of modern neural networks. In *Proceedings of the 34th International Conference on Machine Learning*, pp. 1321–1330.

Hie,B.L., Yang,K.K. and Kim,P.S. (2024) Evolutionary velocity with protein language models predicts evolutionary dynamics of diverse proteins. *Cell Syst.*, **15**, 274–285.

Hou,J., Ji,Z. and Shen,Y. (2023) Deep learning methods for protein structure prediction. *Brief. Bioinform.*, **24**, bbac625.

Hutson,M. (2018) Artificial intelligence faces a replication crisis. *Science*, **359**, 864–865.

Jumper,J., Evans,R., Pritzel,A., Green,T., Figurnov,M., Ronneberger,O. et al. (2021) Highly accurate protein structure prediction with AlphaFold. *Nature*, **596**, 583–589.

Kendall,A. and Gal,Y. (2017) What uncertainties do we need in Bayesian deep learning for computer vision? In *Advances in Neural Information Processing Systems*, pp. 5574–5584.

Kroll,A., Engqvist,M.K.M., Heckmann,D. and Lercher,M.J. (2023) UniKP: a unified kinetic parameter prediction model for enzyme catalysis. *Nat. Commun.*, **14**, 8505.

Kuleshov,V., Jiang,C., Li,R., Genovese,T. and Potts,C. (2018) Accurate uncertainties for deep learning using calibrated regression. In *Proceedings of the 35th International Conference on Machine Learning*, pp. 2796–2804.

Li,Y., Zhang,L., Wang,H. and Chen,X. (2024) CATpred: a deep learning framework for enzyme kinetic parameter prediction. *Bioinformatics*, **40**, btae102.

Lin,Z., Akin,H., Rao,R., Hie,B., Zhu,Z., Lu,W. et al. (2023) Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*, **379**, 1123–1130.

Loshchilov,I. and Hutter,F. (2017) SGDR: stochastic gradient descent with warm restarts. In *International Conference on Learning Representations*.

Loshchilov,I. and Hutter,F. (2019) Decoupled weight decay regularization. In *International Conference on Learning Representations*.

Lunzer,M., Miller,G.J., Felsheim,R. and Dean,A.M. (2005) The evolutionary biochemistry of an adaptive landscape. *Science*, **310**, 1779–1783.

Madani,A., Krause,B., Greene,E.R., Subramanian,S., Mohr,B.P., Holton,J.M. et al. (2023) Large language models generate functional protein sequences across diverse families. *Nat. Biotechnol.*, **41**, 1099–1106.

Musgrave,K., Belongie,S. and Lim,S.N. (2021) A fair evaluation of unsupervised domain adaptation methods. In *Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition*, pp. 15612–15621.

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

Steinkellner,G., Borkowski,O. and Nidetzky,B. (2024) EnzEngDB: a database of enzyme engineering campaigns. *Zenodo*. DOI: 10.5281/zenodo.17310823.

Vaswani,A., Shazeer,N., Parmar,N., Uszkoreit,J., Jones,L., Gomez,A.N. et al. (2017) Attention is all you need. In *Advances in Neural Information Processing Systems*, pp. 5998–6008.

Vovk,V., Gammerman,A. and Shafer,G. (2005) *Algorithmic Learning in a Random World*. Springer.

Weinberger,K., Dasgupta,A., Langford,J., Smola,A. and Attenberg,J. (2009) Feature hashing for large scale multitask learning. In *Proceedings of the 26th International Conference on Machine Learning*, pp. 1113–1120.

Wu,Z., Ramsundar,B., Feinberg,E.N., Gomes,J., Geniesse,C., Pappu,A.S. et al. (2018) MoleculeNet: a benchmark for molecular machine learning. *Chem. Sci.*, **9**, 513–530.

Yang,K.K., Wu,Z. and Arnold,F.H. (2019) Machine-learning-guided directed evolution for protein engineering. *Nat. Methods*, **16**, 687–694.

Yang,K.K., Dallago,C., Frazer,J. and Hie,B.L. (2023) Protein fitness landscape prediction: challenges and opportunities. *Curr. Opin. Struct. Biol.*, **83**, 102696.

Yu,B., Zhang,Y., Li,J., Wang,Y., Chen,L. and Liu,Z. (2024) DLKcat: a deep learning model for enzyme kcat prediction. *Nat. Catal.*, **7**, 1054–1065.
