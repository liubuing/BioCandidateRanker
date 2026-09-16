# Manuscript Reference Audit

Date: 2026-09-16. Scope: `manuscript/Bioinformatics_Main_Final_EN.md` reference list
(~40 entries) and in-text citations. The CN manuscript mirrors the EN list and inherits
every correction below. Verified against publisher pages / PubMed / PMLR on 2026-09-16.

## Verdict summary

- **Fabricated or fatally misattributed (must replace): 5** — entries 11, 17, 21, 36, 43
  (UniKP, CATpred, Yang 2023, Hou 2023, and the DLKcat entry).
- **Wrong metadata (authors/year/volume; fix in place): 4** — entries 15 (Hie), 24
  (Kuleshov), 37 (EnzEngDB), 42 (DLKcat in-text form).
- **Content mismatch (citation does not support the claim): 1** — entry 27 (Musgrave).
- **Spot-verified correct: 12**; remaining classics are consistent with known records but
  should get a mechanical DOI check before submission.

## Must-replace entries

### 11. "Kroll et al. 2023 = UniKP" — wrong paper under this title

Cited as: Kroll,A., Engqvist,M.K.M., Heckmann,D. and Lercher,M.J. (2023) UniKP: a unified
kinetic parameter prediction model for enzyme catalysis. Nat. Commun., 14, 8505.

Facts: UniKP is **Yu, H. et al. (2023)** "UniKP: a unified framework for the prediction of
enzyme kinetic parameters", *Nat. Commun.*, **14**, 8505 (doi:10.1038/s41467-023-44113-1).
Kroll's real 2023 Nat. Commun. paper is **TurNuP**: Kroll,A., Ranjan,S., Engqvist,M.K.M.
and Lercher,M.J. (2023) "Turnover number predictions for kinetically uncharacterized
enzymes using machine and deep learning", *Nat. Commun.*, **14**, 8442
(doi:10.1038/s41467-023-41847-8). The dataset citation in Methods 2.2 should therefore be
(DLKcat: Li et al. 2022; UniKP: Yu et al. 2023), optionally + TurNuP where the curated
kcat set is meant.

### 21. "Li et al. 2024 CATpred, Bioinformatics btae102" — fabricated

Cited as: Li,Y., Zhang,L., Wang,H. and Chen,X. (2024) CATpred: a deep learning framework
for enzyme kinetic parameter prediction. Bioinformatics, 40, btae102.

Facts: doi 10.1093/bioinformatics/btae102 belongs to Flexiplex (Cheng et al. 2024, a
demultiplexer, unrelated). The real work is **Boorla,V.S., Hölzer,M., Maranas,C.D. et al.
(2025)** "CatPred: a comprehensive framework for deep learning in vitro enzyme kinetic
parameters", *Nat. Commun.* (doi:10.1038/s41467-025-57215-9; code
github.com/maranasgroup/CatPred). CatPred is a direct competitor (kcat/Km/Ki) and must
also appear in Related Work with a comparison position, not just the reference list.

### 42. "Yu, B. et al. 2024 DLKcat" — wrong authors/year/volume/pages

Cited as: Yu,B., Zhang,Y., Li,J., Wang,Y., Chen,L. and Liu,Z. (2024) DLKcat: a deep
learning model for enzyme kcat prediction. Nat. Catal., 7, 1054–1065.

Facts: DLKcat is **Li,F., Yuan,L., Lu,H. et al. (2022)** "Deep learning-based kcat
prediction enables improved enzyme-constrained model reconstruction", *Nat. Catal.*,
**5**, 662–672 (doi:10.1038/s41929-022-00798-z).

### 36. "Hou, Ji, Shen 2023, Brief. Bioinform. bbac625" — fabricated

Cited as: Hou,J., Ji,Z. and Shen,Y. (2023) Deep learning methods for protein structure
prediction. Brief. Bioinform., 24, bbac625.

Facts: doi 10.1093/bib/bbac625 is scDCCA (single-cell RNA-seq clustering; Wang et al.).
Replace with a real protein-structure/deep-learning review of choice, or drop the
citation — the 30% identity threshold is better anchored by Rost (1999), which is already
present.

### 43. "Yang, K.K., Dallago,C., Frazer,J. and Hie,B.L. (2023) ... 102696" — fabricated

Cited as: Yang,K.K., Dallago,C., Frazer,J. and Hie,B.L. (2023) Protein fitness landscape
prediction: challenges and opportunities. Curr. Opin. Struct. Biol., 83, 102696.

Facts: Curr. Opin. Struct. Biol. 83:102696 is a cryo-EM time-resolved methods paper
(Klebl et al.). The likely intended review is **Kosonocky,C.W., Alamdari,S., Yang,K.K.
and Amini,A.P. (2025)** in *Curr. Opin. Struct. Biol.* — confirm exact title/volume from
the publisher page before inserting.

## Fix-in-place entries

### 15. Hie et al. — year and volume wrong

Cited as 2024, Cell Syst., 15, 274–285. Correct: **2022, Cell Syst., 13(4), 274–285.e9**
(doi:10.1016/j.cels.2022.03.001). Pages were right; authors (Hie, Yang, Kim) right.

### 24. Kuleshov et al. 2018 — author list wrong

Cited as Kuleshov,V., Jiang,C., Li,R., Genovese,T. and Potts,C. Correct: **Kuleshov,V.,
Fenner,N. and Ermon,S.** (ICML 2018, PMLR 80). Title/venue correct.

### 37. EnzEngDB — wrong authors, year, and venue

Cited as Steinkellner,G., Borkowski,O. and Nidetzky,B. (2024), Zenodo. Facts: the
database is real and the project's Zenodo record (10.5281/zenodo.17310823) is genuine,
but the formal publication is the **NAR Database issue** paper by the Mora group
("Enzyme Engineering Database (EnzEngDB)..." *Nucleic Acids Res.*, **54**, D564,
gkaf1142). Cite the NAR paper with its real authors and keep the Zenodo DOI as the
dataset citation.

### 27. Musgrave et al. 2021 — content mismatch

The paper is a fair-evaluation study for unsupervised domain adaptation (CVPR 2021); it
is cited to support multi-seed evaluation best practice. Keep or replace with a
seed-protocol reference (e.g., Bouthillier et al. 2021, "Accounting for variance in
machine learning benchmarks", MLSys) that actually supports the claim.

## Spot-verified correct (no action)

Angelopoulos & Bates 2023; Baltrušaitis et al. 2019; Bar-Even et al. 2011; Beltagy et
al. 2020; Burges et al. 2005; Cho et al. 2014; Elnaggar et al. 2022 (TPAMI 44, 7112–7127);
Gilmer et al. 2017; Guo et al. 2017; Jumper et al. 2021 (Nature 596, 583–589); Kendall &
Gal 2017; Lin et al. 2023 (Science 379, 1123–1130); Lunzer et al. 2005 (Science 310,
1779–1783); Madani et al. 2023 (Nat. Biotechnol. 41, 1099–1106); Rives et al. 2021; Rogers
& Hahn 2010; Rost 1999; Schomburg et al. 2004; Starr & Thornton 2016; Steinegger &
Söding 2017; Vaswani et al. 2017; Vovk et al. 2005; Weinberger et al. 2009; Wu et al.
2018; Yang, Wu & Arnold 2019 (Nat. Methods 16, 687–694); Nix & Weigend 1994; Orth et al.
2010; Hutson 2018; Pineau et al. 2021; Loshchilov & Hutter 2017/2019; Gasteiger et al.
2020; Schütt et al. 2018.

## Numeric errors found during the same pass

- Table 1 and Supplementary Table S2: three-seed Pearson SD is **0.0222** (full-precision
  recompute from `artifacts/esm2-t6-opt-seed{42,123,7}/test_metrics.json`), not 0.0218.
  Fixed in the EN/CN markdown sources on 2026-09-16; the docx/pdf exports must be
  regenerated from them.
- `artifacts/esm2-t12-seed42/best.pt` exists (trained, no test metrics). Limitations
  4.3 says larger variants "were not evaluated" — literally true (no test evaluation),
  but state explicitly that one t12 run was trained without evaluation to preempt a
  reviewer question.

## In-text citation repairs (EN/CN)

- Abstract/Intro/Methods "(Yu et al., 2024; Kroll et al., 2023; Li et al., 2024)" →
  "(Yu et al., 2023; Li et al., 2022; Kroll et al., 2023)" with corrected entries, plus
  "(Boorla et al., 2025)" where recent work is invoked.
- "(Hou et al., 2023)" (homology-aware benchmark practice) → remove or replace per entry
  36 above.
- "(Steinkellner et al., 2024)" → corrected EnzEngDB entry everywhere (Methods 2.6,
  Discussion 4.2).
- "(Hie et al., 2024; ...)" → "(Hie et al., 2022; ...)" everywhere (Intro, 3.5, Discussion).
