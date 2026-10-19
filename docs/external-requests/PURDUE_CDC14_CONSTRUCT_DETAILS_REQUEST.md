# Purdue Cdc14 follow-up: exact fusion constructs and enzyme concentration

To: Mark C. Hall `<mchall@purdue.edu>`

Subject: Follow-up on Cdc14 dataset 10.4231/429k-qe94: exact GST fusion sequences and assay enzyme concentration

Dear Dr. Hall,

Thank you for depositing the Cdc14 dataset under CC0. We have now downloaded the full
43-file bundle from PURR and are auditing the human Cdc14A/Cdc14B DiFMUP steady-state
kinetics for a frozen, homology-cold external evaluation.

The deposited Table S1 and the `Fig5E_S3D_hCdc14_DiFMUP_Normalization.xlsx` workbook give
us the fitted kcat and KM values, and the Methods describe the assay buffer and temperature.
Four details are still needed before any row can be admitted under our prespecified
protocol; they are not inferable from the deposited files without guessing:

1. **Exact N-terminal fusion sequences.** The Methods state that hCdc14A (residues 1-413)
   and hCdc14B (residues 1-411) were expressed from pDEST15 with an N-terminal GST tag
   (plasmids pHLP674, pHLP730, pHLP731, pHLP678, pHLP732, pHLP733). Could you provide the
   exact amino-acid sequence of the expressed fusion for each plasmid, including the
   GST-linker-attB1 junction and any vector-derived residues, or the plasmid maps/GenBank
   files? Please also state whether the GST tag was cleaved before the DiFMUP assays.
2. **hCdc14B 3A positions.** Table S1 names hCdc14B Q390A and hCdc14B 3A, but the local
   article text only spells out the hCdc14A 3A substitutions (L399A/R400A/K403A). Could you
   confirm the three hCdc14B 3A positions?
3. **Per-assay enzyme concentration.** The `Fig5E...Normalization.xlsx` rates are already in
   s-1, but the workbook does not record the enzyme concentration used for each
   fluorescence-to-rate conversion. Could you provide the active enzyme concentration (nM)
   for each hCdc14A/B construct/assay, or the Prism analysis parameter that was used?
4. **Attempted-measurement completeness.** If any hCdc14A/B DiFMUP concentrations, trials,
   or constructs were measured but excluded from Table S1, please indicate them so the
   roster is complete rather than selected on outcome.

We will preserve attribution, keep your files unchanged, exclude unresolved rows rather
than infer missing values, and share the machine-readable audit before any confirmatory
model evaluation. These rows would not be used for model selection.

Best regards,

`[SENDER NAME]`

`[ORGANIZATION]`

`[CONTACT]`

## Internal audit

- Dataset: `https://doi.org/10.4231/429k-qe94` (CC0 1.0), 43 files, 50 MB bundle.
- Local acquisition receipt:
  `artifacts/external/temporal-absolute-kinetics/purdue-429k-qe94/acquisition-receipt-20261019.json`.
- Local staged candidates:
  `artifacts/external/temporal-absolute-kinetics/purdue-429k-qe94/candidate-kcat-km-20261019.csv`.
- MMseqs2 zero-hit evidence:
  `artifacts/external/temporal-absolute-kinetics/purdue-429k-qe94/cdc14-vs-unikp.tsv`.
- The follow-up is limited to construct identity, mutation positions, enzyme concentration,
  and roster completeness; it does not request predictions or selected positive results.
