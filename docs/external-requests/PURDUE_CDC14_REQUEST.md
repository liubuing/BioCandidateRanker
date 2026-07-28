# Purdue Cdc14 Data Request

To: Mark C. Hall `<mchall@purdue.edu>`

Repository fallback: Purdue Libraries research-data form at
`https://purdue.libanswers.com/form?queue_id=2673`

Subject: Request for files and construct metadata for dataset 10.4231/429k-qe94

Dear Dr. Hall,

We are building a prospectively governed external evaluation of enzyme kinetics models and
are auditing the dataset associated with DOI `10.4231/429k-qe94` and the JBC article DOI
`10.1016/j.jbc.2024.107644`.

The DataCite record indicates that the Purdue deposit contains 43 files, but the repository
file endpoints currently fail from our environment with TLS unexpected-EOF errors. Could
you provide the deposited files directly or help restore access to the repository file
inventory and downloads?

For scientific admission, we specifically need:

- the original row-level kinetic workbooks and source-cell definitions;
- complete attempted measurements, including inactive, censored, or failed entries;
- exact final assayed human Cdc14 amino-acid chains, including tags, linkers, truncations,
  mutations, cleavage scars, and processed termini;
- unambiguous DiFMUP and phosphopeptide substrate identities;
- original values, uncertainties, units, replicates, pH, temperature, buffer, enzyme
  concentration, substrate series, fixed components, and fitting method;
- confirmation that the deposited CC0 terms cover the supplied measurement files.

We will preserve attribution and source identities, keep the original files unchanged, and
exclude unresolved rows rather than infer missing values. We can share the resulting
machine-readable audit before any confirmatory model evaluation. We will not use the data
for model selection if it enters the frozen external benchmark.

Preferred delivery: original workbook files, CSV/TSV exports where available, exact FASTA,
a short data dictionary, and checksums or a stable repository link.

Thank you for your help.

Best regards,

`[SENDER NAME]`

`[ORGANIZATION]`

`[CONTACT]`

## Internal Evidence

- Contact source: DataCite marks Mark C. Hall as `ContactPerson`.
- Article correspondence source: `https://pmc.ncbi.nlm.nih.gov/articles/PMC11407943/`.
- Current blocker: `artifacts/external/temporal-absolute-kinetics/purdue-429k-qe94/blocker-evidence.json`.
- Acquisition log: `artifacts/external/temporal-absolute-kinetics/purdue-429k-qe94/acquisition-attempts.json`.
