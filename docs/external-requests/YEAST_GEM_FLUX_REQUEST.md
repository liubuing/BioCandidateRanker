# yeast-GEM Experimental Flux Metadata Request

To: Eduard Kerkhoven `<eduardk@chalmers.se>`

Public repository route:
`https://github.com/SysBioChalmers/yeast-GEM/issues/new?template=help.md`

Subject: Data dictionary and provenance for yeast-GEM v9.1.0 anaerobic flux table

Dear Dr. Kerkhoven and yeast-GEM maintainers,

We are auditing
`data/physiology/flux_data_anaerobic.tsv` from yeast-GEM release `v9.1.0` for use as
experimental flux labels mapped to a versioned yeast genome-scale model.

The table contains 126 candidate rows and several explicit `r_*` reaction identifiers, but
the file is headerless and we could not locate enough metadata to interpret it without
guessing. Could you provide or point us to the authoritative data dictionary and row-level
source mapping?

We need:

- column names and whether values are absolute or glucose-normalized;
- original units, sign convention, direction, and biomass normalization;
- strain/genotype, replicate, medium/feed, oxygenation, temperature, growth phase, and
  sampling time;
- measurement method, uncertainty, and row-level citation/source-table location;
- definitions for the `Celton2012`, `Jouhten2008`, `Nissen01`-`Nissen04`, and
  `Wasylenko2014` groups;
- mappings for rows without `r_*` IDs and clarification of `r_4783`, which is absent from
  our frozen Yeast-MetaTwin model;
- confirmation of which rows are experimental measurements rather than model-derived
  values.

We will retain the CC BY 4.0 attribution, keep FBA simulations separate from experimental
labels, and exclude unresolved mappings rather than match reaction names heuristically. We
can share our audit receipt and proposed normalized schema.

Thank you for your guidance.

Best regards,

`[SENDER NAME]`

`[ORGANIZATION]`

`[CONTACT]`

## Internal Evidence

- Source file: `https://github.com/SysBioChalmers/yeast-GEM/blob/v9.1.0/data/physiology/flux_data_anaerobic.tsv`.
- Current audit: `artifacts/experimental-yeast-flux/blocker.json`.
- Reproducible auditor: `scripts/audit_experimental_yeast_flux.py`.
