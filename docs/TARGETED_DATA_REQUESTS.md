# Targeted Data Requests

The software implementation is closed. New scientific work starts from targeted data
requests or an experimental collaboration, not broad literature discovery. Request only
data that can unblock a frozen gate. Do not request model predictions as labels.

## Shared Delivery Contract

Every delivery should include:

- a stable dataset DOI or repository URL, version, retrieval date, and license covering the
  measurement table;
- original source files plus SHA256 checksums; CSV/TSV/JSON is preferred over screenshots;
- one stable row ID per measurement and a DOI or citation for the primary experiment;
- the exact final assayed amino-acid chain, including tags, linkers, cleavage scars,
  mutations, processed termini, and component stoichiometry;
- substrate names and PubChem CIDs or exact isomeric structures;
- original values, uncertainties, units, replicate IDs, censoring, and aggregation rules;
- pH, temperature, buffer, enzyme concentration, substrate concentration series, fixed
  cosubstrates, cofactors, metals, oxygen conditions, and coupling-enzyme capacities;
- explicit experimental, curated, inferred, predicted, or simulated evidence semantics;
- permission to use the delivered rows for model development and scientific reporting.

Incomplete deliveries are audited once. A bounded correction is attempted; unresolved
fields are recorded and the affected rows are skipped without blocking other work. HTTP
404 endpoints are recorded in `unfinished-acquisitions.json` and are not retried during the
same cycle.

## Absolute Kinetics Expansion

Target: at least 108 additional eligible records and five new global MMseqs families. The
54-substrate threshold is already satisfied.

Priority recipients:

1. Purdue Cdc14 repository custodians for DOI `10.4231/429k-qe94`: request the advertised
   43 files, workbook cell provenance, exact human constructs, DiFMUP/phosphopeptide
   identities, and complete assay conditions.
2. Authors of already audited zero-hit families: request only the missing exact construct,
   saturation, oxygen, or coupling-capacity evidence identified in each source blocker.
3. Experimental collaborators: request a prospective panel spanning at least five unrelated
   families, with complete attempted measurements and at most 20 admitted rows per family.

The temporal protocol at `configs/temporal_absolute_kinetics_protocol.json` is the binding
acceptance specification. Model predictions remain prohibited until all gates pass and the
final record list is frozen.

## Governed Km

The local UniKP pickle is insufficient; see
`artifacts/unikp-km-training-blocker-audit.json`. Request from the Kroll/BRENDA/UniKP data
custodians:

- original positive linear Km value and physical unit for every row;
- BRENDA/SABIO-RK/source-record accession and primary publication DOI;
- exact protein accession and sequence, substrate structure, organism, EC, pH, temperature,
  buffer, and assay notes;
- replicate and geometric-aggregation membership;
- the original train/test assignment and permission terms for each underlying record;
- a separate, unused SABIO-RK or other source-level external set where available.

The delivery will be converted to strict normalized JSONL only after row-level provenance,
units, permission, duplicate policy, and family leakage checks pass.

## Absolute Activity

Do not request generic fitness or relative activity as an absolute endpoint. A usable
delivery must declare one compatible physical dimension, preferably
`umol product min^-1 mg protein^-1`, and include:

- the protein-mass denominator and whether it represents total or active enzyme;
- exact constructs, substrate identities, complete assay conditions, calibration curves,
  instrument settings, and blank corrections;
- source values rather than normalized percentages or campaign percentiles;
- enough unrelated families to create family-aware train, validation, and test partitions.

EnzEngDB fitness remains a campaign-relative ranking endpoint. It cannot be converted to
absolute activity. The current blocker audit is
`artifacts/absolute-activity-training-blocker-audit.json`.

## Experimental Flux

Prefer chemostat or isotope-resolved `13C-MFA` measurements. Each row must include:

- strain/genotype, biological replicate, medium, feed rates, oxygenation, temperature,
  growth phase, and sampling time;
- measured reaction/metabolite flux, original unit, uncertainty, method, and source row;
- reaction direction and a mapping to a named genome-scale model reaction, including the
  model file/version and mapping confidence;
- exchange and intracellular flux semantics, biomass normalization, and sign convention;
- explicit `experimental` evidence type; FBA solutions remain `simulated`.

The target reaction may not appear in its FBA input vector. Conditions, strains, and model
identities must be split without leakage before training.

## Prospective Ranking Collaboration

An independent custodian must provide complete attempted-candidate rosters before labels
are visible to the modeling team. The roster includes active, inactive, and censored
candidates. Required files and lifecycle are defined in
`docs/PROSPECTIVE_RANKING_BENCHMARK.md`.

## Request Template

Subject: Request for row-level enzyme kinetics data and exact assay constructs

> We are preparing a prospectively governed enzyme-model evaluation. We are not asking for
> model predictions or selected positive results. Could you provide the original row-level
> measurement table, exact final assayed amino-acid sequences including tags/cleavage scars,
> substrate concentration series and fixed-component conditions, uncertainties/replicates,
> and the license or permission covering reuse? A preferred delivery is CSV/TSV plus FASTA,
> a data dictionary, and checksums. We will preserve source attribution, exclude unresolved
> rows rather than infer missing values, and share the resulting audit with you before any
> confirmatory evaluation.
