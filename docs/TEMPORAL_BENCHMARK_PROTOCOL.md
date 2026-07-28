# Temporal Absolute-Kinetics Benchmark Protocol

Frozen: 2026-07-20

## Objective

Build a prospective external benchmark for absolute enzyme kinetics that is temporally
newer than the UniKP development snapshot, legally reusable, record-level traceable, and
protein-homology-cold relative to every UniKP development sequence.

## Source Gates

- The primary measurement article must be published on or after 2023-01-01.
- The measurement table must be explicitly covered by CC0 1.0 or CC BY 4.0.
- A later repository deposit does not make an older experiment temporally eligible.
- BRENDA, SABIO-RK, or other third-party compilations are excluded unless source-level
  redistribution rights and record-level provenance are both established.
- Preprints require a stable DOI and an explicit dataset license.

## Record Gates

Each record must contain direct experimental `kcat`, `Km`, or `kcat/Km`; an exact protein
sequence or an unambiguous accession plus variant; a single defined substrate that can be
mapped to a structure identifier; units; a measurement-article DOI; dataset DOI or stable
record URL; and source file plus row/table-cell provenance. Relative activity, fitness,
predictions, pseudo-labels, and bulk environmental apparent kinetics are excluded.

Derived `kcat` is accepted only when enzyme concentration is reported and the conversion
from `Vmax` is explicit and reproducible. Conditions such as pH and temperature are
retained when reported but are not silently imputed.

## Leakage Control

Temporal eligibility does not establish independence. Every sequence is searched against
the full UniKP development corpus with MMseqs2 at 30% identity, 80% coverage, coverage
mode 0. Any hit excludes the sequence from the confirmatory pool. Exact
sequence-substrate overlap must be zero. Exact substrate and Murcko-scaffold overlap are
reported as separate applicability strata.

UniKP lacks record-level citations, so publication overlap can be excluded only when
known and otherwise remains explicitly unavailable.

## Readiness Gate

The benchmark is not frozen for model evaluation until it contains at least:

- 300 accepted absolute-kinetics records;
- 30 MMseqs2-homology-cold protein families;
- 50 unique substrates;
- no more than 20 records from one family.

Selection and family caps are label-independent. A family needs at least three records
for a family-level metric. Failure of any gate means no model evaluation.

## No-Peeking Policy

Candidate metadata, licenses, sequence mappings, substrate mappings, and homology may be
audited during curation. Model predictions are prohibited until the final record list,
primary checkpoint, comparators, and metrics are frozen. The final benchmark receives one
evaluation only and cannot subsequently be used for architecture or hyperparameter
selection.
