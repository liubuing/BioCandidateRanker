# Experimental Collaboration Specification

## Absolute Kinetics Panel

Design the next panel to fill the frozen deficits directly:

- at least five protein families with no development-corpus MMseqs hit at 30% identity and
  80% coverage, coverage mode 0;
- enough measurements to contribute at least 108 eligible rows after exclusions and the
  20-record-per-family cap;
- defined small-molecule substrates with exact structures;
- complete attempted measurements, including failed, inactive, and censored outcomes in a
  separate campaign roster;
- initial-rate substrate titrations extending to at least five times fitted Km when feasible;
- demonstrated saturation of fixed substrates, cofactors, metals, oxygen, and coupling
  capacity, or reciprocal/global kinetic fits that estimate the saturation limit directly;
- active-enzyme concentration or a documented formula for every derived kcat.

Sequence and non-label metadata may be screened before measurements are exposed. The final
record cap and family selection cannot use kinetic values. After the pool reaches all gates,
freeze the final rows and model checkpoints before one evaluation.

## Activity Panel

Use one physical endpoint and one denominator throughout the collaboration. If
mass-specific activity is selected, report `umol product min^-1 mg protein^-1`, exact
protein-mass basis, calibration, and active fraction. Do not mix fluorescence units,
fitness, percent-of-WT, kcat, or mass-specific activity in one regression target.

## Flux Panel

Prefer chemostat steady states with matched `13C-MFA`, extracellular rates, growth rate,
biomass composition, and uncertainty. Freeze strain and condition splits before training.
Retain the exact model-mapping table as a versioned artifact. Experimental measurements and
FBA simulations must remain separate targets or evidence strata.

## Candidate Ranking Panel

The experimental custodian owns labels until the blind prediction deposit is accepted.
Candidate selection is label-independent, all attempted candidates remain in the roster,
and the final evaluation is run once. EnzEngDB and IMDH are closed and cannot be reused for
selection.
