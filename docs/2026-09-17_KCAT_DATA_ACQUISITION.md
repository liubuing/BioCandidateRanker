# Same-endpoint kcat acquisition plan (2026-09-17)

## Current gate

The frozen temporal absolute-kinetics protocol requires at least 300 eligible records, 30 global MMseqs families, and 50 unique substrates. The curated pool currently has 192 records, 25 families, and 54 substrates. The unresolved minimum is **108 records and five new global families**. These are admission gaps, not a target for indiscriminate row collection. The current pool remains unscored.

## Request queue

| Priority | Source and recipient | Exact missing evidence | Upper contribution before final checks | Draft | Status |
|---|---|---|---|---|---|
| 1 | Purdue Cdc14 deposit `10.4231/429k-qe94`; Mark C. Hall | Accessible original workbooks, source cells, exact constructs, substrate identities, assay conditions | Unknown; likely one family, subject to global clustering and row caps | [Purdue request](external-requests/PURDUE_CDC14_REQUEST.md) | Draft, unsent |
| 2 | Multicopper oxidase `10.1021/acs.biochem.6c00183`; Tina R. Tuveng | O2 saturation in the same steady-state donor assays | 27 reported evidence rows, **at most 20 admitted per family** | [Oxygen request](external-requests/MULTICOPPER_OXYGEN_REQUEST.md) | Draft, unsent |
| 3 | Phage enzymes `10.1002/advs.202512937`; corresponding authors | Exact final assayed His6-MBP/HMT and His6-TF fusion sequences and tag-cleavage status | **At most four** of eleven reported rows pass the frozen substrate-range gate | [Fusion request](external-requests/PHAGE_FUSION_CONSTRUCT_REQUEST.md) | Draft, unsent |
| 4 | VaPAL2 `10.1016/j.jbc.2026.111301`; corresponding authors | Active-core termini and exact GN10-SL structure | Unquantified pending row-level review; no admission assumed | [VaPAL2 request](external-requests/VAPAL2_ACTIVE_CORE_REQUEST.md) | Draft, unsent |
| 5 | Independent assay partner, recipient to identify | Prospective, complete multi-family kcat panel under frozen measurement contract | Recruitment target: at least five unrelated families and at least 20 attempts per family; accepted yield unknown | [Panel invitation](external-requests/ABSOLUTE_KCAT_PANEL_INVITATION.md) | Draft, unsent |

The first four sources cannot be assumed to close the gap: Purdue and VaPAL2 yields are unknown, and the known numerical maximum from the multicopper and phage sources is 24 before all remaining gates. Even if each source adds a different family, the four named sources offer at most four candidate families. Global clustering could reduce that count. The independent panel is therefore a necessary acquisition track, not a contingency deferred until all author replies arrive.

## Execution order

1. Verify the recipient and sender fields in each draft. External messages require explicit authorization; no messages have been sent. Do not represent draft status as a submitted request.
2. Send the four source-specific requests after authorization. Record send date, sender, recipient, source DOI, message ID, and response location in the acquisition ledger. Use a single respectful follow-up only if warranted; do not repeatedly poll inaccessible repository endpoints.
3. Identify an independent assay partner and agree on a candidate roster and assay conditions before label disclosure. The 5 × 20 attempt design yields only 100 attempted records and cannot by itself guarantee 108 *eligible* additions; plan spare families or rows based on preflight attrition, while preserving the 20-per-family admission cap.
4. For each delivery, archive original files and checksums; preserve source-cell and row identities; audit license, publication date, direct measurement, exact construct, substrate structure, assay conditions, saturation, overlap, and global MMseqs family assignment. Record failures without inferring missing evidence.
5. Recompute eligible counts and families. Freeze the final record list and model checkpoints only if every gate passes. Run a single final evaluation. If the gate remains blocked, retain the narrower manuscript claim.

## Decision boundary

The public articles and local blocker audits support *requests* for missing evidence, not provisional acceptance of their measurements. Author replies may clarify, confirm a failure, or supply new primary data. None of those outcomes changes the frozen protocol. The independent custodian must keep labels hidden until record selection and model checkpoints are locked.

Binding specification: `configs/temporal_absolute_kinetics_protocol.json`. Per-source evidence: `artifacts/external/temporal-absolute-kinetics/<source>/blocker-evidence.json`.
