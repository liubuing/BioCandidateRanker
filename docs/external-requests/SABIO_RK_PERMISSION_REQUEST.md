# SABIO-RK Permission And Sequence Request

To: SABIO-RK team `<sabiork@h-its.org>`

Institutional escalation: HITS `<info@h-its.org>`

Scientific contact if needed: Dr. Ulrike Wittig `<Ulrike.Wittig@h-its.org>`

Subject: Permission and exact-assay-sequence request for governed Km model training

Dear SABIO-RK team,

We are preparing a governed enzyme Km model using direct experimental records with stable
row provenance. A targeted audit confirmed that SABIO-RK exposes the core fields needed for
many records, including stable kinetic-law and reaction IDs, original linear Km values and
units, normalized values, substrates, PubMed citations, UniProt accessions, and some assay
conditions.

Before using any records, we need written clarification and permission for the following
scope:

- internal model training on an attributed SABIO-RK export;
- storage of normalized row-level derived data with stable SABIO-RK IDs;
- generation and storage of model checkpoints, aggregate metrics, and non-row-level
  scientific reports;
- the intended sharing or publication boundaries for normalized data and checkpoints;
- any required attribution, access controls, or non-commercial-use declarations.

We also need a defensible way to identify the exact final assayed amino-acid chain for each
selected entry. A canonical UniProt accession alone cannot resolve recombinant tags,
truncations, mutations, processing, isoforms, or cleavage scars. Does SABIO-RK expose a
construct-level field or source-backed mapping that proves when the canonical UniProt chain
was assayed unchanged? If not, could you recommend a supported export or workflow for
recovering exact construct evidence from the linked source records?

We will not infer missing units or sequences, will not use model predictions as labels, and
will exclude records that cannot be mapped exactly. We can provide our data contract and
audit receipt on request.

Thank you for clarifying the permitted workflow.

Best regards,

`[SENDER NAME]`

`[ORGANIZATION]`

`[CONTACT]`

## Internal Evidence

- Terms: `https://sabiork.h-its.org/ui/terms`.
- Current audit: `artifacts/external/governed-km-sabio-rk/blocker.json`.
- Reproducible acquisition audit: `scripts/acquire_sabio_rk_km.py`.
