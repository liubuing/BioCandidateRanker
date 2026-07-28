# Prospective Campaign Collaboration Invitation

Recipient: `[USER MUST SUPPLY A REAL EXPERIMENTAL COLLABORATOR]`

Subject: Invitation to run an independently custodied blind enzyme candidate campaign

Dear `[COLLABORATOR NAME]`,

We would like to collaborate on a prospective, one-time blind evaluation of an enzyme
candidate-ranking model. The campaign must be genuinely unobserved by the modeling team and
must retain all attempted candidates, including inactive and censored outcomes.

The proposed separation is:

- your team supplies a label-free candidate roster and exact model inputs;
- an independent named custodian keeps numeric assay labels inaccessible to us;
- we freeze the model checkpoint and verify family nonoverlap before labels are released;
- we deposit a complete prediction for every frozen candidate;
- the custodian verifies the deposit and releases labels once;
- the final campaign-level evaluation is executed once and cannot be used for tuning.

The default protocol requires at least three campaigns, three distinct eligible families,
and 20 attempted candidates per campaign. Candidate inclusion cannot depend on assay
outcomes. Exact sequences, substrate identities, complete attempt accounting, and family
nonoverlap evidence are required.

The complete handoff package is available under
`artifacts/external/prospective-candidate-ranking/handoff`. It contains the protocol, roster
template, model-input schema, campaign registry, family audit, model freeze, custodian
agreement, and executable commands. Blank templates are not evidence; the actual named
custodian must complete and approve them.

If this design fits an upcoming experimental campaign, please identify the campaign lead,
independent label custodian, candidate families, approximate roster sizes, and expected
timeline before sharing any numeric labels.

Best regards,

`[SENDER NAME]`

`[ORGANIZATION]`

`[CONTACT]`
