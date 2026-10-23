"""Review 26 metadata groups / 106 rows; preserve the historical pool and gates.

Reads the frozen pool-quality register (``pool-quality-register-20261005``),
re-verifies every recorded evidence hash and the protocol hash, groups the 106
pending-campaign rows into 26 metadata groups, and applies the hand-recorded
per-source dispositions (hold_* or conditions_supported_candidate) with evidence
locators and unresolved gaps. Writes ``review.json`` and ``record-overlay.json``
under ``artifacts/external/ranking-campaign-review-20261006``.

Governance role: overlay-only review -- no numeric response or model score
enters the dispositions, ``final_evaluation_ready_groups`` and
``accepted_additions`` stay 0, absolute record statuses and labels are
unchanged, and no predictions are generated.
"""
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/external/temporal-absolute-kinetics'
PREVIOUS = ROOT / 'artifacts/external/pool-quality-register-20261005'
OUT = ROOT / 'artifacts/external/ranking-campaign-review-20261006'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def main():
    """Verify the frozen register, apply per-source dispositions, and write the overlay."""
    previous = read(PREVIOUS / 'summary.json')
    for file, digest in previous['evidence_sha256'].items():
        assert sha(ROOT / file) == digest, file
    protocol = ROOT / 'configs/temporal_absolute_kinetics_protocol.json'
    assert sha(protocol) == previous['protocol_sha256']
    all_rows = read(PREVIOUS / 'records.json')
    selected = [r for r in all_rows if r['ranking_status'] == 'pending_campaign_or_analysis_protocol']
    assert len(selected) == 106
    buckets = defaultdict(list)
    for r in selected:
        buckets[r['ranking_group_key']].append(r)
    # No numeric response or model score enters these source/method dispositions.
    decisions = {
      'biostudies-S-EPMC10795190': ('hold_complex_comparison', 'raw/PMC10795190-fullText.xml',
        'Activity Assays at 25 and 5 C; source CSV component_stoichiometry and fixed-substrate columns',
        'Matched temperature/metal does not make this a simple variant comparison: HisGS2 versus HisGS4:HisZ4, with different fixed cosubstrate concentrations.',
        ['Predefine complex representation and catalytic-site normalization for every comparator.', 'Treat regulator addition as a separate complex/context task; do not concatenate chains or count regulator as a variant.']),
      'biostudies-S-EPMC11760959': ('conditions_supported_candidate', 'raw/PMC11760959-fullText.xml',
        'Figure 5 substrate-scope initial-rate caption; existing assay provenance',
        'Same-substrate WT/S58A complex pair; source supports initial rates below 10% conversion and three independent experiments, fit SE.',
        ['Two candidates only; precision and minimum-size rules remain unfrozen.', 'Preserve separate required A/B chains and common complex representation.']),
      'biostudies-S-EPMC12645435': ('hold_source_uncertainty', 'raw/cb5c00120_si_001.pdf',
        'SI p7 and prior detail-audit main Figure 3; three non-inosine substrate pairs',
        'Each substrate pair has a common chemical-condition description, but source-level replicate/error interpretation remains unresolved.',
        ['Do not propagate the D62N numerical discrepancy to these six rows: D62N is not in this selection.', 'Resolve duplicate-versus-triplicate scope and avoid interpreting the SD column name as verified SD.']),
      'biostudies-S-EPMC12667214': ('conditions_supported_candidate', 'raw/PMC12667214-fullText.xml',
        'Enzyme assay; accepted CMPK reciprocal-substrate rows and PPase rows',
        'CMP and ATP series use the same CMPK context but remain separate substrate tasks; PPase has its own temperature and endpoint assay.',
        ['CMP/ATP reuse variants and are statistically dependent.', 'Verify endpoint linearity for the five-minute PPase assay and coupled-reporter nonlimitation for CMPK before scoring.']),
      'biostudies-S-EPMC12751057': ('hold_assay_comparability', 'raw/PMC12751057-fullText.xml',
        'Spectrophotometric assay under steady-state conditions',
        'Shared pH7/25C buffer and initial-rate definition confirmed; enzyme doses, time windows and replicate numbers differ.',
        ['WT has four independent measurements, variants two.', 'Check substrate depletion and initial-rate window for high-enzyme E55Q; different enzyme dose alone is not rejection.']),
      'biostudies-S-EPMC12781114': ('hold_assay_comparability', 'raw/PMC12781114-fullText.xml',
        'Steady-State Kinetics; Table 1 footnotes',
        'WT/P78G share buffer, temperature, 1% ethanol and enzyme concentration; two-substrate fit and WT NADH inhibition require endpoint-comparability review.',
        ['The warning on unreliable NADH parameters is not automatic proof that kcat is wrong.', 'Predefine handling of model-dependent extrapolated turnover and unequal substrate ranges.']),
      'edmond-S0HJ48': ('hold_assay_comparability', 'provenance.json',
        'Existing assay provenance; candidate rows 017-020',
        'Two reciprocal-substrate groups share the same A.comosus WT/deletion pair; assay provenance supports common buffer but fixed-cosubstrate comparability not closed here.',
        ['Inspect deposited kinetic workbooks and fixed-cosubstrate conditions.', 'Do not count reciprocal curves as two independent studies; this source uses inherited evidence this round.']),
      'europepmc-PMC10520331': ('conditions_supported_candidate', 'raw/PMC10520331-fullText.xml',
        'CM, CDT, and coupled CM-CDT in vitro kinetic assays',
        'Separate accepted CM/CDT tasks: four-time-point initial rates, background correction and explicit enzyme normalization; SD from two independently prepared biological replicates documented.',
        ['Keep isolated activity endpoints separate; no sequential coupled CM+CDT rows added.', 'Preserve mono/bifunctional and split-construct dependence; two replicates do not establish high-precision ranks.']),
      'europepmc-PMC11656708': ('conditions_supported_candidate', 'raw/bi4c00350_si_001.pdf',
        'Previously extracted SI assay methods pages 28-29; prior campaign receipt',
        'Four substrate groups retain common CHES/pH9/22C assay evidence from prior supplemental review.',
        ['Retain dose, concentration range, Zn loading and fit limitations.', 'Do not restore rows excluded by the frozen family cap; SI evidence is carried forward, not newly refitted.']),
      'europepmc-PMC11659886': ('conditions_supported_candidate', 'raw/PMC11659886-fullText.xml',
        'Michaelis-Menten Kinetics; prior 29-record detail audit',
        'Both substrates have documented continuous initial rates, same-buffer calibration and enzyme-normalized fits; 21+3 records share one study.',
        ['Keep prior figure-caption range flags and undefined +/- errors open.', 'Preserve three parent blocks for sulfated substrate; nonsulfated WT comparisons are dependent on the same study.']),
      'europepmc-PMC12284513': ('conditions_supported_candidate', 'raw/PMC12284513-fullText.xml',
        'Steady-state kinetic parameters paragraph',
        'After ethanol splitting, each group has common enzyme dose, peroxide, temperature, buffer and C16 substrate; GC-MS initial product rates measured over 0.1-2 minutes.',
        ['2.5% ethanol has two candidates; 10% has three; minimum size remains undecided.', 'Do not pool solvent conditions or replace rates with one-hour conversion yields.']),
      'europepmc-PMC12362431': ('hold_readout_scope', 'raw/PMC12362431-fullText.xml',
        'Wild-Type and Variant Forms steady-state characterization; Table 1',
        'NADPH-consumption turnover is not automatically a product-specific GDP-L-fucose endpoint; one variant changes product composition.',
        ['Define a readout-compatible question before inclusion.', 'Resolve fixed-cofactor/occupancy comparability; do not multiply turnover by reported product fractions to invent labels.']),
      'europepmc-PMC12838360': ('conditions_supported_candidate', 'raw/PMC12838360-fullText.xml',
        'Racemase activity assays; Kinetic analyses; Table 3',
        'L/D substrate directions remain separate, with common 30-minute assay and reported linearity support; reverse detection uses the appropriate oxidase.',
        ['Do not import additional variants beyond the retained pool.', 'Resolve Table3 fixed-5mM footnote versus saturation-series description and variant-specific initial-rate evidence; current selected C193S/WT/A79C subset is kept explicit.'])
    }
    groups = []
    records = []
    for key, members in buckets.items():
        source = members[0]['source_id']
        status, evidence, locator, finding, gaps = decisions[source]
        ep = BASE / source / evidence
        assert ep.exists()
        group = dict(group_key=key, source_id=source, substrate_id=members[0]['substrate_id'],
            candidate_ids=[r['candidate_id'] for r in members], records=len(members),
            status=status, finding=finding, unresolved=gaps, evidence_file=str(ep.relative_to(ROOT)),
            evidence_sha256=sha(ep), evidence_locator=locator,
            distinct_sequence_signatures=len({tuple(r['sequence_signature']) for r in members}),
            source_dependency_block=source, evaluation_ready=False)
        groups.append(group)
        for r in members:
            records.append(dict(source_id=source, candidate_id=r['candidate_id'], group_key=key,
                                updated_campaign_status=status, absolute_status_unchanged=r['absolute_status'],
                                prediction_eligible=False))
    assert len(groups) == 26 and len(records) == 106
    assert {(r['source_id'],r['candidate_id']) for r in records} == {(r['source_id'],r['candidate_id']) for r in selected}
    counts = {}
    for status in sorted({g['status'] for g in groups}):
        subset = [g for g in groups if g['status'] == status]
        counts[status] = dict(groups=len(subset), records=sum(g['records'] for g in subset), sources=len({g['source_id'] for g in subset}))
    result = dict(date='2026-10-06', reviewed_records=106, reviewed_metadata_groups=26,
        reviewed_sources=len(decisions), disposition_counts=counts,
        final_evaluation_ready_groups=0, accepted_additions=0, labels_modified=False,
        absolute_record_statuses_changed=False, model_predictions_generated=False,
        explanation='Conditions-supported means common chemical/measurement context has evidence; initial-rate quality, uncertainty, independence and statistical readiness may still be unresolved. No final pass is assigned.',
        previous_register_sha256=sha(PREVIOUS / 'records.json'), previous_summary_sha256=sha(PREVIOUS/'summary.json'),
        protocol_sha256=sha(protocol), groups=groups)
    OUT.mkdir(parents=True, exist_ok=True)
    for name, obj in [('review.json',result), ('record-overlay.json',records)]:
        (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['groups','explanation']},indent=2))


if __name__ == '__main__':
    main()
