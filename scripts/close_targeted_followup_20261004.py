"""Close public-source follow-up without promoting unresolved measurements.

Compares deposited 9QEV/9QEZ/5JGY construct sequences against the CaCA
reference, re-checks the precomputed homology hits, and records per-source
holds (mutant loop discrepancy, mature N-terminus boundary) plus evidence
paragraphs from the source XMLs. Writes
targeted-followup-20261004/source-evidence.json and closeout.json with
new_accepted_records=0; the frozen protocol and prior global audit hashes are
re-verified and no recluster or predictions are run.
"""
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/external/temporal-absolute-kinetics'
OUT = BASE / 'targeted-followup-20261004'


def sha(path):
    """Return the SHA256 hex digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    """Verify constructs, record source dispositions, and write the closeout."""
    refs = json.loads((BASE / 'public-batch2-20261004/sequence-register.json').read_text())
    native = next(r['sequence'] for r in refs if r['name'] == 'CaCA')
    sequences = {}
    for code in ['9QEV', '9QEZ', '5JGY']:
        cif = MMCIF2Dict(str(OUT / (code + '.cif')))
        sequences[code] = ''.join(cif['_entity_poly.pdbx_seq_one_letter_code_can'][0].split())
        if code.startswith('9Q'):
            assert '10.1002/pro.70396' in cif['_citation.pdbx_database_id_DOI']
    assert sequences['9QEV'] == native + 'HHHHHH'
    assert len(sequences['9QEZ']) == 180
    hits = (OUT / 'homology/homology_hits.tsv').read_text().splitlines()
    assert hits and all(line.startswith('5JGY\t') for line in hits)
    assert not (OUT / 'additional-homology/homology_hits.tsv').read_text().strip()
    diff = json.loads((OUT / 'sequence-discrepancies.json').read_text())
    assert diff['PSA_reference_matches_supplement']
    assert diff['CaCA_loop_narrative'] != diff['CaCA_loop_deposited']
    evidence = {}
    for pmc, terms in [('PMC12648627', ['Glu85‐Thr', '4.46', '9QEZ']),
                       ('PMC11372007', ['RH32919', 'signal sequence']),
                       ('PMC12302338', ['5JGY'])]:
        path = BASE / 'public-batch2-20261004' / (pmc + '.xml')
        tree = ET.parse(path)
        paragraphs = [''.join(p.itertext()) for p in tree.iter('p')]
        evidence[pmc] = dict(source_sha256=sha(path),
                             paragraphs=[t for t in paragraphs if any(term in t for term in terms)])
    (OUT / 'source-evidence.json').write_text(json.dumps(evidence, indent=2, ensure_ascii=False), encoding='utf-8')
    protocol = ROOT / 'configs/temporal_absolute_kinetics_protocol.json'
    assert sha(protocol) == 'a9ebdc022d77925ebdf32abd016e4ab09fae6a1ab4421dde7bb675107742e36f'
    previous = ROOT / 'artifacts/external/temporal-global-family-audit-20261004-targeted/readiness-audit.json'
    audit = json.loads(previous.read_text())
    previous_receipt = json.loads((BASE / 'targeted-curation-20261004/closeout.json').read_text())
    assert sha(previous) == previous_receipt['global_audit_sha256']
    receipt = dict(date='2026-10-04', new_accepted_records=0, current_counts=audit['counts'],
                   remaining_records=91, remaining_families=2, model_predictions_run=False,
                   readiness_gate_passes=False, global_recluster_run=False,
                   reason_no_recluster='No source admitted; accepted pool unchanged.',
                   protocol_sha256=sha(protocol), previous_global_audit_sha256=sha(previous),
                   source_dispositions=[
                       dict(pmc='PMC12302338', status='deposited_construct_homology_excluded',
                            tested_construct='PDB 5JGY', qualifying_hits=len(hits),
                            evidence=hits, scope='Tested AKR4C13 construct, not every protein mentioned in article'),
                       dict(pmc='PMC12648627', status='held_construct_discrepancy', pending_measurements=1,
                            full_deposited_mutant_hits=0, tag_free_mutant_hits=0,
                            reason='Deposited mutant loop differs from explicit narrative donor-loop sequence; assay construct cannot yet be assigned unambiguously.'),
                       dict(pmc='PMC11372007', status='held_mature_sequence_boundary', pending_measurements=2,
                            supplement_translation_matches_reference=True, precursor_hits=0,
                            reason='Original-signal construct RH32919 identified; exact secreted mature N terminus remains unresolved. No-hit precursor alone does not resolve actual assay sequence.')],
                   files_sha256={str(p.relative_to(OUT)): sha(p) for p in OUT.rglob('*')
                                 if p.is_file() and p.name != 'closeout.json'})
    (OUT / 'closeout.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'accepted_additions': 0, 'counts': audit['counts'], 'protected_artifacts_unchanged': True}))


if __name__ == '__main__':
    main()
