"""Record source-level decisions without admitting records or running predictions.

Closes the 2026-10-07 priority-source follow-up into
priority-source-followup-20261007/closeout.json: an FtsQ/FLS reference
homology screen (deposition-derived reference plus in-silico variants, not
confirmed constructs) and a polymer-substrate scope exclusion. The frozen
protocol hash is asserted at import time and every precomputed MMseqs hit line
is re-checked against the pinned thresholds; accepted_additions stays 0.
"""
from pathlib import Path
import hashlib
import json
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/external/priority-source-followup-20261007'

def sha(path):
    """Return the SHA256 hex digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()

protocol = ROOT / 'configs/temporal_absolute_kinetics_protocol.json'
assert sha(protocol) == 'a9ebdc022d77925ebdf32abd016e4ab09fae6a1ab4421dde7bb675107742e36f'
hits = defaultdict(list)
for line in (OUT / 'fls-reference-hits.tsv').read_text().splitlines():
    q, target, identity, length, qc, tc, ev, bits = line.split('\t')
    assert float(identity) >= .3 and float(qc) >= .8 and float(tc) >= .8
    hits[q].append(dict(target=target, identity=float(identity), query_coverage=float(qc), target_coverage=float(tc)))
assert len(hits) == 5
result = {
    'date': '2026-10-07', 'accepted_additions': 0,
    'model_predictions_run': False, 'protocol_sha256': sha(protocol),
    'reference_sha256': sha(ROOT / 'artifacts/homology-final/unique_proteins.fasta'),
    'fls_reference_screen': {
        'sequence_source': 'https://www.rcsb.org/fasta/entry/4QQ8/display',
        'scope': 'PDB reference including deposited termini and four in-silico variants, NOT confirmed experimental constructs',
        'mutation_residue_checks': {'I28': True, 'L482': True, 'F484': True},
        'mmseqs_thresholds': {'identity': .3, 'coverage': .8, 'cov_mode': 0},
        'qualifying_hits': dict(hits),
        'query_sha256': sha(OUT / 'fls-reference-variants.fasta'),
        'hits_sha256': sha(OUT / 'fls-reference-hits.tsv'),
    },
    'source_decisions': [
        {'pmcid': 'PMC11357561', 'candidate_rows': 5,
         'status': 'not_admitted_deprioritized_reference_homology',
         'reason': 'All five reference queries have qualifying development hits. Exact experimental constructs remain unconfirmed; no blanket final construct exclusion claimed.',
         'unresolved': ['experimental construct sequence', 'kinetics assay temperature and initial-rate evidence', 'supplement retrieval'],
         'parent_construct_reference_doi': '10.1039/C7GC02898A'},
        {'pmcid': 'PMC11969528', 'candidate_rows': 6,
         'status': 'excluded_from_current_absolute_protocol_substrate_scope',
         'enzyme': 'AxyAgu115A GH115 glycosidase',
         'substrate': 'beechwood xylan',
         'reason': 'Polymeric substrate does not supply a single defined molecular structure required by the frozen protocol; do not substitute released 4-O-methyl glucuronic acid as substrate.',
         'optional_ranking_use': 'pending separate explicitly frozen polymer-compatible design, construct and assay review',
         'license_caveat': 'Article CC-BY-4.0; Figshare supplement landing page states CC-BY-NC-4.0. Do not apply article license automatically to supplemental material.',
         'supplement_landing_url': 'https://acs.figshare.com/articles/journal_contribution/28655777'},
    ],
    'pool': {'historical_records': 209, 'carry_forward_records': 201, 'held_records': 8,
             'carry_forward_families': 28, 'carry_forward_substrates': 65,
             'record_gap': 99, 'family_gap': 2},
    'formal_run_status': 'not_ready',
    'download_caveat': 'fls-pmc.zip is HTML challenge content, not a usable ZIP; no supplement recovered',
    'next_screen_order': ['defined small-molecule substrate', 'primary measurement date and license',
                          'reference homology prescreen', 'exact constructs and variants',
                          'measurement and condition extraction', 'global family and readiness audit'],
}
(OUT / 'closeout.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps({'accepted_additions': 0, 'reference_queries_with_hits': len(hits), 'total_hits': sum(map(len, hits.values()))}))
