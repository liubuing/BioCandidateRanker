"""Reproducible sequence/provenance audit; no automatic admission or predictions."""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import xml.etree.ElementTree as ET
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/external/small-molecule-followup-20261007'
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def cds(name):
    r = SeqIO.read(OUT / (name + '.gb'), 'genbank')
    return next(f.qualifiers['translation'][0] for f in r.features if f.type == 'CDS')

protected = {
    'configs/temporal_absolute_kinetics_protocol.json': 'a9ebdc022d77925ebdf32abd016e4ab09fae6a1ab4421dde7bb675107742e36f',
    'artifacts/external/temporal-global-family-audit-20261004-targeted/readiness-audit.json': '487782f6b33c349ea1016eb43ea7119611a902d9e3f3b64be2c872538894d5de',
    'artifacts/homology-final/unique_proteins.fasta': '69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9',
}
for rel, expected in protected.items():
    assert sha(ROOT / rel) == expected, rel
queries = {}
hits = defaultdict(list)
for prefix in ['', 'second-']:
    for r in SeqIO.parse(OUT / (prefix + 'queries.fasta'), 'fasta'):
        queries[r.id] = str(r.seq)
    for line in (OUT / (prefix + 'hits.tsv')).read_text().splitlines():
        fields = line.split('\t')
        assert float(fields[2]) >= .3 and float(fields[4]) >= .8 and float(fields[5]) >= .8
        hits[fields[0]].append(fields)
native = cds('PQ112773')
pdb = str(SeqIO.read(OUT / '9F1T.fasta', 'fasta').seq)
assert len(native) == len(pdb) == 488
differences = [{'position': i, 'GenBank': a, 'PDB': b} for i, (a, b) in enumerate(zip(native, pdb), 1) if a != b]
assert len(differences) == 4
assert pdb == str(SeqIO.read(OUT / '9F3Z.fasta', 'fasta').seq)
lac = cds('MZ817083')
assert lac[87] == 'G' and lac[510] == 'M'

tree = ET.parse(OUT / 'PMC11312515.xml')
table = next(x for x in tree.findall('.//table-wrap') if x.findtext('label') == 'Table 1')
table_text = ' '.join(table.itertext())
measurements = []
for name, value, error, ph, concentration in [
    ('ABTS', 11.20, .38, 2.5, 1.85),
    ('2,6-dimethoxyphenol', .01, .00, 8.0, 94.4),
    ('potassium ferrocyanide', 10.27, 1.07, 4.0, 18.88),
]:
    assert f'{value:.2f}' in table_text
    measurements.append(dict(
        candidate_id='lcoe229-followup-' + str(len(measurements) + 1),
        doi='10.3390/ijms25158521', pmcid='PMC11312515',
        source_location='Table 1; Methods 4.5', substrate=name,
        kcat_s_1=value, reported_plusminus_s_1=error,
        error_type='not confirmed from Table 1; do not infer SD from other figures',
        assay_temperature_C=24, assay_pH=ph,
        enzyme_concentration_ug_mL=concentration,
        buffer='50 mM McIlvaine', assay_volume_uL=200, minimum_replicates=3,
        kinetic_caveat='sigmoidal fit; reported concentration parameter is K0.5, not Km' if name == 'ABTS' else '',
        status='pending_assayed_construct_sequence_conflict',
        ready_for_predictions=False,
        source_sha256=sha(OUT / 'PMC11312515.xml'),
    ))
(OUT / 'staged-measurements.json').write_text(json.dumps(measurements, indent=2) + '\n')
reviews = [
    dict(pmcid='PMC11312515', status='priority_held_sequence_conflict', staged_measurements=3,
         native_pdb_differences=differences,
         other_gaps=['Article reports 1641 bp; deposited CDS is 1467 bp.',
                     'Tagged and tag-free PDB entries have identical 488-aa entity sequences; tag retention and assay construct assignment unresolved.',
                     'Substrate structure identifiers and final global audit pending.']),
    dict(pmcid='PMC12116972', status='held_parent_numbering_mismatch',
         expected={'88': 'I', '511': 'D'}, observed={'88': lac[87], '511': lac[510]},
         reference_hit_count=0, variants_generated=False,
         reason='No unsupported sequence offset or mutation applied.'),
    dict(pmcid='PMC12938175', status='deprioritized_reference_hit_and_parent_unresolved',
         reference='A4XGA6', reference_hit_count=len(hits['A4XGA6_reference']),
         reason='Name-based UniProt candidate is not confirmed assayed CsCE; R17/S142/L184 residues mismatch K17/L142/N184. No variants generated or final construct exclusion inferred.'),
    dict(pmcid='PMC12899640', status='deprioritized_component_hit_and_source_consistency_hold',
         reference='BAA12912.1', reference_hit_count=len(hits['BAA12912_reference']),
         reason='Unprocessed reference differs from reported numbering; exact sfGFP fusion and processed enzyme unconfirmed. Source gives parental kcat=12 and Km=0.5 (ratio 24), but later parental efficiency=12; do not silently reconcile. No mutant kcat inferred from efficiency.'),
    dict(pmcid='PMC11433769', status='pending_exact_construct',
         reason='Recovered Methods: 40 C, pH 7, NADH 0.25 mM, WT 0.5 uM/M12 0.1 uM, 1-750 mM 4-hydroxy-2-butanone, triplicates. No exact reference sequence verified; not screened.'),
    dict(pmcid='PMC12483075', status='fulltext_retry_http500',
         reason='Fetch error is not evidence of biological ineligibility.'),
]
receipt = dict(date='2026-10-07', source_reviews=reviews, accepted_additions=0,
               new_fulltexts=3, additional_staged_measurements=3,
               reference_queries=len(queries), queries_with_hits=sum(bool(hits[q]) for q in queries),
               query_results={q: dict(length=len(s), qualifying_hits=len(hits[q]), hits=hits[q]) for q,s in queries.items()},
               mmseqs_version='5d152c612b6ad2a56f657b7a02c127eceaea2a75',
               thresholds={'min_identity': .3, 'min_both_coverage': .8, 'cov_mode': 0},
               no_hit_scope='Four query entries represent two enzyme sources; no-hit is not full construct or admission approval.',
               protected_hashes=protected, model_predictions_run=False,
               pool={'carry_forward': 201, 'held': 8, 'families': 28, 'substrates': 65},
               artifacts={p.name: sha(p) for p in OUT.iterdir() if p.is_file() and p.name != 'audit.json'})
(OUT / 'audit.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps({k: receipt[k] for k in ['accepted_additions','new_fulltexts','additional_staged_measurements','reference_queries','queries_with_hits']}))
