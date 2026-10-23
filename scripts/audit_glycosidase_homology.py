"""Audit reference and article-supported construct interpretations, without prediction.

Reconstructs AnGH3 (deposited 737-residue translation vs 19-residue signal
removed) and BGA3P (deposited reference vs article I2V) interpretations, then
screens all four query sequences against the development reference with a
pinned MMseqs2 build run through WSL. Writes discovery-20261003/{
audit-queries.fasta, homology-audit.json, excluded-measurements.csv,
source-decisions.json}: all seven angh3/bgap measurements are re-labelled
excluded_development_homology, no record is accepted, and the frozen protocol
and prior audit hashes must be unchanged.
"""
import hashlib
import csv
import copy
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/external/temporal-absolute-kinetics/discovery-20261003'
REF = ROOT / 'artifacts/homology-final/unique_proteins.fasta'
VERSION = '5d152c612b6ad2a56f657b7a02c127eceaea2a75'

def sha(p):
    """Return the SHA256 hex digest of a file."""
    return hashlib.sha256(p.read_bytes()).hexdigest()

def linux(p):
    """Map a Windows path to its WSL /mnt/d/ equivalent for the MMseqs invocation."""
    return '/mnt/d/' + p.relative_to(Path('D:/')).as_posix()

def main():
    """Screen both construct interpretations per enzyme and record the exclusion decision."""
    assert sha(REF) == '69cbc91b88b69f7d4e8fe97f55d7ab157d7c438795989cd6a18fdb0ce28ab2f9'
    protected = [ROOT/'configs/temporal_absolute_kinetics_protocol.json', ROOT/'artifacts/external/temporal-global-family-audit/readiness-audit.json']
    hashes = {str(p.relative_to(ROOT)):sha(p) for p in protected}
    gb = (OUT/'XM_655340.1.gb').read_text()
    a = ''.join(re.search(r'/translation="([^"]+)"', gb)[1].split())
    bpath = ROOT/'artifacts/external/temporal-absolute-kinetics/discovery-20260929/XED97098.1-reference-only.fasta'
    b = ''.join(bpath.read_text().splitlines()[1:])
    assert len(a)==737 and a[19:].startswith('TGQVL')
    assert len(b)==675 and b.startswith('MISNKLPKMF')
    seqs = {'angh3_reference_737':a, 'angh3_mature_718':a[19:],
            'bgap_reference_675':b, 'bgap_article_I2V_675':b[:1]+'V'+b[2:]}
    query = OUT/'audit-queries.fasta'
    query.write_text(''.join(f'>{k}\n{v}\n' for k,v in seqs.items()), encoding='ascii')
    exe = ['wsl','-d','Ubuntu-24.04-D','--','/home/liubuing/.local/bin/mmseqs']
    version = subprocess.check_output(exe+['version']).decode().strip()
    assert version == VERSION, version
    hits = OUT/'homology-hits.tsv'
    cmd = exe + ['easy-search',linux(query),linux(REF),linux(hits),linux(OUT/'search-tmp'),
                 '--min-seq-id','0.3','-c','0.8','--cov-mode','0','--threads','4',
                 '--format-output','query,target,fident,alnlen,qcov,tcov,evalue,bits']
    with (OUT/'mmseqs.log').open('wb') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    result = {}
    for name in seqs:
        rows = [x.split('\t') for x in hits.read_text().splitlines() if x.split('\t')[0]==name]
        qualified = [r for r in rows if float(r[2])>=.3 and float(r[4])>=.8 and float(r[5])>=.8]
        result[name] = dict(length=len(seqs[name]), qualifying_hits=len(qualified),
                            best_hit=max(qualified,key=lambda r:float(r[7])) if qualified else None)
    assert hashes == {str(p.relative_to(ROOT)):sha(p) for p in protected}
    receipt = dict(date='2026-10-03', mmseqs_version=version, command=cmd,
                   reference_sha256=sha(REF), query_sha256=sha(query), hits_sha256=sha(hits),
                   protected_hashes=hashes, results=result, predictions_generated=False,
                   raw_sources={p.name:sha(p) for p in [OUT/'PMC12095705.xml', OUT/'PMC11286919.xml', OUT/'XM_655340.1.gb', bpath]},
                   sequence_interpretation={
                       'angh3':'Accession translation has 737 residues. Removing reported 19-residue signal yields 718 and TGQVL start; article length wording is inconsistent. Both forms screened.',
                       'bgap':'Article forward coding sequence supports I2V relative to deposited MIS... reference; stop codon argues against C-terminal vector fusion. Both interpretations screened; no independent construct deposition verified.'})
    (OUT/'homology-audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    old = ROOT/'artifacts/external/temporal-absolute-kinetics/discovery-20260929/candidate-measurements.csv'
    with old.open(encoding='utf-8', newline='') as f:
        measurements = list(csv.DictReader(f))
    cell = copy.deepcopy(measurements[0])
    cell.update(candidate_id='20261003-007', substrate_name='cellobiose',
                kcat_s_1='73.5', km_mM='2.71',
                source_location='Table 2 / cellobiose / Michaelis-Menten / Kcat (s-1)')
    measurements.append(cell)
    for row in measurements:
        prefix = 'angh3' if row['source_id']=='angh3-2024' else 'bgap'
        keys = [k for k in result if k.startswith(prefix)]
        assert all(result[k]['qualifying_hits'] > 0 for k in keys)
        row['status'] = 'excluded_development_homology'
        row['decision_evidence'] = 'homology-audit.json; ' + ','.join(keys)
    with (OUT/'excluded-measurements.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(measurements[0]))
        w.writeheader()
        w.writerows(measurements)
    decision = dict(date='2026-10-03', sources_screened=2, measurement_conditions_screened=7,
                    distinct_enzyme_substrate_pairs=4, excluded_measurements=7,
                    accepted_additions=0, exclusion='development_homology',
                    fit_policy='Michaelis-Menten rows retained for extraction; Hanes-Woolf alternative fits are not independent measurements.',
                    stopping_rule='A failed homology gate suffices for exclusion; substrate identifiers and exact overlap need not be completed for acceptance.',
                    sequence_caveat='Construct interpretations remain article-derived; exclusion is supported for both reference and stated construct interpretations.',
                    artifacts={p.name:sha(p) for p in [OUT/'homology-audit.json', OUT/'excluded-measurements.csv']},
                    previous_snapshot='discovery-20260929 retained unchanged; superseded for these two source decisions',
                    predictions_generated=False)
    (OUT/'source-decisions.json').write_text(json.dumps(decision,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
