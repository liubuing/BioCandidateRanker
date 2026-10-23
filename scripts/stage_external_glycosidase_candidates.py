"""Reproduce a pending-only literature extraction; never add accepted records.

Reads the two recorded glycosidase sources (AnGH3, BgaP) and the protected
global-audit counters, writes six directly reported kcat/Km measurements to
``candidate-measurements.csv`` under
``artifacts/external/temporal-absolute-kinetics/discovery-20260929`` (all with
status ``pending_not_accepted``), fetches the XED97098.1 reference-only FASTA
with pinned length, and records ``curation-receipt.json`` with the pending
gates, protected-file hashes, and artifact hashes.

Governance role: fail-closed staging -- no accepted additions, no predictions,
mature-construct and substrate-structure gates explicitly unresolved, and
protected config/audit hashes re-verified before writing the receipt.
"""
import csv
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/external/temporal-absolute-kinetics/discovery-20260929'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    """Stage pending measurements and the reference FASTA; write the receipt."""
    protected = [ROOT / 'configs/temporal_absolute_kinetics_protocol.json',
                 ROOT / 'artifacts/external/temporal-global-family-audit/readiness-audit.json']
    before = {str(p.relative_to(ROOT)): digest(p) for p in protected}
    OUT.mkdir(parents=True, exist_ok=True)
    sources = [
        dict(source_id='angh3-2024', doi='10.1186/s40643-024-00784-2',
             url='https://link.springer.com/article/10.1186/s40643-024-00784-2',
             published='2024-07-29', license='CC-BY-4.0', accession='XM_655340.1',
             organism='Aspergillus nidulans', enzyme='AnGH3',
             locator='Results / Kinetic parameters; Table 2 caption',
             sequence_status='pending_mature_assayed_construct',
             note='mRNA accession is not a verified mature protein sequence; cellobiose kcat not transcribed or inferred from rounded quantities.'),
        dict(source_id='bgap-2025', doi='10.1038/s41598-025-02561-3',
             url='https://www.nature.com/articles/s41598-025-02561-3',
             published='2025-05-21', license='CC-BY-4.0', accession='XED97098.1',
             organism='Paenibacillus sp. 8', enzyme='BgaP',
             locator='Results / Properties of recombinant BgaP, kinetic paragraph; Methods / Kinetic studies',
             sequence_status='reference_retrieved_construct_pending',
             note='Native reference accession alone does not establish recombinant N terminus. Methods describe start/stop codons and Sanger verification; reconcile actual construct before acceptance. Temperature rows share the same enzyme-substrate pairs.')]
    rows = []
    values = [(0, 'p-nitrophenyl-beta-D-glucopyranoside', 65, 6.0, 275, .0607),
              (1, 'o-nitrophenyl-beta-D-galactopyranoside', 10, 6.6, 570.41, 3.00),
              (1, 'o-nitrophenyl-beta-D-galactopyranoside', 25, 6.6, 1641.65, 3.62),
              (1, 'o-nitrophenyl-beta-D-galactopyranoside', 40, 6.6, 3136.76, 4.51),
              (1, 'lactose', 10, 6.6, .262, 51.76),
              (1, 'lactose', 30, 6.6, 1.607, 55.50)]
    for i, (s, substrate, temp, ph, kcat, km) in enumerate(values, 1):
        source = sources[s]
        rows.append(dict(candidate_id=f'20260929-{i:03d}', source_id=source['source_id'],
                         measurement_article_doi=source['doi'], source_url=source['url'],
                         source_location=source['locator'], organism=source['organism'],
                         enzyme=source['enzyme'], reference_accession=source['accession'],
                         substrate_name=substrate, temperature_C=temp, pH=ph,
                         kcat_s_1=kcat, km_mM=km, kcat_origin='directly_reported',
                         substrate_structure_id='', status='pending_not_accepted'))
    with (OUT / 'candidate-measurements.csv').open('w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=XED97098&rettype=fasta&retmode=text'
    fasta_path = OUT / 'XED97098.1-reference-only.fasta'
    if not fasta_path.exists():
        raw = urllib.request.urlopen(url, timeout=30).read()
        sequence = ''.join(raw.decode().splitlines()[1:])
        assert raw.startswith(b'>XED97098.1 ') and len(sequence) == 675
        fasta_path.write_bytes(raw)
    fasta = fasta_path.read_text()
    assert fasta.startswith('>XED97098.1 ') and len(''.join(fasta.splitlines()[1:])) == 675
    receipt = dict(schema_version=1, checked_on='2026-09-29', sources=sources,
                   candidate_measurements=len(rows), enzyme_substrate_pairs=3,
                   accepted_additions=0, predictions_generated=False,
                   gates_pending=['assayed_construct_identity', 'substrate_structure_identifiers',
                                  'exact_overlap', 'frozen_mmseqs_homology',
                                  'global_family_mapping_and_label_independent_cap'],
                   sequence_download=dict(url=url, reference_only=True, sha256=digest(fasta_path)),
                   existing_global_audit=json.loads(protected[1].read_text())['counts'],
                   protected_hashes=before,
                   artifacts={p.name: digest(p) for p in [OUT / 'candidate-measurements.csv', fasta_path]},
                   claim_boundary='Pending literature measurements, not a benchmark and not additional independent validation evidence.')
    assert before == {str(p.relative_to(ROOT)): digest(p) for p in protected}
    (OUT / 'curation-receipt.json').write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(candidate_measurements=6, enzyme_substrate_pairs=3, accepted_additions=0, protected_files_unchanged=True)))


if __name__ == '__main__':
    main()
