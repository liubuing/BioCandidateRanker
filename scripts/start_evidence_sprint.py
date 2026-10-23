"""Metadata/source screening only; no training or candidate scoring.

Writes the EstSIT01 closeout (superseding the measurement-consistency hold),
runs two pinned Europe PMC search queries, screens the returned metadata
against prior local material, downloads and license-screens four selected full
texts (one becomes ``priority_sequence_prescreen``, non-eligible licenses are
recorded as such), and writes ``discovery-register.json`` plus per-article
table extracts under ``artifacts/external/temporal-absolute-kinetics/
evidence-sprint-20261003``.

Governance role: discovery bookkeeping -- every record stays
``metadata_only_not_qualified`` unless explicitly dispositioned;
``accepted_additions`` stays 0, batch leads with incompatible licenses are
excluded from direct import, and protected config/audit hashes are re-verified
before the register is written.
"""
import hashlib
import json
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/external/temporal-absolute-kinetics'
OUT=BASE/'evidence-sprint-20261003'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fetch(url,path):
    """Return cached bytes for a URL, downloading to path on first use."""
    if not path.exists():
        path.write_bytes(urllib.request.urlopen(url,timeout=30).read())
    return path.read_bytes()

def main():
    """Write the closeout, run pinned searches, screen full texts, and register results."""
    OUT.mkdir(exist_ok=True)
    protected=[ROOT/'configs/temporal_absolute_kinetics_protocol.json',ROOT/'artifacts/external/temporal-global-family-audit/readiness-audit.json']
    hashes={str(p.relative_to(ROOT)):sha(p) for p in protected}
    review=BASE/'estsit01-curation-20261003/reviewed-measurements.csv'
    closeout=dict(date='2026-10-03',supersedes='group-wide measurement consistency hold in construct-review.json',
        sources=[dict(doi='10.1186/s40643-024-00776-2',rows=1,
            status='pending_structure_crosscheck_and_source_finalization',
            measurement_consistency_flag=False,reason='Later-paper Table 3 discrepancy is not evidence against this earlier measurement.'),
            dict(doi='10.3390/molecules29245839',rows=4,status='pending_measurement_clarification',
            measurement_consistency_flag=True,reason='14.9%-27.6% conversion discrepancy under same-preparation assumptions; not proven error, no label rescaling.')],
        public_review='One bounded public-material review completed; do not repeat unchanged source investigation without new evidence.',
        draft_questions='docs/ESTSIT01_SOURCE_REVIEW_20261003_CN.md',author_contact_sent=False,
        accepted_additions=0,original_labels_preserved=True,prior_review_sha256=sha(review))
    (OUT/'estsit01-closeout.json').write_text(json.dumps(closeout,indent=2)+'\n')
    queries=[
      '(TITLE_ABS:"kinetic characterization" AND (TITLE_ABS:orthologs OR TITLE_ABS:homologs OR TITLE_ABS:panel)) AND FIRST_PDATE:[2023-01-01 TO 2026-10-03] AND OPEN_ACCESS:y',
      '(TITLE_ABS:"high-throughput" AND TITLE_ABS:"enzyme kinetics") AND FIRST_PDATE:[2023-01-01 TO 2026-10-03] AND OPEN_ACCESS:y']
    entries={}; searches=[]
    oldfiles=[p for p in BASE.rglob('*') if p.suffix in {'.json','.md'} and OUT not in p.parents]
    oldtext='\n'.join(p.read_text(encoding='utf-8',errors='replace').lower() for p in oldfiles)
    for i,q in enumerate(queries,1):
        url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode(dict(query=q,format='json',pageSize=100))
        raw=OUT/f'query-{i}.json';data=json.loads(fetch(url,raw))
        searches.append(dict(query=q,url=url,total_hits=data['hitCount'],returned=len(data['resultList']['result']),sha256=sha(raw)))
        for x in data['resultList']['result']:
            key=x.get('doi') or x['id']
            entries[key]=dict(doi=x.get('doi'),pmcid=x.get('pmcid'),title=x['title'],
                              prior_local_mention=key.lower() in oldtext,status='metadata_only_not_qualified')
    selected=['PMC12292525','PMC10882141','PMC12044592','PMC10629569']
    for pmc in selected:
        url=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML'
        file=OUT/(pmc+'.xml')
        try:
            tree=ET.fromstring(fetch(url,file))
            doi=next(n.text for n in tree.iter('article-id') if n.get('pub-id-type')=='doi')
            licenses=[''.join(n.itertext()) for n in tree.iter('license')]
            tables=[dict(id=n.get('id'),text=''.join(n.itertext())) for n in tree.iter('table-wrap') if 'cat' in ''.join(n.itertext()).lower()]
            entries[doi].update(status='fulltext_initial_screen',fulltext_url=url,fulltext_sha256=sha(file),
                                license_evidence=licenses,kinetics_table_mentions=len(tables))
            if pmc=='PMC12292525':
                entries[doi].update(disposition='priority_sequence_prescreen',direct_kcat_rows_verified=2,
                    reason='Table 2: two plant CCD enzymes; article CC BY 4.0. Sequences, constructs and homology not yet audited; not two established global families.')
            else:
                entries[doi].update(disposition='article_license_not_eligible',
                    reason='Article license is non-commercial/no-derivatives; separate compatible measurement dataset license not verified. No direct import.')
            (OUT/(pmc+'-tables.json')).write_text(json.dumps(tables,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        except Exception as e:
            for x in entries.values():
                if x.get('pmcid')==pmc:x.update(status='acquisition_blocked',error=str(e))
    audit=dict(date='2026-10-03',searches=searches,unique_metadata_records=len(entries),entries=list(entries.values()),
        batch_leads=[dict(name='EnzyExtractDB',url='https://github.com/ChemBioHTP/EnzyExtract',
             status='excluded_direct_import_license',reason='README declares data CC BY-NC 4.0; MIT applies to code. May use original-paper references as discovery leads only.'),
             dict(name='Fordyce data portal',url='https://www.fordycelab.com/data',status='priority_source_index_pending',
             reason='Portal advertises 23765 master rows; not a count of eligible kcat measurements. Resolve primary dates, licenses and endpoints per study.')],
        source_size_policy='>=10 rows is prioritization only for this new sprint, not a new acceptance gate; prior rounds preserved.',
        actual_accepted_additions=0,predictions_generated=False,protected_hashes=hashes)
    assert hashes=={str(p.relative_to(ROOT)):sha(p) for p in protected}
    (OUT/'discovery-register.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(metadata_records=len(entries),fulltexts=sum(x['status']=='fulltext_initial_screen' for x in entries.values()),accepted=0)))
    for x in entries.values():
        if x['status']=='fulltext_initial_screen':print(x['pmcid'],x['doi'],x['kinetics_table_mentions'],x['license_evidence'][0][:150])

if __name__=='__main__':main()
