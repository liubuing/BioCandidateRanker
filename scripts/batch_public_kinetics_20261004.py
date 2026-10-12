"""Reproducible public-source discovery; never modifies accepted data or runs models."""
import concurrent.futures as cf
import hashlib
import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'artifacts/external/temporal-absolute-kinetics'
OUT = BASE / 'public-batch-20261004'
DIRECTIONS = {
    'esterase': 'esterase',
    'plant_transferase': '(plant AND (glycosyltransferase OR methyltransferase))',
    'sugar_isomerase': '(isomerase AND (sugar OR xylose OR arabinose OR ribose))',
    'oxidoreductase': '(dehydrogenase AND (plant OR yeast OR thermophilic))',
    'carbohydrate': '(xylanase OR amylase OR cellulase)',
    'lyase': '(aldolase OR decarboxylase)',
    'phosphorylase': '(phosphorylase OR phosphomutase)',
    'industrial_hydrolase': '(cutinase OR lipase OR lactonase)',
}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def fetch(url, path):
    if not path.exists():
        req = urllib.request.Request(url, headers={'User-Agent': 'PublicKineticsEvidenceAudit/1.0'})
        data = urllib.request.urlopen(req, timeout=45).read()
        path.write_bytes(data)
    return path.read_bytes()

def search(item):
    name, terms = item
    q = f'{terms} AND (kcat OR "turnover number" OR "kinetic characterization") AND FIRST_PDATE:[2023-01-01 TO 2026-10-04] AND OPEN_ACCESS:y'
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode(dict(query=q, format='json', pageSize=100))
    p = OUT / f'search-{name}.json'
    try:
        data = json.loads(fetch(url, p))
        return dict(direction=name, query=q, url=url, hit_count=data['hitCount'], returned=len(data['resultList']['result']), sha256=sha(p)), data['resultList']['result']
    except Exception as e:
        return dict(direction=name, query=q, url=url, error=str(e)), []

def fulltext(entry):
    pmc = entry['pmcid']
    p = OUT / f'{pmc}.xml'
    url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML'
    try:
        tree = ET.fromstring(fetch(url, p))
        licenses = [dict(text=' '.join(n.itertext()), attributes=n.attrib,
                         links=[v for x in n.iter() for k,v in x.attrib.items() if k.endswith('href')]) for n in tree.iter('license')]
        urls = [u.lower() for x in licenses for u in x['links']]
        urls += [u.lower() for x in licenses for u in re.findall(r'https?://[^\s<>]+', x['text'])]
        eligible = any(re.search(r'creativecommons.org/licenses/by/4\.0', u) or 'creativecommons.org/publicdomain/zero/1.0' in u for u in urls)
        restricted = any('by-nc' in u or 'by-nd' in u for u in urls)
        tables = [dict(id=n.get('id'), text=' '.join(n.itertext())) for n in tree.iter('table-wrap') if re.search(r'k\s*cat|turnover', ' '.join(n.itertext()), re.I)]
        paragraphs = [' '.join(n.itertext()) for n in tree.iter('p') if re.search(r'accession|genbank|uniprot|deposited|k\s*cat|turnover', ' '.join(n.itertext()), re.I)]
        evidence = dict(licenses=licenses, tables=tables, paragraphs=paragraphs)
        ep = OUT / f'{pmc}-evidence.json'
        ep.write_text(json.dumps(evidence, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        entry.update(fulltext_url=url, fulltext_sha256=sha(p), evidence_sha256=sha(ep),
                     license_screen='compatible_article_license' if eligible and not restricted else 'restricted_or_unresolved_license',
                     kinetics_table_mentions=len(tables), status='fulltext_screened_pending_manual_verification')
    except Exception as e:
        entry.update(status='fulltext_fetch_failed', error=str(e))
    return entry

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    protected = [ROOT/'configs/temporal_absolute_kinetics_protocol.json', ROOT/'artifacts/external/temporal-global-family-audit/readiness-audit.json']
    before = {str(p.relative_to(ROOT)): sha(p) for p in protected}
    snapshot = OUT/'prior-doi-mentions.json'
    if not snapshot.exists():
        prior = set()
        for p in BASE.rglob('*'):
            if p.suffix in {'.json', '.md', '.csv'} and OUT not in p.parents:
                prior.update(x.rstrip('.,;') for x in re.findall(r'10\.\d{4,9}/[-._;()/:a-z0-9]+', p.read_text(encoding='utf-8', errors='replace').lower()))
        snapshot.write_text(json.dumps(sorted(prior), indent=2)+'\n')
    prior = set(json.loads(snapshot.read_text()))
    entries = {}; searches = []
    with cf.ThreadPoolExecutor(max_workers=3) as pool:
        for search_info, records in pool.map(search, DIRECTIONS.items()):
            searches.append(search_info)
            for r in records:
                doi = r.get('doi', '').lower(); key = doi or r.get('pmcid') or r['id']
                if key not in entries:
                    entries[key] = dict(doi=doi, pmcid=r.get('pmcid'), title=r.get('title'), publication_date=r.get('firstPublicationDate'),
                                        prior_local_mention=doi in prior if doi else False, directions=[], status='metadata_only')
                entries[key]['directions'].append(search_info['direction'])
    # First pass prioritizes new benign industrial/plant sources across directions.
    # Title filter only affects review priority, not scientific inclusion criteria.
    deprioritize = re.compile(r'review|pathogen|virulen|toxin|toxic|antibiotic|antimicrobial|cancer|tumor|virus|infection', re.I)
    selected = []
    for direction in DIRECTIONS:
        candidates = [e for e in entries.values() if direction in e['directions'] and e['pmcid'] and not e['prior_local_mention'] and not deprioritize.search(e['title'])]
        for e in candidates[:6]:
            if e not in selected: selected.append(e)
    with cf.ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(fulltext, selected))
    assert before == {str(p.relative_to(ROOT)):sha(p) for p in protected}
    audit = dict(date='2026-10-04', searches=searches, unique_metadata_records=len(entries),
                 prior_local_doi_mentions=sum(e['prior_local_mention'] for e in entries.values()),
                 fulltexts_acquired=sum('fulltext_sha256' in e for e in entries.values()),
                 compatible_article_licenses=sum(e.get('license_screen')=='compatible_article_license' for e in entries.values()),
                 accepted_additions=0, predictions_generated=False, protected_hashes=before,
                 note='Article license and table mentions are screening signals only. Manual primary measurement, sequence, substrate, homology and global cap audits required. Search pages capped at 100 each; not exhaustive.',
                 entries=list(entries.values()))
    (OUT/'batch-register.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k not in {'searches','entries','protected_hashes','note'}}, ensure_ascii=False))
    for e in selected:
        print(e['pmcid'],e['doi'],e.get('license_screen',e['status']),e.get('kinetics_table_mentions'),e['title'])

if __name__ == '__main__': main()
