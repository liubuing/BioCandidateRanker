"""Reproducible public-source discovery; never modifies accepted data or runs models.

Runs the eight predeclared Europe PMC direction queries (first publication
2023-01-01 through 2026-10-04, open access only), snapshots search pages and
full-text XML with license and kinetics-table evidence, and records everything
in public-batch-20261004/batch-register.json with accepted_additions=0. The
CC-BY/CC0 license screen and table mentions are screening signals only; the
frozen protocol and pool hashes must be identical before and after.

中文：本脚本做可复现的公开来源发现，绝不修改已接受数据、绝不运行模型：执行八条预先
声明的 Europe PMC 方向查询（首发 2023-01-01 至 2026-10-04，仅开放获取），对检索页与
全文 XML 连同许可及动力学表格证据做快照，并把一切记入
public-batch-20261004/batch-register.json，其中 accepted_additions=0。CC-BY/CC0 许可
筛查与表格提及只是筛查信号；冻结协议与池哈希在前后必须一致。治理角色：发现批次审计
收据（fail-closed，保护产物不可变动）。
"""
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
    """Return the SHA256 hex digest of a file.

    中文：返回文件的 SHA256 十六进制摘要。
    """
    return hashlib.sha256(p.read_bytes()).hexdigest()

def fetch(url, path):
    """Download url into path once, then always serve the cached local bytes.

    中文：把 url 下载到 path 一次，之后始终提供本地缓存字节。
    """
    if not path.exists():
        req = urllib.request.Request(url, headers={'User-Agent': 'PublicKineticsEvidenceAudit/1.0'})
        data = urllib.request.urlopen(req, timeout=45).read()
        path.write_bytes(data)
    return path.read_bytes()

def search(item):
    """Run one direction's Europe PMC query and return its receipt and raw results.

    中文：执行一个方向的 Europe PMC 查询并返回其收据与原始结果。
    """
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
    """Acquire one article's full text, screen its license, and extract evidence.

    中文：获取一篇文章的全文，筛查其许可并提取证据。
    """
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
    """Discover, screen, and register the public batch without accepting anything.

    中文：发现、筛查并登记公开批次，不接受任何记录。
    """
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
