"""Second, explicitly selected batch of 20 benign enzymology primary-source leads."""
import concurrent.futures as cf
import json
import batch_public_kinetics_20261004 as discovery

BASE=discovery.BASE
OUT=BASE/'public-batch2-20261004'
SELECTED=['PMC13396203','PMC12770603','PMC12896513','PMC12960341','PMC12734981',
          'PMC12648627','PMC11118615','PMC13269172','PMC12983273','PMC11975207',
          'PMC12009067','PMC12994016','PMC12384609','PMC12302338','PMC10339645',
          'PMC11199699','PMC13358373','PMC12817654','PMC10572201','PMC11372007']

def main():
    OUT.mkdir(exist_ok=True)
    prior=BASE/'public-batch-20261004/batch-register.json'
    original=json.loads(prior.read_text(encoding='utf-8'))
    by_id={e.get('pmcid'):e for e in original['entries']}
    entries=[dict(by_id[x]) for x in SELECTED]
    assert len(entries)==20 and all(e['status']=='metadata_only' for e in entries)
    discovery.OUT=OUT
    with cf.ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(discovery.fulltext,entries))
    audit=dict(date='2026-10-04',selection='20 distinct previously unreviewed title leads from prior 637 records; emphasis on defined small substrates and published sequence provenance, not predicted performance.',
               prior_register_sha256=discovery.sha(prior),entries=entries,
               fulltexts=sum('fulltext_sha256' in x for x in entries),
               compatible_article_licenses=sum(x.get('license_screen')=='compatible_article_license' for x in entries),
               accepted_additions=0,predictions_generated=False)
    (OUT/'batch-register.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    for e in entries:print(e['pmcid'],e['doi'],e.get('license_screen'),e.get('kinetics_table_mentions'),e['title'])
if __name__=='__main__':main()
