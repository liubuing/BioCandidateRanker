"""Stage published benign enzyme components; homology is a screen, not acceptance.

Reads the staged article XMLs under
``artifacts/external/temporal-absolute-kinetics/public-batch-20261004``, fetches
the listed nucleotide accessions (and the UniProt Q9HJS7 reference), validates
each GenBank CDS translation against its deposited protein, and stages the
published catalytic components -- plus the documented LbXyn43B_T274A
substitution -- into ``component-queries.fasta`` and
``sequence-register.json``.

Governance role: fail-closed staging only -- every entry carries scope
"published catalytic component; full assayed construct unresolved"; length
discrepancies are recorded, never repaired; fetch/validation failures are
captured per row instead of aborting the batch; ``accepted_additions`` stays 0.

中文：暂存已发表的良性酶组件；同源性只是筛查而非接受：取全文清单列出的
核苷酸 accession，校验并暂存组件，写登记册。
"""
import hashlib
import json
import urllib.request
from pathlib import Path
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/external/temporal-absolute-kinetics/public-batch-20261004'
SOURCES=[
 ('AaCSE1','PV962689','PMC12774947',324),('MeMAGL3','PV962700','PMC12774947',344),
 ('LbXyn43B','KX174294','PMC12432433',None),
 ('CiCQE1','PX093951','PMC12405362',None),('CiCQE3','PX093953','PMC12405362',None),
 ('LlCCR','PQ867195','PMC13249509',None),('LlDBR1','PQ867197','PMC13249509',None),
 ('LlDBR2','PQ867198','PMC13249509',None),('LlALDH1','PX215327','PMC13249509',None),
 ('LlALDH2','PX215326','PMC13249509',None),('LlP2GT1','PQ867199','PMC13249509',None),
 ('LlP2GT2','PQ867200','PMC13249509',None),('LlP4GT','PQ867201','PMC13249509',None)]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    """Fetch, validate, and stage each published component; write the register.
中文：取数、校验并暂存每个已发表组件；写登记册。
"""
    rows=[];seqs={}
    for name,accession,pmc,length in SOURCES:
        article=OUT/(pmc+'.xml'); assert accession in article.read_text(encoding='utf-8')
        row=dict(name=name,accession=accession,pmcid=pmc,source_sha256=sha(article),scope='published catalytic component; full assayed construct unresolved')
        p=OUT/(accession+'.gb')
        url=f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=gb&retmode=text'
        try:
            if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=35).read())
            record=SeqIO.read(p,'genbank')
            coding=[f for f in record.features if f.type=='CDS' and 'translation' in f.qualifiers]
            assert len(coding)==1, 'Expected one unambiguous translated CDS'
            seq=coding[0].qualifiers['translation'][0]
            translated=str(coding[0].extract(record.seq).translate()).rstrip('*')
            assert translated==seq, 'CDS translation disagrees with deposited protein'
            if length is not None and len(seq)!=length:
                row.update(paper_reported_length=length,length_discrepancy=len(seq)-length,
                           note='Paper length includes one extra residue versus deposited CDS translation; retain deposited sequence for prescreen only, do not silently repair.')
            assert set(seq)<=set('ACDEFGHIKLMNPQRSTVWY')
            seqs[name]=seq
            row.update(status='sequence_staged',record_id=record.id,protein_id=coding[0].qualifiers.get('protein_id'),length=len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),record_sha256=sha(p),url=url)
        except Exception as e:row.update(status='sequence_fetch_or_validation_failed',error=repr(e))
        rows.append(row)
    # Independent UniProt reference for a thermostable archaeal esterase.
    p=OUT/'Q9HJS7.fasta';url='https://rest.uniprot.org/uniprotkb/Q9HJS7.fasta'
    row=dict(name='Ta0887',accession='Q9HJS7',pmcid='PMC13498779',scope='reference component; synthesized insert and tag require crosscheck',url=url)
    try:
        assert 'Q9HJS7' in (OUT/'PMC13498779.xml').read_text(encoding='utf-8')
        if not p.exists():p.write_bytes(urllib.request.urlopen(url,timeout=35).read())
        seq=str(SeqIO.read(p,'fasta').seq);assert len(seq)==235
        seqs['Ta0887']=seq;row.update(status='sequence_staged',length=len(seq),record_sha256=sha(p),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())
    except Exception as e:row.update(status='sequence_fetch_or_validation_failed',error=str(e))
    rows.append(row)
    if 'LbXyn43B' in seqs:
        seq=seqs['LbXyn43B'];assert seq[273]=='T'
        seqs['LbXyn43B_T274A']=seq[:273]+'A'+seq[274:]
        rows.append(dict(name='LbXyn43B_T274A',pmcid='PMC12432433',status='sequence_staged',scope='published Table 4 substitution on deposited component; not full construct',length=len(seq),parent='LbXyn43B'))
    (OUT/'component-queries.fasta').write_text(''.join(f'>{k}\n{v}\n' for k,v in seqs.items()))
    (OUT/'sequence-register.json').write_text(json.dumps(dict(entries=rows,query_sha256=sha(OUT/'component-queries.fasta'),accepted_additions=0),indent=2)+'\n')
    for r in rows:print(r['name'],r['status'],r.get('length'),r.get('error',''))
if __name__=='__main__':main()
