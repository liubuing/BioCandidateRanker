"""Export a deterministic development-only CatPred plumbing check."""
from pathlib import Path
import csv, gzip, hashlib, json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/catpred-mode-b/smoke-data'

def main():
    OUT.mkdir(exist_ok=True)
    records={}; receipt={'purpose':'Development plumbing check only; no test rows, no performance claim', 'partitions':{}}
    for partition, count in [('train',32),('validation',8)]:
        source=ROOT/f'artifacts/external/sota-homology-cold/data/{partition}.csv'
        with source.open(encoding='utf-8',newline='') as f:
            rows=list(csv.DictReader(f))[:count]
        target=OUT/f'{partition}.csv'
        with target.open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f);w.writerow(['pair_id','reactant_smiles','log10kcat_max','pdbpath','sequence'])
            for r in rows:
                name='seq_'+hashlib.sha256(r['sequence'].encode('ascii')).hexdigest()[:24]+'.pdb'
                records[name]={'name':name,'seq':r['sequence']}
                w.writerow([r['pair_id'],r['smiles'],r['log10_kcat'],name,r['sequence']])
        receipt['partitions'][partition]={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'rows':len(rows),'selection':'first rows in frozen order','sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
    with gzip.open(OUT/'protein_records.json.gz','wt',encoding='utf-8') as f:
        json.dump(records,f)
    receipt['unique_sequences']=len(records)
    (OUT/'receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
