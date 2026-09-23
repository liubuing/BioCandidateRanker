"""Export only frozen training/validation data for CatPred and a length stress pilot."""
from pathlib import Path
import csv, gzip, hashlib, json

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/catpred-mode-b'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    source_rows={}
    for partition in ['train','validation']:
        source=ROOT/f'artifacts/external/sota-homology-cold/data/{partition}.csv'
        with source.open(encoding='utf-8',newline='') as f:source_rows[partition]=list(csv.DictReader(f))
    for name in ['development-data','length-pilot-data']:
        out=BASE/name
        out.mkdir(exist_ok=True)
        if (out/'receipt.json').exists():raise FileExistsError(out/'receipt.json')
        records={}; receipt={'purpose':name,'partitions':{}}
        for partition,all_rows in source_rows.items():
            selected=all_rows if name=='development-data' else sorted(all_rows,key=lambda r:(-len(r['sequence']),r['pair_id']))[:32 if partition=='train' else 16]
            target=out/f'{partition}.csv'
            with target.open('w',encoding='utf-8',newline='') as f:
                w=csv.writer(f);w.writerow(['pair_id','reactant_smiles','log10kcat_max','pdbpath','sequence'])
                for r in selected:
                    key='seq_'+hashlib.sha256(r['sequence'].encode('ascii')).hexdigest()[:24]+'.pdb'
                    records[key]={'name':key,'seq':r['sequence']}
                    w.writerow([r['pair_id'],r['smiles'],r['log10_kcat'],key,r['sequence']])
            receipt['partitions'][partition]={'rows':len(selected),'sha256':sha(target),'source_sha256':sha(ROOT/f'artifacts/external/sota-homology-cold/data/{partition}.csv'),'max_length':max(len(r['sequence']) for r in selected),'rows_over_2047':sum(len(r['sequence'])>2047 for r in selected)}
        with gzip.open(out/'protein_records.json.gz','wt',encoding='utf-8') as f:json.dump(records,f)
        receipt['records_sha256']=sha(out/'protein_records.json.gz')
        receipt['unique_sequences']=len(records)
        receipt['estimated_float32_esm_cache_bytes']=sum(min(len(r['seq']),2047)*1280*4 for r in records.values())
        (out/'receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
        print(json.dumps(receipt))

if __name__=='__main__':main()
