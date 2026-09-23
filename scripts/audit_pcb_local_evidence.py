"""Read-only dataset audit; never predicts or opens temporal-pool labels."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'artifacts/pcb-readiness'
CORPUS = Path(r'C:\biological\Metabolic model prediction\Integrated_Yeast_MetaTwin_Deployment\04_prediction_plugins\UniKP\datasets\Kcat_combination_0918_wildtype_mutant.json')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source_manifest=read('configs/unikp_source_manifest.json')
    assert digest(CORPUS)==source_manifest['identity']['sha256'], 'source identity changed'
    raw=json.loads(CORPUS.read_text(encoding='utf-8'))
    keys=sorted(set().union(*(r.keys() for r in raw)))
    coverage={key:sum(r.get(key) not in (None,'') for r in raw) for key in keys}
    units=Counter(str(r.get('Unit','<absent>')) for r in raw)
    # Metadata conflict counts are raw exact-string pairs, not canonicalized pairs.
    groups={}
    for row in raw:
        groups.setdefault((row.get('Sequence'),row.get('Smiles')), []).append(row)
    conflicts={field:sum(len({str(r.get(field,'')) for r in rows})>1 for rows in groups.values())
               for field in ['Organism','ECNumber','Type']}
    manifest=read('artifacts/external/sota-homology-cold/data-manifest.json')
    split_path=ROOT/manifest['split_manifest']['path']
    assert digest(split_path)==manifest['split_manifest']['sha256'], 'split identity changed'
    partitions={}; sets={}
    for name, entry in manifest['partitions'].items():
        path=ROOT/entry['path']; assert digest(path)==entry['sha256'], f'{name} changed'
        with path.open(encoding='utf-8',newline='') as f:
            reader=csv.DictReader(f); rows=list(reader); fields=reader.fieldnames
        assert len(rows)==entry['rows']
        partitions[name]={'rows':len(rows),'sha256':entry['sha256'],'columns':fields,
                          'unique_sequences':len({r['sequence'] for r in rows}),
                          'unique_smiles':len({r['smiles'] for r in rows})}
        sets[name]={'sequence':{r['sequence'] for r in rows},'smiles':{r['smiles'] for r in rows},
                    'pair':{(r['sequence'],r['smiles']) for r in rows}}
    overlaps={}
    for a,b in [('train','validation'),('train','test'),('validation','test')]:
        overlaps[f'{a}__{b}']={k:len(sets[a][k]&sets[b][k]) for k in sets[a]}
    temporal=read('artifacts/external/temporal-global-family-audit/readiness-audit.json')
    selection=read('artifacts/external/enzengdb-v1/selection-homology-cold.json')
    package_versions={}
    for package in ['torch','transformers','rdkit','numpy','pandas','scikit-learn','sentencepiece','fair-esm','torch-scatter','descriptastorus','progres']:
        try: package_versions[package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError: package_versions[package]=None
    gpu=subprocess.run(['nvidia-smi','--query-gpu=name,memory.total,driver_version','--format=csv,noheader'],capture_output=True,text=True)
    result={'generated_at':datetime.now(timezone.utc).isoformat(),'predictions_generated':False,
            'corpus':{'sha256':digest(CORPUS),'raw_rows':len(raw),'field_nonempty_counts':coverage,
                      'raw_unit_counts':dict(units),'raw_pair_metadata_conflict_groups':conflicts,
                      'conflict_grouping':'raw exact Sequence/Smiles strings; includes rows later rejected',
                      'license_status':source_manifest['license_status'],
                      'provenance':source_manifest['upstream_record_provenance']},
            'partitions':partitions,'exact_unique_overlaps':overlaps,
            'context_contract':{'inputs':['organism','ec','enzyme_type','reaction'],
                                'unikp_reaction':'default empty; adapter does not populate',
                                'pH_temperature':'not modeled by this adapter/collator',
                                'aggregation':'median log10 target; metadata retained from first source row'},
            'temporal_readiness':{'counts':temporal['counts'],'blockers':temporal['blockers']},
            'enzengdb_selected':{'campaigns':len(selection['campaign_counts']),
                                'records':sum(selection['campaign_counts'].values())},
            'runtime':{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
                       'packages':package_versions,'gpu':gpu.stdout.strip()}}
    (OUT/'local-evidence-audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    entries=[
      ['development_kcat','16443 aggregated pairs','log10 kcat','s^-1 before log10','internal already observed','corpus only','redistribution not verified','absolute regression; metadata diagnostic'],
      ['EnzEngDB cold','6423 rows / 51 campaigns','campaign fitness','campaign-specific','already observed','dataset DOI and campaign records; primary assay review needed','CC-BY-4.0 per existing audit','within-campaign ranking only'],
      ['IMDH','512 genotypes / 1024 cofactor observations','derived ln kcat','absolute source-unit constant unresolved','already observed; homologous','Dryad 10.5061/dryad.7nd70','CC0 per existing audit','within-cofactor ranking only'],
      ['temporal absolute pool','192 rows / 25 families / 54 substrates','absolute kinetics; primary kcat','governed canonical units','unscored; gate failed','source-specific governed evidence','governed per-source gate','curation only until frozen gate passes'],
      ['prospective ranking','no independent roster admitted','campaign-specific ranking','must be preregistered','pending external data','custodian roster required','pending','not yet evaluable']]
    with (OUT/'endpoint-register.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['dataset','size','endpoint','unit','evaluation_status','provenance','license','permitted_use']);w.writerows(entries)
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
