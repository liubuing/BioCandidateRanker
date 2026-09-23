"""Guard the new baseline against row-order, cache, and partition contamination."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('catapro_runner',Path(__file__).parents[1]/'scripts/run_catapro_mode_b.py')
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)

def fixture_cache(tmp_path,monkeypatch):
    monkeypatch.setattr(runner,'ROOT',tmp_path)
    monkeypatch.setattr(runner,'OUT',tmp_path)
    protocol=tmp_path/'protocol.json';protocol.write_text('{}')
    monkeypatch.setattr(runner,'PROTOCOL',protocol)
    parts={'train':[{'pair_id':'a','log10_kcat':'1'}],'validation':[{'pair_id':'b','log10_kcat':'2'}]}
    monkeypatch.setattr(runner,'load_partitions',lambda limit:parts)
    dest=tmp_path/'smoke';dest.mkdir()
    manifest={'limit':1,'identity':{'protocol_sha256':runner.sha(protocol)},'partitions':{}}
    for name,rows in parts.items():
        file=dest/f'{name}.npy';np.save(file,np.ones((1,1959),dtype=np.float32))
        manifest['partitions'][name]={'row_ids':[r['pair_id'] for r in rows],
            'path':str(file.relative_to(tmp_path)),'sha256':runner.sha(file)}
    runner.save_json(dest/'feature-manifest.json',manifest)
    return manifest,dest

def test_changed_row_identity_rejected(tmp_path,monkeypatch):
    manifest,dest=fixture_cache(tmp_path,monkeypatch)
    manifest['partitions']['train']['row_ids']=['different']
    runner.save_json(dest/'feature-manifest.json',manifest)
    with pytest.raises(ValueError,match='ordering'):runner.load_feature_data('smoke')

def test_changed_feature_bytes_rejected(tmp_path,monkeypatch):
    _,dest=fixture_cache(tmp_path,monkeypatch)
    np.save(dest/'train.npy',np.zeros((1,1959),dtype=np.float32))
    with pytest.raises(ValueError,match='Features changed'):runner.load_feature_data('smoke')

def test_protocol_change_rejected(tmp_path,monkeypatch):
    fixture_cache(tmp_path,monkeypatch)
    runner.PROTOCOL.write_text('{"changed":true}')
    with pytest.raises(ValueError,match='Protocol changed'):runner.load_feature_data('smoke')

def test_only_train_validation_requested(tmp_path,monkeypatch):
    monkeypatch.setattr(runner,'ROOT',tmp_path)
    manifest={'partitions':{'test':{'path':'DO_NOT_OPEN_TEST.csv'}}}
    for name in ['train','validation']:
        p=tmp_path/f'{name}.csv';p.write_text('pair_id,log10_kcat\n'+name+',1\n')
        manifest['partitions'][name]={'path':p.name,'sha256':runner.sha(p),'rows':1}
    path=tmp_path/'manifest.json';runner.save_json(path,manifest)
    protocol=tmp_path/'protocol.json'
    runner.save_json(protocol,{'data_manifest':path.name,'data_manifest_sha256':runner.sha(path)})
    monkeypatch.setattr(runner,'PROTOCOL',protocol)
    assert set(runner.load_partitions())=={'train','validation'}

def test_sequence_truncation_and_rare_residues():
    seq='A'*500+'M'*50+'U'*500
    prepared=runner.sequence_text(seq).split()
    assert prepared==list('A'*500+'X'*500)
