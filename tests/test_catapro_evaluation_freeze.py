"""Changing a sealed checkpoint must fail before test prediction."""
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
import evaluate_catapro_once as evaluator

def test_sealed_checkpoint_change_rejected(tmp_path):
    p=tmp_path/'checkpoint.pt';p.write_bytes(b'original')
    entries=[{'path':p.name,'sha256':evaluator.run.sha(p)}]
    evaluator.verify_files(entries,tmp_path)
    p.write_bytes(b'changed')
    with pytest.raises(ValueError,match='Frozen file changed'):
        evaluator.verify_files(entries,tmp_path)

def test_missing_sealed_file_rejected(tmp_path):
    with pytest.raises(FileNotFoundError):
        evaluator.verify_files([{'path':'missing.pt','sha256':'0'*64}],tmp_path)
