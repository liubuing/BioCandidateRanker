"""Download pinned official MolT5 files into this workspace, recording hashes."""
from pathlib import Path
import hashlib,json,os
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/catapro-mode-b/models/molt5-base-smiles2caption'
os.environ['HF_HOME']=str(ROOT/'artifacts/catapro-mode-b/hf-cache')
os.environ['HF_HUB_DISABLE_XET']='1'
from huggingface_hub import snapshot_download

revision='7b7d4b0ab8b66b351e669b1f66272418ba15c3d9'
snapshot_download('laituan245/molt5-base-smiles2caption',revision=revision,local_dir=OUT,
                  allow_patterns=['README.md','config.json','spiece.model',
                                  'special_tokens_map.json','tokenizer_config.json','tokenizer.json'],max_workers=2)
weight=OUT/'pytorch_model.bin'
expected='111fc277b175d7b0586905424448e20f1dce1dc1b1e912c6e447921b7d8621b1'
if not weight.exists():
    partial=OUT/'pytorch_model.bin.partial'
    url=f'https://huggingface.co/laituan245/molt5-base-smiles2caption/resolve/{revision}/pytorch_model.bin'
    with requests.get(url,stream=True,timeout=(15,60)) as response:
        response.raise_for_status();count=0
        with partial.open('wb') as f:
            for chunk in response.iter_content(8*1024*1024):
                f.write(chunk);count+=len(chunk)
                if count%(128*1024*1024)==0:print(f'Downloaded {count//(1024*1024)} MiB',flush=True)
    h=hashlib.sha256()
    with partial.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    if h.hexdigest()!=expected:raise ValueError('Official LFS weight SHA256 mismatch')
    partial.replace(weight)
files={}
for p in OUT.iterdir():
    if not p.is_file() or p.name=='download-receipt.json' or p.name.endswith('.partial'):continue
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    files[p.name]={'sha256':h.hexdigest(),'bytes':p.stat().st_size}
if files['pytorch_model.bin']['sha256']!=expected:raise ValueError('Existing weight SHA256 mismatch')
receipt={'repository':'laituan245/molt5-base-smiles2caption','revision':revision,'files':files}
(OUT/'download-receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt,indent=2))
