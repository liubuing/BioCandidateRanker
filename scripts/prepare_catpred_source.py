"""Fetch pinned upstream archive; extract only reviewed-scope code/config files."""
from pathlib import Path,PurePosixPath
import hashlib,json,urllib.request,zipfile
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'artifacts/catpred-mode-b';OUT.mkdir(exist_ok=True)
commit='b314a28a84d388755237de90d9f336f951b291b3'
archive=OUT/'upstream.zip'
if not archive.exists():
    with urllib.request.urlopen(f'https://codeload.github.com/maranasgroup/CatPred/zip/{commit}',timeout=45) as r, archive.with_suffix('.partial').open('wb') as f:
        while chunk:=r.read(1024*1024):f.write(chunk)
    archive.with_suffix('.partial').replace(archive)
tree=json.loads((ROOT/'artifacts/pcb-readiness/CatPred/tree.json').read_text())
assert tree['sha']==commit and not tree.get('truncated',False)
index={r['path']:r for r in tree['tree'] if r['type']=='blob'}
receipt={'commit':commit,'files':[]}
with zipfile.ZipFile(archive) as z:
    for entry in z.infolist():
        if entry.is_dir():continue
        parts=PurePosixPath(entry.filename).parts[1:]
        if '..' in parts:raise ValueError('Unsafe archive path')
        path='/'.join(parts)
        selected=(path.startswith(('catpred/','scripts/')) and path.endswith(('.py','.sh'))) or path in ['train.py','predict.py','setup.py','setup.cfg','LICENSE','LICENSE.txt','README.md','environment.yml']
        if not selected:continue
        data=z.read(entry)
        assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==index[path]['sha'],path
        dest=OUT/'upstream'/path;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
        receipt['files'].append({'path':path,'sha256':hashlib.sha256(data).hexdigest()})
(OUT/'source-receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print('Verified source files:',len(receipt['files']))
