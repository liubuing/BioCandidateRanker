"""Snapshot upstream metadata and small source files without executing upstream code."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import urllib.request
import base64

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/pcb-readiness'

def fetch(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'BioCandidateRanker-source-audit'})
    with urllib.request.urlopen(req, timeout=10) as response:
        return response.read()

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    receipt = {'retrieved_at': datetime.now(timezone.utc).isoformat(), 'repositories': []}
    for repo, branch in [('maranasgroup/CatPred', 'main'), ('zchwang/CataPro', 'master')]:
        item = {'repository': repo, 'branch': branch, 'files': [], 'errors': []}
        receipt['repositories'].append(item)
        try:
            commit = json.loads(fetch(f'https://api.github.com/repos/{repo}/commits/{branch}'))
            sha = commit['sha']; item['commit'] = sha
            tree = json.loads(fetch(f'https://api.github.com/repos/{repo}/git/trees/{sha}?recursive=1'))
            dest = OUT / repo.split('/')[1]; dest.mkdir(exist_ok=True)
            (dest/'tree.json').write_text(json.dumps(tree, indent=2), encoding='utf-8')
            for entry in tree['tree']:
                path = entry['path']; name = Path(path).name.lower()
                include = (('/' not in path and name.startswith(('readme','license','requirements','environment'))) or
                           path in ['setup.py', 'predict.py', 'train.py', 'reproduce_training.sh',
                                    'inference/predict.py', 'inference/act_model.py', 'inference/run_catapro.sh'] or
                           (repo.endswith('CataPro') and path.endswith('.py')))
                if entry['type'] != 'blob' or not include or entry.get('size', 0) > 150000:
                    continue
                url = f'https://raw.githubusercontent.com/{repo}/{sha}/{path}'
                try:
                    blob = json.loads(fetch(entry['url']))
                    data = base64.b64decode(blob['content'])
                    target = dest/'source'/path; target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(data)
                    item['files'].append({'path': path, 'url': url, 'bytes': len(data),
                                          'sha256': hashlib.sha256(data).hexdigest()})
                except Exception as exc:
                    item['errors'].append({'path': path, 'error': str(exc)})
                (OUT/'source-receipt.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
                print(repo, path, flush=True)
        except Exception as exc:
            item['errors'].append({'error': str(exc)})
        (OUT/'source-receipt.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
    print(json.dumps(receipt, indent=2))

if __name__ == '__main__':
    main()
