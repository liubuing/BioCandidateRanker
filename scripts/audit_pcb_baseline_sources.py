"""Snapshot upstream metadata and small source files without executing upstream code.

Records the default-branch commit and full Git tree for maranasgroup/CatPred
and zchwang/CataPro, downloads only allow-listed text files up to 150 KB
(README/license/requirements, entry points, CataPro *.py), and stores them
under artifacts/pcb-readiness with every URL and SHA256 in
source-receipt.json. Upstream code is never executed; errors are recorded
per file instead of aborting the snapshot.

中文：本脚本为上游元数据与小型源文件拍快照而不执行上游代码：记录 maranasgroup/CatPred
与 zchwang/CataPro 默认分支的 commit 及完整 Git tree，只下载允许清单内不超过 150 KB 的
文本文件（README/license/requirements、入口脚本、CataPro 的全部 *.py），存入
artifacts/pcb-readiness，并把每个 URL 与 SHA256 记入 source-receipt.json。上游代码绝不
执行；错误按文件记录而非中止快照。治理角色：来源快照审计收据（只读取证，不跑上游）。
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import urllib.request
import base64

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/pcb-readiness'

def fetch(url):
    """Fetch a URL with the audit user agent and return the raw response bytes.

    中文：以审计用 User-Agent 抓取 URL 并返回原始响应字节。
    """
    req = urllib.request.Request(url, headers={'User-Agent': 'BioCandidateRanker-source-audit'})
    with urllib.request.urlopen(req, timeout=10) as response:
        return response.read()

def main():
    """Snapshot both upstream repositories and incrementally write source-receipt.json.

    中文：为两个上游仓库拍快照并增量写出 source-receipt.json。
    """
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
