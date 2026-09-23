"""Post-hoc paired cluster bootstrap, with fixed seeds averaged within each draw."""
from pathlib import Path
import csv
import hashlib
import json
import argparse
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/paired-family-analysis'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def rows(path):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))

def aligned(path, ids, labels, prediction_key, label_key):
    data = rows(path)
    mapping = {r['pair_id']: r for r in data}
    if len(mapping) != len(data) or set(mapping) != set(ids):
        raise ValueError(f'Pair IDs differ: {path}')
    y = np.array([float(mapping[i][label_key]) for i in ids])
    if not np.allclose(y, labels, rtol=0, atol=1e-6):
        raise ValueError(f'Labels differ: {path}')
    p = np.array([float(mapping[i][prediction_key]) for i in ids])
    if not np.isfinite(p).all():
        raise ValueError('Nonfinite prediction')
    return p

def metrics_by_draw(predictions, labels, groups, weights):
    """Return draw x metric; each row of weights counts complete families."""
    n_groups = weights.shape[1]
    result = []
    for prediction in predictions:
        row_stats = np.column_stack([np.ones(len(labels)), (prediction-labels)**2,
            abs(prediction-labels), prediction, labels, prediction**2,
            labels**2, prediction*labels])
        stats = np.zeros((n_groups, 8))
        np.add.at(stats, groups, row_stats)
        total = weights @ stats
        n, sq, ae, x, y, xx, yy, xy = total.T
        correlation = (xy-x*y/n) / np.sqrt(np.maximum(0, xx-x*x/n)*np.maximum(0, yy-y*y/n))
        result.append(np.column_stack([np.sqrt(sq/n), ae/n, correlation]))
    return np.mean(result, axis=0)

def main():
    global OUT
    parser = argparse.ArgumentParser()
    parser.add_argument('--current-code', action='store_true')
    args = parser.parse_args()
    if args.current_code:
        OUT = OUT/'current-code'
    protocol_path = ROOT/'configs/paired_family_analysis_protocol.json'
    protocol = json.loads(protocol_path.read_text(encoding='utf-8'))
    if args.current_code:
        amendment_path = ROOT/'configs/current_code_paired_analysis_protocol.json'
        amendment = json.loads(amendment_path.read_text(encoding='utf-8'))
        assert amendment['implementation_sha256'] == sha(__file__)
        for entry in amendment['frozen_files']:
            assert sha(ROOT/entry['path']) == entry['sha256']
        protocol['status'] = amendment['scope']
    for entry in protocol['frozen_files']:
        if sha(ROOT/entry['path']) != entry['sha256']:
            raise ValueError('Frozen input changed: '+entry['path'])
    data = rows(ROOT/'artifacts/external/sota-homology-cold/data/test.csv')
    ids = [r['pair_id'] for r in data]
    assert len(set(ids)) == len(ids)
    labels = np.array([float(r['log10_kcat']) for r in data])
    split = json.loads((ROOT/'artifacts/homology-final/homology_split.json').read_text())
    mapping = {}
    for r in split['rows']:
        if r['split'] == 'test':
            if r['sequence_id'] in mapping:
                assert mapping[r['sequence_id']] == r['cluster']
            mapping[r['sequence_id']] = r['cluster']
    families = [mapping['seq_'+hashlib.sha256(r['sequence'].encode('ascii')).hexdigest()[:24]] for r in data]
    unique, groups = np.unique(families, return_inverse=True)
    g = len(unique)
    rng = np.random.default_rng(protocol['bootstrap_seed'])
    weights = rng.multinomial(g, np.full(g, 1/g), size=protocol['bootstrap_replicates'])
    all_weights = np.vstack([np.ones(g), weights])
    metrics = {}; inputs = []
    for model in ['reference', 'CataPro', 'UniKP']:
        predictions = []
        for seed in protocol['seeds']:
            if model == 'reference':
                path = OUT/f'reference-seed{seed}.csv'; pk, lk = 'prediction', 'label'
                receipt = json.loads((OUT/f'reference-seed{seed}-receipt.json').read_text())
                assert receipt['predictions_sha256'] == sha(path)
            elif model == 'CataPro':
                path = ROOT/f'artifacts/catapro-mode-b/internal-test/seed{seed}-predictions.csv'
                pk, lk = 'prediction_log10_kcat', 'label_log10_kcat'
            else:
                path = ROOT/f'artifacts/external/sota-homology-cold/unikp-mode-b/seed{seed}/predictions.csv'
                pk, lk = 'log10_kcat_predicted', 'log10_kcat_observed'
            predictions.append(aligned(path, ids, labels, pk, lk))
            inputs.append({'path': str(path.relative_to(ROOT)), 'sha256': sha(path)})
        metrics[model] = metrics_by_draw(predictions, labels, groups, all_weights)
        assert np.isfinite(metrics[model]).all()
    contrasts = {}
    for comparator in ['CataPro', 'UniKP']:
        delta = metrics['reference'] - metrics[comparator]
        contrasts[comparator] = {name: {'difference': float(delta[0,j]),
            'ci95': np.quantile(delta[1:,j], [.025,.975]).tolist()}
            for j, name in enumerate(['rmse','mae','pearson'])}
    report = {'scope': protocol['status'], 'protocol_sha256': sha(protocol_path),
        'implementation_sha256': sha(__file__), 'rows': len(data), 'families': g,
        'largest_family_rows': int(np.bincount(groups).max()), 'replicates': len(weights),
        'metrics': {k: dict(zip(['rmse','mae','pearson'], v[0].tolist())) for k,v in metrics.items()},
        'reference_minus_comparator': contrasts, 'inputs': inputs,
        'limitations': 'Conditional on fixed checkpoints/seeds and observed internal split; post-hoc, no multiplicity adjustment; CI crossing zero does not establish equivalence.'}
    if args.current_code:
        report['amendment_sha256'] = sha(amendment_path)
    target = OUT/'paired-bootstrap.json'
    if target.exists():
        raise FileExistsError('Do not overwrite a completed analysis')
    target.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='inputs'}, indent=2))

if __name__ == '__main__':
    main()
