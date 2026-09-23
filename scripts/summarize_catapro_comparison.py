"""Audit saved predictions and append CataPro to a new, non-overwriting comparison."""
from pathlib import Path
import csv,json,hashlib,statistics,sys
import torch

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from biocandidate.evaluation import regression_metrics

def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    out=ROOT/'artifacts/catapro-mode-b/internal-test'
    result=read(out/'summary.json')
    freeze=read(ROOT/'configs/catapro_internal_test_freeze.json')
    assert sha(ROOT/'configs/catapro_internal_test_freeze.json')==result['freeze_sha256']
    for entry in freeze['files']:assert sha(ROOT/entry['path'])==entry['sha256'],entry['path']
    with (ROOT/'artifacts/external/sota-homology-cold/data/test.csv').open(encoding='utf-8',newline='') as f:reference=list(csv.DictReader(f))
    for entry in result['prediction_files']:
        p=ROOT/entry['path'];assert sha(p)==entry['sha256']
        seed=p.name.split('-')[0].removeprefix('seed')
        with p.open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
        assert [r['pair_id'] for r in rows]==[r['pair_id'] for r in reference]
        labels=torch.tensor([float(r['label_log10_kcat']) for r in rows])
        assert torch.equal(labels,torch.tensor([float(r['log10_kcat']) for r in reference]))
        metrics=regression_metrics(torch.tensor([float(r['prediction_log10_kcat']) for r in rows]),labels)
        for key in ['rmse','mae','pearson']:assert abs(metrics[key]-result['seeds'][seed][key])<1e-6
    old_path=ROOT/'artifacts/external/sota-homology-cold/summary.json';old=read(old_path)
    rows=old['rows']+[{'model':'CataPro fixed-split retraining','seeds':[7,42,123],**result['mean_sd']}]
    combined={'generated_on':'2026-09-18','old_summary_sha256':sha(old_path),'catapro_summary_sha256':sha(out/'summary.json'),
              'test_rows':1646,'rows':rows,'limitations':['Same data split; model capacities, input modalities and selection rules differ.',
              'Seed SD is not a confidence interval for test-sample uncertainty.','Internal previously observed test only.']}
    (out/'comparison.json').write_text(json.dumps(combined,indent=2),encoding='utf-8')
    lines=['# CataPro 现代基线实测结果','', '2026-09-18。完整特征缓存、三种子固定划分训练及单次内部测试已完成。','',
           '| 方法 | RMSE | MAE | Pearson |','|---|---|---|---|']
    for row in rows:
        lines.append('| '+row['model']+' | '+' | '.join(f"{row[k]['mean']:.4f} ± {row[k]['sd']:.4f}" for k in ['rmse','mae','pearson'])+' |')
    lines += ['', '**判断：CataPro 的平均 RMSE 和 MAE 数值略低、Pearson 数值更高；本项目尚未证明优于该现代基线。** RMSE 仅相差约 0.0040，不能据三种子标准差认定统计显著。所有结果为原有 1,646 条内部测试集，均不是独立外部验证。', '',
              '训练集 13,157 条、验证集 1,640 条。CataPro 采用官方预测头与默认优化参数，并用固定验证集早停，替代官方十折训练；UniKP 与 DLKcat 按其既有冻结协议训练。模型容量、输入模态及模型选择方式不同，不能把这个比较解释为单独检验融合结构优劣。','',
              'CataPro 三种子训练分别运行 28、23、31 轮，最佳 epoch（从 0 开始）为 7、2、10。总计约 72 秒为缓存特征后的预测头训练时间，不含编码器下载或从零蛋白表示提取。','',
              '复用 ProtT5 缓存经过权重、预处理和短／中／长序列数值检查；MolT5 与原版 CPU 提取器及 MACCS 也通过检查。7 项针对缓存、数据分区和冻结文件的测试通过。保存预测的行顺序、标签、哈希和重算指标均核验一致。','',
              '下一步应补 CatPred 并进行按家族的配对差值分析；在没有新证据前，不因看到测试结果而调整 CataPro 或本项目的超参数。论文重点仍需落在可独立复验的生物学问题与泛化边界。','',
              '记录：artifacts/catapro-mode-b/internal-test/summary.json、comparison.json、seed*-predictions.csv；冻结协议：configs/catapro_internal_test_freeze.json；运行说明：docs/CATAPRO_MODE_B_RUNBOOK_CN.md。']
    (ROOT/'docs/CATAPRO_INTERNAL_COMPARISON_CN.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('\n'.join(lines[:11]));print('Three saved prediction files verified; no inference performed.')

if __name__=='__main__':main()
