"""Install bounded-memory inputs and atomic recovery without editing upstream.

Patches catpred.train.run_training in memory (each source anchor must match
exactly once, or the run refuses to start) to call the atomic restore/save
hooks, then starts upstream training. Files on disk in the upstream tree are
never modified.

中文：本脚本在不改动上游源码的前提下安装有界内存输入与原子恢复：在内存中改接
catpred.train.run_training（每个源码锚点必须恰好匹配一次，否则运行拒绝启动）以调用原子
restore/save 钩子，然后启动上游训练。上游树中的磁盘文件从不被修改。治理角色：上游
训练入口适配器（注入确定性恢复钩子，配合 catpred_recovery 协议使用）。
"""
import importlib,inspect,textwrap
from catpred_memory_entry import install_adapters
import catpred_recovery

def install_recovery():
    """Rewrite run_training in memory to insert epoch-boundary recovery calls.

    中文：在内存中改写 run_training，插入 epoch 边界的恢复调用。
    """
    module=importlib.import_module('catpred.train.run_training')
    source=textwrap.dedent(inspect.getsource(module.run_training))
    old='        for epoch in trange(args.epochs):'
    new='        start_epoch, best_score, best_epoch, n_iter = _recover_restore(save_dir, model, optimizer, scheduler, train_data_loader)\n        for epoch in trange(start_epoch, args.epochs):'
    assert source.count(old)==1;source=source.replace(old,new)
    old='        # Evaluate on test set using model with best validation score'
    new='            _recover_save(save_dir, model, optimizer, scheduler, train_data_loader, epoch, best_score, best_epoch, n_iter)\n\n'+old
    assert source.count(old)==1;source=source.replace(old,new)
    module.__dict__.update(_recover_restore=catpred_recovery.restore,_recover_save=catpred_recovery.save)
    namespace={};exec(compile(source,__file__,'exec'),module.__dict__,namespace)
    module.run_training=namespace['run_training']
    cross=importlib.import_module('catpred.train.cross_validate');cross.run_training=module.run_training
    import catpred.train
    catpred.train.run_training=module.run_training

if __name__=='__main__':
    install_adapters();install_recovery()
    from catpred.train import catpred_train
    catpred_train()
