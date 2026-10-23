"""Install bounded-memory inputs and atomic recovery without editing upstream.

Patches catpred.train.run_training in memory (each source anchor must match
exactly once, or the run refuses to start) to call the atomic restore/save
hooks, then starts upstream training. Files on disk in the upstream tree are
never modified.
"""
import importlib,inspect,textwrap
from catpred_memory_entry import install_adapters
import catpred_recovery

def install_recovery():
    """Rewrite run_training in memory to insert epoch-boundary recovery calls."""
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
