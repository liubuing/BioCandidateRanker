"""CatPred adapter: noncanonical residues use existing zero embedding index 20.

ESM features and raw sequence remain unchanged. Canonical residues are unchanged.
The upstream source tree is not modified. This adaptation must be disclosed.
"""
import inspect
import textwrap

def patch_source(source):
    original='[letter_to_num[a] for a in seq]'
    replacement='[letter_to_num.get(a, 20) for a in seq]'
    if source.count(original)!=1:
        raise ValueError('Unexpected upstream tokenizer; refusing patch')
    return source.replace(original,replacement)

def install_adapter():
    import catpred.models.model as module
    source=textwrap.dedent(inspect.getsource(module.MoleculeModel.forward))
    namespace={}
    exec(compile(patch_source(source),__file__,'exec'),module.__dict__,namespace)
    module.MoleculeModel.forward=namespace['forward']

if __name__=='__main__':
    install_adapter()
    from catpred.train import catpred_train
    catpred_train()
