"""CatPred adapter: noncanonical residues use existing zero embedding index 20.

ESM features and raw sequence remain unchanged. Canonical residues are unchanged.
The upstream source tree is not modified. This adaptation must be disclosed.

中文：本脚本是 CatPred 适配器：非规范残基复用既有的零嵌入索引 20。ESM 特征与原始序列
保持不变；规范残基不受影响；上游源码树不被修改。治理角色：上游训练入口适配器——
tokenizer 语义适配，此适配必须披露。
"""
import inspect
import textwrap

def patch_source(source):
    """Rewrite the tokenizer to fall back to index 20, refusing any other source text.

    中文：把 tokenizer 改写为回退到索引 20；源文本与此不符即拒绝改写。
    """
    original='[letter_to_num[a] for a in seq]'
    replacement='[letter_to_num.get(a, 20) for a in seq]'
    if source.count(original)!=1:
        raise ValueError('Unexpected upstream tokenizer; refusing patch')
    return source.replace(original,replacement)

def install_adapter():
    """Recompile MoleculeModel.forward with the patched tokenizer on the imported module.

    中文：在已导入的模块上用改写后的 tokenizer 重新编译 MoleculeModel.forward。
    """
    import catpred.models.model as module
    source=textwrap.dedent(inspect.getsource(module.MoleculeModel.forward))
    namespace={}
    exec(compile(patch_source(source),__file__,'exec'),module.__dict__,namespace)
    module.MoleculeModel.forward=namespace['forward']

if __name__=='__main__':
    install_adapter()
    from catpred.train import catpred_train
    catpred_train()
