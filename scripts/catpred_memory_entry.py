"""Read verified CPU ESM cache tensors through mmap; preserve model numerics.

Installs two runtime adapters before upstream training: the noncanonical
residue tokenizer adapter and an mmap + weights_only loader that only accepts
ESM cache entries resolved under the verified cache root. The upstream source
tree is never edited; behavior differs only in how tensors are loaded.

中文：本脚本在上游训练前安装两个运行时适配器：非规范残基 tokenizer 适配器，以及仅接受
解析到已验证缓存根之下 ESM 缓存条目的 mmap + weights_only 加载器，从而以只读 mmap 方式
读取已验证的 CPU ESM 缓存张量并保持模型数值不变。上游源码树从不被修改；行为差异仅在
张量的加载方式。治理角色：上游训练入口适配器（内存受控、数值保持、需披露）。
"""
import catpred_noncanonical_entry

def install_memory_adapter():
    """Patch catpred's torch loader for mmap CPU loads scoped to the ESM cache root.

    中文：改接 catpred 的 torch 加载器，把 mmap CPU 加载限定在 ESM 缓存根内。
    """
    import torch
    import catpred.data.cache_utils as cache
    from catpred.security.deserialization import ensure_trusted_path
    original=cache.load_torch_artifact
    root=(cache.CACHE_PATH/'esm/proteins').resolve()

    def load(path, *, purpose, map_location=None, roots=None, allow_unsafe=None):
        if purpose!='esm cache entry' or map_location!='cpu':
            return original(path,purpose=purpose,map_location=map_location,roots=roots,allow_unsafe=allow_unsafe)
        resolved=ensure_trusted_path(path,purpose=purpose,roots=roots)
        resolved.relative_to(root)
        return torch.load(str(resolved),map_location='cpu',weights_only=True,mmap=True)

    cache.load_torch_artifact=load

def install_adapters():
    """Install both the noncanonical-tokenizer and bounded-memory adapters.

    中文：同时安装非规范 tokenizer 与有界内存两个适配器。
    """
    catpred_noncanonical_entry.install_adapter()
    install_memory_adapter()

if __name__=='__main__':
    install_adapters()
    from catpred.train import catpred_train
    catpred_train()
