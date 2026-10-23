"""Read verified CPU ESM cache tensors through mmap; preserve model numerics.

Installs two runtime adapters before upstream training: the noncanonical
residue tokenizer adapter and an mmap + weights_only loader that only accepts
ESM cache entries resolved under the verified cache root. The upstream source
tree is never edited; behavior differs only in how tensors are loaded.
"""
import catpred_noncanonical_entry

def install_memory_adapter():
    """Patch catpred's torch loader for mmap CPU loads scoped to the ESM cache root."""
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
    """Install both the noncanonical-tokenizer and bounded-memory adapters."""
    catpred_noncanonical_entry.install_adapter()
    install_memory_adapter()

if __name__=='__main__':
    install_adapters()
    from catpred.train import catpred_train
    catpred_train()
