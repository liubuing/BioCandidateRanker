"""Read verified CPU ESM cache tensors through mmap; preserve model numerics."""
from pathlib import Path
import catpred_noncanonical_entry

def install_memory_adapter():
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
    catpred_noncanonical_entry.install_adapter()
    install_memory_adapter()

if __name__=='__main__':
    install_adapters()
    from catpred.train import catpred_train
    catpred_train()
