"""BioCandidateRanker public API.

Stable import surface for downstream consumers: the multimodal ranker and its
configuration dataclass. Training, evaluation, and the governance/fail-closed
pipelines are reached through :mod:`biocandidate.cli` and the dedicated audit
modules rather than through the package root.
"""

from .config import ModelConfig
from .model.ranker import BioCandidateRanker

__all__ = ["BioCandidateRanker", "ModelConfig"]
__version__ = "0.1.0"
