"""BioCandidateRanker public API.

Stable import surface for downstream consumers: the multimodal ranker and its
configuration dataclass. Training, evaluation, and the governance/fail-closed
pipelines are reached through :mod:`biocandidate.cli` and the dedicated audit
modules rather than through the package root.

中文：BioCandidateRanker 的公开 API。为下游使用者提供稳定的导入入口：多模态排序器
及其配置 dataclass。训练、评估以及治理/fail-closed 流水线需通过
:mod:`biocandidate.cli` 和各专用审计模块访问，而非包根目录。
"""

from .config import ModelConfig
from .model.ranker import BioCandidateRanker

__all__ = ["BioCandidateRanker", "ModelConfig"]
__version__ = "0.1.0"
