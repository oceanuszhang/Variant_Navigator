"""ESM2 log-likelihood scoring of a variant, executed on Modal.

STUB — owned by the ESM2/Modal workstream. The contract below is what the rest
of the pipeline calls; keep the dataclass fields stable and the internals are
free to change.

Interpretation: a strongly negative `delta_ll` means the variant is unlikely
under the language model, i.e. structurally disruptive. Resistance mutations are
often mildly negative (tolerated but perturbing) rather than catastrophic, so
treat this as one feature, not a verdict.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field


@dataclass
class VariantScore:
    variant: str
    wt_log_likelihood: float
    mt_log_likelihood: float
    delta_ll: float                       # mt - wt; negative = disfavoured
    per_position: dict[int, float] = field(default_factory=dict)
    model: str = "esm2_t33_650M_UR50D"
    backend: str = "modal"

    @property
    def disruptive(self) -> bool:
        return self.delta_ll < -1.0


def score_variant(
    wildtype: str,
    mutant: str,
    variant: str,
    *,
    model: str = "esm2_t33_650M_UR50D",
) -> VariantScore:
    """Score one variant. Raises until the Modal endpoint is wired up."""
    raise NotImplementedError(
        "ESM2 scoring not wired yet. Set MODAL_TOKEN_ID/MODAL_TOKEN_SECRET and "
        "implement the Modal call here; see pipeline/esm2_score.py docstring."
    )


def available() -> bool:
    """True when credentials are present, so main.py can degrade gracefully."""
    return bool(os.getenv("MODAL_TOKEN_ID") and os.getenv("MODAL_TOKEN_SECRET"))
