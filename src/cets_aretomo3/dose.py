"""Acquisition order and dose bookkeeping for a raw (tilt-sorted) stack.

CETS stores only the exclusive pre-exposure per image; the per-image exposure and the acquisition order are
kept in the companion manifest. Nothing here assumes a constant dose or repeats an increment: the values are
taken from the source (``_TLT.txt`` / mdoc) or given explicitly.
"""

from typing import List, Optional, Sequence


def pre_exposure(acq_index_1b: Sequence[int], exposure: Sequence[float], convention: str = "exclusive") -> List[float]:
    """Accumulated dose before (``exclusive``) or including (``inclusive``) each raw image, from the per-image
    exposures and the 1-based acquisition order (per raw row)."""
    if convention not in ("exclusive", "inclusive"):
        raise ValueError(f"dose convention must be 'exclusive' or 'inclusive', got {convention!r}")
    n = len(acq_index_1b)
    if len(exposure) != n:
        raise ValueError(f"{len(exposure)} exposures for {n} rows")
    if sorted(acq_index_1b) != list(range(1, n + 1)):
        raise ValueError(f"acquisition indices are not a permutation of 1..{n}: {list(acq_index_1b)}")
    order = sorted(range(n), key=lambda i: acq_index_1b[i])
    out = [0.0] * n
    running = 0.0
    for i in order:
        if convention == "inclusive":
            running += float(exposure[i])
            out[i] = running
        else:
            out[i] = running
            running += float(exposure[i])
    return out


def constant_exposure(n: int, dose_per_tilt: Optional[float]) -> Optional[List[float]]:
    return None if dose_per_tilt is None else [float(dose_per_tilt)] * n
