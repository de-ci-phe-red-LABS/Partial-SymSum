"""
Degree attainment: the test vector only works if the chosen derivative-variable
set drives the permutation to its maximum algebraic degree d = 2**rounds on the
observed bits.

Checks provided:

  * ``derivative_nonzero`` - the d-fold derivative is not all zeros,
    i.e. the set *reaches* degree d.
  * ``reaches_max_degree`` - d-fold derivative != 0 AND (d+1)-fold == 0,
    the clean two-sided certificate of "exactly degree d".
  * ``find_degree_attaining_set`` - search symmetric directions for a size-d set
    that reaches degree d.

If a set does NOT reach degree d, the d-fold derivative collapses to all zeros
(a trivial ZeroSum) and no SymSum / SymSum_part behaviour can be observed.
"""

from __future__ import annotations
import random

from core.derivatives import kfold_derivative
from core.symmetry import symmetric_direction, is_zero


def direction_pool(cipher):
    """All symmetric directions (one per mirrored bit-pair of every lane)."""
    return [symmetric_direction(cipher, lane, z)
            for lane in range(cipher.n_lanes)
            for z in range(cipher.half)]


def word_pool(cipher, lane):
    """Symmetric directions confined to a single lane/word.  For slowly
    diffusing primitives (Bash-f) a degree-attaining set must be drawn from one
    word or one vertical plane -- a set drawn from the whole state almost never
    reaches the full degree."""
    return [symmetric_direction(cipher, lane, z) for z in range(cipher.half)]


def derivative_nonzero(cipher, base, directions, rounds):
    return not is_zero(kfold_derivative(cipher, base, directions, rounds))


def reaches_max_degree(cipher, base, directions, rounds, pool=None, rng=None):
    """True iff the len(directions)-fold derivative is nonzero and adding any
    one more direction makes it vanish (degree is exactly len(directions))."""
    if is_zero(kfold_derivative(cipher, base, directions, rounds)):
        return False
    rng = rng or random
    pool = pool or direction_pool(cipher)
    used = set(directions)
    extra = next(d for d in pool if d not in used)
    return is_zero(kfold_derivative(cipher, base, directions + [extra], rounds))


def build_max_degree_set(cipher, base, rounds, cap, pool=None, rng=None):
    """Greedily grow a symmetric direction set, always keeping the running
    derivative non-vanishing, until no direction extends it or ``cap`` is
    reached.  ``len(result)`` is the attained algebraic degree on the symmetric
    subspace (which may be below the bound 2**rounds, as for low-round Bash)."""
    rng = rng or random.Random(0)
    pool = pool or direction_pool(cipher)
    candidates = list(pool)
    rng.shuffle(candidates)
    directions = []
    while len(directions) < cap:
        chosen = None
        for cand in candidates:
            if cand in directions:
                continue
            if derivative_nonzero(cipher, base, directions + [cand], rounds):
                chosen = cand
                break
        if chosen is None:
            break
        directions.append(chosen)
    return directions


def find_degree_attaining_set(cipher, base, rounds, d, pool=None, rng=None):
    """Build a size-``d`` degree-attaining set, or raise reporting the true
    maximum reachable degree."""
    dirs = build_max_degree_set(cipher, base, rounds, d, pool, rng)
    if len(dirs) < d:
        raise RuntimeError(
            f"{cipher.name} at {rounds} rounds reaches degree {len(dirs)} "
            f"only (wanted {d}); use more rounds or a different base.")
    return dirs
