#!/usr/bin/env python3
"""
Bash-f[1536] example: ZeroSum / SymSum / SymSum_part at 4 rounds.

Bash-f reaches the full degree d* = 2^4 = 16 at 4 rounds -- BUT only with a
carefully chosen variable set.  Bash-f diffuses slowly (L3 mixes only inside a
vertical plane), so a derivative-variable set drawn from the *whole state*
almost never reaches degree 16; the derivative then collapses to all zeros.
Confining the variables to a single 64-bit word makes the set degree-attaining.
This is the "reach the highest degree, or you get all zeros" caveat in its
sharpest form, and matches the paper's own observation.

delta = 2 (constant added after S3), so SymSum holds at d and d-1.

Runtime: ~15-20 s (the ZeroSum term is a 2^17-point sum).
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import ciphers, symmetry as sym, degree
from core.derivatives import kfold_derivative
from examples import _demo

BASH = ciphers.BASH


def _caveat_whole_state(rounds=4, target=16, trials=3, seed=0):
    """Show that random whole-state variable sets do NOT reach degree 16 -> the
    16-fold derivative is a trivial all-zero ZeroSum.  (A single-word set, used
    below, does reach it.)"""
    base = sym.self_symmetric_state(BASH, random.Random(seed))
    pool = degree.direction_pool(BASH)
    print("-" * 70)
    print(" CAVEAT: random whole-state derivative-variable sets for Bash-f")
    print("-" * 70)
    for t in range(trials):
        dirs = random.Random(seed + 100 + t).sample(pool, target)
        out = kfold_derivative(BASH, base, dirs, rounds)
        print(f"   random whole-state {target}-set #{t + 1}: "
              f"{target}-fold is {'ALL ZEROS' if sym.is_zero(out) else 'non-zero'}")
    print(f" => random whole-state sets do not reach degree {target}, so the")
    print(f"    {target}-fold derivative collapses to a trivial all-zero ZeroSum.")
    print(f"    A degree-attaining set must be confined to a single word (below).")
    print()


# A degree-16 variable set confined to word (lane) 0, verified for the base
# state used here (random seed 0).  Finding such a set by search is slow, so the
# verified test vector is recorded directly -- the run still checks it reaches
# degree 16 (non-zero 16-fold) and that the 17-fold is a ZeroSum.
_BASH_WORD0_Z = [0, 1, 3, 5, 7, 10, 14, 17, 19, 20, 21, 23, 26, 27, 28, 29]

if __name__ == "__main__":
    _caveat_whole_state()
    dirs = [sym.symmetric_direction(BASH, 0, z) for z in _BASH_WORD0_Z]
    _demo.run(BASH, rounds=4, attained_degree=16,
              pool=degree.word_pool(BASH, lane=0), directions=dirs,
              pool_note="variables confined to a single 64-bit word (lane 0) -- "
                        "required for Bash-f to reach the full degree")
