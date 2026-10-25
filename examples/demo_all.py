#!/usr/bin/env python3
"""
One consolidated demonstration for all three primitives.

At a FIXED round (4 rounds, maximum degree d = 2^4 = 16) each cipher uses a
fixed self-symmetric input message and a fixed derivative-variable set that
reaches the highest degree, and we sweep the derivative order to exhibit, in a
single run, the three regimes and the decay between them:
    k = d+1  -> ZeroSum        (XOR-sum is all zeros)
    k = d    -> SymSum         (fully symmetric; also k = d-1 when delta = 2)
    k < d    -> SymSum_part    (partial symmetry, decaying toward random)

The fixed message is the deterministic self-symmetric state produced by
``self_symmetric_state(cipher, random.Random(0))`` (printed for each cipher);
the derivative-variable sets below are verified to attain degree 16 at 4 rounds.
Runs on plain Python 3, no dependencies.  Runtime: ~40 s.
"""

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import ciphers, symmetry as sym
from core.symmetry import symmetric_direction as SD
from core.derivatives import kfold_derivative

ROUNDS = 4
D_MAX = 16  # 2^ROUNDS

# Fixed degree-16 derivative-variable sets (lane, z), verified for the
# seed-0 self-symmetric base state at 4 rounds.
VARS = {
    "keccak": [(11, 7), (9, 24), (10, 16), (12, 12), (5, 20), (4, 2), (2, 20),
               (6, 9), (12, 0), (3, 25), (3, 7), (4, 19), (7, 2), (11, 8),
               (8, 20), (13, 24)],
    "xoodoo": [(6, 15), (2, 0), (1, 11), (6, 8), (0, 10), (10, 15), (1, 7),
               (10, 5), (2, 11), (2, 1), (4, 7), (0, 12), (10, 12), (1, 2),
               (10, 1), (6, 9)],
    # Bash-f diffuses slowly: the variables must lie in a single word (lane 0).
    "bash": [(0, z) for z in (0, 1, 3, 5, 7, 10, 14, 17, 19, 20, 21, 23, 26,
                              27, 28, 29)],
}


def fixed_message(cipher):
    return sym.self_symmetric_state(cipher, random.Random(0))


def lanes_hex(cipher, state, limit=6):
    w = cipher.lane_bits // 4
    xs = [f"{lane:0{w}x}" for lane in state[:limit]]
    return " ".join(xs) + (" ..." if cipher.n_lanes > limit else "")


def run_fixed_round(key, cipher):
    base = fixed_message(cipher)
    dirs = [SD(cipher, lane, z) for (lane, z) in VARS[key]]
    extra_z = next(z for z in range(cipher.half)
                   if (VARS[key][0][0], z) not in VARS[key])
    extra = SD(cipher, VARS[key][0][0], extra_z)

    bar = "=" * 66
    print(bar)
    print(f" {cipher.name}   rounds = {ROUNDS}   d = {D_MAX}   delta = {cipher.delta}"
          f"  (SymSum at d{', d-1' if cipher.delta == 2 else ''})")
    print(bar)
    print(f" fixed message (self-symmetric, first lanes): {lanes_hex(cipher, base)}")
    if key == "bash":
        print(" derivative variables confined to a single word (Bash-f diffuses slowly)")
    print()
    print(f"  {'order k':>7} | {'PS %':>7} | regime")
    print("  " + "-" * 48)
    for k in range(D_MAX + 1, D_MAX - 4, -1):
        out = kfold_derivative(cipher, base, (dirs + [extra])[:k], ROUNDS)
        label, ps = sym.classify(cipher, out)
        tag = {D_MAX + 1: "  <- d+1", D_MAX: "  <- d"}.get(k, "")
        print(f"  {k:>7} | {ps:>7.2f} | {label}{tag}")
        if k == D_MAX + 1:
            assert sym.is_zero(out), "expected ZeroSum at d+1"
        if k == D_MAX:
            assert not sym.is_zero(out), "expected non-zero SymSum at d"
    print()


if __name__ == "__main__":
    print("\n##### ZeroSum / SymSum / SymSum_part (fixed message, 4 rounds) #####\n")
    for key, cipher in (("keccak", ciphers.KECCAK),
                        ("xoodoo", ciphers.XOODOO),
                        ("bash", ciphers.BASH)):
        run_fixed_round(key, cipher)
