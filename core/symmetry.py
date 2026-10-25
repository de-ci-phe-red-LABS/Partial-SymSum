"""
Symmetry measurement: the fold map phi, the asymmetry count A = wt(phi),
percentage symmetry, and the ZeroSum / SymSum / SymSum_part classifier.

tau is the within-lane half swap (z <-> z + half).  A state is *self-symmetric*
when every lane equals its own half swap, i.e. the low and high halves of each
lane coincide.  The fold

    phi(state) = XOR of (low half, high half) of each lane         (N/2 bits)

marks exactly the mirrored bit-pairs that differ; A = wt(phi) and

    PS% = 100 * (1 - A / (N/2)).

A fully symmetric state has A = 0 (PS = 100%); a random state has PS ~ 50%.
"""

from __future__ import annotations
import random


def half_mask(cipher):
    return (1 << cipher.half) - 1


def self_symmetric_state(cipher, rng=None):
    """Random self-symmetric state: each lane's two halves are equal."""
    rng = rng or random
    m = half_mask(cipher)
    out = []
    for _ in range(cipher.n_lanes):
        h = rng.getrandbits(cipher.half) & m
        out.append(h | (h << cipher.half))
    return tuple(out)


def is_self_symmetric(cipher, state):
    m = half_mask(cipher)
    return all((lane & m) == ((lane >> cipher.half) & m) for lane in state)


def symmetric_direction(cipher, lane_idx, z):
    """A derivative direction that flips the mirrored pair (z, z+half) of one
    lane, so adding it preserves self-symmetry."""
    assert 0 <= z < cipher.half
    bit = (1 << z) | (1 << (z + cipher.half))
    return tuple(bit if i == lane_idx else 0 for i in range(cipher.n_lanes))


def fold(cipher, state):
    """phi(state) as a single integer of N/2 bits (lanes concatenated)."""
    m = half_mask(cipher)
    f = 0
    for i, lane in enumerate(state):
        lf = (lane & m) ^ ((lane >> cipher.half) & m)
        f |= lf << (cipher.half * i)
    return f


def asymmetry_count(cipher, state):
    return bin(fold(cipher, state)).count("1")


def percentage_symmetry(cipher, state):
    n = cipher.half * cipher.n_lanes
    return 100.0 * (1.0 - asymmetry_count(cipher, state) / n)


def is_zero(state):
    return all(lane == 0 for lane in state)


def classify(cipher, state):
    """Return (label, PS%). ZeroSum => all zeros; SymSum => nonzero but fully
    symmetric; SymSum_part => partially symmetric; else random-like."""
    if is_zero(state):
        return "ZeroSum (all zeros)", 100.0
    a = asymmetry_count(cipher, state)
    n = cipher.half * cipher.n_lanes
    ps = 100.0 * (1.0 - a / n)
    if a == 0:
        return "SymSum (fully symmetric)", ps
    # 3-sigma binomial threshold for the random null (mean 50%).
    sigma = 100.0 * (0.5 / n ** 0.5)
    if ps >= 50.0 + 3 * sigma:
        return "SymSum_part (partial symmetry)", ps
    return "random-like", ps
