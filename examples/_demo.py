"""
Shared driver for the per-cipher examples.

For a chosen cipher and (small) round count it:

  1. builds a self-symmetric input state (the base),
  2. greedily finds a degree-attaining symmetric derivative-variable set -- the
     "test vector that reaches the highest degree" -- and reports the attained
     degree d* (which for low-round Bash is below the bound 2**rounds),
  3. sweeps the derivative order k around d* and labels each regime:
        k = d*+1  -> ZeroSum        (XOR-sum is all zeros)
        k = d*    -> SymSum         (fully symmetric; also k=d*-1 when delta=2)
        k = d*-1  -> SymSum/Part    (depends on delta)
        k = d*-2, d*-3 -> SymSum_part (partial symmetry, above the 3-sigma null)

It also makes the "did not reach the highest degree => all zeros" point
explicit whenever d* < 2**rounds.
"""

from __future__ import annotations
import random

from core import symmetry as sym, degree
from core.derivatives import kfold_derivative


def _describe_direction(cipher, d):
    for i, lane in enumerate(d):
        if lane:
            z = (lane & -lane).bit_length() - 1
            return (i, z)
    return None


def _lanes_hex(cipher, state, limit=6):
    w = cipher.lane_bits // 4
    xs = [f"{lane:0{w}x}" for lane in state[:limit]]
    tail = " ..." if cipher.n_lanes > limit else ""
    return " ".join(xs) + tail


def run(cipher, rounds, attained_degree, seed=0, depth=5, pool=None,
        pool_note=None, directions=None):
    """Demonstrate ZeroSum / SymSum / SymSum_part for ``cipher`` at ``rounds``.

    ``attained_degree`` is the algebraic degree the derivative-variable set
    reaches at this round count (= 2**rounds here).  It is verified
    empirically: the d*-fold derivative is non-zero and the (d*+1)-fold is all
    zeros.  ``pool`` restricts the directions searched (e.g. a single word for
    Bash); defaults to the whole state.  ``directions`` supplies a precomputed
    degree-attaining set, skipping the search (used for Bash, whose single-word
    search is slow)."""
    rng = random.Random(seed)
    base = sym.self_symmetric_state(cipher, rng)
    if pool is None:
        pool = degree.direction_pool(cipher)
    bound = cipher.max_degree(rounds)
    dstar = attained_degree

    if directions is not None:
        D = list(directions)
    else:
        D = degree.build_max_degree_set(cipher, base, rounds, dstar,
                                        pool, random.Random(seed + 1))
    if len(D) < dstar:
        raise RuntimeError(f"could not reach degree {dstar} for {cipher.name}")
    extra = next(x for x in pool if x not in set(D))

    bar = "=" * 70
    print(bar)
    print(f" {cipher.name}   rounds = {rounds}   delta = {cipher.delta}  "
          f"(SymSum holds at d{', d-1' if cipher.delta == 2 else ''})")
    print(bar)
    print(f" algebraic-degree bound 2^{rounds} = {bound}")
    print(f" attained degree on symmetric subspace d* = {dstar}"
          + ("   (< bound -- see note)" if dstar < bound else "   (= bound)"))
    print()
    print(" test vector (reaches the highest degree):")
    if pool_note:
        print(f"   {pool_note}")
    print(f"   base self-symmetric state (first lanes): {_lanes_hex(cipher, base)}")
    dirs_desc = [_describe_direction(cipher, d) for d in D]
    print(f"   {dstar} symmetric derivative variables (lane, z): {dirs_desc}")
    print()
    print(f"  {'order k':>7} | {'PS %':>7} | regime")
    print("  " + "-" * 52)
    for k in range(dstar + 1, max(dstar - depth, 1) - 1, -1):
        dirs = (D + [extra])[:k]
        out = kfold_derivative(cipher, base, dirs, rounds)
        label, ps = sym.classify(cipher, out)
        note = {dstar + 1: "  <- d*+1", dstar: "  <- d*"}.get(k, "")
        print(f"  {k:>7} | {ps:>7.2f} | {label}{note}")

        if k == dstar + 1:
            assert sym.is_zero(out), "expected ZeroSum at order d*+1"
        if k == dstar:
            assert not sym.is_zero(out), "expected non-zero SymSum at order d*"
    print()
    if dstar < bound:
        print(f" NOTE: d* = {dstar} < 2^{rounds} = {bound}.  Any {bound}-fold derivative")
        print(f"       here is ALL ZEROS (a trivial ZeroSum): the variable set cannot")
        print(f"       reach degree {bound} at {rounds} rounds, so SymSum / SymSum_part")
        print(f"       are observable only at order d* = {dstar} and below.")
        print()
    return dstar
