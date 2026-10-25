"""
Uniform adapter over the three permutations.

A *state* is a tuple of ``n_lanes`` lane integers, each ``lane_bits`` wide.
Every cipher exposes the same interface so the derivative / symmetry code is
written once:

    c = KECCAK            # or XOODOO, BASH
    y = c.perm(x, rounds) # reduced-round permutation on a lane tuple
    d = c.max_degree(r)   # algebraic-degree bound at r rounds (= 2**r)
    c.delta               # SymSum degree gap (2 for Keccak/Bash, 1 for Xoodoo)

The symmetry involution tau is the within-lane half swap (bit z <-> z+half,
half = lane_bits//2); it commutes with every step of each round function except
the round-constant addition, which is exactly what SymSum exploits.
"""

from __future__ import annotations
from primitives import keccak, xoodoo, bash


class Cipher:
    def __init__(self, name, lane_bits, n_lanes, delta, _perm):
        self.name = name
        self.lane_bits = lane_bits
        self.n_lanes = n_lanes
        self.delta = delta
        self.half = lane_bits // 2
        self.n_bits = lane_bits * n_lanes
        self._perm = _perm

    def perm(self, state, rounds):
        return self._perm(state, rounds)

    @staticmethod
    def max_degree(rounds):
        return 2 ** rounds


# --- Keccak-p[1600]: lanes indexed i = x + 5y ---------------------------------

def _keccak_perm(state, rounds):
    A = [[state[x + 5 * y] for y in range(5)] for x in range(5)]
    B = keccak.keccak_p(A, rounds)
    return tuple(B[i % 5][i // 5] for i in range(25))


KECCAK = Cipher("Keccak-p[1600]", 64, 25, 2, _keccak_perm)


# --- Xoodoo[384]: lanes indexed i = x + 4y ------------------------------------

def _xoodoo_perm(state, rounds):
    A = [[state[x + 4 * y] for x in range(4)] for y in range(3)]
    B = xoodoo.xoodoo(A, rounds)
    return tuple(B[i // 4][i % 4] for i in range(12))


XOODOO = Cipher("Xoodoo[384]", 32, 12, 1, _xoodoo_perm)


# --- Bash-f[1536]: 24 lanes, already flat -------------------------------------

def _bash_perm(state, rounds):
    return tuple(bash.bash_f(list(state), rounds))


BASH = Cipher("Bash-f[1536]", 64, 24, 2, _bash_perm)


ALL = {"keccak": KECCAK, "xoodoo": XOODOO, "bash": BASH}
