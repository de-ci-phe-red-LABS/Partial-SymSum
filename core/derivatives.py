"""
Higher-order derivatives as XOR-sums over an affine subspace.

A k-fold derivative with respect to k directions d_1..d_k is

    D = XOR over all 2^k subsets T of { directions }  of  F(base XOR sum(T)).

Each direction here is a *symmetric* direction (flips a mirrored bit-pair), so
every point ``base XOR sum(T)`` stays self-symmetric and the output XOR-sum can
be meaningfully folded.  This is precisely the k-fold derivative of
Definitions 1-2 of the paper, restricted to the symmetric affine subspace.
"""

from __future__ import annotations


def xor_states(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def kfold_derivative(cipher, base, directions, rounds):
    """XOR-sum of the reduced-round permutation over the 2^k-point affine
    subspace ``base + span(directions)``."""
    k = len(directions)
    n = cipher.n_lanes
    acc = [0] * n
    point = list(base)
    prev = 0
    for mask in range(1 << k):
        # Gray-code update: flip the single direction that changed.
        gray = mask ^ (mask >> 1)
        changed = gray ^ prev
        if changed:
            j = changed.bit_length() - 1
            d = directions[j]
            for i in range(n):
                point[i] ^= d[i]
        prev = gray
        y = cipher.perm(tuple(point), rounds)
        for i in range(n):
            acc[i] ^= y[i]
    return tuple(acc)
