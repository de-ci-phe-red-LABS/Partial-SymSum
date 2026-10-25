#!/usr/bin/env python3
"""
Correctness tests for the three primitives and the derivative/symmetry pipeline.

Run:  python3 tests/test_vectors.py

Covers:
  * Keccak / SHA-3  -> NIST FIPS-202 known-answer tests (via hashlib) + chi^-1.
  * Xoodoo          -> full forward/inverse round-trip.
  * Bash-f          -> full forward/inverse round-trip, and a bit-for-bit
                       cross-check against the lab MATLAB reference.
  * derivatives/degree/symmetry -> ZeroSum at d*+1, SymSum at d* (small configs).
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from primitives import keccak, xoodoo, bash
from core import ciphers, symmetry as sym, degree
from core.derivatives import kfold_derivative


def test_sha3_kats():
    for name in ("sha3_224", "sha3_256", "sha3_384", "sha3_512"):
        fn = getattr(keccak, name)
        for msg in (b"", b"abc", b"The quick brown fox jumps over the lazy dog"):
            assert fn(msg).hex() == getattr(hashlib, name)(msg).hexdigest(), name
    print("[ok] SHA-3 known-answer tests (224/256/384/512)")


def test_keccak_chi_inv():
    rng = random.Random(0)
    A = [[rng.getrandbits(64) for _ in range(5)] for _ in range(5)]
    assert keccak.chi_inv(keccak.chi(A)) == A
    print("[ok] Keccak chi / chi^-1 round-trip")


def test_xoodoo_roundtrip():
    rng = random.Random(1)
    for nr in (1, 3, 6, 12):
        A = [[rng.getrandbits(32) for _ in range(4)] for _ in range(3)]
        assert xoodoo.xoodoo_inv(xoodoo.xoodoo(A, nr), nr) == A, nr
    print("[ok] Xoodoo forward/inverse round-trip (1,3,6,12 rounds)")


def test_bash_roundtrip():
    rng = random.Random(2)
    for nr in (1, 2, 4, 8, 24):
        S = [rng.getrandbits(64) for _ in range(24)]
        assert bash.bash_f_inv(bash.bash_f(S, nr), nr) == S, nr
    print("[ok] Bash-f forward/inverse round-trip (1,2,4,8,24 rounds)")


# Reference output of the lab MATLAB Bash_f.m on the state S{i}=i-1 (i=1..24),
# 24 rounds.  MATLAB uses the (buggy) initial constant 0x3BF5080AC8BA9400 and the
# spec tap 0xDC2BE1997FE0D8AE; reproduced here bit-for-bit to validate the port.
_BASH_MATLAB_REF = [
    0x7CBCBF806CC1729D, 0xC03B8BD1CC72E2CD, 0x3A304145431859FF, 0x4EDCA53E4D12BCAC,
    0x1F9073C0C193F25A, 0xC3C94E56C529D4EF, 0x20A0D5435BDB4E1D, 0xC6D20C85CCBB761E,
    0x131A6393D17AEA79, 0xA8842A3B5D7869CA, 0x316FD90A147B5869, 0xCB2DE0B3CC92D2B5,
    0x2E89A9E7EADCC4ED, 0x1F2A8A0D056E7512, 0xCC64CDDE8BB77D09, 0x7429C34070D2DC7C,
    0xD9E50B712642AB30, 0x50FDA984D9811F28, 0x6878792BE874B3F6, 0x1F7D730572E58872,
    0x4AC691FD382FD63F, 0x4319D51D7FE2A456, 0x7132FAF27A897CD3, 0x04527A5D499AA405,
]


def test_bash_matlab_reference():
    out = bash.bash_f(list(range(24)), 24,
                      C0=0x3BF5080AC8BA9400, tap=0xDC2BE1997FE0D8AE)
    assert out == _BASH_MATLAB_REF, "Bash-f disagrees with MATLAB reference"
    # And the spec-correct constant must give a different value.
    assert bash.bash_f(list(range(24)), 24) != _BASH_MATLAB_REF
    print("[ok] Bash-f matches lab MATLAB reference bit-for-bit")


def test_degree_and_symmetry():
    # Small, fast configs: verify ZeroSum at d*+1 and SymSum at d*.
    # Keccak/Xoodoo reach the degree bound with whole-state variables; Bash
    # diffuses slowly, so its variables are confined to a single word.
    configs = [
        (ciphers.KECCAK, 3, degree.direction_pool(ciphers.KECCAK), 8),
        (ciphers.XOODOO, 3, degree.direction_pool(ciphers.XOODOO), 8),
        (ciphers.BASH, 3, degree.word_pool(ciphers.BASH, 0), None),
    ]
    for c, rounds, pool, expect in configs:
        base = sym.self_symmetric_state(c, random.Random(0))
        cap = c.max_degree(rounds)
        D = degree.build_max_degree_set(c, base, rounds, cap, pool,
                                        random.Random(0))
        dstar = len(D)
        if expect is not None:
            assert dstar == expect, (c.name, dstar)
        extra = next(x for x in pool if x not in set(D))
        zs = kfold_derivative(c, base, D + [extra], rounds)
        assert sym.is_zero(zs), f"{c.name}: expected ZeroSum at d*+1"
        ss = kfold_derivative(c, base, D, rounds)
        assert not sym.is_zero(ss), f"{c.name}: expected non-zero at d*"
        assert sym.asymmetry_count(c, ss) == 0, f"{c.name}: expected SymSum at d*"
        print(f"[ok] {c.name}: degree d*={dstar}, "
              f"ZeroSum at k={dstar + 1}, SymSum at k={dstar}")


if __name__ == "__main__":
    test_sha3_kats()
    test_keccak_chi_inv()
    test_xoodoo_roundtrip()
    test_bash_roundtrip()
    test_bash_matlab_reference()
    test_degree_and_symmetry()
    print("\nALL TESTS PASSED")
