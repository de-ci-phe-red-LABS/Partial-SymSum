#!/usr/bin/env python3
"""
Keccak-p[1600] example: ZeroSum / SymSum / SymSum_part at 4 rounds.

4-round Keccak reaches the full algebraic degree d* = 2^4 = 16 on the symmetric
subspace.  delta = 2 (the round constant is added after chi), so full symmetry
survives at orders d and d-1; partial symmetry appears just below.

Runtime: ~15-20 s (the ZeroSum term is a 2^17-point sum).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import ciphers
from examples import _demo

if __name__ == "__main__":
    _demo.run(ciphers.KECCAK, rounds=4, attained_degree=16)
