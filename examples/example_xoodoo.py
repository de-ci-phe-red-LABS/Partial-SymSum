#!/usr/bin/env python3
"""
Xoodoo[384] example: ZeroSum / SymSum / SymSum_part at 3 rounds.

3-round Xoodoo reaches the full algebraic degree d* = 2^3 = 8.  delta = 1 (the
round constant is added *before* chi), so -- unlike Keccak/Bash -- full symmetry
survives only at order d; partial symmetry sets in immediately below.

Runtime: < 1 s.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import ciphers
from examples import _demo

if __name__ == "__main__":
    _demo.run(ciphers.XOODOO, rounds=3, attained_degree=8)
