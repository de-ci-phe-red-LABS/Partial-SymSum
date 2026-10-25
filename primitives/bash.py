"""
Bash-f[1536] permutation (Belarusian standard STB 34.101.77 / ePrint 2016/587).

State model
-----------
24 lanes of 64 bits, ``S[0..23]``, arranged as 3 rows x 8 columns
(row = i // 8, col = i % 8).  Lane ``S[i]`` is an ``int`` with bit
``z = (S[i] >> z) & 1`` (z in 0..63).

Each round applies, per column v in 0..7, the linear map ``L3`` then the
non-linear ``S3`` to the vertical triple ``(S[v], S[v+8], S[v+16])``; then the
lane permutation ``P``; then the round constant is XORed into ``S[23]``.
Because the constant is added *after* the non-linear ``S3``, the SymSum degree
gap is delta = 2 (SymSum holds at orders d and d-1), as for Keccak.

Ported from the lab MATLAB reference (Bash_f.m, S3.m, L3.m, P.m).  Note: the
MATLAB round-constant literals lose their low bits to IEEE-754 ``double``
rounding (e.g. the initial constant becomes ...9400 instead of ...94B1); this
port uses the spec-correct octet-reversed values.
"""

from __future__ import annotations

MASK64 = (1 << 64) - 1

# Octet-reversed spec constants (Bash convention).
#   paper:  C0  = B194BAC80A08F53B      ->  0x3BF5080AC8BA94B1
#   paper:  tap = AED8E07F99E12BDC      ->  0xDC2BE1997FE0D8AE
C0_BASH = 0x3BF5080AC8BA94B1
TAP_BASH = 0xDC2BE1997FE0D8AE


def rotl64(v: int, d: int) -> int:
    """RotHi: cyclic shift toward high bits (left rotation)."""
    d &= 63
    if d == 0:
        return v & MASK64
    return ((v << d) | (v >> (64 - d))) & MASK64


def L3(W, m1, n1, m2, n2):
    W0, W1, W2 = W
    t0 = W0 ^ W1 ^ W2
    t1 = W1 ^ rotl64(W0, m1) ^ rotl64(t0, n1)
    inner = W1 ^ rotl64(t0, n1)
    t2 = W2 ^ rotl64(W2, m2) ^ rotl64(inner, n2)
    return (t0, t1, t2)


def S3(W):
    W0, W1, W2 = W
    o0 = W0 ^ (W1 | (~W2 & MASK64))
    o1 = W1 ^ (W0 | W2)
    o2 = W2 ^ (W0 & W1)
    return (o0 & MASK64, o1 & MASK64, o2 & MASK64)


def _pi0(v): return (v + 2 * (v % 2) + 7) % 8
def _pi1(v): return v + 1 - 2 * (v % 2)
def _pi2(v): return (5 * v + 6) % 8


def P(S):
    out = [0] * 24
    for u in range(24):
        if u < 8:
            tgt = _pi0(u) + 8
        elif u < 16:
            tgt = _pi1(u - 8) + 16
        else:
            tgt = _pi2(u - 16)
        out[tgt] = S[u]
    return out


def _next_C(C, tap=TAP_BASH):
    if C & 1:
        return (C >> 1) ^ tap
    return C >> 1


def _column_constants():
    """The 8 (m1,n1,m2,n2) tuples used within every round (reset each round)."""
    m1, n1, m2, n2 = 8, 53, 14, 1
    tuples = []
    for _ in range(8):
        tuples.append((m1, n1, m2, n2))
        m1 = (7 * m1) % 64
        n1 = (7 * n1) % 64
        m2 = (7 * m2) % 64
        n2 = (7 * n2) % 64
    return tuples


_COL_CONST = _column_constants()


def bash_f(S, nrounds=24, C0=C0_BASH, tap=TAP_BASH):
    """Apply ``nrounds`` Bash-f rounds to the 24-lane state ``S``.
    Returns a new list; ``S`` is not modified.  ``C0``/``tap`` are exposed so
    the MATLAB reference (which uses double-rounded constants) can be matched
    bit-for-bit in the validation."""
    S = list(S)
    C = C0
    for _ in range(nrounds):
        for v in range(8):
            m1, n1, m2, n2 = _COL_CONST[v]
            idx = (v, v + 8, v + 16)
            w = L3((S[idx[0]], S[idx[1]], S[idx[2]]), m1, n1, m2, n2)
            w = S3(w)
            S[idx[0]], S[idx[1]], S[idx[2]] = w
        S = P(S)
        S[23] ^= C
        C = _next_C(C, tap)
    return S


# ---------------------------------------------------------------------------
# Inverse (used only for the round-trip self-test)
# ---------------------------------------------------------------------------

def _build_s3_inverse():
    inv = [0] * 8
    for v in range(8):
        a0, a1, a2 = (v >> 0) & 1, (v >> 1) & 1, (v >> 2) & 1
        o0 = a0 ^ (a1 | (1 - a2))
        o1 = a1 ^ (a0 | a2)
        o2 = a2 ^ (a0 & a1)
        inv[(o0 & 1) | ((o1 & 1) << 1) | ((o2 & 1) << 2)] = v
    return inv


_S3_INV = _build_s3_inverse()


def S3_inv(W):
    W0, W1, W2 = W
    o0 = o1 = o2 = 0
    for z in range(64):
        ov = ((W0 >> z) & 1) | (((W1 >> z) & 1) << 1) | (((W2 >> z) & 1) << 2)
        v = _S3_INV[ov]
        o0 |= (v & 1) << z
        o1 |= ((v >> 1) & 1) << z
        o2 |= ((v >> 2) & 1) << z
    return (o0, o1, o2)


def _invert_gf2(colimages, n):
    rows = []
    for i in range(n):
        m = 0
        for j in range(n):
            m |= (((colimages[j] >> i) & 1) << j)
        rows.append(m | (1 << (n + i)))
    for col in range(n):
        piv = next(r for r in range(col, n) if (rows[r] >> col) & 1)
        rows[col], rows[piv] = rows[piv], rows[col]
        for r in range(n):
            if r != col and (rows[r] >> col) & 1:
                rows[r] ^= rows[col]
    inv_cols = [0] * n
    for i in range(n):
        high = rows[i] >> n
        for j in range(n):
            if (high >> j) & 1:
                inv_cols[j] |= (1 << i)
    return inv_cols


def _triple_to_int(w):
    return w[0] | (w[1] << 64) | (w[2] << 128)


def _int_to_triple(v):
    return (v & MASK64, (v >> 64) & MASK64, (v >> 128) & MASK64)


def _build_L3_inverse(consts):
    L_inv = {}
    for (m1, n1, m2, n2) in consts:
        cols = []
        for j in range(192):
            w = _int_to_triple(1 << j)
            cols.append(_triple_to_int(L3(w, m1, n1, m2, n2)))
        L_inv[(m1, n1, m2, n2)] = _invert_gf2(cols, 192)
    return L_inv


_L3_INV = _build_L3_inverse(_COL_CONST)


def L3_inv(W, m1, n1, m2, n2):
    inv_cols = _L3_INV[(m1, n1, m2, n2)]
    v = _triple_to_int(W)
    out = 0
    for j in range(192):
        if (v >> j) & 1:
            out ^= inv_cols[j]
    return _int_to_triple(out)


def _p_target(u):
    if u < 8:
        return _pi0(u) + 8
    if u < 16:
        return _pi1(u - 8) + 16
    return _pi2(u - 16)


_PTGT = [_p_target(u) for u in range(24)]


def P_inv(S):
    # Forward: out[target(u)] = S[u]; hence inverse: out[u] = S[target(u)].
    return [S[_PTGT[u]] for u in range(24)]


def bash_f_inv(S, nrounds=24, C0=C0_BASH):
    """Inverse of ``bash_f`` (for validation)."""
    S = list(S)
    # Reconstruct the constant sequence used in the forward direction.
    seq = []
    C = C0
    for _ in range(nrounds):
        seq.append(C)
        C = _next_C(C)
    for C in reversed(seq):
        S[23] ^= C
        S = P_inv(S)
        for v in range(8):
            m1, n1, m2, n2 = _COL_CONST[v]
            idx = (v, v + 8, v + 16)
            w = S3_inv((S[idx[0]], S[idx[1]], S[idx[2]]))
            w = L3_inv(w, m1, n1, m2, n2)
            S[idx[0]], S[idx[1]], S[idx[2]] = w
    return S
