"""
Xoodoo[384] permutation.

State model
-----------
3 planes ``y in {0,1,2}``, each of 4 lanes ``x in {0,1,2,3}`` of 32 bits.
``A[y][x]`` is a 32-bit ``int`` with bit ``z = (A[y][x] >> z) & 1`` (z in 0..31).

Round = theta, rho_west, iota, chi, rho_east  (Daemen et al., "Xoodoo cookbook").
The round constant is added to lane ``A[0][0]`` (iota).  Because iota is applied
*before* chi, the constant is injected before the non-linear layer -> the
SymSum degree gap is delta = 1 (SymSum holds at order d only), unlike
Keccak/Bash where the constant comes after chi (delta = 2).

A full inverse permutation is provided so the implementation can be
round-trip validated.
"""

from __future__ import annotations

MASK32 = (1 << 32) - 1

# Round constants, in order of application for the full 12-round Xoodoo.
RC_XOODOO = [
    0x00000058, 0x00000038, 0x000003C0, 0x000000D0,
    0x00000120, 0x00000014, 0x00000060, 0x0000002C,
    0x00000380, 0x000000F0, 0x000001A0, 0x00000012,
]


def rotl32(v: int, n: int) -> int:
    n &= 31
    if n == 0:
        return v & MASK32
    return ((v << n) | (v >> (32 - n))) & MASK32


def _shift_plane(plane, t, v):
    """Cyclic shift of a 4-lane plane by ``t`` lanes (x) and ``v`` bits (z)."""
    return [rotl32(plane[(x - t) % 4], v) for x in range(4)]


# ---------------------------------------------------------------------------
# Forward steps
# ---------------------------------------------------------------------------

def theta(A):
    P = [A[0][x] ^ A[1][x] ^ A[2][x] for x in range(4)]
    E = [_shift_plane(P, 1, 5)[x] ^ _shift_plane(P, 1, 14)[x] for x in range(4)]
    return [[A[y][x] ^ E[x] for x in range(4)] for y in range(3)]


def rho_west(A):
    return [A[0][:],
            _shift_plane(A[1], 1, 0),
            _shift_plane(A[2], 0, 11)]


def iota(A, rc):
    B = [A[0][:], A[1][:], A[2][:]]
    B[0][0] ^= rc
    return B


def chi(A):
    B = [[0] * 4 for _ in range(3)]
    for x in range(4):
        B[0][x] = A[0][x] ^ ((~A[1][x]) & A[2][x] & MASK32)
        B[1][x] = A[1][x] ^ ((~A[2][x]) & A[0][x] & MASK32)
        B[2][x] = A[2][x] ^ ((~A[0][x]) & A[1][x] & MASK32)
    return B


def rho_east(A):
    return [A[0][:],
            _shift_plane(A[1], 0, 1),
            _shift_plane(A[2], 2, 8)]


def xoodoo_round(A, rc):
    return rho_east(chi(iota(rho_west(theta(A)), rc)))


def xoodoo(A, nrounds=12, start=0):
    """Apply ``nrounds`` Xoodoo rounds starting at constant index ``start``."""
    A = [row[:] for row in A]
    for r in range(start, start + nrounds):
        A = xoodoo_round(A, RC_XOODOO[r % len(RC_XOODOO)])
    return A


# ---------------------------------------------------------------------------
# chi inverse (3-bit S-box across the 3 planes at each (x, z))
# ---------------------------------------------------------------------------

def _build_chi3_inverse():
    inv = [0] * 8
    for v in range(8):
        a = [(v >> y) & 1 for y in range(3)]
        b = [a[y] ^ ((1 - a[(y + 1) % 3]) & a[(y + 2) % 3]) for y in range(3)]
        bv = sum(b[y] << y for y in range(3))
        inv[bv] = v
    return inv


_CHI3_INV = _build_chi3_inverse()


def chi_inv(A):
    out = [[0] * 4 for _ in range(3)]
    for x in range(4):
        for z in range(32):
            bv = sum(((A[y][x] >> z) & 1) << y for y in range(3))
            v = _CHI3_INV[bv]
            for y in range(3):
                out[y][x] |= ((v >> y) & 1) << z
    return out


# ---------------------------------------------------------------------------
# theta inverse (invert the linear parity map on a 128-bit plane, once)
# ---------------------------------------------------------------------------

def _plane_to_int(p):
    return p[0] | (p[1] << 32) | (p[2] << 64) | (p[3] << 96)


def _int_to_plane(v):
    return [(v >> (32 * x)) & MASK32 for x in range(4)]


def _theta_L(pint):
    """P -> P ^ (P<<(1,5)) ^ (P<<(1,14)) on a plane, as a 128-bit int."""
    p = _int_to_plane(pint)
    s1 = _shift_plane(p, 1, 5)
    s2 = _shift_plane(p, 1, 14)
    return _plane_to_int([p[x] ^ s1[x] ^ s2[x] for x in range(4)])


def _invert_gf2(colimages):
    """Given L as 128 column images (L(e_j)), return L^{-1} column images."""
    n = 128
    # rows of augmented [M | I], each a 256-bit int; M[i][j] = bit i of col j.
    rows = []
    for i in range(n):
        m = 0
        for j in range(n):
            m |= (((colimages[j] >> i) & 1) << j)
        rows.append(m | (1 << (n + i)))
    # Gaussian elimination.
    for col in range(n):
        piv = next(r for r in range(col, n) if (rows[r] >> col) & 1)
        rows[col], rows[piv] = rows[piv], rows[col]
        for r in range(n):
            if r != col and (rows[r] >> col) & 1:
                rows[r] ^= rows[col]
    # Extract inverse columns: inv_col[j] = bits of column j of M^{-1}.
    inv_cols = [0] * n
    for i in range(n):
        high = rows[i] >> n
        for j in range(n):
            if (high >> j) & 1:
                inv_cols[j] |= (1 << i)
    return inv_cols


_THETA_LINV = _invert_gf2([_theta_L(1 << j) for j in range(128)])


def theta_inv(A):
    Pp = _plane_to_int([A[0][x] ^ A[1][x] ^ A[2][x] for x in range(4)])
    # Solve L(P) = P' -> P = sum of inverse columns selected by P' bits.
    P = 0
    for j in range(128):
        if (Pp >> j) & 1:
            P ^= _THETA_LINV[j]
    p = _int_to_plane(P)
    s1 = _shift_plane(p, 1, 5)
    s2 = _shift_plane(p, 1, 14)
    E = [s1[x] ^ s2[x] for x in range(4)]
    return [[A[y][x] ^ E[x] for x in range(4)] for y in range(3)]


def rho_west_inv(A):
    return [A[0][:],
            _shift_plane(A[1], -1, 0),
            _shift_plane(A[2], 0, -11)]


def rho_east_inv(A):
    return [A[0][:],
            _shift_plane(A[1], 0, -1),
            _shift_plane(A[2], -2, -8)]


def xoodoo_round_inv(A, rc):
    return theta_inv(rho_west_inv(iota(chi_inv(rho_east_inv(A)), rc)))


def xoodoo_inv(A, nrounds=12, start=0):
    A = [row[:] for row in A]
    for r in reversed(range(start, start + nrounds)):
        A = xoodoo_round_inv(A, RC_XOODOO[r % len(RC_XOODOO)])
    return A
