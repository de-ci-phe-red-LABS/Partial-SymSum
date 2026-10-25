# Partial-SymSum

**Reference implementation and reproducible examples for the paper**

> **The Decaying Symmetry: Investigating Partial Symmetry in SymSum**
> *de.ci.phe.red Lab, Indian Institute of Technology Bhilai* — INDOCRYPT 2026 (paper #215)

## Authors

| Author | GitHub | Affiliation |
|---|---|---|
| **Sahiba Suryawanshi** | [@suryasahiba](https://github.com/suryasahiba) | de.ci.phe.red Lab, IIT Bhilai |
| **Ratnesh Shrivastava** | [@ratneshshrivastava31](https://github.com/ratneshshrivastava31) | de.ci.phe.red Lab, IIT Bhilai |
| **Dhiman Saha** | [@c-i-p-h-e-r](https://github.com/c-i-p-h-e-r) | de.ci.phe.red Lab, IIT Bhilai |

## About

**Partial-SymSum** is a small, dependency-free Python artifact that lets a reader
of the paper **see the three regimes for themselves** — **ZeroSum**, **SymSum**,
and **SymSum<sub>part</sub>** (partial symmetry) — on the three primitives studied
in the paper: **Keccak-*p* / SHA-3**, **Xoodoo / Xoodyak**, and **Bash-f / Bash**.

It shows how the SymSum distinguisher *decays* into partial symmetry as the
derivative order drops below the algebraic degree, and makes explicit the caveat
that makes or breaks the experiment: the derivative-variable set must actually
**reach** the highest degree, or the result collapses to a trivial all-zero sum.

Everything runs on a stock CPython 3.8+ interpreter; there is nothing to install.

---

## 1. The one idea you need: algebraic degree

Each output bit of an SPN permutation is a polynomial (ANF) over 𝔽₂ in the input
bits. The **algebraic degree** is the largest monomial degree. The non-linear
layer (χ / S-box) has degree 2, so after `r` rounds the degree is bounded by

```
d = 2^r            (the "highest degree")
```

A **k-th order derivative** is the **XOR-sum of the permutation over a
k-dimensional affine subspace** (the 2ᵏ points spanned by `k` input directions).
Its defining property:

```
deg( D^k F ) = deg(F) − k
```

So, at a round count whose maximum degree is `d`:

| derivative order `k` | result | regime |
|---|---|---|
| `k ≥ d+1` | degree < 0 ⇒ **the XOR-sum is all zeros** | **ZeroSum** |
| `k = d` (and `d−1` if δ=2) | top-degree part only, and it is constant-independent | **SymSum** (fully symmetric) |
| `k = d−2, d−3, …` | some round-constant-dependent monomials survive | **SymSum<sub>part</sub>** (partial symmetry) |
| `k` much smaller | — | random (PS ≈ 50%) |

**SymSum** works because the highest-degree monomials that depend on the round
constant (Rcon) do **not** reach the top degree; the `d`-fold derivative isolates
only the Rcon-independent top part, so a self-symmetric input yields a symmetric
output. As the order drops, Rcon-dependent monomials re-enter and symmetry
**decays** — that is `SymSum_part`, the subject of the paper.

### δ: where the constant is added matters

| primitive | Rcon added | δ | SymSum holds at |
|---|---|---|---|
| Keccak-*p* / SHA-3 | **after** χ | **2** | orders `d` **and** `d−1` |
| Bash-f | **after** S3 | **2** | orders `d` **and** `d−1` |
| Xoodoo | **before** χ | **1** | order `d` only |

You can see this directly in the examples: Xoodoo's symmetry starts decaying one
order earlier than Keccak's/Bash's.

---

## 2. The caveat this artifact makes explicit: *reach the highest degree, or you get all zeros*

The `d`-fold derivative is non-zero **only if the chosen derivative-variable set
actually drives the permutation to degree `d` on the observed bits.** If the set
does **not** reach the highest degree, the `d`-fold derivative has nothing left
to isolate and **collapses to all zeros** — a *trivial* ZeroSum that looks like a
result but shows nothing, and SymSum / SymSum_part cannot be observed at all.

That is why the shipped **test vectors are degree-attaining**, and the code
verifies this two ways for every example:

* the `d`-fold derivative is **non-zero** (the set *reaches* degree `d`), and
* the `(d+1)`-fold derivative is **all zeros** (degree is exactly `d`).

**A live example of the caveat.** Bash-f reaches the full degree `2^r`, but it
diffuses slowly (`L3` mixes only inside a vertical plane), so a
derivative-variable set drawn from the **whole state almost never** reaches it —
the `d`-fold derivative is then **all zeros**. At 4 rounds, confining the
variables to a **single 64-bit word** (or a single vertical plane) is what makes
the set reach degree 16. `example_bash.py` shows both sides: three random
whole-state 16-sets each give an all-zero 16-fold, while a single-word set
reaches degree 16 and exhibits the full SymSum → SymSum_part decay (fully
symmetric down to order `d−2 = 14`, matching the paper).

---

## 3. Layout

```
Partial-SymSum/
├── primitives/          # the three permutations (reduced-round capable)
│   ├── keccak.py        #   Keccak-p[1600] + SHA-3-224/256/384/512
│   ├── xoodoo.py        #   Xoodoo[384] (+ full inverse)
│   └── bash.py          #   Bash-f[1536] (+ full inverse)
├── core/
│   ├── ciphers.py       # uniform adapter over the three primitives
│   ├── derivatives.py   # k-fold derivative = XOR-sum over an affine subspace
│   ├── symmetry.py      # fold φ, asymmetry 𝒜, PS%, ZeroSum/SymSum/Part classifier
│   └── degree.py        # degree-attainment check + degree-attaining set search
├── examples/
│   ├── demo_all.py        # all 3 ciphers: fixed message, ZeroSum/SymSum/SymSum_part
│   ├── example_keccak.py
│   ├── example_xoodoo.py
│   └── example_bash.py
├── tests/
│   └── test_vectors.py  # SHA-3 KATs, permutation round-trips, degree/symmetry checks
├── README.md
├── requirements.txt
└── LICENSE
```

---

## 4. Running it

```bash
# all three ciphers in one run: a fixed message per cipher and the
# ZeroSum / SymSum / SymSum_part table at 4 rounds
python3 examples/demo_all.py           # ~40 s

# or one worked example per cipher
python3 examples/example_xoodoo.py     # < 1 s
python3 examples/example_bash.py       # a few seconds
python3 examples/example_keccak.py     # ~15-20 s  (ZeroSum term is a 2^17 sum)

# correctness tests (SHA-3 known-answer tests, round-trips, degree checks)
python3 tests/test_vectors.py
```

### Sample output (Xoodoo, 3 rounds, d\* = 8, δ = 1)

```
  order k |    PS % | regime
        9 |  100.00 | ZeroSum (all zeros)      <- d*+1
        8 |  100.00 | SymSum (fully symmetric) <- d*
        7 |  100.00 | SymSum (fully symmetric)
        6 |  100.00 | SymSum (fully symmetric)
        5 |   86.98 | SymSum_part (partial symmetry)
        4 |   59.90 | random-like
        3 |   44.27 | random-like
```

A symmetric state is split into two half-lanes σ₁‖σ₂; the **fold**
φ = σ₁ ⊕ σ₂ marks the differing (asymmetric) bit-pairs, the **asymmetry count**
is 𝒜 = wt(φ), and **PS% = 100·(1 − 𝒜/n)**. A fully symmetric output has PS = 100 %;
a random one sits at 50 % ± 3σ. "SymSum_part" is flagged when PS clears the 3σ
threshold of the binomial null.

---

## 5. How this maps to the paper

| paper concept | here |
|---|---|
| higher-order derivative (Defs 1–2) | `core/derivatives.py` |
| structurally / self-symmetric states | `core/symmetry.py` |
| fold φ, asymmetry 𝒜, PS% (Def 6) | `core/symmetry.py` |
| degree-separation / reaching degree `d` | `core/degree.py` |
| ZeroSum / SymSum / SymSum<sub>part</sub> | the regime table printed by each example |
| δ = 2 (Keccak, Bash) vs δ = 1 (Xoodoo) | `core/ciphers.py`, visible in the examples |

The examples use **small round counts** so they run in seconds; the paper's
headline distinguishers use more rounds (e.g. 5-round SHA-3-384 with the χ⁻¹
trick), where the same decay is deeper. The qualitative behaviour is identical.

---

## 6. Correctness / validation

* **Keccak / SHA-3** — validated against the NIST FIPS-202 known-answer tests
  (`sha3_224/256/384/512`), plus a χ⁻¹ round-trip.
* **Xoodoo** — validated by a full forward/inverse round-trip at 1–12 rounds.
* **Bash-f** — validated by a full forward/inverse round-trip, and cross-checked
  **bit-for-bit against the lab's MATLAB reference** implementation (the reference
  output is baked into `tests/test_vectors.py`, so the check runs without MATLAB).
  *Note:* the MATLAB initial round-constant literal is `…9400`, i.e. it is missing
  the low byte of the spec value `B194BAC80A08F53B` (octet-reversed
  `0x3BF5080AC8BA94B1`). This port uses the **spec-correct** constants
  `C0 = 0x3BF5080AC8BA94B1`, `tap = 0xDC2BE1997FE0D8AE`; the MATLAB value is
  available via the `C0` argument of `bash_f` purely to reproduce the reference
  bit-for-bit (with `C0 = 0x3BF5080AC8BA9400`).

---

## 7. License

MIT — see [LICENSE](LICENSE). © 2026 de.ci.phe.red Lab, IIT Bhilai.
