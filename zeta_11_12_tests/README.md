# Numerical checks for the 11/12 guide

These scripts carry out the numerical checks described in Section 9 of
[*Sixth powers in a mean square: an accessible guide to the 11/12 argument for ζ(s)*](../zeta_11_12_guide.pdf)
("the guide"). They test the identities and the inequalities between parameters in Sections 1–8 of the guide on
finite samples, numbered displays and unnumbered statements alike. Equation numbers below, such as (7.5), refer to
that PDF.

## What these checks do and do not show

- **Identities** are tested on finite samples, either exactly or to rounding error. Examples are the Möbius–Gauss
  identity (6.7), the exact Poisson formula (6.4), the cube inversion (8.1) and the local factors (7.3)–(7.5).
- **Inequalities between scale parameters** in Section 8, (8.2)–(8.11), are tested on hundreds of thousands of
  random admissible parameter sets.
- **Analytic estimates are not tested as such.** The mean-square bound (3.8), the large sieve (7.7), the completed
  estimate (7.8), the cube-removal estimate (8.4), the transfer estimate (8.10) and the final estimate (8.12) are
  bounds with unspecified constants. For three of them, (3.8), (7.7) and (7.8), the scripts print the size of
  the left side relative to the right side, as a sanity check, without any pass/fail criterion.
- **External inputs** are used, not proved. These are the transformation law of Kubota's cubic theta function
  (Patterson; Dunn–Radziwiłł) and the quadratic large sieve over number fields (Goldmakher–Louvel). The transcribed
  theta coefficients are checked numerically, both against the automorphy of θ (test 7.5) and end to end through the
  reflection formula at single cusps (test 7.6).

Passing tests are evidence that the formulas in the guide are stated with the right normalizations, signs and
conventions. They are not a proof.

## Requirements and running

Python 3.8 or later, with numpy, scipy, sympy and mpmath. The last run used Python 3.12, numpy 2.2,
scipy 1.16, sympy 1.14 and mpmath 1.3.

```sh
pip install -r requirements.txt
python3 run_all.py            # every script, about 2 minutes; outputs go to outputs/*.txt
python3 run_all.py --quick    # skips the 75-second theta reflection test
python3 test_sec6.py          # any single script can also be run on its own
```

`run_all.py` prints one summary line per script and exits with status 1 if any check fails or any script crashes.
The outputs of the last run are in [`outputs/`](outputs/), where all 74 checks passed.

## Reading the output

Each check prints one line:

```
[PASS] 6.5 [guide (6.7)] Moebius-Gauss identity (6.7), 70 squarefree n: max error 4.71e-15 (tol 1e-12)
```

- `6.5` is the test id used in the scripts and in the table below.
- `[guide (6.7)]` is the display of the guide being checked. The map from test ids to displays is the dictionary
  `GUIDE` in `common.py`.
- `max error` is usually the largest absolute or relative difference between the two sides over all cases
  tried. For combinatorial checks it is a count of failures, with tolerance 0, and for a few tests it is a
  normalized ratio: see the table.
- Indented lines without `[PASS]`/`[FAIL]` are sanity checks of bounds, or printed controls that are *supposed*
  to fail.

## Conventions

These are the guide's conventions.

| Symbol | Meaning |
|---|---|
| ω, O | ω = e^{2πi/3}, O = ℤ[ω]; the element a + bω is stored as the integer pair `(a, b)` |
| N(a + bω) | a² − ab + b² |
| primary | n ≡ 1 mod 3, i.e. a ≡ 1 and b ≡ 0 mod 3 |
| λ | 1 + 2ω = i√3 (`LAM = (1, 2)`) |
| χ_n(x) | the sextic residue symbol (x/n)₆: at a prime π, the sixth root of unity ≡ x^{(N(π)−1)/6} mod π, or 0 if π divides x |
| e(z) | exp(4πi Im z/√3), the additive character of the guide's Sections 6–7 |
| ce(z) | exp(2πi(z + z̄)), the additive character of Kubota, Patterson and Dunn–Radziwiłł, used for the theta function |
| α(n) | n/\|n\| |
| γ_j(n) | N(n)^{−1/2} Σ_{x mod n} χ_n(x)^j e(x/n) |

Arithmetic in O is exact, using Python integers. Residue symbols are computed by Euler's criterion modulo each
prime, and Gauss sums by brute force over complete residue systems, except where `fastgauss.py` is used. That
module is itself checked against brute force in test 8.0. Complex sums are in double precision. The
Mellin integrals and the Dirichlet series use 20–30-digit arithmetic (mpmath).

## Files

| File | Contents |
|---|---|
| `zw.py` | exact arithmetic in ℤ[ω]: products, norms, divisibility, reduction mod m, residue systems, sextic symbols, Möbius function, e(z) |
| `fastgauss.py` | fast lists of primary primes, cubic residue symbols, and cubic Gauss sums g(c) via twisted multiplicativity |
| `theta.py` | Kubota's cubic theta function from Patterson's coefficients τ(ℓ), in the normalization of Dunn–Radziwiłł (5.6)–(5.8), and the action of SL₂(O) on upper half-space |
| `common.py` | shared helpers (cached factorization, χ_n, μ_K, lists of primary n, γ_j, G) and the `report` function with the test-id → guide map |
| `test_sec1_2.py` | Sections 1–2 |
| `test_sec3.py` | Section 3 |
| `test_sec4_5.py` | Sections 4–5 |
| `test_sec4_8_extra.py` | statements of Sections 4–8 not covered elsewhere, including the j = 0 correction in Section 7 |
| `test_sec6.py` | Section 6 |
| `test_sec7_8.py` | Sections 7–8 |
| `test_sec7_completed_sanity.py` | sanity check of the size of the completed estimate (7.8) |
| `test_sec7_reflection_j1.py` | end-to-end check of the theta transformation at single cusps (test 7.6) |
| `run_all.py` | runs everything, saves the outputs and summarizes them |

Every script starts with a docstring that lists its tests, the displays they check, and the run time.

## The checks

The results are from the last run, saved in `outputs/`. "Exact" means the two sides agreed identically, or to
rounding error, in every case tried.

### Section 1: the integer coefficients (`test_sec1_2.py`)

| Test | Guide | What is checked | Result |
|---|---|---|---|
| 1.1 | (1.1) | χ₋₃ as defined equals the Kronecker symbol (−3/m) for m ≤ 3000 and is completely multiplicative (2000 random pairs) | exact |
| 1.2a | (1.2) | coefficientwise c ∗ (1 ∗ χ₋₃) = δ for all m ≤ 10⁵, i.e. the series of c times ζ(s)L(s, χ₋₃) is 1 | exact |
| 1.2b | (1.2) | Σ_{m≤10⁵} c(m)m^{−s} against 1/(ζ(s)L(s, χ₋₃)) at s = 3, 3 + 4i, 4 | 1.4e−13 (truncation) |
| 1.3a | (1.3) | the table for all 9690 prime powers p^j ≤ 10⁵ with p ≠ 3, and c(3) = −1, c(3^j) = 0 for j ≥ 2 | exact |
| 1.3b | after (1.3) | c(m) = ∏ c(p^e) for all m ≤ 10⁵ | exact |
| 1.3c | after (1.3) | c(7) = −2, c(49) = 1, c(25) = −1 | exact |
| 1.5a | (1.5) | coefficientwise for all m ≤ 10⁵ | exact |
| 1.5b | after (1.5) | the local factors at 2 and 3 are 1 − 4^{−s} and 1 − 3^{−s} | exact |
| 1.5c | (1.5) | numerically at s = 3 + 2i | 7e−14 (truncation) |

c(m) is computed from its definition (1.2) with sympy's ordinary Möbius function.

### Section 2: Eisenstein integers (`test_sec1_2.py`)

| Test | Guide | What is checked | Result |
|---|---|---|---|
| 2.0, 2.1 | (2.1) | ω² + ω + 1 = 0; \|a + bω\|² = a² − ab + b² for \|a\|, \|b\| ≤ 40 | exact |
| 2.1b | after (2.1) | there are exactly six units, ±1, ±ω, ±ω² | exact |
| 2.2 | (2.2) | n ≡ 1 mod 3 in ℤ[ω] (tested as 3 \| n − 1) iff a ≡ 1, b ≡ 0 mod 3; each n coprime to 3 has exactly one primary associate, and a multiple of λ has none (14 640 elements) | exact |
| 2.3a | before (2.3) | for every prime p < 3000: p ≡ 1 mod 3 has exactly two primary elements of norm p, which are conjugate with product p; p ≡ 2 mod 3 has no element of norm p and −p is primary; 3 = −λ² | exact |
| 2.3b | — | μ_K computed from the factorization of N(n) in ℤ agrees with μ_K from Eisenstein factorizations | exact |
| 2.3c | before (2.3) | #{primary n : N(n) = m} = Σ_{d\|m} χ₋₃(d) for all m ≤ 10⁵ with 3 ∤ m (a consequence of unique factorization) | exact |
| 2.3 | (2.3) | Σ_{N(n)=m} μ_K(n) = c(m) for all 33 333 values m ≤ 10⁵ coprime to 6 | exact |
| 2.3d | after (2.3) | the local factors (1 − p^{−s})² for split p and 1 − p^{−2s} for inert p | exact |
| 2.4 | (1.4) = (2.4) | A₁(D) computed both ways, with a complex smooth W, for D = 300, 3000, 40 000 | 1e−14 |

Here μ_K is computed from N(n), independently of c, so the two sides of (2.3) are independent computations.

### Section 3: the average (`test_sec3.py`)

| Test | Guide | What is checked | Result |
|---|---|---|---|
| 3.0a | before (3.1) | 6 \| N(π) − 1 for the first 3000 primary primes π ∤ 6 | exact |
| 3.0b, c | after (3.1) | χ_n is multiplicative in the numerator and in the denominator (random samples) | exact |
| 3.1 | (3.4) | "primary and coprime to 6" ⇔ a ≡ 1, b ≡ 0 mod 3, with a, b not both even (\|a\|, \|b\| ≤ 40) | exact |
| 3.2a | (3.4) | the integer-coordinate formula equals M(D, H) for D = 60, H = 90, complex W | exact |
| 3.2b | (3.5), (3.6) | the double sum over columns, with S_H computed directly, equals M(D, H) | 5e−16 |
| 3.2c | (3.7) | A_u(D) = Σ_{(m,6)=1} c_u(m) W(m/D) for 40 rows | 1e−15 |
| 3.2d | after (3.7) | printed: c_u(p) depends on the residue symbols (u = 2 and u = 5 + ω at p = 67) | — |
| 3.3 | (3.8) | sanity: M(D, D^{1.1})/(D·H) = 0.35, 0.32, 0.31 for D = 100, 200, 400 | bounded |

### Sections 4–5: special rows and the zero-free region (`test_sec4_5.py`, `test_sec4_8_extra.py`)

| Test | Guide | What is checked | Result |
|---|---|---|---|
| 4.1 | (4.1), (4.2) | lattice counts S_H(n, n) for four n and H ≤ 10⁴; the value shown is \|count − main term\|/(τ(n)(√H + 1)) | ≤ 0.23 |
| 4.2 | (4.3) | χ_n(v⁶) = 1_{(v,n)=1} (random v, n) | exact |
| 4.11 | (4.3) | A_{v⁶}(D) = Σ_{(n,v)=1} μ_K(n)W(N(n)/D) for composite v | 1e−14 |
| 4.3 | (4.5) | A_{π⁶} − A₁ = −Σ_{π\|n} μ_K(n)W(N(n)/D) for eight π, D = 1500; the size relative to D/N(π) is printed | 8e−17 |
| 4.4 | after (4.5) | J log Y/Y = 0.50, 0.52, 0.51 for Y = 10³, 10⁴, 10⁵ (printed); the rows π⁶ are distinct | exact |
| 4.10 | (4.6) | the extraction inequality on actual values (D = 1500, Y = 200, J = 21); the error term is 0.036·JD²/Y² | holds |
| 4.5 | (4.7) | χ_n(v³) = (v/n)₂ and χ_n(v²) = (v/n)₃, with both symbols computed independently by Euler's criterion | exact |
| 4.6 | (4.8) | (v/n)₂ = (v/N(n)) for \|v\| ≤ 30 and 150 values of n | exact |
| 4.7a | (2.3), (3.7) | c₁(m) = c(m) for m ≤ 2000 | exact |
| 4.7b | (4.9) | c_{v³}(m) = c(m)(v/m) for m ≤ 2000 and v ∈ {2, −1, 3, 5, −3, 7, −6, 10} | exact |
| 4.8 | (4.10) | the Euler product has coefficients c(m)(v/m), m ≤ 2000 | exact |
| 4.9 | after (4.10) | printed: the row counts H^{1/2} and H^{1/12} | — |
| 5.1 | (5.1), (5.2) | the exponents 11/6 + 5ϑ/6 and 5/3 − ϑ/3 (sympy) | symbolic |
| 5.2 | (1.5), (5.4) | Σ_{(m,6)=1} c(m)m^{−s} = 1/F(s) at s = 3 and 2.5 + 4i | within truncation error |
| 5.3 | (5.5) | ∫ W(m/D)D^{−s} dD/D = m^{−s}Ŵ(s) (30-digit quadrature) | 1e−35 |

### Section 6: Poisson summation (`test_sec6.py`, `test_sec4_8_extra.py`)

| Test | Guide | What is checked | Result |
|---|---|---|---|
| 6.3a | (6.1) | Φ̂ for Φ(t) = e^{−πt} by two-dimensional quadrature, against (2/√3)e^{−4πt/3} | 8e−22 |
| 6.1 | (6.2) | the factorization of χ̄_{n₁}χ_{n₂} | exact |
| 6.2 | after (6.2) | ψ is primitive: each of 128 local components is nontrivial | exact |
| 6.3b | (6.4) | the exact Poisson formula with Φ(t) = e^{−πt}, five cases including g ≠ 1 and z₁ = 1 | 5e−15 |
| 6.10 | (3.5) + (6.4) | end to end: the whole smoothed mean square computed directly equals its Poisson-dual form (D = 40, H = 60) | 2e−15 |
| 6.11 | (6.5) | the frequency-range algebra (sympy) | symbolic |
| 6.5 | (6.7) | the Möbius–Gauss identity for 70 squarefree n, including composite n | 5e−15 |
| 6.6 | after (6.7) | G(n) depends only on n mod 4 | exact |
| 6.7 | (6.8) | the paired identity for 22 coprime pairs; a control with the conjugation on the wrong factor fails (error 1.3) | 4e−15 |
| 6.8 | before (6.10) | a_ξ(vm) = a_ξ(v)a_ξ(m)χ_m(v)⁴ (22 pairs) | 4e−15 |
| 6.9 | (6.10) | the character identity behind he⁵v⁴ = kf⁴ | 6e−15 |

### Section 7: theta and the large sieve (`test_sec7_8.py`, `test_sec4_8_extra.py`, `test_sec7_reflection_j1.py`, `test_sec7_completed_sanity.py`, `test_sec6.py`)

| Test | Guide | What is checked | Result |
|---|---|---|---|
| 7.3 | before (7.1) | among indices λ^{−3}x with x primary, the nonzero theta coefficients are exactly x = nb³, n squarefree (908 indices) | exact |
| 7.1 | (7.1) | d₀(λ^{−3}nb³) = 3^{5/2}\|b\|·conj((λ/n)₃)γ₂(n) against the coefficients in `theta.py` (104 indices) | 4e−15 |
| 8.0 | (7.1) | the fast cubic Gauss sums of `fastgauss.py` against brute force | 3e−15 |
| 7.2 | after (7.1) | (λ/n)₃ depends only on n mod 9 | exact |
| 7.5 | input | automorphy of θ under z ↦ z + 1 and under the inversion | 2e−16 |
| 7.6 | Section 7 | the transformation with one z̄-derivative, end to end at single cusps, with the kernel ∏± Γ(1 + t ± 1/6)/Γ(1 − t ± 1/6), for 39 cusp configurations (all cusp types, all values of the multiplier). Controls with the sign flipped, with the kernel of two derivatives, or with the wrong power of conj(α(c)) fail at size 1 | 5.5e−9 |
| 7.4a | (7.3), (7.4) | the additive expansions of χ_π^j (1 ≤ j ≤ 5) and of 1_{π∤x}, eight primes | 7e−15 |
| 7.4b | (7.3), (7.5) | the reflected local sums Σ_{h≠0} χ_π(h)^{−j−2}e(εh^{−1}x/π): a Gauss sum times χ_π(εx)^{−j−2} for j ≠ 4, a Ramanujan sum for j = 4; \|g\| = √N(π) | 3e−14 |
| 7.12a | (7.4) | the constant (inactive) term 1 − 1/N(π) of the j = 0 expansion | 9e−16 |
| 7.12b | (7.5) | the active j = 0 branch summed directly equals −γ₂(π)χ_π(ε)^{−2}N(π)^{−1/2}χ_π(x)^{−2}: modulus N(π)^{−1/2} if π ∤ x, and 0 if π \| x (643 triples (π, ε, x)) | 5e−15 |
| 7.11 | (7.6) | B_{π,1}(x) = χ_π(x)^{−3} = χ_π(x)³ = (x/π)₂ | 3e−15 |
| 7.9 | after (7.6) | λ⁴ = 9; the least nonzero norm on λ^{−4}O is 1/81, attained exactly at the six units divided by 9 | exact |
| 7.8 | (7.7) | sanity: K = U = 400 (109 squarefree rows, 138 squarefree columns, even columns included); LHS/((K + U)Σ\|β\|²) ≤ 0.254 for six coefficient choices | bounded |
| — | (7.8) | sanity (`test_sec7_completed_sanity.py`): Σ_k \|T(X; k, f)\|² ≈ (0.2–0.35)ℋ for X = 1500, ℋ ≤ 1000, f ∈ {1, 7 + 9ω} | bounded |
| 7.10 | after (7.8) | χ_n(v²)χ_n(f)⁴ = χ_n(fv²)⁴ | 4e−15 |
| 6.14 | §7.3, quadratic twists | γ₃(n) depends only on n mod 4 (90 squarefree n, 12 classes), and χ_n(4) only on n mod 2 | exact |
| 6.12 | (7.10) | cubic twists: μ_K(n)γ₂(n) = χ_n(−1)G(n)α(n)γ₋₁(n) for 70 squarefree n, including composite n | 4e−15 |
| 6.13 | (7.11) | γ₁(n) = conj(χ_n(4))γ₂(n)²γ₃(n) and γ₂(n)³ = μ_K(n)α(n) for 70 squarefree n | 6e−15 |

### Section 8: cube removal and the descent (`test_sec7_8.py`)

| Test | Guide | What is checked | Result |
|---|---|---|---|
| 8.1 | (8.1) | the cube inversion for three pairs (k, f) at X = 4000, with six nontrivial cubes | 2e−15 |
| 8.2 | (8.2)–(8.11) | the inequalities between scales on 200 000 random admissible parameter sets, and the retained-support chain (8.8)–(8.9) on a further 200 000 sets with a retained r. A control that keeps only F′ ≤ ℋ′ instead of (8.8) violates X′/N(j) ≥ ΣF/ℋ in about 9000 sets | no violation |
| 8.4 | (8.7) | the local cancellation identity | exact |
| 8.3 | (6.11), (8.13) | the initial scales (D²/(HB²))/(D/B) = D^{−ϑ}/B (sympy) | symbolic |

## References

- OpenAI, *The Quasi-Riemann Hypothesis*, manuscript dated October 5, 2026,
  [github.com/openai/math](https://github.com/openai/math/tree/main/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026).
- A. Dunn and M. Radziwiłł, *Bias in cubic Gauss sums: Patterson's conjecture*, Ann. of Math. (2) 200 (2024),
  967–1057, [arXiv:2109.07463](https://arxiv.org/abs/2109.07463). Section 5.1 and Appendix A give the theta
  function and its cusp coefficients.
- L. Goldmakher and B. Louvel, *A quadratic large sieve inequality over number fields*, Math. Proc. Cambridge
  Philos. Soc. 154 (2013), 193–212, [arXiv:1112.1642](https://arxiv.org/abs/1112.1642).
