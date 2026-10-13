# t0: the numbered displays of the Introduction, (0.1)–(0.10)

`t0_intro.py` tests the ten displays of the Introduction; `rows.json` has one row per display. Run `./run_all.sh`
(about 10 seconds; needs PARI/GP). The helpers of `../../part_A_lines_0001_0549` (Eisenstein-integer arithmetic and
sextic symbols) are imported read-only.

| Display | Label | Kind | Test |
|---|---|---|---|
| (0.1) | `eq:intro-curve` | D | discriminant −432d² (sympy) and j = 0 (PARI), 40 values of d |
| (0.2) | `eq:intro-cm-map` | D | (ωx)³ + d = x³ + d (sympy); order 3 |
| (0.3) | `eq:intro-ap` | E | brute-force count = −Legendre sum = PARI ellap, 1366 pairs (d, p) |
| (0.4) | `eq:intro-L` | D | PARI ellan = product of the local factors (0.5), m ≤ 2500, 10 values of d |
| (0.5) | `eq:intro-local` | C | all 1200 bad primes for 0 < \|d\| ≤ 200 are additive; a_p = 0 there; d = 16 is good at 2 |
| (0.6) | `eq:intro-fe` | C | Λ(s) = εΛ(2 − s) with PARI lfun, 8 curves, 4 complex points |
| (0.7) | `eq:intro-fe-u` | E | Λᵘ(s) = εΛᵘ(1 − s), same curves and points |
| (0.8) | `eq:intro-psi` | D | conj((4d/n)₆) = ((4d)⁵/n)₆ exactly, 2258 pairs |
| (0.9) | `eq:intro-hecke` | E | a_m = Σ b(m₁) Σ_{N n = m₂} n conj((4d/n)₆), exactly in ℤ[ω], m ≤ 2500 |
| (0.10) | `eq:intro-inverse` | E | Σ μ_K(n) n conj((4d/n)₆) = coefficient of ∏(1 − a_p X + p X²), exactly |

The checks with a wrong variant (0.3, 0.6, 0.7, 0.9, 0.10) require that variant to fail. Last run (2026-10-13):
10 checks, 0 failures.
