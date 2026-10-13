# t6: displays B.1-B.38 (Appendix `app:s1-general`, the CM transformation formula at full generality)

Guide: `cm-eleven-twelfths-guide.tex` (2026-10-13 numbering), lines 3793-4380. Inventory: `sources/inventory/inventory_t6.json`
(38 displays). One row per display is in `rows.json`. Script paths in `rows.json` are relative to `tests_cm/`.

Run everything with `./run_all.sh`. It runs the re-runs of the existing suites and the new scripts in parallel and takes
**107 s wall time** (about 620 s CPU) on this machine. All checks pass and the overall exit status is 0. Outputs are in
`out/`. Every script exits nonzero on any failure, and the new scripts honour `CM_TESTS_FORCE_FAIL=1` (self-tested).

## Counts by kind

| kind | count | displays |
|---|---|---|
| E (exact identity) | 18 | B.8-B.10, B.12-B.18, B.20, B.21, B.23, B.24, B.29-B.32 |
| I (explicit inequality) | 2 | B.19 (majorant `eq:d132`), B.22 (`eq:d134`) |
| A (asymptotic, unspecified constants) | 8 | B.3, B.11, B.25, B.27, B.28, B.34, B.35, B.38 |
| D (definition) | 10 | B.1, B.2, B.4-B.7, B.26, B.33, B.36, B.37 |
| C, R | 0 | (the OAI78 statements are restated in CM form and compared line by line with `paper.tex`; see the rows) |

**All 20 E/I displays were tested, and all pass.** Six of the A displays have a finite size check (B.3, B.11, B.25, B.27,
B.28, B.38's exponent identity). Two (B.34, B.35) are marked "not testable by computation". The D displays were checked
for consistency with OAI78 and with the implementations, with a numerical check where one exists (B.2, B.4, B.7, B.26,
B.36, B.37).

**Failures: none. Displays I think are wrong: none.**

## Re-runs of existing tests (`rerun_existing.sh`, outputs `out/rerun_*.txt`)

| script | covers | result now |
|---|---|---|
| `phase2_2026-10-12/s1/t_s1_endtoend_active.py` gauss/bessel x pi 0/1 | B.9, B.10, B.24 (ungrouped), B.29, B.30 | PASS: rel. err 1.7e-12 (Gaussian) and 5.2e-8 (Bessel); marked 8.1e-13 / 1.8e-8; all 6 controls fail (0.43-4.4) |
| `phase2_2026-10-12/s1/t_s1_kernel.py` | B.1, B.7, B.11, B.15/B.23 constants, B.18 mode, B.21 growth, B.25 | PASS (mode 4e-31, constant 4e-31, line independence 3.5e-8) |
| `phase2_2026-10-12/s1/t_s1_local_phases.py` | local data in B.10, B.30 | PASS (6.3e-14; falsifiers 1.73) |
| `proofs_2026-10-12/b3/b3_symbolic_checks.py` | B.23 (constant), Gamma shifts of B.21 | PASS (26 checks) |
| `proofs_2026-10-12/b3/t_b3_prop_endtoend.py` | B.9 with `P` empty, cusps with lambda dividing c | PASS (1.1e-12) |

## New scripts (this folder)

- `t6_lib.py` is the shared setup. It holds the theta coefficient tables d_0, d_+, d_-, the 252 translates and the cusp
  data of the s1 end-to-end setting (L = 18, base twist A, one prime of norm 7). The code is copied from
  `t_s1_endtoend_active.py`. `fastgauss.py` and `zw.py` are imported read-only from `phase2_2026-10-12/s1/`.
- `t6_theta_translates.py [0|1]` (40 s each, both primes of norm 7) checks the following pointwise:
  - T1 `eq:d129`: Fourier inversion at all 2268 classes mod 18p, for j = 0..5 (error 1.3e-13).
  - T2 `eq:d130`: the cusp expansions for all 11 matrices H, from the identity expansion evaluated at H(z', v') (1.8e-15).
  - T3 `eq:s1-automorphy`: Kubota automorphy at all 252 translates and at 8 extra cusps with lambda dividing c (2.3e-13).
    The unconjugated kappa fails.
  - T4 `eq:s1-cusp` and `eq:d131`: the second z-bar derivative at 260 cusps (3.6e-13). The falsifiers (cv)^-2 and
    conj(alpha)^2 fail.
  - T5 `eq:d132`: the size bound on all 60770 table entries, at most 6 units per index, and the partial sums.
  - T6a `eq:s1-fe` and `eq:s1-dual`: J_2(s) by quadrature of f(v), using the identity expansion above v0 and the 252
    cusp expansions below it, against the dual Dirichlet series at s = -2 and -2.5+1.5i (1.9e-13).
  - T6b `eq:s1-direct`: the same quadrature at s = 4 and 3.5+2i, against D_2 L_2 computed as separate series (2e-12).
  - T6c: the fast gamma_2 against its brute-force definition.
- `t6_mellin_rep.py` (20 s) checks the Mellin representation `eq:d127` for two non-compact tests, two values of X,
  j = 0..5 and two lines (6e-13). It also checks the Mellin form of `eq:s1-central` with Psi = nu chi_.(m), for rows m with
  a square factor, a factor at lambda and 2, and a unit (1.3e-12). Finally it checks that the inner integral of
  `eq:d144` with `eq:d145` (mark D = p) equals the marked sum with V_G (9e-12).
- `t6_kernel_bookkeeping.py` (16 s) runs:
  - K1 `eq:s1-mellin-tail`, through its integration-by-parts identity and the bound with an explicit constant.
  - K2 M_{A,B}(V_G).
  - K3 `eq:d128`, normalised by M_{A,ceil(4A)+j+2}.
  - K4 the `eq:s1-stirling` sup.
  - K5 `eq:d134`.
  - K6 sympy bookkeeping for `eq:d132`, `eq:d133`, the substitution in `eq:s1-fe`, the exponent inequality of
    `eq:s1-cm-completed` and the exponent identity of `eq:d146`.
  - K7 `eq:d140` and `eq:d141`, exactly in Z[omega].
  - K8 |C| <= 1/81 and C = -i conj(alpha(c))^2 zeta, for `eq:s1-scalar` and `eq:d139`.
  - K9 `eq:d122`: phi computed from sextic symbols is 18-periodic and equals the closed forms used in the tests.
- `t6_moment_sanity.py` (3 s) is a size check of `eq:d137` (O-dyads) and `eq:s1-cm-completed`. It sums |C^(2)_V(Z;m)|^2
  over all rows N(m) <= 8Z for Z = 100 and 300, with V = exp(-y-1/y). Every ratio to the claimed bound is at most 0.021.

## Remarks (no display is wrong)

1. **Notation, not an error.** The bullet before `eq:s1-automorphy` (guide line 3952) says G is in Gamma_1(3). This
   follows OAI78 (paper.tex L1891), where Gamma_1(3) = {g = I mod 3}. Elsewhere the guide calls this group Gamma(3)
   (line 1432, and line 3514 in Appendix A), and the usual meaning of Gamma_1(3) is larger. The matrices G of the tests
   satisfy G = I mod 3 in all four entries. The Notation paragraph of Appendix B could say that Gamma_1(3) is OAI78's
   name for Gamma(3).
2. **Coverage limits.** The marked formula (B.29, B.30) is tested end to end only with one marked prime and an empty base
   set (Q = {p}, P_0 = {}). Larger Q need tables to N(c) of order 400, which costs far more. The inclusion-exclusion
   step is finite algebra, and the local factors of the active mark are tested at 15 primes in `t_s1_local_phases.py` (e).
   In the active-prime end-to-end test every translate has (lambda, c) = 1. The cusps with lambda dividing c are covered
   by T3 and T4 here and by `t_b3_prop_endtoend.py` with P empty.
3. The figures quoted in the guide's "Numerical checks" paragraph agree with the re-runs: 1.7e-12 and 5.2e-8, the spread
   2.4e-12, the marked error 8e-13, and control errors between 0.43 and 4.4 over 120 control evaluations.
