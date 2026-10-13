# Display tests, agent t7: Appendices C and D of the CM guide

Scope: the 30 numbered displays of `sources/inventory/inventory_t7.json`, that is Appendix C (Hecke
L-functions of infinity type zbar/|z|: C.1-C.18) and Appendix D (the zero detector for a CM target: D.1-D.12) of
`cm-eleven-twelfths-guide.tex` (version of 2026-10-13). One row per display is in `rows.json`
(label, printed number, kind, what was done, script, result).

Run everything with `./run_all.sh` (python3 with numpy, scipy, sympy, mpmath; PARI/GP 2.17 `gp`). It runs the five new
scripts of this folder and then `rerun_existing.sh`, which re-runs, read-only, the existing suites
`tests_cm/phase2_2026-10-12/s2` (t1-t4) and `.../s3` (t1-t5), writing their output to `outputs/rerun_*.txt`.
`T7_SKIP_RERUN=1 ./run_all.sh` runs only the new scripts (about 20 s). Every script exits 1 if a check fails or a
falsifier (a deliberately wrong variant) does not fail; with `CM_TESTS_FORCE_FAIL=1` the first check counts as failed
(self-test of the exit status; checked: exit 1). Helpers: `t7lib.py` (imports `zw.py` read-only from
`tests_cm/part_B_lines_0550_0997`); `c4_fe_pari.py` and `d1_members.py` also import `s3lib.py` read-only from
`phase2_2026-10-12/s3`. No other suite was modified.

## Counts by kind

| kind | count | displays |
|---|---|---|
| E (exact identity) | 19 | C.1, C.3-C.10, C.12-C.15, C.17, C.18, D.6, D.7, D.9, D.10 |
| I (explicit inequality) | 1 | D.8 |
| A (asymptotic, unspecified constants) | 1 | C.16 (size check; proved in Lemma s2-logcontrol) |
| D (definition) | 8 | C.2, C.11, D.1-D.5, D.12 |
| C (quoted) | 0 | (C.3 and C.7 are quoted from [QF, Lemma 17, Example 20], but are E and were tested; both agree with the source) |
| R (restatement) | 1 | D.11 = eq:d158 |

All 20 E/I displays were tested, and all pass. No display failed, and none is suspected wrong.

## New tests (this folder)

The test characters (`t7lib.test_chars`) are eight primitive phi with phi(zeta) = zeta, Q = 7, 9, 21, 28, 75, 91,
117, 225. Their lambda-exponents are 0, 1, 2, and they involve the primes over 2, inert 5, and split 7 and 13. A-D are
the characters of s2, built independently here, and their root numbers agree with the s2 README.

| id | displays | what | result |
|---|---|---|---|
| T7-C1a/b | C.1 | psi((a)) = phi(a) conj alpha(a): independent of generator, multiplicative, psi((b))/psi((c)) = conj alpha(b/c) for b = c mod f; alphabar: phi = unit mod 3 | PASS 1.4e-15; falsifier phi(a)alpha(a): 1.73 |
| T7-C2 | C.2 | gamma(alphabar) = -1 exactly (in Z[zeta_18]); W independent of the generator of (f) | PASS |
| T7-C3 | C.3 | twisted Gauss sum for every h mod f (incl. h not coprime), exactly in Z[zeta_M], M = lcm(6, Q); abs(gamma) = 1 exactly | PASS (exact); falsifiers phi(h) and an imprimitive phi fail |
| T7-C7 | C.7 | Gaussian transform: sympy closed form plus 1-D quadratures | PASS (2e-25); falsifier e(z) = exp(2 pi i Im z) |
| T7-C8a/b | C.8 | Wirtinger step (sympy); transform of conj(w) e^{-2 pi v abs(w)^2/sqrt3} by quadrature | PASS (2e-25) |
| T7-C9a | C.9 | Poisson step for every residue t mod f, 3 values of v | PASS 1.5e-13; falsifiers Jacobian, e(-ty/f) |
| T7-C9b | C.4, C.9 | both lines of eq:d153 = direct Theta(v) = right side of eq:d148, 5 values of v | PASS 1.9e-15; 4 falsifiers fail |
| T7-C5 | C.5, C.10 | Mellin integral (trapezoid in log v on direct theta sums) = element sum = 6 x ideal sum, Re s > 1, psi and conj psi | PASS 1.6e-12 |
| T7-C11 | C.11 | B, B^- by quadrature = incomplete-gamma series used by s2/t1 | PASS 4e-14 |
| T7-C12 | C.12 | three expressions of eq:d156 by separate quadratures/series | PASS 6.5e-15 |
| T7-C13 | C.13 | split identity and its conjugate, also at Re s < 1; constants (kappa, W/Q) symbolically | PASS 1.1e-11 |
| T7-C6 | C.6 | Q^{s/2} I_psi(s) = W Q^{(1-s)/2} I_conjpsi(1-s) with both sides from direct theta quadratures (independent of eq:d148) | PASS 7e-15 |
| T7-C6p | C.6 | PARI lfuncheckfeq with eps = W = -gamma sqrt(Q)/f supplied (not lfunrootres); coefficients and W built in gp at 60 digits from exact data; PARI L(2.5+it) = Dirichlet series | PASS (-253 bits; falsifiers >= -5 bits) |
| T7-C14 | C.14 | eq:s2-fe-unc with PARI values of L(s,psi), L(1-s,conj psi); sympy equivalence with eq:d150 | PASS 9e-16 |
| T7-C15 | C.15 | left-line identity at 6 heights; max abs(L)/(Q^{3/5}(3+abs(t))^{6/5}) = 0.30 | PASS 4e-16 |
| T7-C16 | C.16 | size check (A): winding number 0 on radius 2-a-2e, max(abs(L)+abs(1/L)) <= 8.4 on radius 2-a-6e (a = 0.6, e = 1e-4, t = 0, 10, 40) | size check |
| T7-D6 | D.6 | direct sum = Euler product, times Moebius series = 1, 4 members (finite/CM, sigma = +-1) | PASS 4e-10 |
| T7-D9 | D.9 | both expressions of eq:d167, eps_a, generator independence, infinity type on x = 1 mod m, agreement with CM member | PASS 7e-15 |
| T7-D10 | D.10 | 4 CM members of X_u (sigma = +-1, Q = 63 ... 14400): lfuncheckfeq -254 bits, abs(eps) = 1, identity at 4 points | PASS 6e-16 |
| T7-D8 | D.8 | same members, a in {0.51, 0.75, 1}, e in {1e-4, 9.9e-4}, 7 heights: ratio <= 0.96 (C_0 = 1 works), equal to (2 pi)^{1-2a-12e} x gamma quotient | PASS; falsifiers fail |

## Existing tests re-run (2026-10-13, all exit 0)

s2/t1_theta_fe (C.4, C.6), s2/t2_growth (C.15, Lemma s2-gamma), s2/t3_reflection (C.18), s2/t4_conductors (C.1(iii),
C.17), s3/t1_labels (D.3-D.6, D.12), s3/t2_conductors (C.17, D.7), s3/t3_fe (C.14, C.17, D.7, D.10), s3/t4_gamma (D.8),
s3/t5_citations. Outputs: `outputs/rerun_*.txt`; summary `outputs/rerun_summary.txt`.

## Notes and caveats

- The definitions D.1 (Theta_0) and D.12 (witnesses) were compared with OAI78: L3333-3335 (Theta = <eta, T^>) and
  eqs. old-eq:1.1a, old-eq:1.1 at L706-714 (M_psi, S_psi), with psi = alpha^sigma theta chi(u)^sigma at base U. The
  zero convention L636-639 and the left-line bound L1446-1453 cited near D.4 and C.15 are as quoted.
- C.3 = [QF, Lemma 17] verbatim; C.7 = [QF, Example 20] after removing the factor 2/sqrt3 that [QF] puts into
  Phi-hat; the Poisson formula used in C.9 is [QF, Theorem 19].
- Characters F and G (inert 5) have real theta series and W = -1, so Theta(Q^{-1/2}) = 0 and L(1/2, psi) = 0. This is
  forced by the functional equation and is not an error; T7-C9b avoids v = Q^{-1/2} for relative errors.
- The conductors of X_u members come from s3lib's randomized congruence search (50 samples per candidate modulus);
  they are checks, not certified conductor computations (as noted in the s3 README).
- PARI values are computed through PARI's own use of the functional equation. The independent evidence is
  lfuncheckfeq, which starts from the Euler coefficients, and T7-C6, which uses only direct theta sums.
- C.16 (and D.11) are analytic; the size check is not a proof. They are proved in Lemma s2-logcontrol, following
  [OAI78, Lemma logarithmic-control, L1554-1600].

## Runtime

Last full run (2026-10-13, `./run_all.sh`, exit 0): 689 s in total. The new scripts take 20 s: c1 1 s, c2 1 s,
c3 12 s, c4 4 s, d1 2 s. The re-runs take 669 s: s2/t1 33 s, s2/t2 181 s, s2/t3 384 s, s2/t4 0 s, s3/t1 4 s,
s3/t2 22 s, s3/t3 21 s, s3/t4 23 s, s3/t5 1 s. See `outputs/run_all_summary.txt`.
