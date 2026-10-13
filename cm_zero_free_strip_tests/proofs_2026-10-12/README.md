# Checks of the Part I proofs added to the CM guide on 2026-10-12

On 2026-10-12 the Part I proofs of `cm-eleven-twelfths-guide.tex` were completed. Steps identical to the zeta case
are cited from OpenAI's paper; steps in which the angular factor alpha changes anything are written out. The new
material was written in four blocks:

| block | guide location | content |
|---|---|---|
| b1 | Sections 3, 5, 6 | proof of the point count and the Euler products; the extraction in detail (uniform in d); holomorphy and continuation of 1/L |
| b2 | Section 7.3 | the CM paired identities (forward and reverse), twisted multiplicativity, the Poisson reduction to the dual mean square |
| b3 | Section 8, Appendix A | the CM theta transformation (Proposition, second derivative), uniformity, decay of the weight |
| b4 | Sections 9.1, 10.5, 10.6 | the completed estimate (with the bridge to the large sieve of Goldmakher-Louvel), cube removal, the transfer estimate |

Each block folder has a README listing every new identity or inequality, its test, and whether an earlier test of
the suite already covers it.

```sh
./run_all.sh          # about 2 minutes; outputs in outputs/; needs python3 (numpy, sympy, mpmath, scipy) and PARI/GP
```

**Last run (2026-10-12): every script passed, and every falsifier failed as required.** Main results:
- b1: every step of the point-count proof holds exactly in Z[omega] for 4412 pairs (d, pi); the wrong variants
  (opposite sign of pi, unconjugated symbol, chi_pi(d) for chi_pi(4d)) fail. The local factors at p | 6d agree with PARI
  for 20 918 pairs. The counting inequalities of the extraction hold.
- b2: the forward and reverse CM paired identities hold in class form (error about 2e-15; the zeta coefficients fail).
  End to end, the smoothed mean square equals the zero-frequency term plus the final dual form to 3e-15 at D = 120,
  H = 40 and 170; four wrong variants fail.
- b3: T(X; Psi) computed directly from its definition equals the CM transformation formula to 1.1e-12, for three
  twists covering all three cusp types; the variants conj(alpha(c))^2, the zeta prefactor -i/81 and the
  one-derivative kernel fail. The Mellin, chain-rule and constant computations hold symbolically.
- b4: the paired identities with a twist (432 characters), the first and second Poisson identities of the transfer,
  the preimage regrouping tau(r) 1_{f'|k'}, cube removal with conj(alpha)^6, the Hecke-family facts behind the use of
  the large sieve, and the weight decay all pass; the zeta-form variants fail.

**Exit status (added 2026-10-12).** Every script exits 1 if a check fails or a falsifier does not fail, and
`run_all.sh` exits 1 if any script exits nonzero or prints a failure line. Most of these scripts print values rather
than verdicts, so the criteria were made explicit in a final `EXIT STATUS` line (no printed computation was changed):
- `b1/t_points.py`: every failure counter is 0, the falsifier J = +pi never holds, every other falsifier misses some
  case; `b1/t_extract.py`: every True/False printed is True and P(7) = 2 (the largest-Y line of T7 is a sanity
  number); `b1/t_bd.gp`: no failure, the falsifier fires at every good prime, and the GP block ran to the end.
- `b2/*.py`: identity errors <= 1e-10, falsifier errors > 1e-6 (at D = 40 the k f^2 falsifier does not fail, see
  the b2 README, so use D >= 100 as `run_all.sh` does).
- `b3/b3_symbolic_checks.py`: exact results equal their targets; mpmath errors (40 digits) <= 1e-20 (1e-30 for the
  closed form of (4)); numeric ratios of (3) within 1e-12 of 1; fitted growth exponents of (6) and (8) within 1e-3.
  `b3/t_b3_prop_endtoend.py`: the PASS criterion it prints, plus 0 failures of the supplement check.
- `b4/b4_checks.py`: every [FAIL] line and every falsifier that DID NOT FAIL counts; `b4/b4_weight_decay.py`: the
  C9 line, plus the zeta-quotient falsifier, which must differ from V# by more than 1e-3 at some x (one new line).
With `CM_TESTS_FORCE_FAIL=1` every script counts its first check as failed (self-test of the wiring).

Not testable by computation: the analytic inequalities themselves (the large sieve, the completed estimate, cube
removal, the transfer estimate), beyond the size checks listed in the block READMEs.
