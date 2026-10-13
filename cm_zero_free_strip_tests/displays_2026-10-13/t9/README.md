# t9: displays of Appendix F (the plain fourth moment with an angular factor)

Agent t9 covers the 26 numbered displays F.1-F.26 of `cm-eleven-twelfths-guide.tex`, Appendix `app:s5-plain`:
the setting, the primitive CM row character, the reflection of a plain sum with an angular factor
(Jacobi-Anger and the Hankel/J_1 transform), the lattice sums, the coefficient class, and the proof. There is one
row per display in `rows.json`. Run `./run_all.sh`. It runs the six new scripts `d1`-`d6`, then re-runs the six
scripts of `tests_cm/phase2_2026-10-12/s5` read-only, writing their output to `outputs/rerun_s5_*.txt`. It exits 1
if any script exits nonzero or prints a failure line. Every new script exits 1 on a failed check or on a
falsifier that does not fail. They use the s5 reporting helpers through `t9lib.py`, and
`CM_TESTS_FORCE_FAIL=1` gives exit status 1.

**Runtime (2026-10-13):** 53 s in total. The new scripts take 33 s and the s5 re-runs take 20 s; the machine load
average was about 12.

## Counts by kind

| kind | count | displays |
|---|---|---|
| E (exact identity, checked on finite data) | 17 | F.1, F.4, F.5, F.8-F.15, F.17, F.20, F.21, F.24-F.26 |
| I (explicit inequality) | 0 | none |
| A (asymptotic bound; size check only) | 4 | F.3, F.16, F.19, F.22 |
| D (definition) | 5 | F.2, F.6, F.7, F.18, F.23 |
| C (quoted from a source) | 0 | none |
| R (restatement) | 0 | none |

**All 17 E displays were tested, and all pass.** The new scripts make 78 checks, all PASS, and 44 falsifiers, all
of which fail as required. The s5 re-runs make 204 checks (t1-t6), with 0 failures.

- **F.8 (eq:s5-reflection)** was already tested end to end by s5/t3 and s5/t5. Both were re-run and pass (max err
  3.3e-11, which comes from the lognormal profile via Mellin-Barnes; the closed-form profiles reach 1.8e-15).
- **F.13 (eq:s5-hankel) and F.17 (eq:d182).** The existing tests (s5/t2, s5/t1 1a) only covered X = 1 and x = 1. d3
  and d4 now test X = 0.37, 4.2 and x = 0.3, 2.5. For F.17 the integral is computed directly in y, at Re s = -0.45
  and 0.7, near both ends of the strip.
- **The intermediate steps had no direct test before.** Earlier tests checked them only implicitly, through the
  final reflection identity. They now have their own tests:
  - F.4, F.5, F.24-F.26: `d2_rowchar.py`;
  - F.9-F.12, F.14, F.15: `d3_poisson.py`;
  - F.20, F.21: `d5_lattice.py`;
  - F.1: `d1_alpha_mult.py`.

  In F.10 and F.11 the planar transform comes from its definition, by 2-d polar quadrature, so the Poisson steps
  are tested independently of F.13. F.4 is checked pointwise on all 1510 ideals of norm at most 2500. This
  includes rows where lambda moves from the conductor of the finite part into r_k (the guide's example chi_pi^5,
  N(pi) = 7, e = 1: N(f) = 7, N(g) = 63), and rows where it moves the other way (N(f) = 468, N(g) = 52).

## Size checks for the A displays

| display | check | result | proof |
|---|---|---|---|
| F.3 | Case (1) on 451 rows (all ideals coprime to 6 with N(k) <= 1500, including k = 1). The mean of the fourth moment stays bounded for X1 X2 from 1e2 to 9e8, far outside the large-sieve range. | Mean at most 9.1e-4, no growth. With e = 0, the principal row gives \|S(X)\|^2 ~ 135 at X = 3e4. | Proved in the guide's sec:s5-proof. |
| F.16 | The proof's explicit bound \|(y d/dy)^j W#(y)\| <= y^{1/2} int \|W_j\| x^{1/2} dx, at 81 points, j <= 2. The decay of y^A \|(y d/dy)^j W#\| per decade, up to y = 1e6. | Holds. | Lemma s5-sharp(c). |
| F.19 | Re-run of s5/t4. | Unmasked sums at most 2.5, decaying to X = 3e5. Masked sums at most 3.3. With e = 0 there is a volume term of about 0.067 X. | Lemma s5-lattice. |
| F.22 | X sum (1 + sqrt X \|w\|/12)^{-3}, and the transform bound for the bump at t <= 64. | The sum is at most 523 for 1 <= X <= 1e6. With J0 = 2 it diverges. The transform bound holds. | Lemma s5-lattice. |

## Definitions

- **F.2.** Matches OAI78 old-eq:1.1 (L704-708) with alpha(l)^e inserted, and is implemented in d6.
- **F.6.** Implemented in s5 `hecke.gauss`; \|gamma\| = 1 for 12 characters, and eps(alpha) = 1.
- **F.7.** Implemented in s5/t1, s5/t3 and d4; the closed forms are checked.
- **F.18.** Implemented in d5 and s5/t4.
- **F.23.** Compared with OAI78 eq:centered-coefficient-form (L13051-13064): it is the same coefficient with
  alpha(l1 l2)^e, and q_u = N(u). The allocation scalar is checked in d1.

## Failures and suspected errors in the guide

**None.** No display failed, and none appears to be wrong.

Two development runs did fail, both because of bugs in the tests themselves, not in the guide. Both were fixed
before the final run:

- F.10 failed because the lattice box for y was too small for residues t with large \|t\|.
- F.17 failed because mpmath's quadosc, run with Bessel-zero spacing in y, did not converge. The fix subtracts the
  leading term near 0 and integrates in v = sqrt y.

## Scripts

| script | displays | checks / falsifiers |
|---|---|---|
| `d1_alpha_mult.py` | F.1, F.23 (allocation scalar) | 4 / 3 |
| `d2_rowchar.py` | F.4, F.5, F.24, F.25, F.26 | 37 / 10 |
| `d3_poisson.py` | F.9-F.15 | 22 / 23 |
| `d4_mellin.py` | F.17, F.16 (size check) | 5 / 2 |
| `d5_lattice.py` | F.20, F.21, F.22 (size checks) | 8 / 5 |
| `d6_fourth_moment.py` | F.3 (size check) | 2 / 1 |
| s5 `t1`-`t6` (re-run) | F.7, F.8, F.13, F.16, F.17, F.19 | 204 lines, 0 failures |

Helpers are reused read-only from `tests_cm/phase2_2026-10-12/s5`: `s5lib.py` and `hecke.py`. No bytecode is
written into that folder.
