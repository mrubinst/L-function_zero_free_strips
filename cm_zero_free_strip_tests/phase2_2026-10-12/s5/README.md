# s5 tests: the plain fourth moment with an angular factor

These tests cover Appendix `app:s5-plain` in `../inserts.tex`. Run `./run_all.sh`; it takes about 20 s.
Outputs are in `outputs/`. Each script prints `RESULT|id|PASS/FAIL|...` lines, one line per falsifier, and
a final `SUMMARY`. A falsifier is a deliberately wrong variant. It "fails as required" when its error
exceeds the stated threshold.

Results (2026-10-12): **106 checks PASS and 98 falsifiers fail as required, with 0 failures overall.**

| script | claim (inserts.tex) | checks / falsifiers | worst error |
|---|---|---|---|
| `t1_dual_profile.py` | Lemma s5-sharp: Mellin transform of J_1(2 sqrt y) is Gamma(s+1/2)/Gamma(3/2-s) (control: J_0 gives Gamma(s)/Gamma(1-s)); closed forms x^{1/2}e^{-x} -> y^{1/2}e^{-y} and x^{3/2}e^{-x} -> y^{1/2}e^{-y}(2-y); Mellin inversion of MW(1-s)G(s) on Re s = -0.3, 0.4, 1.0 reproduces W^sharp for a C_c^oo bump; W^sharp(y) ~ MW(3/2) y^{1/2}; annular coefficients O(2^-j); Stirling exponent 2 sigma - 1 | 8 / 5 | 2.5e-31 (Bessel-Mellin), 3.7e-9 (inversion) |
| `t2_hankel.py` | eq. s5-hankel: planar transform of (z/\|z\|)^e W(\|z\|^2) under e(-zw), and ghat(0) = 0 (bump and x^{1/2}e^{-x}, e = +1 and -1) | 8 / 13 | 7.8e-14 |
| `t3_reflection.py` | Lemma s5-reflection: X^{-1/2} sum Xi W = eps Y^{-1/2} sum conj(Xi) W^sharp, with Y = 3N(f)/(4 pi^2 X) and eps = e gamma(phi)(f/\|f\|)^e, for 12 primitive characters (alpha, conj alpha, twists by sextic symbols mod primes of norm 7, 13, 19, a product, the cubic character mod 2) and three profiles; \|eps\| = 1; the conductor changes only at lambda | 51 / 68 | 3.3e-11 |
| `t4_lattice.py` | Lemma s5-lattice: unmasked L(X,t) bounded and tending to 0 for X up to 3e5 (4 characters theta, t = 0 and 3.7, e = +1 and -1); calibration of the e = 0 volume term pi/(6 sqrt3) X \|I(t)\|; mask inclusion-exclusion; the centered and uncentered products | 21 / 3 | 1.6e-13 (identity) |
| `t5_row_reflection.py` | Lemma s5-rowchar and Step 3: OAI78's reflection with deleted and restored factors (L12736-12744) in CM form on 5 rows with redundant zeros; the bound C_k^CM N(r_k) <= 27 C_k N(R_0k); mask deletion (2.1a) with alpha | 15 / 8 | 1.3e-15 |
| `t6_kappa.py` | Step 6: the kappa-dependent inequalities of OAI78 ((3.6), L12947-12950, (3.11), comparison margin) at kappa = 3/4 and 1, by random sampling; s_kappa = 7/8 | 3 / 1 | 0 violations |

The falsifiers include:
- the finite-order gamma quotient Gamma(s)/Gamma(1-s) and the J_0 kernel (t1, t2, t3, t5);
- the prefactor -i and the unconjugated angle (w/\|w\|)^e (t2);
- the conductor scale of the finite-order part, i.e. no change at lambda (t3), and OAI78's (C_k, R_0k)
  in place of (C_k^CM, r_k) (t5);
- root numbers without the factor e or without (f/\|f\|)^e (t3);
- no angular factor (e = 0), which has a volume main term of size ~0.067 X at t = 0 (t4);
- the claim O(2^-j) for OAI78's J_0 kernel, and the exponent 2 sigma in Stirling (t1);
- kappa = 0.70 (t6).

Remarks:
- `t5` uses row characters of the form tau(n) prod_p chi_p(n)^{v_p}, i.e. residue symbols *modulo* the moving
  primes, which can be evaluated at any element. OAI78's psi_k(n) = tau(n) chi_n(k) differs from such a
  character by a fixed-ray factor (sextic reciprocity), so it lies in the same coefficient class. Primitive
  conductors are computed exactly (`hecke.py`: the finite part on (O/m)^x, then reduction to the conductor).
- The analytic test profiles x^{1/2}e^{-x}, x^{3/2}e^{-x} and the log-normal one are not annular. Lemma
  s5-reflection only needs (z/\|z\|)^e W(\|z\|^2) to be Schwartz, which holds for them. Annular profiles are
  tested in t1 and t2; their W^sharp decays too slowly (like exp(-c y^{1/4})) for a direct
  high-precision test of the full reflection sum.
- Helpers: `s5lib.py` (Eisenstein arithmetic, sextic symbols, e(z), reporting) and `hecke.py` (Hecke
  characters of infinity type (a/\|a\|)^e, conductors, Gauss sums).
