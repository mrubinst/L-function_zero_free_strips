# s4 tests: the CM inverse moment (Appendix "The inverse moment for co-oriented coefficients")

```sh
python3 s4_checks.py          # about 20 s; output in s4_checks_output.txt
```
The script uses the helpers of `tests_cm/part_B_lines_0550_0997` (`cmlib.py`, `common.py`, `zw.py`,
`fastgauss.py`). The conventions are OAI78's (= the guide's): primary = 1 mod 3; chi_n(x) = (x/n)_6, zero on
nonunits; e(z) = exp(4 pi i Im z/sqrt 3); gamma_j(n) = N(n)^{-1/2} sum chi_n(x)^j e(x/n);
G(n) = conj(chi_n(4)) gamma_3(n); alpha(n) = n/|n|.

**Last run (2026-10-12): 12 of 12 checks PASS, and all 12 falsifiers fail as required (SUMMARY 24/24).**

| id | identity in inserts.tex | test | result | falsifier (must fail) |
|---|---|---|---|---|
| S4.1 | (s4-S1)-(s4-S4), the four switching identities, pointwise | all 366 squarefree primary u coprime to 6 with N(u) <= 1300; 157 of them are composite, with up to 3 prime factors, and inert primes are included | PASS, max err 2.5e-15 | zeta cubic class conj(alpha)gamma_2 on the left; side 2 without conjugation; G in place of conj(G); zeta Moebius class mu. All have err 2 |
| S4.2 | (s4-S5), the anti-oriented collapse mu alpha gamma_{-1} = gamma_2 chi(-1) conj(G) | same 366 u | PASS, 1.8e-15 | claim that the anti-oriented class is conj(alpha)^2 gamma_2: err 2 |
| S4.3 | twisted multiplicativity a(vw) = a(v)a(w)chi_w(v)^4, a = conj(alpha)^2 gamma_2 nu', nu' = chi_.(lambda) | 159 coprime pairs, N(vw) <= 1300 | PASS, 3.2e-15 | chi_w(v)^2: err 1.73 |
| S4.4 | ray-quotient forms of the switches: a(z1)conj(a(z2)) gamma(chi_z1 conj chi_z2) = m(z1)conj(m(z2)) G([z1][z2]^-1), and m(z1)conj(m(z2)) gamma(conj chi_z1 chi_z2) = a(z1)conj(a(z2)) G([z2][z1]^-1), with m = mu conj(alpha), a = conj(alpha)^2 gamma_2. These are OAI78's ray factors at L10232 (P_T = 1) and L10784-10785 | 80 coprime pairs, including the pairs with z1 = 1 or z2 = 1; G is evaluated on classes mod 72 (also checked: G(n) is a class function for 150 squarefree n) | PASS, 2.4e-15 and 2.1e-15 | zeta cubic class: err 2; wrong quotient G([z1][z2]^-1) in the second: err 1.73 |
| S4.5 | (s4-inversion), the CM marked cube inversion with conj(alpha(h))^6 and d_h(A) = d(h^3 A) | exact finite sums with W a complex bump on (1/2, 2); nu = chi_.(lambda), puncture at a prime of norm 37, label f of norm 43, row k; cases X = 900 (unmarked, two rows) and X = 3000 with a two-slot product mark {P7, P19} x {P13, P31} | PASS, relative 7.1e-16 | conj(alpha(h))^3 (the zeta power): rel err 0.44; unshifted mark d(A): rel err 0.46 |
| S4.6 | (s4-overlap), the overlap identity of the initialization, with alpha(j)^eps and slot coefficients mu(p) conj(alpha(p)^eps) a_i(p) | direct M_u Q_u (two slots, two primes each, complex coefficients, nonconstant slot weights, nu = chi_.(lambda)) against the sum over assigned subsets; eps = +1 and -1; three rows, one sharing a slot prime | PASS, relative 1.6e-15 | slot coefficients without conj(alpha(p)^eps): rel err 4.45 |
| S4.7 | (s4-sixth), the CM sixth-power column identity with alpha(d)^eps | D = 400, eps = +1 and -1, three pairs (u, a), one with u sharing a prime with a | PASS, relative 6.5e-16 | without alpha(d)^eps: rel err 0.60 |

Earlier tests that cover related identities (not rerun here): `tests_cm/proofs_2026-10-12/b4` C1-C4 (the paired
identities with twists, and the Part I first and second Poisson identities end to end), and C6 (the unmarked cube
inversion with conj(alpha)^6). S4.5 adds the puncture, the character nu and a product-form mark.

## Not tested

- The analytic statements: Lemma s4-canonical, Proposition s4-marked, and Hypothesis (R). Hypothesis (R), the CM
  reflected energy, is item 1 of Section 16 and is not proved in this insert.
- The full first canonical Poisson transformation of OAI78 with cube variables (its ray factor
  G(u1 u2^{-1}) R(u1 u2^{-1}, P_T) is described in OAI78 only "up to a fixed-ray factor"). The CM change there is
  the column twist of Lemma s4-transfer, which is exact by linearity; its u-part is S4.1/S4.4.
