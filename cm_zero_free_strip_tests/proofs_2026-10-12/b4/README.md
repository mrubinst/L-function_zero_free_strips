# b4 tests: new identities and inequalities in b4.tex

The scripts are in `drafts_cm_fill/compile_b4/` and use the helpers of
`tests_cm/part_B_lines_0550_0997` (`cmlib.py`, `common.py`, `zw.py`, `fastgauss.py`).

```sh
python3 compile_b4/b4_checks.py        # C1-C8, C10, C11; about 20 s; output b4_checks_output.txt
python3 compile_b4/b4_weight_decay.py  # C9; about 15 s; output b4_weight_decay_output.txt
```

**Last run (2026-10-12): all checks PASS. Every falsifier failed as required.**

## Exact identities (class E)

| id | statement in b4.tex | test | result | existing test |
|---|---|---|---|---|
| C1 | (b2-forward), which b4 uses (and (b4-paired-fw) in the standalone file): forward paired identity with a twist xi_1. The literal sum over the 432 characters of frak G, with c_{xi'}(xi_1) from xi_1(t)G(t^{-1}). | 60 coprime squarefree pairs, xi_1 a nontrivial character. Falsifiers: zeta coefficients conj(alpha)gamma_2; xi_1 conjugated on the wrong column. | PASS, 2.2e-15. Falsifiers err 2. | B27 (`t3_gauss.py`) only for xi_1 = 1 |
| C2 | (b2-reverse), which b4 uses (and (b4-paired-rv) in the standalone file): reverse identity, output mu conj(alpha) xi_1 with l_{xi_1}(xi) from (xi G)(t). | Same 60 pairs, xi nontrivial. Falsifier: output mu xi_1 (zeta form). | PASS, 2.6e-15. Falsifier err 2. | B73 (`t3_gauss.py`) |
| C3 | Lemma b4-second-id: the second Poisson identity, end to end. | Sum over y in Z[omega] of Phi(N(y)/Y) abs(P_{C,e0,t}(y;U))^2, with the self-dual Gaussian Phi(t) = exp(-2 pi t/sqrt3) (so Phihat = Phi) and complex U, compared with Z' + the sum over g, e, w, h_2, xi'. Three configurations: (C, t, e0) = (1+3w, 4+3w, 1+3w), (1+3w, 1, 1), ((1+3w)(-5), 4+3w, -5); ell = 30-40, Y = 250-500, xi_1 = chi_n(-1)(omega/n)_3. Falsifiers: Moebius coefficient mu xi_1 (no conj(alpha)); cubic coefficient conj(alpha)gamma_2. | PASS, relative 2.4e-16. Falsifiers err 7.8e-4 (the zero-frequency diagonal dominates both sides). | none |
| C4 | (b4-first-id): the first Poisson identity of the transfer, end to end. | A(W), with the same Gaussian Phi, compared with Z + the sum over xi_1, C, t, e0, f, h, x1, x2, with p_{hf^2} containing conj(alpha). Two configurations (L, F, H) = (20, 7, 25), (25, 1, 30); xi = (omega/n)_3. Falsifier: p_y without conj(alpha). | PASS, relative 1.2e-15. Falsifier err 0.24. | none |
| C5 | Regrouping in Lemma b4-second: the sum of mu(e)mu(w) over the preimages of (r, f', k') equals tau(r) 1_{f' divides k'}. | 400 random (r, f', k') built from actual Eisenstein primes; brute force over (t, C, e, w, e0) with the conditions of Lemma b4-second-id. | PASS, 0 mismatches | S10.9 checks the local identity (10.9) only |
| C6 | (b4-long) with beta_0(b,f): the long part of the cube inversion. | X = 3000, H_c = 2, two (k, f), xi = (omega/n)_3, compared with the h-sum of the T(X/N(h)^3). Falsifier: conj(alpha)^3 in the inversion and in T. | PASS, 8.5e-19 (abs(P_long) about 5e-4). Falsifier err 1e-3. | S10.1/B68 test (10.1) only |
| C7a-d | The Hecke-family facts behind the use of [GL] (INSERT 1, "The large sieve"). (a) (x/k)_2 kappa_lambda(x)^{e_k} is trivial on the 6 units. (b) (n/k)_2(k/n)_2 depends only on k, n mod 4. (c) For coprime k1 = k2 mod 4, (u/k1k2)_2 = 1 for every unit u. (d) (-2n'/k)_2 = (-2/k)_2(n'/k)_2. | 170 squarefree k; 3000 random pairs (144 class pairs); 156 pairs; 500 pairs. | PASS, 0 violations | none |
| C10 | Scale identities in Lemmas b4-first and b4-second. (i) (1/(LF)) H/(N(e0)N(t)ell) = w_{C,e0}. (ii) X0' = L/(N(r)N(f')). (iii) w Y N(g)/(N(e)ell) = c_* Sigma N(r)/L^2. (iv) Y N(h2)N(g)^2/(N(e)ell^2) = c_* Sigma F N(k')N(r)^2/(HL). (v) w Y ell = c_* Sigma/(N(C)^2 N(t)). (vi) (Sigma R/L^2) R (L/R)^2 = Sigma. | sympy | PASS | S10.5-8 test related scale relations |

## Explicit inequalities (class I)

| id | statement | test | result |
|---|---|---|---|
| C8 | Lemma b4-sqfree, local costs: cost_pi <= N(pi)^2 if pi divides s_0, and cost_pi <= N(pi)^{v_pi(g)} otherwise. These give (b4-aY), a^2 Y << H^2 N(g)/X. | Enumeration of v_pi(s_0) in {0,1} and v_pi(g) <= 12. | PASS |
| C11 | Lemma b4-first: the map (f, h) to h f^2 has multiplicity <= #{f squarefree : f^2 divides y} <= tau(y). | All squarefree f with N(f) <= 200 and h in a box. | PASS (max multiplicity 3) |
| (dyadic) | Lemma b4-sqfree: the sum over dyadic U, B of (1 + UB^3/Z)^{-2} is << log^2(2+Z), and the same sum weighted by sqrt U is << sqrt Z, for every Z > 0. | Direct summation for Z = 10^-3 to 10^6 (inline in the session; trivial to add). | ratios <= 0.4 and <= 3.5 |
| (cube) | Lemma b4-cube: H + 2H^2 F H_c^3/X << Sigma and X/N(h)^3 >= 1 for N(h) <= H_c. | Already tested | S10.2-3 / B70 PASS |
| (child) | Lemma b4-second: X'/N(j) >= Sigma F/H >= 1, H' >= 1, and 1 <= X'/X0' < 4. | Already tested | S10.10, S10.12-13 PASS |

## Size checks (class A: bounds with unspecified constants)

| id | statement | test | result |
|---|---|---|---|
| C9 | Lemma b3-decay, which b4 uses (and Lemma b4-weight in the standalone file): abs(V^sharp(x)) << min(x^{1/4}, x^{-A}) for the CM quotient Gamma(4/3+t)Gamma(5/3+t)/(Gamma(4/3-t)Gamma(5/3-t)). | mpmath. Moving the contour from Re t = 0 to -1/4, -1 and 3 leaves V^sharp unchanged (diff <= 1e-14), which confirms holomorphy in Re t > -4/3. abs(V^sharp)/x^{1/4} stays bounded as x -> 0, and x^3 abs(V^sharp) -> 0 as x grows. The test weight is log-normal, V(y) = exp(-(log y)^2), with Vhat(s) = sqrt(pi) e^{s^2/4}: entire, with Gaussian decay on vertical lines, which is all the proof uses. A compactly supported bump makes the vertical integrals ill-conditioned at Re t = 3. Control: the zeta quotient gives a different function. | PASS (size check only) |

## Not tested (analytic statements with unspecified constants)

These are proved in b4.tex but cannot be tested numerically in any meaningful way:
- Lemma b4-sqfree and Proposition b4-completed (= (9.2));
- Lemma b4-cube (= (10.4));
- Lemmas b4-first and b4-second as inequalities (their exact identities are C3-C5);
- Proposition 2, the transfer estimate.

The large sieve (9.1) is cited, as before.

## Suggested additions for the numerics section of the guide

C1, C3, C4 and C6 are the CM-specific checks. Each has a zeta-form falsifier that fails, so they fit the
convention of Section 18 ("the variant is the corresponding zeta normalization ... and the variant was required
to fail"). Suggested bullet:

> In Sections 9.1, 10.5 and 10.6, the following hold:
> - the paired identities with twists of Lemma b2-paired, summed over all 432 characters;
> - the two Poisson identities of the transfer, Lemmas b4-first and b4-second-id, end to end, with the self-dual Gaussian row
>   weight, to relative error 2.4e-16;
> - the long-cube formula (b4-long) with conj(alpha)^6.
>
> The zeta forms fail in each case.

(Labels refer to b4.tex. In the combined build `compile_b4/cm-eleven-twelfths-guide.pdf` they are (10.24), Lemma 1 (b2), Lemmas 4 and 5.)
