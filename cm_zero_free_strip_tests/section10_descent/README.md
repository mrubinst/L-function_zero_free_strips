# Section 10 (revised 2026-10-12): numerical and symbolic checks

Target: Section 10 of `../../cm-eleven-twelfths-guide.tex`, "Remove the cubes, then reduce the number of rows", as
revised on 2026-10-12. The revision ports the presentation of Section 8 of the zeta guide: the precise transfer
estimate, the explicit induction, the retained-support argument, and exclusion removal. Test ids refer to the display
numbers (10.1)-(10.23) of the revised guide.

```sh
python3 test_sec10_descent.py      # about 5 seconds; output in outputs/test_sec10_descent.txt
```
It uses the helpers of `../part_B_lines_0550_0997` (`cmlib.py`, `common.py`, `zw.py`, `fastgauss.py`). The citation
check reads OpenAI's manuscript source at
`sources/openai/qrh-october-5-2026-paper2.tex` (a copy of OpenAI's `build/paper2.tex`).

**Last run (2026-10-12): 16 PASS (15 computational, 1 citation check), 0 FAIL. Every falsifier failed as required.**
Falsifiers exist for S10.1, S10.2-3, S10.9, S10.10, S10.11a and S10.11; the other tests are symbolic, scale-algebra or
inequality checks with no natural wrong variant.

| id | display | claim | class | test | result |
|---|---|---|---|---|---|
| S10.1 | (10.1) | cube inversion with mu_K(h) conj(alpha(h))^6 Psi(h)^3/N(h) | E | X=2000, two (k,f,xi); falsifier conj(alpha(h))^3 | PASS, 3.7e-16 |
| S10.2-3 | (10.2), (10.3) | small cubes: H + H^2 F N(h)^3/X <= 2 Sigma, X/N(h)^3 >= 1 | I | 79 506 exact rational cases; falsifier H_c^3 = min(X, X^2/H) | PASS |
| S10.4 | after (10.4) | H^2 <= X empties the supremum; nonempty supremum has N(b)^3 > X^2/H^2 | I | 464 023 exact cases | PASS |
| S10.5-7 | (10.5)-(10.7) | Y_{C,e0} -> Y_enlarged; Y_natural = L^2F^2/H; ratio Sigma/(LF) | E | sympy | PASS |
| S10.8 | (10.8) | H' = H/(N(b)^3F^2R^2); without enlargement H/(F^2R^2) (= H at F=R=1) | E | sympy | PASS |
| S10.9 | (10.9) | local cancellation | E | both cases; falsifier | PASS |
| S10.10 | (10.10) | F' <= N(f') <= N(k') <= HL/(Sigma F N(r)^2) <= H' | I | 100 000 retained triples with actual f' \| k' in Z[omega]; falsifier N(r) < R | PASS |
| S10.11a | text before (10.11) | a(vn) = a(v)a(n)chi_n(v)^4, a = conj(alpha)^2 gamma_2 xi | E | 300 coprime squarefree pairs, two xi; falsifier chi_n(v)^2 | PASS, 5.1e-15 |
| S10.11 | (10.11) | exclusion-removal identity with the CM coefficients | E | 8 cases, X'=2000, two xi, r with 1-3 primes; falsifiers: j dropped from the twist; no sign mu_K(j) | PASS, 3.0e-15 |
| S10.11b | after (10.11) | E_r <= tau(r) sum_{j\|r} E(H, X/N(j), F N(j)) | I | both sides computed from (7.4), 3 cases | PASS (ratio <= 0.31) |
| S10.12-13 | (10.12), (10.13) | the replacement preserves H', Sigma' = L/R, H'/Sigma'; X'/N(j) >= Sigma F/H >= 1 | E/I | 200 000 random retained sets | PASS |
| S10.15 | (10.15) | the three scale inequalities | I | sympy (H'/Sigma' = H/(Sigma F R)) + 100 000 samples | PASS |
| S10.16-20 | (10.16), (10.19), (10.20) | standing hypotheses, contraction, surviving hypotheses, chain length <= ceil(C0/kappa) | I | 100 000 random chains (half extremal) | PASS |
| S10.18 | (10.18) | base case: six units; abs(a_xi(n)) <= 1; counting bound | E/I | computed CM averages, 4 (X,F) | PASS |
| S10.23 | (10.23) | H/Sigma = D^-theta/B; kappa = theta/2 for large D | E/I | sympy + 20 000 samples | PASS |
| S10.cite | Section 10 | Lemma 4.4, Lemma 5.3, Proposition 5.4 (C^{4m+12}, [u/16,4v]), Proposition 5.1, Section 5.5, Section 7, Lemma 7.3 of OAI1112 | C | numbering recomputed from the LaTeX source | PASS (citation check) |
| S10.22 | (10.22) | E << D^eps Sigma | A | sanity only: E/(Sigma ‖W‖^2) <= 0.17 on a grid; control a(n) = 1 | size check |

Caveat (reviewer, 2026-10-12): the 100 000 random chains of S10.16-20 test the scale bookkeeping of the induction
under the stated transition formulas (contraction, surviving hypotheses, chain length). They do not test the
analytic transfer inequality (10.14) or cube removal (10.4), which are cited.

Exit status: the script exits 1 if any check fails or a falsifier does not fail (through `allok` and the atexit hook
of `cmlib.report`/`cmlib.falsifier`); `CM_TESTS_FORCE_FAIL=1` makes the first check fail, as a self-test.

Not tested (analytic, cited): (10.4) cube removal (OAI1112 Lemma 5.3, CM form); (10.14) the transfer estimate
(Proposition 2, CM form of OAI1112 Proposition 5.4); (10.17) and (10.21), which are the inductive assertion and its
consequence and so depend on (10.4) and (10.14).
