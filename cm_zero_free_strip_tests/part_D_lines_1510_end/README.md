# Part D audit: lines 1510 to 1789 of `cm-eleven-twelfths-guide.tex`

This folder audits the section "CM newforms of higher weight: why the method should extend" (lines 1510 to 1666) and
the appendix "The reflected series and explicit theta coefficients" (lines 1668 to 1754), plus the bibliography. The
paper was not modified.

## Running the tests

```
python3 run_all.py          # all 57 test lines; about 6 minutes (t08 alone takes about 3.5 minutes)
python3 run_all.py --fast   # skips t08
```
Exit status (2026-10-12): `pd_common.report` records failures and an atexit hook makes each script exit 1 if a check
failed; `run_all.py` exits 1 if any script exits nonzero. `CM_TESTS_FORCE_FAIL=1` makes the first check of each
script fail, as a self-test.

Each script prints one `[PASS]` or `[FAIL]` line per claim. Each line includes a **falsifier**: a deliberately wrong
variant that must fail. A test passes only if the claim holds **and** its falsifier fails. Outputs are saved in
`outputs/`.

**Last run:** 57 PASS, 0 FAIL.

| file | content |
|---|---|
| `pd_common.py` | reporting, PARI/GP wrapper, Eisenstein helpers |
| `zw.py`, `fastgauss.py` | unmodified copies of `tests_guide/` helpers |
| `appendix_theta.py` | `tau`, `tau_1`, `tau_2`, `t_sigma`, `d_sigma` implemented literally from the appendix tables (brute-force `g(r,c)`) |
| `t01_eta_quotients.py` | the eta quotients, the Part I rows, and the functional equations (PARI/GP plus an independent q-product) |
| `t02_chain_rule.py` | the j-th derivative chain rules, symbolically (sympy) |
| `t03_mellin.py` | the direct and reflected Mellin transforms for general j |
| `t04_kernels.py` | `R_j`, its poles and growth, `V#`, `R_j(1)`, the Hecke gamma quotient, and the Hankel kernel |
| `t05_misc.py` | the weights, the Poisson transform, the co-oriented classes, the angular integral, kappa, the vertical shift, the term count, `(b^3/k)_2`, and the cusp-scalar phases |
| `t06_sympowers.py` | the factorization of the symmetric powers (local Euler polynomials, PARI `ellap`) |
| `t07_theta_appendix.py` | the theta coefficients: automorphy, the expansions at the other two cusps, support and bound, the comparison with DR, and convergence |
| `t08_reflection_j.py` | an end-to-end single-cusp theta reflection with j = 1, 2, 3 derivatives (B = 80000, 41 cusp configurations) |

## Claim, class, test and result

Classes:

- **E**: an exact identity that can be checked finitely.
- **I**: an explicit inequality between numbers or exponents.
- **A**: an asymptotic or analytic bound, or an expectation.
- **C**: a cited result.

"(partial)" means that only a finite or representative part of the claim was tested.

| # | line | label | claim | class | test | range | result |
|---|---|---|---|---|---|---|---|
| 1 | 1512-1513 | | psi_k = conj(alpha)^k eta are the unitary Hecke characters with infinity type (zbar/\|z\|)^k | C | none (convention; consistent with L1033 and L1057) | | untested |
| 2 | 1513-1516 | | L(f, s+k/2) = L(s, psi_k or its conjugate) for CM newforms of weight k+1 | C | D01.7, D01.12 (spot checks for k = 3, 2) | | PASS (spot) |
| 3 | 1518-1522 | Expected statement | L(s, conj(alpha)^k eta) != 0 for Re s > 7/8, equivalently L(f,s) != 0 for Re s > k/2+7/8 | A | none (a conjecture, not claimed); the arithmetic of the equivalence is in D05.1 | | untested |
| 4 | 1523 | | critical strip k/2 < Re s < k/2+1 | I | D05.1 | k = 1..29 | PASS |
| 5 | 1524 | | 11/12 analogue k/2 + 11/12 | I | D05.1 | | PASS |
| 6 | 1527-1529 | | same length exponents for each j; constants depend on j | A | none | | untested |
| 7 | 1532 | | q-expansions of the two eta products to q^3000 | E | D01.1 (PARI = Euler/Jacobi product) | n <= 3000 | PASS |
| 8 | 1537 | | eta(3z)^8 is a newform of weight 4 and level 9 with CM by K (k = 3) | E | D01.2 | | PASS |
| 9 | 1537-1538 | | 207 split primes below 3000, 414 generator checks | E | D01.4 | | PASS |
| 10 | 1539 | | a_p = 2 Re(pi^3) for both primary generators | E | D01.5 | 207 primes p | PASS |
| 11 | 1539 | | p = 2 mod 3: a_p = 0 and a_{p^2} = -p^3 | E | D01.6 | a_p for p <= 40000; a_{p^2} for p <= 197 | PASS |
| 12 | 1540 | | L(f, s+3/2) = L(s, alpha^3), as the primitive character of conductor lambda | E | D01.7: a_m = sum over primary g with N(g) = m of g^3, exact | m <= 40000 | PASS |
| 13 | 1540-1541 | | the target is conj(alpha)^3 and eta is trivial | E | D01.7 (real coefficients) | | PASS |
| 14 | 1541-1542 | | expected strip 19/8 < Re s < 5/2, centre 2 | I | D05.1 | | PASS |
| 15 | 1543 | | eta(2z)^3 eta(6z)^3 is a newform of weight 3 and level 12 (k = 2) | E | D01.3 (also nebentypus chi_{-3}, CM by -3) | | PASS |
| 16 | 1543-1545 | | a_p = 2 Re(pi^2 (2/pi)_3) = 2 Re(pi^2 chi_pi(4)), 414 checks | E | D01.8 | 207 primes p | PASS |
| 17 | 1545-1546 | | for p not dividing 12 with p = 2 mod 3: a_p = 0 and a_{p^2} = p^2 (nebentypus chi_{-3}) | E | D01.9, D01.3 | p <= 197 | PASS |
| 18 | 1546-1547 | | a_4 = 0, not 4 | E | D01.10 | | PASS |
| 19 | 1547 | | L(f, s+1) = L(s, alpha^2 chi_.(4)), primitive: conductor 2, value -1 at lambda | E | D01.12, exact identity | m <= 40000 | PASS |
| 20 | 1548 | | target conj(alpha)^2 conj(chi_.(4)); strip 15/8 < Re s < 2, centre 3/2 | I | D05.1 | | PASS |
| 21 | 1549 | | a_3 = -3, so the factor at 3 is (1+3^{-s})^{-1}; the factor at 2 is 1 | E | D01.11 (a_{3^k} = (-3)^k; a_{2^k} = 0) | 3^k, 2^k <= 40000 | PASS |
| 22 | 1551-1553 | display | L(f, s+1) = (1+3^{-s})^{-1} L^(6)(s, alpha^2 chi_.(4)) | E | D01.13 (chi_pi(4) by Euler's criterion) | n <= 3000 | PASS |
| 23 | 1557-1560 | display | gamma_2(c) conj(alpha(cn^3)) psi(cn^3) = gamma_2(c) conj(alpha(cn^3))^{k+1} eta(cn^3) | E | D05.2 | 60 random (c, n, k) | PASS |
| 24 | 1561-1562 | prop:cm-rows | the row formula for general k ("formal in psi(p)") | C | none (prop:cm-rows lies in another chunk) | | untested |
| 25 | 1562-1563 | | numerators have finite order; denominators have infinity type (zbar/\|z\|)^k; slot amplitudes are unchanged | A | none (structural) | | untested |
| 26 | 1568-1569 | | at z = 0, z' is affine in zbar with slope -(cv)^{-2}, and zbar', v' do not depend on zbar | E | D02.4 | | PASS |
| 27 | 1570-1573 | display | d_zbar^j F(z', zbar', v') = (-(cv)^{-2})^j (d_{z'}^j F)(-delta/c, 1/(N(c)v)) | E | D02.1 (symbolic, complete monomial basis; a proof for each j) | j <= 6 | PASS |
| 28 | 1574 | | every derivative of order j >= 1 kills the constant term | E | D02.3 | j <= 6 | PASS |
| 29 | 1576, 1581-1582 | | absolute convergence exactly for Re s > 1/6 - j/2; continuation elsewhere | A | D03.2 (endpoint exponent, numerical divergence and convergence) | j = 1, 2, 3 | PASS (partial) |
| 30 | 1577-1581 | display | the direct Mellin transform with j derivatives | E | D03.1 (40-digit quadrature) | j <= 6, three (l, s) each | PASS |
| 31 | 1583 | | j = 1, 2 give the two table entries | E | D03.3 | | PASS |
| 32 | 1583-1584 | | "checked numerically for j <= 5, error below 1e-25" | E | D03.1 | | PASS (2e-40) |
| 33 | 1585-1586 | | reflected side: prod Gamma((j+2)/2 - s +- 1/6), cusp factor conj(alpha(c))^{2j}, dual factor alpha(mu)^j | E | D03.4 (closed form, including constants); D08.1-3 (end to end) | j <= 4; j <= 3 | PASS |
| 34 | 1587-1591 | display | R_j(t); R_1 and R_2 are the table kernels | E | D04.1 | | PASS |
| 35 | 1591-1592 | | R_j is holomorphic for Re t > 1/6 - (j+1)/2, which lies left of -1/4 | I | D04.2 (exact pole list; genuine simple pole) | j <= 40 | PASS |
| 36 | 1593 | | \|R_j(sigma+iu)\| << (1+\|u\|)^{4 sigma} on -1/4 <= sigma <= A | A | D04.3 (finite sup; Stirling limit 1) | j <= 6, A = 3, u <= 1e7 | PASS (partial) |
| 37 | 1593-1594 | | V# << min(x^{1/4}, x^{-A}) | A | D04.4 (Gaussian V, A = 2) | x in [1e-8, 1e8], j <= 3 | PASS (partial) |
| 38 | 1595-1598 | display | R_j(1) = (3j-4)(3j-2)(3j+2)(3j+4)/1296 ~ j^4/16 | E | D04.5 (exact, sympy) | j <= 60 | PASS |
| 39 | 1600 | | dual lengths are unchanged: the argument N(mu)X/N(c)^2 does not depend on j | E | D03.5 (symbolic plus numeric); D08 | | PASS |
| 40 | 1601 | | alpha(mu)^j is a unimodular coefficient | E | D03.4 | | PASS |
| 41 | 1602-1603 | | Prop cm-completed and Z^{1/4+e}, Z^{3/16+e} should hold for every k | A | none (an expectation) | | untested |
| 42 | 1605 | | the functional equation has Gamma(s+k/2) | C | D01.15 (PARI lfuncheckfeq for both examples) | | PASS (spot) |
| 43 | 1606-1607 | | the exponent of \|Gamma(1-s+k/2)/Gamma(s+k/2)\| is 1 - 2 sigma, independent of k; cm-growth keeps its exponents | A | D04.6 | k <= 8, six sigma | PASS (partial) |
| 44 | 1607-1608 | | the Poisson transform of abar^k e^{-cv\|a\|^2} is a multiple of xibar^k e^{-c'\|xi\|^2/v} v^{-1-k} | E | D05.3 (constant i^{-k}) | k <= 5 | PASS |
| 45 | 1609 | | it vanishes at xi = 0, so L(s, psi) is entire | A | D05.3 (vanishing only) | | PASS (partial) |
| 46 | 1611-1613 | | co-oriented classes mu conj(alpha)^k <-> conj(alpha)^{k+1} gamma_2; cube factor conj(alpha)^{3(k+1)} | E | D05.4 | 40 squarefree n, k <= 4 | PASS |
| 47 | 1613-1614 | | the theta steps use d_zbar^{k+1} theta | E | D08 | | PASS |
| 48 | 1614-1615 | | Gamma(s+k/2)/Gamma(1-s+k/2) is a Hankel kernel of order k with the same vertical growth | E | D04.7 (Mellin transform of J_k), D04.6 | k <= 4 | PASS |
| 49 | 1615-1616 | | int alpha(z)^{+-k} f(\|z\|^2) dz = 0 | E | D05.5 (grid refinement) | k <= 6 | PASS |
| 50 | 1618-1621, 1623 | | label sets, slot normalization and exponents carry over verbatim | A | none | | untested |
| 51 | 1622 | | kappa = 3/4 | I | D05.6 | | PASS |
| 52 | 1625-1626 | | Part I rows: u_0 = 1 for eta(3z)^8 and u_0 = 4 for eta(2z)^3 eta(6z)^3 | E | D01.14a, b (Dirichlet inverse) | n <= 1500 | PASS |
| 53 | 1630 | | k = 0: j = 1 and the weight is gamma_2 conj(alpha) eta | E | D05.2 | | PASS |
| 54 | 1631-1632 | | the recipe covers all unitary characters with purely angular infinity type; negative k reduces by conjugation | C | none (standard) | | untested |
| 55 | 1633-1637 | display | L(s, psi N^{it}) = L(s - it, psi) | E | D05.7 | | PASS |
| 56 | 1638-1640 | | the conjugate target leaves alpha^{k-1}; for k = 1 the constant term is not killed | E | D05.2; D02.3 (j = 0 falsifier) | | PASS |
| 57 | 1640-1641 | | alpha^{k-1} is produced by d_z^{k-1} | E | D02.5 | n <= 6 | PASS |
| 58 | 1641-1642 | | d_z^j chain rule (-(cbar v)^{-2})^j d_{zbar'}^j, "checked for j <= 3" | E | D02.2 | j <= 6 | PASS |
| 59 | 1643 | | for k >= 2 either orientation should work | A | none | | untested |
| 60 | 1644-1645 | | L_an(Sym^m E_d, s) = L(Sym^m E_d, s + m/2) | E | D06.3 | | PASS |
| 61 | 1646-1654 | display | Sym^m factorization into L(s, (psi_d^*)^{m-2i}), times L(s, chi_{-3}^{m/2}) for even m | E | D06.1 (local Euler polynomials, PARI a_p) | p <= 600, m <= 8, 7 values of d (5960 cases) | PASS |
| 62 | 1655-1656 | | last factor zeta if m = 0 mod 4, L(s, chi_{-3}) if m = 2 mod 4 | E | D06.2, D06.1 | | PASS |
| 63 | 1657 | | this is the CM factorization of [MW, Sec. 3.5] | C | MW source (antsVII.tex, Sec. 3.5, eq. (cmfact)) read: identical after s -> s + m/2 | | PASS (inspection) |
| 64 | 1657-1659 | | the remaining factor is covered by OpenAI's theorem; zero-free for Re s > 7/8 | A | none (conditional on the expectation) | | untested |
| 65 | 1659-1661 | | pole at s = 1 when 4 \| m; not primitive for m >= 2 | E | D06.1, D06.2 | | PASS (nonvanishing of other factors at s = 1 not tested) |
| 66 | 1662-1665 | | six restatements are needed for each k | A | none | | untested |
| 67 | 1673-1677 | | theta inputs come from [DR, Sec. 5.1 and Appendix]; DR's GRH-conditional asymptotic is not an input | C | DR source read: Sec. 5.1 is the cusp section, Appendix rows 1/10/19; DR's main theorem assumes GRH | | PASS (inspection) |
| 68 | 1679-1686 | eq:reflection | T(X; Psi) = finite sum of C sum d_s(l) alpha(l)^2 N(l)^{-1/2} psi Prod B V#(N(l)X/N(c)^2) | E | D08.2 (single cusp, j = 2, kernel R_2, appendix coefficients; 41 cusps; error 2e-8) | | PASS (partial) |
| 69 | 1687-1691 | | sigma in {0,+,-}; A; c = c_0 prod pi; c_0 \| 8 lambda^2; psi(x) = e(beta x/(lambda^3 c_0)) | C | D05.8 (c_0 \| 24 for L = 72) | | consistent |
| 70 | 1692 | | at most 5184 * 2^{\|P\|} terms | C | D05.8 (5184 = N(72) = number of h_0 mod L in OAI78's count) | | consistent (L = 72 is implicit) |
| 71 | 1692-1693 | | \|C\| <= 1/81 | C | none (OAI78's \|zeta\| <= 1/81) | | untested |
| 72 | 1693-1695 | | constants as in OAI78's completed reflection with two derivatives | C | D08.2; D05.10 (phases -i, -1 of the table's cusp scalars, L1226) | | PASS (partial) |
| 73 | 1696 | | absolutely convergent series | A | D07.8 (dyadic tails of the majorant decay like 2^{-1.5}) | N(9l) <= 3000 | PASS (partial) |
| 74 | 1697-1698 | | V# and B_{pi,j} are defined in eq:weight and eq:local | E | inspection (lines 771, 788) | | PASS |
| 75 | 1700-1707 | | g(r,c); (x/(-2))_3 = x mod 2; \|g(r,c)\| = \|c\| (used for the bound) | E | D07.5 | squarefree c with N(c) <= 300, 4 twists | PASS |
| 76 | 1708-1722 | table | tau(l) | E | D07.1 (automorphy under T and E), D07.6, D07.7 (same as DR (5.7)) | | PASS |
| 77 | 1723-1735 | table | tau_1 and tau_2 | E | D07.2 (direct computation of theta at gamma_10, gamma_19), D07.7 (DR (5.13)-(5.14)) | | PASS |
| 78 | 1737-1742 | display | t_0, t_+, t_- and d_sigma | E | D07.2, D07.3, D08 | | PASS |
| 79 | 1743-1744 | | these are precisely the coefficients of the conjugated theta functions at the three cusps | E | D07.3 | | PASS |
| 80 | 1744-1747 | | support is in u lambda^m n b^3 with m >= -4, n squarefree, n and b primary | E | D07.4 | N(9l) <= 3000, 7532 coefficients | PASS |
| 81 | 1748 | display | \|d_s(u lambda^m n b^3)\| <= 27 * 3^{m/6} \|b\| | I | D07.4 (sharp: equality for row 4) | | PASS |
| 82 | 1750-1751 | | (b^3/k)_2 = (b/k)_2 | E | D05.9 | 400 random pairs | PASS |
| 83 | 1752-1754 | | fixing u, m, b leaves squarefree columns n | A | D07.4 (n squarefree) | | PASS (structural) |
| 84 | 1757-1788 | bibliography | bibliographic data | C | MW Sec. 3.5 and DR Sec. 5.1 / Appendix / GRH confirmed in sources; volumes and pages not checked | | untested |
| 85 | 1585 (and 1226) | | the general-j cusp factor agrees with the table's -(i/81) conj(alpha(c))^2 and -(1/81) conj(alpha(c))^4 | E | D05.10 | | PASS |

### Counts

| class | claims | PASS | FAIL | untested |
|---|---|---|---|---|
| E | 51 | 51 | 0 | 0 |
| I | 7 | 7 | 0 | 0 |
| A | 15 | 7 (partial) | 0 | 8 |
| C | 12 | 7 (spot, inspection or consistency) | 0 | 5 |
| **total** | **85** | **72** | **0** | **13** |

## Discrepancies

**No genuine errors were found in lines 1510 to 1789.** Every exact identity, count and numerical value tested
correctly, and every falsifier failed.

### Convention ambiguities and notational remarks

None of these is an error.

1. **Infinity type (L1513).** The phrase "infinity type (zbar/|z|)^k" uses the convention psi((a)) = psi_inf(a) for
   a = 1 mod f. Under the other common convention, psi((a)) = psi_inf(a)^{-1}, the same character has infinity type
   (z/|z|)^k. The paper uses its convention consistently (L1033, L1057).
2. **Reuse of psi (L1691).** In eq. (reflection), psi(x) = e(beta x/(lambda^3 c_0)) is an additive character. Elsewhere
   psi is the Hecke target (L1557). The meaning is clear from context.
3. **The count 5184 (L1692).** The bound 5184 * 2^{|P|} equals N(L) * 2^{|P|} with L = 72. That is OAI78's count of terms
   (classes h_0 mod L, times subsets of active primes). The paper never names L, and the reader must accept that L = 72
   satisfies OAI78's periodicity requirement for the CM twist Psi_0. The statement c_0 | 8 lambda^2 is consistent with
   this, since lambda^2 h_0/72 = -h_0/24.
4. **Evaluation point (L1572).** The point (-delta/c, 1/(N(c)v)) omits the zbar' coordinate (-deltabar/cbar), as eq.
   (chain) does. The tests use the full point.
5. **Cusp scalars (L1585 and L1226).** The paper states only the factor conj(alpha(c))^{2j}. The full j-dependence of the
   cusp scalar is (-1)^j conj(alpha(lambda^-3))^{-j} conj(alpha(c))^{2j} = (-i)^j conj(alpha(c))^{2j}. This gives the
   table's -i/81 for j = 1 and -1/81 for j = 2, and D05.10 confirms it.

## Not tested, and why

- **Expectations, by design not theorems (rows 3, 6, 25, 41, 50, 59, 64, 66).** These are statements that estimates and
  lemmas "should hold" or "carry over" for every k. Testing them would require the analytic arguments themselves.
- **Rows 1 and 54.** These are definitions and conventions, or standard facts about Hecke characters.
- **Row 24 (prop:cm-rows for general k).** The proposition belongs to another chunk. The extension to general k is
  formal in psi(p).
- **Row 71 (|C| <= 1/81).** This is OAI78's bound |zeta| <= 1/81. The phases are checked (D05.10), but the full scalar
  requires OAI78's Fourier expansion of Psi.
- **Row 84 (bibliographic volumes and pages).** These are not mathematical claims, and the arXiv abstract pages carry no
  journal references to compare against.
- **Partial items.**
  - Row 68: eq. (reflection) was tested end to end at a single cusp, with no moving primes. The finite additive
    expansion of Psi and the local factors B_{pi,j} (eq. (local), another chunk) are not exercised.
  - Rows 36, 37, 43, 73: the analytic bounds were tested on finite ranges and with one test function.
