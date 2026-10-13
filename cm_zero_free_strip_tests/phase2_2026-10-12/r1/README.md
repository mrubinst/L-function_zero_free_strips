# r1 tests: abstract, introduction and Sections 1-7 of the CM guide

Guide: `sources/cm_zero_free_strip.tex` (version of 2026-10-12, 3499 lines). Equation numbers are
those of the compiled guide (`sources/cm_zero_free_strip.aux`).

All three scripts use `r1lib.py`, a new implementation of the arithmetic of Z[omega] written for this review. It does
not import the helpers of `tests_cm/`, so its results are an independent check of the conventions: primary = 1 mod 3,
chi_n(x) = (x/n)_6 by Euler's criterion with zero values kept, alpha(n) = n/|n|, e(z) = exp(4 pi i Im z/sqrt 3),
Gauss sums by brute force over a complete residue system (Hermite form of n Z[omega]), primary generators.

```sh
./run_all.sh                         # all three, about 10 s; or individually:
python3 -I test_classgroup.py        # about 3 s
python3 -I test_extraction_rows.py   # about 1 s
python3 -I test_step45_bookkeeping.py   # about 5 s
```
Outputs are in `outputs/`. **Last run (2026-10-12): 30 checks (18 + 5 + 7, falsifier rows included), all PASS; every falsifier fails as required.**

## What is new, and what was already covered

Most identities of Sections 2-7 are already tested in `tests_cm/part_A_lines_0001_0549` (A..), `tests_cm/part_B_lines_0550_0997`
(B..) and `tests_cm/proofs_2026-10-12/b1, b2` (b1.., b2-T..). I re-read those tests against the current text and reuse
them where they cover the claim (see the r1 report, section "What was verified", returned as text to the main agent). The scripts here cover
claims that were untested, only implicitly tested, or listed as "suggested, not run".

| id | guide location | claim | test (falsifier) | result |
|---|---|---|---|---|
| C1 | Sec. 7.1, Step 2 of Prop. 3; edit E7 | (O/72)^x has 2592 elements; frak G = {x = 1 mod 3} has 432, is a group, and maps onto the 12 classes of (O/4)^x with fibres of 36 | exact enumeration | PASS |
| C1b | Sec. 7.1 | every class of frak G contains a squarefree primary n coprime to 6 (N(n) <= 9000) | enumeration | PASS (432/432) |
| C2 | Sec. 7.1, 7.3 ("G depends only on n mod 4", OAI A.1) | \|G\| = 1, G(n) depends only on n mod 4, and equals OAI's closed form conj(chi_n(4)) (1+i^-b+i^a+i^(b-a))/2 | 367 squarefree n, N <= 1300 (fals.: n mod 2 does not determine G; dropping conj(chi_n(4)) changes G) | PASS, err 2e-14 |
| C3 | Sec. 7.3 "The finite group"; edit E9 | R(a,b) = chi_b(a)/chi_a(b) is +-1, depends only on the classes mod 4, symmetric | 5861 coprime pairs (primes and composites), all 144 class pairs (fals.: classes mod 2 do not determine R) | PASS |
| C3b | same, OAI (4.1) | on classes mod 4: R is a bicharacter, R(c,c) = chi_c(-1), G(xy) = G(x)G(y)R(x,y) | exact on the 12 classes | PASS |
| C4 | Sec. 7.3 | chi_x(-1) = (-1)^((N(x)-1)/6), a function of x mod 4 | all primary x, N <= 3000 (fals.: chi_x(-1) = 1, as for the cubic symbol) | PASS |
| C5 | (7.7) | chi_x(-1) conj(G(x)) G(y) R(x,y) = G(y x^-1) on frak G | all 144 class pairs mod 4 (fals.: drop chi_x(-1); G(x y^-1)) | PASS, err 4e-15 |
| C6 | (7.3), (7.9), (7.10) | forward and reverse paired identities, class form, also with xi_1 = chi_n(omega) and with a character mod 4 | 70 random coprime squarefree pairs, N(z1 z2) <= 2600, brute-force Gauss sums mod z1 z2 (fals.: zeta coefficient conj(alpha) gamma_2; G(t) for G(t^-1); reverse identity without conj(alpha(z1)) alpha(z2)) | PASS, err 6e-14 |
| C7 | (7.3) literally | with the 12 characters of (O/4)^x: \|c_xi\| <= 1, sum c_xi = G(1) = 1, and sum_xi c_xi a_xi(z1) conj(a_xi(z2)) equals the left side | 40 pairs | PASS, err 6e-14 |
| E1 | Sec. 5, report item R1-6 (optional edits E10-E13) | for pi \| d (pi coprime to 6), chi_n(u_d pi^6) = chi_n(u_d) for every column, so A_{u_d pi^6}(D) = A_{u_d}(D) exactly | d = 7*13*19, all 6 primes above d, 1696 columns (fals.: pi = -5 does not divide d: the rows differ by 0.033) | PASS |
| E2 | Lemma 2(iii) with its explicit constant (b1.18 "not run" concerned (iv); (iii) is the substantive part) | \|A_u0 - A_u0pi^6\| <= (1+16b)\|\|W\|\| D/Y for pi in P_Y, and #{n : pi \| n, N(n) <= bD} <= 1 + 16bD/Y | 216 cases: d in {1, 2, -3, 7}, D in {1500, 6000}, Y in {30, 120, 600}; the difference is also recomputed directly from both rows (sensitivity: the count is not <= bD/(10Y)) | PASS (max ratios 0.022 and 0.042) |
| E3 | Lemma 2(ii) | rows u_0 pi^6 distinct, N = N_0 N(pi)^6 <= N_0 Y^6 | same cases | PASS |
| S1 | Prop. 3, Step 4 | N(r0)N(f) <= b0 D and H N(k) <= C_Phi N(f)^2 N(m1)N(m2) <= C_Phi b0^2 D^2/N(r0)^2 on the support | 3777 random support points from actual squarefree elements, D = 6000, H = D^1.1 | PASS |
| S2 | Prop. 3, Steps 4-5 (b2-T15, "suggested, not run") | N(k) <= calH = C1 D^2/(H B^2); rho in [1,4); 0 < kappa <= C1; x_i in [a0/4, b0]; W_0(N(r0 f m_i)/D) = W_0(rho x_i); hat Phi argument = kappa/(x1 x2) | same (fals.: rho in [1,2); support [a0/2, b0]; N(f) for N(f)^2 in kappa) | PASS |
| S3 | Prop. 3, Step 5 | (X/N(j), F N(j)) satisfy (7.15) with the same calH | every prime j \| r0 in the sample | PASS |
| S4 | Prop. 3, Step 5 | sum_{B <= N(r0) < 2B} (H N(r0)/D^2)(D/B)^2 << H uniformly in B | exact sums, all dyadic B <= 2D | PASS (max ratio 1.000, attained at B = 1) |

Not testable by computation, and not tested here: (MS) itself; the analytic inequalities of Lemma B.2 of OAI1112; the
existence of the constants c_1, Y_1 (prime number theorem); the cited theorems (Deuring, Hecke, Hasse, Kubota-Patterson,
Goldmakher-Louvel).
