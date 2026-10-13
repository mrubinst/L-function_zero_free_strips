# Displays (10.1)-(10.41) of the CM guide (agent t3)

Scope: Section 10 of `../../../cm-eleven-twelfths-guide.tex`, "Remove the cubes, then reduce the number of rows":
cube inversion, small cubes, cube removal, the two-Poisson transfer, exclusion removal, the transfer proposition,
the induction, and the proofs of cube removal (Section 10.5) and of the transfer estimate (Section 10.6). The
41 displays are listed in `sources/inventory/inventory_t3.json`. One row per display is in `rows.json`.

```sh
python3 run_t3.py              # everything; about 2 minutes wall clock (parallel), 4 minutes CPU
python3 run_t3.py --no-table   # skips t2-t5 of part B (the constituents of the table (10.27)); about 20 s
python3 t3_sec10_new.py        # only the new checks; about 6 s
```

`run_t3.py` re-runs these scripts read-only, with `PYTHONDONTWRITEBYTECODE=1`, so nothing is written into the other
suites. Outputs go to `outputs/`, with a summary in `outputs/run_t3_summary.txt`.

- `section10_descent/test_sec10_descent.py` (S10.*)
- `proofs_2026-10-12/b4/b4_checks.py` (C1-C11)
- `part_B_lines_0550_0997/t7_descent.py` (B68-B87)
- `part_B_lines_0550_0997/t2_poisson.py`, `t3_gauss.py`, `t4_theta.py`, `t5_reflection_j2.py` (B17, B27, B46, B50,
  B53, B55; with B68 these are the rows of the table (10.27), assembled as B88 in part B)
- the new `t3_sec10_new.py` (D10.*, S10.32, S10.34)

The runner exits with status 1 in any of these cases: a script exits nonzero; a script prints a failure line (the
same patterns as `part_B/run_all.py`); or a test id that `rows.json` relies on is not reported as PASS. The new
script uses `cmlib.report`/`cmlib.falsifier`, so a failed check, or a falsifier that does not fail, gives exit
status 1. I confirmed this with `CM_TESTS_FORCE_FAIL=1`, which gave exit status 1.

**Last run (2026-10-13):**
- All 8 scripts passed. All 61 required test ids passed. Every falsifier (69 in total) failed as required.
- Wall clock: 127 s. CPU time: 215 s user.
- Per script: t5 126 s, t2 76 s, b4 19 s, t4 13 s, the others 2-6 s.

## Counts by kind (41 displays)

| kind | count | displays |
|---|---|---|
| E (exact identity) | 11 | 10.1, 10.9, 10.12, 10.14, 10.15, 10.27, 10.32, 10.33, 10.36, 10.38, 10.39 |
| I (explicit inequality) | 9 | 10.3, 10.13, 10.16, 10.18, 10.21, 10.22, 10.23, 10.26, 10.35 |
| A (analytic bound) | 8 | 10.4, 10.17, 10.20, 10.24, 10.25, 10.34, 10.40, 10.41 |
| D (definition) | 10 | 10.2, 10.5, 10.6, 10.7, 10.8, 10.10, 10.11, 10.19, 10.30, 10.37 |
| C (cited) | 0 | - |
| R (restatement) | 3 | 10.28 = (10.4), 10.29 = (10.5), 10.31 = (10.7) |

**E and I:** all 20 were tested in this run and all passed.
- 13 are covered by existing tests, re-run now: (10.1), (10.3), (10.9), (10.12), (10.13), (10.14), (10.15),
  (10.18), (10.21), (10.22), (10.23), (10.26), (10.27).
- 3 are end-to-end identities from the b4 suite, re-run now: (10.32), (10.33), (10.36).
- 4 had no direct test, so new tests were written: (10.16), (10.35), (10.38), (10.39).
- New partial checks were also added for (10.21) (sympy) and (10.27) (the zeta column against the sources).

**A:** none can be proved by computation. The guide proves them in these places:
- (10.4) in Lemma 6 (Section 10.5);
- (10.17), (10.40) and (10.41) in Section 10.6 (Lemmas 7-9);
- (10.20), (10.24) and (10.25) by the induction of Section 10.4.

Finite checks for the A displays:
- (10.4): the side claims (empty supremum, exact long-cube identity C6).
- (10.17): the two Poisson identities C3 and C4. Also J_t = 4m+12 and [u/16, 4v] match OAI1112 Proposition 5.4 in
  the source.
- (10.20) and (10.25): scale bookkeeping on 100 000 random chains. The sanity check E/(Sigma ||W||^2) <= 0.17 on
  computed CM averages.
- (10.34): a size check. sup over c of ||K_c||_{C^q}/||W||_{C^q}^2 is at most 0.25 for q <= 2.
- (10.40) and (10.41): the derivative orders C^{2j+4}, C^{2m+4} and C^{4m+12} match OAI1112 Lemma 7.1,
  Lemma 7.3 and Proposition 5.4 in the source. 2(2m+4)+4 = 4m+12 (sympy). E <= A(W) on computed averages.

**D:** each definition was checked against the code that implements it and against OAI1112 where the guide follows
it:
- P, Q, A(W), X_0', a^sharp and K are implemented verbatim in C3 and C4, which pass end to end.
- H_c and L_b in S10.2-3, S10.4 and C6.
- r = tg/e, f' = Cew, k' = e0 e h2: r and f' are squarefree and coprime (D10.10). The definitions are word for word
  those of OAI1112, proof of Lemma 7.3.
- c_* = max(1, 4 C_Phi v^2) gives N(hf^2) <= Y_{C,e0} (D10.30). The falsifier c_* = max(1, C_Phi v^2) fails.

## New tests (`t3_sec10_new.py`)

| id | display | test | result |
|---|---|---|---|
| D10.10 | (10.10) | r = tg/e, f' = Cew squarefree and coprime on 3000 admissible tuples of actual primes; source text of OAI1112 Lemma 7.3 | PASS |
| D10.30 | (10.8), (10.30) | N(hf^2) <= 4 C_Phi v^2 L^2 F^2 N(e0)/(H N(C)^2) <= Y_{C,e0}, exact rationals, largest actual N(h), N(f); 1953 cases | PASS; falsifier c_* = max(1, C_Phi v^2) fails |
| D10.11 | (10.11), (10.37) | N(h2) <= C_Phi (2v)^2 ell^2 N(e)/(Y N(g)^2) gives N(k') <= HL/(Sigma F N(r)^2); 856 exact cases | PASS; falsifier b_* = 4v fails |
| D10.16 | (10.16) | all three inequalities, 50 000 retained (r, f', k') with f' \| k' in Z[omega], every j \| r, exact | PASS; falsifier N(r) for N(r)^2 fails |
| D10.21 | (10.21) | F/(XF)(X\|\|W\|\|)^2 = X\|\|W\|\|^2 <= Sigma\|\|W\|\|^2 (sympy) | PASS |
| D10.35 | (10.35) | 34 blocks (C,t,e0), L=20, F=7, Sigma=150, H=10, CM coefficients mu conj(alpha) xi_1: every hf^2 in N(y) <= Y; multiplicity <= #{f : f^2 \| y} (max 2); sum_{f,h}\|P(hf^2)\|^2 <= m sum_{N(y)<=Y}\|P(y)\|^2 | PASS (worst ratio 0.022); falsifier fails |
| D10.38 | (10.38) | exact on 3000 random (z1, z2, C, e0 \| C, e, h2), shared primes allowed; conj(chi_z(e)) = chi_z(e)^5 | PASS (6e-15); 2 falsifiers fail |
| D10.39 | (10.39) | exact on 4000 random tuples (2486 with nonzero left side) | PASS; falsifier 1_{(n,t)=1} fails |
| D10.41 | (10.41), (10.17), (10.40) | 2(2m+4)+4 = 4m+12; E <= A(W) for Phi = e^{1-t} on 3 computed averages; C^{2j+4}, C^{2m+4}, C^{4m+12}, [u/16,4v] in the OAI1112 source | PASS |
| D10.27 | (10.27) | zeta column in the zeta guide source (sources/zeta_11_12_guide.tex) and OAI1112: mu_K(n), conj(alpha) gamma_2, conj(alpha(b))^3, one zbar-derivative, alpha(ell) | PASS |
| S10.32 | (10.32) | size check: Z/(H Phihat(0) \|\|W\|\|^2) <= 0.084 for L <= 3000, F <= 100, no growth | PASS (size check) |
| S10.34 | (10.34) | size check: sup over the Phihat argument of \|\|K_c\|\|_{C^q(I^2)}/\|\|W\|\|_{C^q}^2: 0.25, 0.06, 0.008 for q = 0, 1, 2 | PASS (size check) |

## Old test ids versus the new display numbers

`section10_descent` uses the display numbers of 2026-10-12. Since then the guide has numbered three more displays
in Section 10: (10.5) eq:d057, (10.6) eq:P and (10.10) eq:d059. It has also numbered all the displays of the
proofs, (10.27)-(10.41). Matched by label and content:

| old id | new display |
|---|---|
| S10.5-7 | (10.7)-(10.9) |
| S10.8 | (10.11) |
| S10.9 | (10.12) |
| S10.10 | (10.13) |
| S10.11 | (10.14) |
| S10.12-13 | (10.15)-(10.16) |
| S10.15 | (10.18) |
| S10.16-20 | (10.19), (10.22), (10.23) |
| S10.18 | (10.21) |
| S10.22 | (10.25) |
| S10.23 | (10.26) |

## Failures and suspected errors

None. No display failed, and I found no display that I think is wrong.

Remarks, none of them errors:
- (10.33) restricts to N(C)N(t) <= vL, while the analogous display in OAI1112 Section 7.1 has <= 2vL. Both are
  correct: for W supported in [u,v], every term with N(C)N(t) > vL vanishes. The range <= 2vL reappears, correctly,
  in Q (10.31), where U is supported in [u/2, 2v].
- S10.10 (for (10.13)) and B79 use floating-point samples with 1e-12 slack. The new D10.16 checks (10.16), and with
  it the step F' <= HL/(Sigma F N(r)^2) of (10.13), in exact rationals.
- (10.27) is not an identity in itself. Each of its rows is tested end to end, with the zeta form required to fail:
  - mu alpha: B17
  - conj(alpha)^2 gamma_2: B27
  - d_zbar^2: B50, B53
  - conj(alpha(b))^6: B46, B68
  - alpha(ell)^2: B55

  The zeta column is checked against the source texts. The zeta guide (sources/zeta_11_12_guide.tex) does not display the
  dual factor alpha(ell) itself, so that entry is checked in OAI1112.
