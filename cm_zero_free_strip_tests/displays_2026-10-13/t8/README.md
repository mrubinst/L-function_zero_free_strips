# Agent t8: displays of Appendix E and Appendix G of the CM guide

Scope: the 35 numbered displays E.1-E.24 and G.1-G.11 of `cm-eleven-twelfths-guide.tex`. Appendix E is "The inverse
moment for co-oriented coefficients" (sec:s4-inverse). Appendix G is "The 11/12 stage for F_CM" (app:s6). They come
from the 31 inventory entries of `sources/inventory/inventory_t8.json`; the align block eq:s4-S1..S5 gives 5
rows. There is one row per display in `rows.json`.

```sh
./run_all.sh        # about 65 s wall time; exit status 1 on any failure; all outputs go to outputs/
```
`run_all.sh` runs the seven new scripts below. It also re-runs, read-only, the existing suites
`tests_cm/phase2_2026-10-12/s4/s4_checks.py` and `tests_cm/phase2_2026-10-12/s6/s6_*.py`, and writes their outputs
to `outputs/rerun_*.txt`. The s6 `run_all.sh` writes into its own folder, so it was not used. The helpers of
`tests_cm/part_B_lines_0550_0997` (cmlib.py, common.py, zw.py, fastgauss.py) are imported unchanged.

**Last run (2026-10-13): every check passes and every falsifier fails as required.** The 13 scripts took 65 s in
total. The forced-failure self-test (`CM_TESTS_FORCE_FAIL=1`) makes the cmlib-based scripts exit with status 1.

## Counts by kind

| kind | count | displays |
|---|---|---|
| E (exact identity) | 16 | E.9-E.13 (s4-S1..S5), E.14 (d173), E.15 (s4-inversion), E.16 (s4-beta), E.17 (d174), E.19 (s4-a0), E.21 (s4-overlap), E.23 (s4-Astar), E.24 (s4-sixth), G.4 (s6-euler), G.7 (d195), G.9 (d196) |
| I (inequality / parameter relation) | 0 | - |
| A (asymptotic bound) | 2 | G.6 (d194), G.8 (s6-bound): size checks only |
| D (definition) | 9 | E.1, E.6, E.7, E.18, E.20, E.22, G.2, G.3, G.5 |
| C (quoted from a source) | 2 | E.4 (d172 = hypotheses of OAI78 lem:marked), G.1 (d190 = OAI78 Props. balanced-low/high) |
| R (restatement) | 6 | E.2, E.3 (eq:d199), E.5 (eq:cm-inverse), E.8 (eq:s1-marked-completion), G.10 (eq:s4-M), G.11 (eq:s1-cm-completed) |

All 16 E displays were tested now, and all pass. Ten of them are covered by re-runs: S4.1, S4.2, S4.5, S4.6 and S4.7 of
s4, and s6-T2a and s6-T5 of s6. Seven were tested by new scripts: d173, s4-beta, d174, s4-a0, s4-Astar, d196, and
s4-overlap pointwise (T8.12, which adds to S4.6). Each E test has a falsifier, and every falsifier fails.

## New scripts (this folder)

| script | tests | what |
|---|---|---|
| `t8_coefficients.py` | T8.1-T8.5, T8.12 | d170 (complete multiplicativity of conj(alpha) on primary elements); s4-a (both forms); d174 (twist on n_a = A_a C t' x_a); s4-R0 (orientation step for both signs, anti-oriented falsifier); d173 (column-twist lemma on random fibrewise identities); s4-overlap pointwise in 4000 cases |
| `t8_signal_init.py` | T8.6, T8.7 | s4-Astar and s4-a0: m(z1)conj m(z2) gamma(conj chi_z1 chi_z2) = sum over theta of Ghat(theta) A_*(z1) conj A_*(z2), with G expanded over the 12 characters of (O/4)^x and 3 twists nu_0; the extraction identity with A_* |
| `t8_beta_block.py` | T8.8 | s4-beta: the reopened completion, split into (h, h') dyads, equals the block with beta^CM = conj(alpha(b))^3 beta, with one beta for 3 rows (so it is independent of k); X = 3000 and 2e5, including b = p^2 |
| `t8_euler_mellin.py` | T8.9-T8.11 | d196: exact Dirichlet coefficients of the Euler product, numerical values of 1/L^S(s, conj psi), and the Mellin side by quadrature of A_1^nu(D) |
| `t8_exponents.py` | EX1-EX4 | exact (sympy): d172 application (m=1, r=1/(1+theta)); d190 (C_I(11/12) = 1/4, omega = Delta_1/2); d198 (max(M,2M-1) at M=1); s6-bound extraction exponents |
| `t8_tex_compare.py` | 22 comparisons | normalized TeX of each display against its source (OAI78, OAI1112, or the guide display it restates), after applying the explicit list of CM edits; for each comparison with edits, the unedited source must differ from the guide |
| `t8_size_checks.py` | SZ1-SZ3 | size checks: d194 for D = 100..6400; s6-bound, A_1^nu(D) at 120 values of D up to 1e6 (mu_K by a sieve, cross-checked); s4-marked / cm-inverse for both signs, with and without a slot; plus the exact inequality (iii) of lem:b1-extract on the twisted rows |

Size-check results (kind A; these bounds are not provable by computation):
- d194: M/D^{2.1} is between 0.25 and 0.29 and M/diagonal is between 0.93 and 0.99, for D = 100..6400, with nu of
  order 12 and with nu = 1. The bound is proved in the proof of Proposition prop:s6-ms.
- s6-bound: |A_1^nu(D)| <= 0.57 D^{1/2}, while the trivial size is about 0.1 D. Also |A_1|/D^{11/12} <= 0.067, and
  it decreases to 0.002 on [1e5, 1e6]. The bound is proved in the proof of Theorem thm:s6-1112.
- s4-marked: the ratio sum |M_u Q_u|^2 / Z^m lies in [0.03, 0.29].

## Failures and suspected errors

**No test failed, and I found no display that I think is wrong.** Remarks:

1. eq:s4-M (E.2) and eq:d171 (E.3) repeat the two lines of eq:d199 (Section sec:78moments) symbol for symbol, once
   psi_u = nu chi_.(u)^eps is written out. eq:s4-T (E.8) is eq:s1-marked-completion under a second name: it writes
   T_d where the other has T^{(2)}_d. eq:d197 (G.10) is eq:s4-M with nu = theta and Z^r = D. These are duplicates,
   not errors.
2. eq:s6-euler (G.4) holds by the definition eq:d192 of nu, so its test is a consistency check.
3. The numerical claims in the guide's own lists match the re-runs:
   - Appendix E says "errors at most 3.3e-15"; the re-run maximum is 3.24e-15.
   - Appendix G says "1.2e-14" for the paired identities; the re-run gives 1.21e-14.
   - Appendix G says "1.5e-14" for the twisted reduction; the re-run gives 1.45e-14.
4. Caveat on `t8_tex_compare.py`. It turns \frac{A}{B} into (A)/(B) and drops the parentheses only around a single
   atom such as q_n or Z^r. It makes one explicit notational edit: OAI78's V_*(q_n q_b^3/X) becomes
   V_*((q_n q_b^3)/X). Everything else is compared string for string.
