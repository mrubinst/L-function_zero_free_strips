# Display tests, agent t1: front matter and Sections 2-6 of Part I

Guide: `sources/cm_zero_free_strip.tex` (equation numbers of `sources/cm_zero_free_strip.aux`,
checked against `sources/inventory/inventory_t1.json`). The 53 displays are (0.11) and (2.1)-(6.5): point counts,
the lattice and the curve's row, the average and (MS), sixth powers and the extraction, why the bound excludes zeros,
and the proofs added there. One row per display is in `rows.json`. The guide was not modified.

```sh
./run_all.sh          # about 100 s: the new tests (t1_*.py, about 10 s), then rerun_existing.sh (about 90 s)
./rerun_existing.sh   # only the re-run of the existing tests
CM_TESTS_FORCE_FAIL=1 ./run_all.sh   # self-test: the first check of each new script must fail (exit 1)
```
Exit status 1 if any script exits nonzero or prints a failure line. Every new check with a falsifier passes only if the
wrong variant fails. Needs python3 (numpy, sympy, mpmath) and PARI/GP. Outputs: `outputs/<script>.txt` for the new
tests and `outputs/existing_<name>.txt` for the re-runs. The existing suites are run read-only, from their own folders,
with `PYTHONDONTWRITEBYTECODE=1`; their own `outputs/` folders are not touched.

**Last run (2026-10-13): total 96 s, every script passed, every falsifier failed as required.**

## Counts

| kind | displays | result |
|---|---|---|
| E (exact identity) | 34 | 34 PASS |
| I (explicit inequality / exponent relation) | 1 | 1 PASS (symbolic) |
| A (asymptotic / analytic, unspecified constants) | 8 | size checks and symbolic bookkeeping only (see below) |
| D (definition / notation) | 7 | consistency checked |
| C (quoted from a source) | 1 | checked against [ZG, (1.2)] |
| R (restatement) | 2 | (0.11) restates eq:normgroup; (5.10) restates (5.8) |
| **total** | **53** | **35 of 35 E/I displays tested and PASS; no FAIL** |

**No display failed, and I found no display that is wrong.**

## New tests (this folder)

The existing tests cover most displays exactly. New scripts were written where no existing test checks the
display itself, or where an earlier README listed a test as "suggested, not run":

| script | id | display | what it checks |
|---|---|---|---|
| `t1_pointcount_steps.py` | T1-P1 | eq:d013 (3.11) | chi^2 is the cubic symbol (t^((p-1)/3) mod pi), chi^3 is the Legendre symbol, kappa(-1) = 1; 188 pairs (p, pi), p < 1200 |
| | T1-P2 | eq:d014 (3.12) | both root counts, exactly in Z[omega], at every residue (b1/t_points.py uses only their product) |
| | T1-P3 | eq:d016 (3.15) | all three equalities and the root count 1 + rho(1-4a) (b1/t_points.py checks only the outer equality) |
| | T1-P4 | eq:b1-S (3.16) | S = -rho(d)kappa(4d)pi and a_p = Tr(rho(d)kappa(4d)pi) directly, 11252 triples (b1/t_points.py checks this only through the chain) |
| `t1_euler_chain.py` | T1-E1 | eq:b1-euler (3.17) | each of the five members evaluated separately (truncated sums and products; PARI `lfun` for B_d/L), 7 values of d including non-minimal d = 320, s = 3 and 2.5+4i; max error 8e-10 against a truncation tolerance of 5e-6 |
| `t1_mellin.py` | T1-M1 | eq:d030 (6.4) | the D-integral of the absolute sum, computed as an integral, equals the product form to 5e-16; also \|q_d(m)\| <= r(m) and the convergence bound |
| | T1-M2 | eq:d031 (6.5) | F_W computed directly from S_d(D) by quadrature for two complex weights; F_W What' = What What' B_d/L = F_W' What to 8e-8 (b1.28, "suggested, not run") |
| `t1_sizes.py` | T1-S1..S3 | eq:d023 (5.8) | N_0^(1/12) = \|4d\|^(5/6) exactly; the Y < 7 remark (first term > D^(1+theta/2)/sqrt 7); the explicit trivial bound of Lemma b1-extract(i); sizes for 0 < \|d\| <= 30, D <= 1e4 (b1.25, "suggested, not run") |

Falsifiers used (all fail as required): chi in place of the Legendre symbol, chi^4 in place of kappa; 1 + rho(V) and
chi in place of kappa in the cube count; discriminant 1 - 2a, unconjugated kappa(4), conj kappa(2); +rho(d)kappa(4d)pi
and rho(d)kappa(d); u_d = 4d, B_d omitted, the unshifted L(E_d, s), conj(alpha) in the Euler product; y^(-sigma) in
place of y^sigma; the crossed pairing F_W What = F_W' What'; |4d|^(5/12); "> D^(1+theta/2)".

## Existing tests re-run (read-only, `rerun_existing.sh`)

| suite | script | ids used | status |
|---|---|---|---|
| part_A_lines_0001_0549 | t0_intro.py, t2_integer.py, t3_lattice.py, t4_average.py, t5_sixth.py | A0.07, A0.08, A2.01-A2.12, A3.01-A3.24, A4.02-A4.09, A5.02-A5.21 | all PASS (SANITY lines as before) |
| part_B_lines_0550_0997 | t1_zeros.py | B02, B04, B06, B07 | all PASS |
| proofs_2026-10-12/b1 | t_points.py, t_extract.py, t_bd.gp | aff, S, Jratio, J, rk, ap; T7-T9 | EXIT STATUS 0 |
| phase2_2026-10-12/r1 | test_extraction_rows.py | E1-E3 | 5 checks, 0 FAIL |
| deuring_identity | hecke_identity.py | identity on 10800 pairs, falsifiers F1-F4 | RESULT: PASS |

Old READMEs cite line numbers of the October 9 version. Displays were matched to tests by content and label. Each
cited test was read to make sure that it checks the display itself, not a neighbouring claim.

## The A displays

None of these can be proved by computation. What was checked:
- eq:MS (4.6): M(D, D^1.1)/D^2.1 is between 0.27 and 0.29 for D <= 1600 (A4.07). Proved by Proposition b2-reduction
  and Sections 7-10.
- eq:d019 (4.8): the error divided by tau(n)(sqrt H + 1) is at most 0.35 (A4.09). The display is identical to the
  zeta guide's eq:diagasym, [ZG, (4.2)].
- eq:approx (5.4): the explicit form of Lemma b1-extract(iii) holds with ratio at most 0.022 (r1 E2).
- eq:extraction (5.5): toy instances with the actual M(D, H) satisfy the proof chain (A5.09).
- eq:Sbound (5.7): the exponent choice was checked (A5.11) and the sizes of S_d(D) (B02).
- eq:d023 (5.8): the exact parts hold (T1-S1, T1-S2, A5.12). The size ratio is at most 4.3e-3.
- eq:d024 (5.9) and eq:d026 (5.11): only the exponent bookkeeping was checked (A5.10; b1 T8-T9). The constants c_1,
  Y_1 come from the prime number theorem and cannot be tested.

## Notes (not errors)

1. eq:d016: replacing 1 - 4a by 1 + 4a does **not** change the sum sum_a kappa(a) rho(1 +- 4a), because a -> -a and
   kappa(-1) = 1. So that variant is a falsifier only for the root count #{t : t^2 - t + a = 0}, and it is used there.
   The sum-level falsifier is 1 - 2a. The guide is correct.
2. eq:d005 is word for word the zeta guide's eq:cm, [ZG, (1.2)] (`sources/zeta_11_12_guide.tex`), apart from the names of the
   summation variables.
3. eq:ud (3.5) and eq:d010 (3.6) state the same identity, \overline{chi_n(4d)} = chi_n((4d)^5), so both are tested by A3.07.
   eq:d025 (5.10) is eq:d023 (5.8) restated as Lemma b1-uniform.

## Rows

| number | label | kind | test (relative to tests_cm/) | result |
|---|---|---|---|---|
| 0.11 | `eq:d001` | R | - | see eq:normgroup |
| 2.1 | `eq:d002` | E | part_A_lines_0001_0549/t2_integer.py A2.01 | PASS (exact) |
| 2.2 | `eq:d003` | E | part_A_lines_0001_0549/t2_integer.py A2.03 | holds (symbolic) |
| 2.3 | `eq:qd` | E | part_A_lines_0001_0549/t0_intro.py A0.07; t2_integer.py A2.05, A2.07 | PASS (exact) |
| 2.4 | `eq:integer-convolution` | E | part_A_lines_0001_0549/t2_integer.py A2.06 | PASS (max err 2.7e-15) |
| 2.5 | `eq:d004` | D | part_B_lines_0550_0997/t1_zeros.py B07; displays_2026-10-13/t1/t1_euler_chain.py T1-E1 | definition |
| 2.6 | `eq:reciprocal` | E | part_A_lines_0001_0549/t2_integer.py A2.07; displays_2026-10-13/t1/t1_euler_chain.py T1-E1 | PASS (exact) |
| 2.7 | `eq:target` | D | part_A_lines_0001_0549/t2_integer.py A2.09; displays_2026-10-13/t1/t1_mellin.py T1-M2 | definition |
| 2.8 | `eq:d005` | C | part_A_lines_0001_0549/t0_intro.py A0.08; t2_integer.py A2.10 | checked against source |
| 2.9 | `eq:d006` | E | part_A_lines_0001_0549/t2_integer.py A2.11, A2.12 | PASS (exact) |
| 3.1 | `eq:d007` | E | part_A_lines_0001_0549/t3_lattice.py A3.01 | PASS (max err 3.6e-12) |
| 3.2 | `eq:d008` | E | part_A_lines_0001_0549/t3_lattice.py A3.03 | PASS (exact) |
| 3.3 | `eq:d009` | E | part_A_lines_0001_0549/t3_lattice.py A3.05 | PASS (exact) |
| 3.4 | `eq:Au` | D | part_A_lines_0001_0549/t4_average.py A4.02, A4.04 | definition |
| 3.5 | `eq:ud` | E | part_A_lines_0001_0549/t3_lattice.py A3.07; deuring_identity/hecke_identity.py | PASS (exact) |
| 3.6 | `eq:d010` | E | part_A_lines_0001_0549/t3_lattice.py A3.07 | PASS (exact) |
| 3.7 | `eq:point-identity` | E | part_A_lines_0001_0549/t3_lattice.py A3.16; proofs_2026-10-12/b1/t_points.py | PASS (exact) |
| 3.8 | `eq:d011` | E | part_A_lines_0001_0549/t3_lattice.py A3.20 | PASS (max err 1.2e-15) |
| 3.9 | `eq:normgroup` | E | part_A_lines_0001_0549/t3_lattice.py A3.08; deuring_identity/hecke_identity.py | PASS (exact) |
| 3.10 | `eq:d012` | E | part_A_lines_0001_0549/t3_lattice.py A3.24, A3.08 | PASS (exact) |
| 3.11 | `eq:d013` | D | displays_2026-10-13/t1/t1_pointcount_steps.py T1-P1 | definition (consistency PASS, exact) |
| 3.12 | `eq:d014` | E | displays_2026-10-13/t1/t1_pointcount_steps.py T1-P2; part_A_lines_0001_0549/t2_integer.py A2.02 | PASS (exact) |
| 3.13 | `eq:b1-count` | E | proofs_2026-10-12/b1/t_points.py (aff) | PASS (exact) |
| 3.14 | `eq:d015` | E | proofs_2026-10-12/b1/t_points.py (S) | PASS (exact) |
| 3.15 | `eq:d016` | E | displays_2026-10-13/t1/t1_pointcount_steps.py T1-P3; proofs_2026-10-12/b1/t_points.py (Jratio) | PASS (exact) |
| 3.16 | `eq:b1-S` | E | displays_2026-10-13/t1/t1_pointcount_steps.py T1-P4; proofs_2026-10-12/b1/t_points.py (J, rk, ap) | PASS (exact) |
| 3.17 | `eq:b1-euler` | E | displays_2026-10-13/t1/t1_euler_chain.py T1-E1; part_A_lines_0001_0549/t2_integer.py A2.07; t3_lattice.py A3.08 | PASS (max err 8e-10) |
| 4.1 | `eq:M` | D | part_A_lines_0001_0549/t4_average.py A4.02, A4.07 | definition |
| 4.2 | `eq:coordinates` | E | part_A_lines_0001_0549/t4_average.py A4.02 | PASS (max rel err 1.4e-15) |
| 4.3 | `eq:double` | E | part_A_lines_0001_0549/t4_average.py A4.02, A4.03 | PASS (max rel err 1.4e-15) |
| 4.4 | `eq:d017` | E | part_A_lines_0001_0549/t4_average.py A4.04 | PASS (max err 1.4e-14) |
| 4.5 | `eq:allrows` | E | part_A_lines_0001_0549/t4_average.py A4.05 | PASS (exact) |
| 4.6 | `eq:MS` | A | part_A_lines_0001_0549/t4_average.py A4.07 | size check: $M/D^{2.1}\in[0.27,0.29]$ |
| 4.7 | `eq:d018` | E | part_A_lines_0001_0549/t4_average.py A4.08 | PASS (exact) |
| 4.8 | `eq:d019` | A | part_A_lines_0001_0549/t4_average.py A4.09 | size check: error$/(\tau(n)(\sqrt H+1))\le0.35$ |
| 5.1 | `eq:d020` | E | part_A_lines_0001_0549/t5_sixth.py A5.02; t3_lattice.py A3.09 | PASS (exact) |
| 5.2 | `eq:d021` | D | proofs_2026-10-12/b1/t_extract.py T7; part_A_lines_0001_0549/t5_sixth.py A5.03 | definition |
| 5.3 | `eq:d022` | E | part_A_lines_0001_0549/t5_sixth.py A5.05 | PASS (max err 6.6e-15) |
| 5.4 | `eq:approx` | A | phase2_2026-10-12/r1/test_extraction_rows.py E2; part_A_lines_0001_0549/t5_sixth.py A5.07 | size check: ratio to explicit bound $\le0.022$ |
| 5.5 | `eq:extraction` | A | part_A_lines_0001_0549/t5_sixth.py A5.09 | size check: explicit chain holds in 3 toy cases |
| 5.6 | `eq:exponents` | I | part_A_lines_0001_0549/t5_sixth.py A5.10 | holds (symbolic) |
| 5.7 | `eq:Sbound` | A | part_A_lines_0001_0549/t5_sixth.py A5.11; part_B_lines_0550_0997/t1_zeros.py B02 | size check: $\|S_d(D)\|/(\\|W\\|D^{11/12})\le0.033$ |
| 5.8 | `eq:d023` | A | displays_2026-10-13/t1/t1_sizes.py T1-S1, T1-S2, T1-S3; part_A_lines_0001_0549/t5_sixth.py A5.12 | size check: ratio $\le4.3\cdot10^{-3}$; exact parts PASS |
| 5.9 | `eq:d024` | A | part_A_lines_0001_0549/t5_sixth.py A5.09, A5.10 | not testable by computation (bookkeeping holds, symbolic) |
| 5.10 | `eq:d025` | R | - | see eq:d023 |
| 5.11 | `eq:d026` | A | proofs_2026-10-12/b1/t_extract.py T8, T9 | not testable by computation (bookkeeping holds, symbolic) |
| 5.12 | `eq:d027` | E | part_A_lines_0001_0549/t5_sixth.py A5.19 | PASS (exact) |
| 5.13 | `eq:quadratic-twists` | E | part_A_lines_0001_0549/t5_sixth.py A5.21 | PASS (exact) |
| 6.1 | `eq:d028` | D | part_B_lines_0550_0997/t1_zeros.py B04; displays_2026-10-13/t1/t1_mellin.py T1-M2 | definition |
| 6.2 | `eq:mellin-zero` | E | part_B_lines_0550_0997/t1_zeros.py B04 | PASS (max err/tol 0.09) |
| 6.3 | `eq:d029` | E | part_B_lines_0550_0997/t1_zeros.py B06 | PASS (max err 6e-39) |
| 6.4 | `eq:d030` | E | displays_2026-10-13/t1/t1_mellin.py T1-M1 | PASS (max rel err 5e-16) |
| 6.5 | `eq:d031` | E | displays_2026-10-13/t1/t1_mellin.py T1-M2 | PASS (max rel err 8e-8) |
