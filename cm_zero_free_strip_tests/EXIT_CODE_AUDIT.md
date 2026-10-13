# Exit-code audit of tests_cm (2026-10-13)

Every test script must exit 0 when all its checks pass and every falsifier fails, and nonzero otherwise. The self-test `CM_TESTS_FORCE_FAIL=1` makes the first check of each script count as failed, so each script must then exit nonzero. Each script below was run once normally and once with `CM_TESTS_FORCE_FAIL=1`, in its suite folder, with the command its runner uses (full ranges; no quick modes), stdin from /dev/null. Times are wall-clock seconds with about 8 jobs in parallel on a shared machine, so they overstate single runs.

**Result: 76 runs (70 distinct scripts; b2_reduction_check, t_s1_endtoend_active and s6_twisted_reduction are run with several arguments, as in the runners). Normal runs: 76 of 76 exit 0. Forced runs: 76 of 76 exit nonzero.** Total script time: 56.1 min normal, 56.1 min forced; wall time of the parallel audit 17.5 min.

| script (command arguments) | normal exit | forced exit | normal s | forced s |
|---|---|---|---|---|
| `part_A_lines_0001_0549/t0_intro.py` | 0 | 1 | 2 | 2 |
| `part_A_lines_0001_0549/t2_integer.py` | 0 | 1 | 94 | 95 |
| `part_A_lines_0001_0549/t3_lattice.py` | 0 | 1 | 16 | 16 |
| `part_A_lines_0001_0549/t4_average.py` | 0 | 1 | 8 | 8 |
| `part_A_lines_0001_0549/t5_sixth.py` | 0 | 1 | 15 | 15 |
| `part_B_lines_0550_0997/t1_zeros.py` | 0 | 1 | 9 | 9 |
| `part_B_lines_0550_0997/t2_poisson.py` | 0 | 1 | 152 | 152 |
| `part_B_lines_0550_0997/t3_gauss.py` | 0 | 1 | 4 | 4 |
| `part_B_lines_0550_0997/t4_theta.py` | 0 | 1 | 26 | 26 |
| `part_B_lines_0550_0997/t5_reflection_j2.py` | 0 | 1 | 232 | 232 |
| `part_B_lines_0550_0997/t6_sieve.py` | 0 | 1 | 5 | 4 |
| `part_B_lines_0550_0997/t7_descent.py` | 0 | 1 | 14 | 14 |
| `part_C_lines_0998_1509/t01_curves.py` | 0 | 1 | 2 | 2 |
| `part_C_lines_0998_1509/t02_gauss.py` | 0 | 1 | 5 | 5 |
| `part_C_lines_0998_1509/t03_reflection.py` | 0 | 1 | 18 | 18 |
| `part_C_lines_0998_1509/t04_growth_fe.py` | 0 | 1 | 144 | 145 |
| `part_C_lines_0998_1509/t05_rows_moments.py` | 0 | 1 | 32 | 32 |
| `part_C_lines_0998_1509/t06_final_inequality.py` | 0 | 1 | 2 | 2 |
| `part_C_lines_0998_1509/t07_source_crossrefs.py` | 0 | 1 | 0 | 0 |
| `part_C_lines_0998_1509/t08_theorem_sanity.py` | 0 | 1 | 9 | 9 |
| `part_C_lines_0998_1509/t09_untested.py` | 0 | 1 | 0 | 0 |
| `part_D_lines_1510_end/t01_eta_quotients.py` | 0 | 1 | 1 | 1 |
| `part_D_lines_1510_end/t02_chain_rule.py` | 0 | 1 | 56 | 56 |
| `part_D_lines_1510_end/t03_mellin.py` | 0 | 1 | 34 | 34 |
| `part_D_lines_1510_end/t04_kernels.py` | 0 | 1 | 54 | 54 |
| `part_D_lines_1510_end/t05_misc.py` | 0 | 1 | 6 | 6 |
| `part_D_lines_1510_end/t06_sympowers.py` | 0 | 1 | 1 | 1 |
| `part_D_lines_1510_end/t07_theta_appendix.py` | 0 | 1 | 62 | 62 |
| `part_D_lines_1510_end/t08_reflection_j.py` | 0 | 1 | 348 | 348 |
| `section10_descent/test_sec10_descent.py` | 0 | 1 | 8 | 8 |
| `proofs_2026-10-12/b1/t_points.py` | 0 | 1 | 21 | 20 |
| `proofs_2026-10-12/b1/t_extract.py` | 0 | 1 | 1 | 1 |
| `proofs_2026-10-12/b1/t_bd.gp` | 0 | 1 | 0 | 0 |
| `proofs_2026-10-12/b2/b2_reduction_check.py 120 170` | 0 | 1 | 2 | 2 |
| `proofs_2026-10-12/b2/b2_reduction_check.py 120 40` | 0 | 1 | 8 | 8 |
| `proofs_2026-10-12/b2/b2_normalization_check.py` | 0 | 1 | 1 | 1 |
| `proofs_2026-10-12/b3/b3_symbolic_checks.py` | 0 | 1 | 5 | 5 |
| `proofs_2026-10-12/b3/t_b3_prop_endtoend.py` | 0 | 1 | 74 | 75 |
| `proofs_2026-10-12/b4/b4_checks.py` | 0 | 1 | 25 | 25 |
| `proofs_2026-10-12/b4/b4_weight_decay.py` | 0 | 1 | 18 | 18 |
| `phase2_2026-10-12/r1/test_classgroup.py` | 0 | 1 | 4 | 4 |
| `phase2_2026-10-12/r1/test_extraction_rows.py` | 0 | 1 | 1 | 1 |
| `phase2_2026-10-12/r1/test_step45_bookkeeping.py` | 0 | 1 | 7 | 7 |
| `phase2_2026-10-12/r2/t1_prop_grouped_endtoend.py` | 0 | 1 | 164 | 164 |
| `phase2_2026-10-12/r2/t2_grouping_cusp_data.py` | 0 | 1 | 14 | 15 |
| `phase2_2026-10-12/r2/t3_gl_family.py` | 0 | 1 | 3 | 3 |
| `phase2_2026-10-12/r3/t1_gauss_and_curves.py` | 0 | 1 | 3 | 3 |
| `phase2_2026-10-12/r3/t2_higher_weight.py` | 0 | 1 | 26 | 27 |
| `phase2_2026-10-12/r3/t3_final_inequality.py` | 0 | 1 | 2 | 2 |
| `phase2_2026-10-12/s1/t_s1_local_phases.py` | 0 | 1 | 47 | 46 |
| `phase2_2026-10-12/s1/t_s1_endtoend_active.py gauss 0` | 0 | 1 | 180 | 181 |
| `phase2_2026-10-12/s1/t_s1_endtoend_active.py gauss 1` | 0 | 1 | 175 | 175 |
| `phase2_2026-10-12/s1/t_s1_endtoend_active.py bessel 0` | 0 | 1 | 174 | 174 |
| `phase2_2026-10-12/s1/t_s1_endtoend_active.py bessel 1` | 0 | 1 | 167 | 167 |
| `phase2_2026-10-12/s1/t_s1_kernel.py` | 0 | 1 | 25 | 25 |
| `phase2_2026-10-12/s2/t1_theta_fe.py` | 0 | 1 | 46 | 46 |
| `phase2_2026-10-12/s2/t2_growth.py` | 0 | 1 | 223 | 223 |
| `phase2_2026-10-12/s2/t3_reflection.py` | 0 | 1 | 424 | 423 |
| `phase2_2026-10-12/s2/t4_conductors.py` | 0 | 1 | 0 | 0 |
| `phase2_2026-10-12/s3/t1_labels.py` | 0 | 1 | 4 | 4 |
| `phase2_2026-10-12/s3/t2_conductors.py` | 0 | 1 | 35 | 35 |
| `phase2_2026-10-12/s3/t3_fe.py` | 0 | 1 | 32 | 32 |
| `phase2_2026-10-12/s3/t4_gamma.py` | 0 | 1 | 36 | 36 |
| `phase2_2026-10-12/s3/t5_citations.py` | 0 | 1 | 0 | 0 |
| `phase2_2026-10-12/s4/s4_checks.py` | 0 | 1 | 25 | 25 |
| `phase2_2026-10-12/s5/t1_dual_profile.py` | 0 | 1 | 20 | 20 |
| `phase2_2026-10-12/s5/t2_hankel.py` | 0 | 1 | 1 | 0 |
| `phase2_2026-10-12/s5/t3_reflection.py` | 0 | 1 | 2 | 2 |
| `phase2_2026-10-12/s5/t4_lattice.py` | 0 | 1 | 1 | 1 |
| `phase2_2026-10-12/s5/t5_row_reflection.py` | 0 | 1 | 2 | 2 |
| `phase2_2026-10-12/s5/t6_kappa.py` | 0 | 1 | 0 | 0 |
| `phase2_2026-10-12/s6/s6_twisted_paired.py` | 0 | 1 | 0 | 0 |
| `phase2_2026-10-12/s6/s6_twisted_reduction.py 120 170` | 0 | 1 | 1 | 1 |
| `phase2_2026-10-12/s6/s6_twisted_reduction.py 120 40` | 0 | 1 | 2 | 2 |
| `phase2_2026-10-12/s6/s6_twisted_reduction.py 200 300` | 0 | 1 | 2 | 2 |
| `phase2_2026-10-12/s6/s6_margins_and_rows.py` | 0 | 1 | 2 | 2 |

## Runners

Each runner was run normally in place (this regenerated the `outputs/` folders) and once with `CM_TESTS_FORCE_FAIL=1` on a copy of `tests_cm` in a scratch directory, so that the forced outputs did not overwrite the real ones. A runner must exit 0 normally and nonzero when forced.

| runner | normal exit | forced exit | normal min | forced min |
|---|---|---|---|---|
| `part_A_lines_0001_0549/run_all.py` | 0 | 1 | 1.1 | 2.0 |
| `part_B_lines_0550_0997/run_all.py` | 0 | 1 | 3.8 | 6.4 |
| `part_C_lines_0998_1509/run_all.py` | 0 | 1 | 2.0 | 3.5 |
| `part_D_lines_1510_end/run_all.py` | 0 | 1 | 5.7 | 8.5 |
| `section10_descent/test_sec10_descent.py` | 0 | 1 | 0.1 | 0.2 |
| `proofs_2026-10-12/run_all.sh` | 0 | 1 | 1.7 | 3.2 |
| `phase2_2026-10-12/run_all.sh (with r1, s3, s5, s6)` | 0 | 1 | 20.6 | 23.8 |
| `phase2_2026-10-12/r1/run_all.sh` | 0 | 1 | 0.1 | 0.1 |
| `phase2_2026-10-12/s3/run_all.sh` | 0 | 1 | 1.2 | 1.2 |
| `phase2_2026-10-12/s5/run_all.sh` | 0 | 1 | 0.4 | 0.4 |
| `phase2_2026-10-12/s6/run_all.sh` | 0 | 1 | 0.1 | 0.1 |

## Notes on the runs

- **Part C `t07_source_crossrefs.py`.** Its first normal run exited 1 with an uncaught `StopIteration`, and the
  unmodified script from before these changes failed the same way. Its check C105 counts the items of the "Six
  statements must still be written out" list in the guide. The live guide, revised on 2026-10-12/13, no longer
  has that list. Part C tests the 2026-10-09 version, so `lib/cmreport.py` now reads the frozen copy
  `cm-eleven-twelfths-guide_2026-10-09.tex`. With it, t07 prints the same 38 claim lines as the stored output and
  exits 0. The table shows this rerun.
- **Forced runs.** Every forced run printed the `CM_TESTS_FORCE_FAIL=1` marker. Its only failure lines were the
  forced check and the summary lines that count it ("1 FAIL", "1 failures", "23/24 as required"). So each nonzero
  exit comes from the hook, not from another failure.
- **Normal runs.** Besides exiting 0, no normal run printed a line that the runners' secondary scan flags. Verdict
  lines and falsifier values agree with the stored outputs from before the change; only timings differ, plus the new
  summary lines.
- **Sub-runners.** The first forced run of `phase2_2026-10-12/run_all.sh` showed exit 1 for all 23 suites. The
  s5 sub-runner then printed no line the scan recognized, so failure was detected through its exit status only.
  After that run, the hook's stderr message became `EXIT STATUS 1: ...` and the sub-runners `r1`, `s3`, `s5`, `s6`
  now print `[FAIL] <script>: exit <rc>, failure lines <n>`. They were rerun (forced on the copy, normal in place);
  the table shows the rerun.
- **Other changes.** Another agent (CODEX) was editing tests_cm at the same time. It added `$PYTHON` interpreter
  selection to the runners, some of the final `raise SystemExit` lines, and `check_failure_propagation.py`, which
  checks the helpers and runners with fake interpreters. Its own validation results are in
  `proofs_2026-10-12/outputs/codex_rerun.txt` and `phase2_2026-10-12/outputs/codex_validation_status.json`.

How to repeat: run each runner normally, then once with `CM_TESTS_FORCE_FAIL=1` on a copy of the folder, since a
forced run overwrites `outputs/`. Every runner must exit 1 when forced and 0 otherwise.

## How failures reach the exit status

**Shared helpers (atexit hook).** Each helper below counts the failed checks it reports. An `atexit` hook then
flushes stdout and stderr, prints `exit status 1: n check(s) failed` on stderr and calls `os._exit(1)` if the count is
positive. So the script exits 1 even if it ends normally or calls `sys.exit(0)`. An uncaught exception still exits 1.
With `CM_TESTS_FORCE_FAIL=1` the first check reported through the helper is counted as failed and printed as a FAIL.

| helper | failure counted when | used by |
|---|---|---|
| `part_A_lines_0001_0549/cmcore.py` `report` | the line is `[FAIL]` (check false, or falsifier not detected) | Part A t0-t5 |
| `part_B_lines_0550_0997/cmlib.py` `report`, `falsifier` | `RESULT|..|FAIL`; a falsifier that DID NOT FAIL | Part B t1-t7, `section10_descent`, `phase2/s6/*` |
| `part_B_lines_0550_0997/common.py` `report` (also the copies in `part_C/lib`, `phase2/s3/lib`) | `[FAIL]` | imported with cmlib |
| `part_C_lines_0998_1509/lib/cmreport.py` `claim` (also the copy in `phase2/s3/lib`) | `[FAIL]` or `[CITED-X]`; forced mode turns the first line of any kind into FAIL | Part C t01-t09 |
| `part_D_lines_1510_end/pd_common.py` `report` | `[FAIL]` | Part D t01-t08 |
| `phase2_2026-10-12/s3/s3lib.py` `report` | `[FAIL]` (including a falsifier that did not fail) | s3 t1-t5 |
| `phase2_2026-10-12/s5/s5lib.py` `report`, `falsifier` | `RESULT|..|FAIL`; a falsifier that DID NOT FAIL | s5 t1-t6 |
| `phase2_2026-10-12/r1/r1lib.py` `Reporter` | already exited via `summary()`; only the forced-failure hook was added | r1 tests |

**Scripts patched individually** (they print their own results). Each now ends with `sys.exit(1 if a check failed
else 0)`, or `quit(status)` in GP, computed from quantities the script already computed and printed. No computation,
tolerance or printed line was changed, except as listed. Where a script printed values with no stated criterion, the
criterion is written in a comment and in one final `EXIT STATUS` line.
- `proofs_2026-10-12/b1/t_points.py`, `b1/t_extract.py`, `b1/t_bd.gp` (which also fails if GP stops with an error
  before the end), `b2/b2_reduction_check.py`, `b2/b2_normalization_check.py`, `b3/b3_symbolic_checks.py`,
  `b3/t_b3_prop_endtoend.py`: these have a final `EXIT STATUS` line. Criteria are listed in `proofs_2026-10-12/README.md`.
- `proofs_2026-10-12/b4/b4_checks.py` (its local `line`/`fals` count failures) and `b4/b4_weight_decay.py`. The
  latter has one new line: the zeta-quotient falsifier must differ from V# by more than 1e-3 at some x. It was
  printed before but never tested.
- `phase2_2026-10-12/r2/t1`, `t2`, `t3`; `r3/t1`, `t2`, `t3`; `s1/*`; `s2/t4`; `s4/s4_checks.py`. Each exits by
  its existing PASS/FAIL verdict and has the forced-failure hook. The other agent added several of the final
  `raise SystemExit` lines in parallel.
- `phase2_2026-10-12/s2/t1_theta_fe.py`, `t2_growth.py`, `t3_reflection.py`: they exit by their verdict. Their
  falsifiers were printed but not tested, so a final `falsifiers:` line now enters the verdict. In t1 each wrong
  variant must have rel. err > 1e-6 for at least one character, since for the real character alphabar three of
  them coincide with the truth, and 'W -> -W' must hold for every character. In t2 and t3 every falsifier must
  have rel. err > 1e-6 for every character.

**Runners.** `part_A`-`part_D` `run_all.py`, `proofs_2026-10-12/run_all.sh`, `phase2_2026-10-12/run_all.sh`
and its sub-runners `r1`, `s3`, `s5`, `s6/run_all.sh` all exit 1 if any script exits nonzero. As a secondary check
they scan each script's output for failure lines: `[FAIL]`, `[CITED-X]`, `RESULT|..|FAIL`, a bare `FAIL` at the start
of a line, `SOME FAIL(ED)`, `DID NOT FAIL`, `FALSIFIER NOT DETECTED`, `EXIT STATUS 1`, `<-- PROBLEM`, and
`n FAIL` or `n failures` with n >= 1. PASS lines and lines with `fails as required` are ignored, so `0 FAIL` and
`0 failures` never count. A scan of every stored output in `tests_cm` found no false positives.
`section10_descent/test_sec10_descent.py` is its own runner. It already exited via `allok`, and the cmlib hook now
also covers it.
