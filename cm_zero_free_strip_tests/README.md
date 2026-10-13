# Numerical and symbolic tests of the CM notes

These scripts test the notes *A zero-free strip for the L-functions of the elliptic curves y² = x³ + d: the 11/12
method and its extension to 7/8* ([`cm_zero_free_strip.pdf`](../cm_zero_free_strip.pdf)). Section 18 of the notes
lists every numbered display with what was done to test it and the result; this folder contains the scripts behind
that table.

## Running the tests

Requirements:
- Python 3 with numpy, scipy, sympy and mpmath (tested with Python 3.12, numpy 2.2, scipy 1.16, sympy 1.14,
  mpmath 1.3);
- PARI/GP, with `gp` on the path (tested with version 2.17.3).

The main command reruns the tests of all 352 numbered displays, including the earlier suites they rely on. It takes
about half an hour:

```sh
cd displays_2026-10-13
./run_all.sh
```

It exits with status 0 if every check passes and 1 otherwise; every runner, however deeply nested, passes a failure
on. `python3 check_runner_propagation.py` in the same folder verifies this: it makes each script in turn fail, on a
scratch copy of this folder, and checks that the failure is reported (about 10 minutes). Set `PYTHON=/path/to/python3`
to choose the interpreter used by the shell runners.

Each test group, and each earlier suite, can also be run on its own; its README says how.

## Layout

| Folder | Contents |
|---|---|
| `displays_2026-10-13/` | The tests of the numbered displays, in ten groups: `t0` the displays (0.1)–(0.10) of the Introduction, `t1` (0.11) and Sections 2–6, `t2` Sections 7–9, `t3` Section 10, `t4` Sections 11–17, `t5` Appendix A, `t6` Appendix B, `t7` Appendices C–D, `t8` Appendices E and G, `t9` Appendix F. In each group, `rows.json` has one row per display: its kind, what was done, the script and the result. These rows make up the table of Section 18. |
| `part_A_lines_0001_0549/` … `part_D_lines_1510_end/` | The claim-by-claim tests of the version of October 9, 2026 (`cm-eleven-twelfths-guide_2026-10-09.tex` in this folder); their line numbers refer to that version. |
| `section10_descent/` | Tests of Section 10 (the descent). |
| `proofs_2026-10-12/` | Tests of the proofs of Part I (Sections 3, 5–7, 9, 10 and Appendix A). |
| `phase2_2026-10-12/` | Tests of the proofs of Part II (Appendices B–G, suites `s1`–`s6`) and independent checks (`r1`–`r3`). |
| `deuring_identity/` | A check of Deuring's identity: the coefficients of 1/L(E_d, s) as Hecke-character sums over ℤ[ω], for 18 values of d and n ≤ 600, against PARI/GP. |
| `sources/` | The TeX files some tests read: the notes and their `.aux` file, OpenAI's two papers, and our guides to the zeta case. See [`sources/README.md`](sources/README.md). |
| `cmsources.py` | The locations of those files; each can be overridden by an environment variable. |
| `check_failure_propagation.py`, `EXIT_CODE_AUDIT.md` | Checks and notes on how failures are reported by the earlier suites. |

## What the tests show

Identities and explicit inequalities are tested on finite samples, exactly where the quantities are algebraic
and otherwise in double or multiprecision arithmetic. Most such tests are paired with a deliberately wrong variant,
typically the zeta normalization, a wrong sign or a wrong power of the angular character, which must fail; this shows
that the test can tell the CM form from the zeta form. Estimates with unspecified constants cannot be verified by
computation; they receive size checks on finite ranges, or none, and rest on their proofs in the notes. Zero counts
come from the argument principle on finitely many boundary samples, and maxima are taken over sampled points; such
results hold on the tested grid only.

The saved outputs in each folder (`outputs/`, `out/`) are from the run of October 13, 2026, in which every check
passed.
