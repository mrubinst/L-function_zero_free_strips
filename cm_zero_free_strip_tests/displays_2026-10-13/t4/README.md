# t4: displays (11.1)-(15.2) and (17.1)-(17.9) of the CM guide

This folder tests the 29 numbered displays in Part II of `cm-eleven-twelfths-guide.tex` (numbering of 2026-10-13):

- Sections 11-15: the family, the probe for a CM target, rows/zeros/slots, counting bad rows, and the final inequality.
- Section 17: higher weight.

The inventory is `sources/inventory/inventory_t4.json`. `rows.json` has one row per display, in the format of
`BRIEF_DISPLAYS.md`. The guide, the sources and the other suites were not modified. The scripts of the
other suites were only re-run read-only (`PYTHONDONTWRITEBYTECODE=1`), and their stdout was captured here.

## Running

```sh
python3 run_all.py          # new tests, then reruns of the existing tests; 108 s wall clock (2026-10-13), exit 1 on any failure
python3 run_all.py --new    # new tests only, about 21 s
```

The run needs `python3` with numpy, sympy and mpmath, and PARI/GP 2.17 (`gp`, for `gcharinit`). Outputs are in
`outputs/`:

- `run_all.txt`: the full log;
- `<script>.txt`: the output of each new test;
- `rerun_<suite>_<script>.txt`: the output of each re-run test.

| file | what it does | time |
|---|---|---|
| `rerun_existing.py` | re-runs the existing tests below; requires exit 0 and status PASS (INFO for C103) for each listed id | 87 s wall, about 4 min CPU |
| `t4_sources.py` | S0: every inventory display is verbatim in the current guide, and its number agrees with the `.aux` file; S1-S3: text comparisons with G78 and OAI78, with falsifier needles | <1 s |
| `t4_local_euler.py` | L1a, L1b: the local content of the row formula (12.2) (new) | 2 s |
| `t4_lfunctions.gp/.py` | F0-F4: PARI Grossencharacters; the functional equation (13.1) for non-self-dual members of the CM family; the growth ratio (13.2); argument-principle sanity checks for (11.1) and (17.1) (new) | 12 s |
| `t4_moments.py` | M1, M2: size checks for (14.2), with a prime slot, and for (14.5) (new) | 7 s |

The re-run existing tests are:

- `part_C_lines_0998_1509`: t01, t02, t04, t05, t06, t08.
- `part_D_lines_1510_end`: t01-t06.
- `phase2_2026-10-12/r3`: t1, t2, t3.

## Counts

| kind | displays | outcome |
|---|---|---|
| E (exact identity) | 13 | 13 tested, 13 PASS |
| I (explicit inequality) | 2 | 2 tested, 2 PASS (exact) |
| A (asymptotic/analytic) | 6 | 6 size or sanity checks, all consistent; each is proved in the guide (place named in `rows.json`) |
| D (definition) | 7 | consistency checked: source text, implementation in a test, or model check |
| C (quoted) | 1 | (13.5) is character for character G78 eq:amplitude |
| R (restatement) | 0 | |

All **15 E/I displays were tested and passed**. Every PASS has a falsifier, a deliberately wrong variant, and every
falsifier failed as required. Thirteen E/I displays are covered by re-runs of existing tests that really check them;
the code of each was read. For two displays the test is new or extended:

- **(12.2) eq:d075, new.** This test checks the local content of the row formula. OAI78's local series
  eq:local-P (L3998-4006) is summed by brute force with the target value at p replaced by three choices:
  - a generic unimodular number, which is not a root of unity;
  - the CM value $\bar\alpha(p)\eta(p)$;
  - $\eta(p)$ itself.

  The sum is compared with the closed form of lem:local-euler, with $a_p=\bar\alpha^3\psi^3$ and
  $b_p=\bar\alpha^2\psi^2\overline{\chi_p(4)}$. The range is all $j=v_p(u)$, all six $\rho$, 8 split and 2 inert
  primes: 1080 identities, max rel. err $9\cdot10^{-16}$. Three falsifiers fail:
  - $\chi_p(4)$ unconjugated;
  - $\psi$ substituted only in the $c$-variable;
  - the closed form at $\bar t$.

  OAI78's evaluations of $g_{\chi_p^r}(p^t,p^{j'})$ and $C_p(t,k,j')$ were first re-derived by brute force (393 cases).
  So the closed form is formal in $\psi(p)$, as Proposition prop:cm-rows says. The derivation of eq:local-P itself
  (unit factors, pair phases) is taken from OAI78. In that series the target enters only through
  $\psi(p)^{e_0}$ and $a_p^l$.
- **(13.1) eq:d077, extended.** The existing C38 checks only the characters $\psi_d^*$ of curves, which are
  self-dual. Here, six members of $\mathcal F_{\rm CM}$ are tested, built as PARI `gchar`s mod 3 and mod 84 with
  conductor norms $Q=9,7,28,112,441,1008$, mostly with non-real root numbers. The checks are:
  - $\Lambda(s)/\overline{\Lambda(1-\bar s)}$ is the same unimodular $\epsilon$ at 5 points, to $8\cdot10^{-38}$.
  - PARI's conductor is $3Q$.
  - `gchar [1,0]` mod 3 is $\bar\alpha$ on primary primes, and $L(5/2+3i)$ agrees with an independent Dirichlet sum.

  The falsifiers $Q$, $\Gamma(s)$ and $\pi^{-s}$ fail.

The A rows are size checks, not proofs:

- **(11.1)** No zeros in $0.88\le\Re s\le2$, $|\Im s|\le30$, for the six characters above. The critical-box winding
  matches `lfunzeros`, which shows that the count detects zeros.
- **(11.2)** C98 was re-run: 5 curves.
- **(13.2)** The measured ratio equals
  $(2\pi)^{2\sigma-1}|\Gamma(\frac32-s)/\Gamma(\frac12+s)|(3+|t|)^{2\sigma-1}$ for every $Q$. Its sup on
  $[-\frac25,\frac75]\times\mathbb R$ is 1955, and the ratio is unbounded near $\sigma=3/2$, consistent with "compact
  $I\subset(-\frac12,\frac32)$".
- **(14.2)** With one prime slot, the ratio is 0.13 at $Z=10^4,4\cdot10^4,1.6\cdot10^5$, with growth exponent at most
  0.05. The existing C103 has no slot.
- **(14.5)** $z=0$: the ratio is 0.29-0.34 for $e=\pm1$. In the control $e=0$, the principal rows $k\in K^{*6}$ give
  252. This is why OAI78 excludes them; with the angular factor they are harmless (guide item (c)).
- **(17.1)** Two eta products and four characters with $k=2,3$ have no zeros right of $k/2+0.88$ up to height 30.

## Failures and suspected errors

**No display failed, and none appears to be wrong.** The following remarks are not errors.

1. (17.4) eq:d090: the evaluation point $(-\delta/c,\,1/(\N(c)v))$ omits the $\bar z'$-coordinate
   $-\bar\delta/\bar c$. Part D noted the same (its remark 4); the tests use the full point.
2. (12.2) eq:d075: (12.2) is a formal consequence once $\mathcal H_{\psi,u}$ is defined as the quotient. Its
   analytic properties are not testable by computation; they use only $|\psi(p)|=1$ (OAI78 lem:local-euler). Those
   properties are convergence in the two regions, $\mathcal H\ll\N(u)^\varepsilon$, and
   $\mathcal H_{\psi,1}=1+O(P_0^{-4/5})$.
3. (11.3) eq:d073: the scope remark of Part C, that the inert statement needs $p\nmid d$, is already fixed in the
   current text ("if $p\nmid6d$ and $p\equiv2\pmod3$").
4. (13.2) eq:d078: its implied constant depends on $I$ and blows up as $\sup I\to3/2$, at the poles of
   $\Gamma(\frac32-s)$. The guide's hypothesis that $I$ be compact is needed.
