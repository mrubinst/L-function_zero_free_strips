# r3 tests: Part II, Sections 15, 17, 18 of the CM guide

All scripts are self-contained (`python3`, sympy, mpmath; PARI/GP `gp` for t1). Outputs of the last run are in
`outputs/`. Every PASS line that tests an identity also checks a falsifier (a deliberately wrong variant that
must fail); the line says which.

```sh
python3 t1_gauss_and_curves.py   # ~2 s (PARI for a_p)
python3 t2_higher_weight.py      # ~15 s
python3 t3_final_inequality.py   # ~10 s
# (the eight sub-audit scripts were working files of the referee and are not included)
```

`r3lib.py` is an independent implementation of the arithmetic of Z[omega] (primary = 1 mod 3, sextic symbols by
Euler's criterion in the residue field, Gauss sums by brute force over a complete residue system, OAI78's
e(z) = exp(4 pi i Im z / sqrt 3)). It does not import the helpers of `tests_cm/`.

## Results (last run: all PASS, every falsifier failed as required)

### t1_gauss_and_curves.py
| id | claim (guide line) | result |
|---|---|---|
| G1 | OAI78 eq:signal: gamma_2^3 = mu alpha, gamma_1 gamma_2 = mu alpha G (75 squarefree c: 45 primes incl. inert, 30 products of 2-3 primes) | err 1.3e-13 |
| G2 | eq:signal-cm (L2257): conj(alpha) gamma_2 gamma_1 = mu G | err 1.3e-13; falsifier alpha: 2.00 |
| G3 | companion (L2469): mu gamma_{-1} = chi(-1) G^{-1} conj(alpha) gamma_2, squarefree n | err 9.1e-14; falsifier alpha: 2.00 |
| G4 | L2262: conj(alpha)^2 gamma_2 eta gamma_1 = mu G psi (eta = chi_.(2+w)) | err 1.3e-13 |
| G5 | L2270: chi_p(4)^3 = 1, so b_p^3 = a_p^2 | 77 primes; falsifier chi_p(2)^3 = 1 fails at 42 |
| C1 | L2121-2122: sqrt(p)(psi*(pi)+psi*(pibar)) = a_p(E_d), p not dividing 6d, both primes (PARI ellap), 12 values of d, p < 1500 | 2756 pairs, err 1.4e-14; falsifier alpha chi(4d): 126 |
| C2 | L2122-2123: psi*(-p) = -1 and a_p = 0 for inert p not dividing 6d | 1461 pairs, exact; falsifier (no angular factor) never equals -1 |
| C3 | L2130: psi*(pibar) = conj(psi*(pi)) | exact |
| C4 | psi_d* is in F_CM: n -> chi_n(4d) is periodic modulo 144\|d\| on primary n prime to 6d (finite-order ray character) | 569 tests, 0 violations; falsifier (shift by 9) differs in 274 |

The conductors N(E_d) printed by t1 (via PARI) are for reference; the relation 3N(f) = N(E_d) was tested in
`tests_cm/part_C_lines_0998_1509` (t04).

### t2_higher_weight.py (Section 17)
| id | claim | result |
|---|---|---|
| H1 | chain rules (L2672-2673, L2742-2743): d_zbar^j and d_z^j at z = zbar = 0, j = 1..5, symbolic | exact; falsifier without the sign fails for odd j |
| H2 | direct Mellin transform with j derivatives (L2677-2682), j = 1..5 | rel. err 3e-31; falsifier (gamma arguments shifted by 1/2) 0.28 |
| H3 | R_j(1) = (3j-4)(3j-2)(3j+2)(3j+4)/1296, j = 1..12; R_1, R_2 are the two kernels of the table (L2324-2325) | exact; falsifier (3j-1) fails |
| H4 | rightmost pole of R_j at 1/6-(j+1)/2 (-5/6, -4/3, ...), left of -1/4; Stirling exponent 4 sigma, j = 1..5 | pass |
| H5 | Gamma(1-s+k/2)/Gamma(s+k/2) has exponent 1-2 sigma for k = 0..5; 6/5 at sigma = -1/10 (Lemma cm-growth); plain-moment profile exponent 2 sigma - 1 | pass |
| H6 | Sym^m factorization (L2748-2755): local Euler polynomials at split and inert primes, m = 1..8 | exact; falsifier (trivial middle character) fails for m = 2 mod 4 |
| H7 | eta(3z)^8: a_p = 2Re(pi^3), a_p = 0 and a_{p^2} = -p^3 for p = 2 mod 3 (incl. a_4 = -8); eta(2z)^3 eta(6z)^3: a_p = 2Re(pi^2 chi_pi(4)), a_p = 0, a_{p^2} = p^2 (p = 2 mod 3, p > 2), a_2 = a_4 = 0, a_3 = -3; to q^3000, 414 generators | 0 mismatches |
| H8 | expected strips 19/8, 5/2, 15/8, 2 and 11/8 | exact |

### t3_final_inequality.py (Section 15 and the r3 Lemma on bins with delta >= 5/6)
| id | claim | result |
|---|---|---|
| E1 | 1-h/4 = 51/64 for h = 13/16 (OAI78 L15513) | exact |
| E2-E4 | -E_* formula, J = j_y/108, completed-square identity eq:balanced-endpoint-certificate | exact; falsifier (49 -> 48) fails |
| E5 | grid minimum of -E_* on 0 <= delta <= 5/6, 0 <= x <= 1/2 is 2.29e-4 >= 49/440640 = 1.11e-4 | pass |
| E6 | -49/440640 + h Delta/4 - Delta = -49/440640 - (51/64) Delta (guide L2550 = OAI78 eq:adaptive-endpoint) | exact |
| E7 | with kappa = 3/4 the bound -49/440640 - Delta is smaller (guide L2552-2553) | exact |
| E8 | z_P(m) = 2(1-2m)/9 at kappa = 3/4; larger than at kappa = 3/4 + 2 Delta | exact |
| E9, E10 | large-delta endpoint E(h) = -1/48 - delta/16 at R = 1-delta, q = delta/2; floor bound -7/1200 | exact |
| E11-E13 | eq:capacity-comparison identity; <= Delta/4 for m >= 1/3, 2q <= 1 (2e4 random points); exactly 0 at kappa = 3/4 | pass |
| E14 | for 5/6 <= delta <= 1 and 1 <= r <= 3/2: 1-alpha+(alpha-delta)r <= 1-delta (Lemma r3-high-bins) | pass; falsifier: fails for delta < 5/6 |
| E15, E16 | E(d) = -1/48 - delta/16 + (d-h)(33/50 - delta/2) at R = 1-delta; slope in [4/25, 73/300]; E(h) <= -7/96 | exact |
| E17 | eq:adaptive-intermediate: root 529/625 = 0.8464; why bins above 5/6 need a separate argument | exact |
| info | the certificate with the amplification slope 5/6 replaced by 7/8, 9/10, 19/20, 1 has negative minimum (-8.0e-4 at 7/8). So 5/6 there is a fixed slope, not a bin ceiling that could be relaxed; bins above 5/6 must be treated by the no-slot count (Lemma r3-high-bins). | informational |

### subaudits/ (checks run by the eight OAI78 sub-audits; outputs in `outputs/subaudits/`)
- `s0001/`: OAI78's own Gauss-sum definitions at 44 primes and 41 composites (eq:signal, eq:signal-cm, companion,
  G = conj(chi(4)) gamma_3, G not multiplicative in 22/41 cases, G depends only on c mod 4 on primes); gamma
  quotients on Re s = -1/10; theta inversion for the conductor-(3) CM character confirming
  Lambda = (3Q)^{s/2}(2 pi)^{-s} Gamma(s+1/2) L (wrong scalings fail).
- `chk1648/`: direct and reflected single-mode Mellin integrals for one and two derivatives (both weights
  v^{2s-1}, v^{2s}), the reflected/direct constants -i 3^{-4} K^{-t} (OAI78 L2397) and -3^{-4} K^{-t} (CM), the
  chain rule for j <= 3, kernel bounds, V# contour independence.
- `scripts_3313/`: OAI78 eq:local-P summed by brute force with a generic unimodular (non-root-of-unity) target
  value against the closed form of lem:local-euler, all j = 0..5 (rel. err 4.5e-16): the local identity is formal
  in psi(p) (Prop. cm-rows); orientation of alpha in eq:signal.
- `t9209/`: both Gauss identities at 31 squarefree composites (4e-15, 3e-15).
- `work_6807/`: OAI78 old-eq:5.13b and eq:selected-local-operation hold identically for any unit target value.
- `scripts_12477/`: the CM Hecke reflection profile M W(1-s) Gamma(s+1/2)/Gamma(3/2-s) (ratio 1.0000-1.0005 for
  e = +-1; OAI78's finite-order profile fails); the CM lattice sums sum alpha(l)^e theta(l) W(N l/X) are ~1e-7 at
  X = 4e4 with and without a mask, versus ~X for e = 0.
- `scripts_14985/`: exact checks of the endgame arithmetic (both lines of eq:common-high-exponent, B0 = 53/32,
  floor -7/1200, small-row margins -79/800 and -63/800, supply margins 23/1443 and 2/185, the certificate
  identity, the capacity comparison, E(h) = -1/48 - delta/16, the root 529/625).
