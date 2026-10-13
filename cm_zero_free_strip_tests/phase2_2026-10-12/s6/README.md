# s6 tests: Appendix B (the 11/12 stage for F_CM)

Run `./run_all.sh` (about 6 s; python3 with numpy). The scripts import the helpers of
`tests_cm/part_B_lines_0550_0997` (cmlib.py, common.py, zw.py, fastgauss.py), unchanged; the path is the optional
first argument. Outputs are in `outputs/`. Conventions are the guide's: primary = 1 mod 3, lambda = 1 + 2 omega,
chi_n(x) = (x/n)_6 with zero values kept, e(z) = exp(4 pi i Im z/sqrt 3), alpha(n) = n/|n|,
a_xi = conj(alpha)^2 gamma_2 xi.

The test twist is nu(n) = exp(2 pi i dlog_2(n mod pi)/12) on primary n prime to S, with pi = 4 + 3 omega
(primary, norm 13) and S = {primes over 2, 3} + {pi}. Its order, 12, does not divide 6. So nu is not a sextic
residue character, and Part I of the guide does not cover psi = conj(alpha) conj(nu) without the twist. The
enlarged class group is G_S = (O/4)^x x (O/pi)^x, of order 144. G, R and chi(-1) factor through O/4, and nu
factors through O/pi.

**Last run (2026-10-12): every check passes and every falsifier fails.**

| id | claim (insert location) | test | result | falsifiers (all fail) |
|---|---|---|---|---|
| s6-T2a | Lemma lem:b2-paired, eq. (7.9) over G_S with xi_1 = conj(nu) (proof of prop:s6-ms, Step 2) | literal sum over all 144 characters of G_S; c_xi(conj nu) by finite Fourier analysis on G_S; brute-force Gauss sums gamma(conj chi_z1 chi_z2); 60 random coprime squarefree pairs (z1, z2) prime to 6 pi | max error 1.2e-14 | xi_1 = nu (2.0); G(t) for G(t^-1) (1.7); untwisted c_xi (2.0) |
| s6-T2b | eq. (7.10) over G_S (used in the transfer, Lemma lem:b4-first) | literal sum over the 144 characters xi_1, for xi = nu and one random xi, 30 pairs each | max error 1.1e-14 | beta = mu alpha xi_1 (2.0) |
| s6-T1 | twisted first Poisson reduction: M_D = Z + sum_xi c_xi(conj nu) S_xi, with S_xi in the final (r0, f, k) form (7.17), every column-type variable prime to S (proof of prop:s6-ms) | end to end. Left side: rows A_u(D) = sum mu alpha nu chi_n(u) W(N(n)/D) over all u, Gaussian Phi. Right side: class form conj(nu(m1)) nu(m2) G(m2 m1^-1) | relative error 2.6e-15 (D,H) = (120,170); 1.5e-14 (120,40); 2.4e-15 (200,300) | orientation nu(m1) conj(nu(m2)) (0.12, 0.18, 0.08); no twist on the dual side (0.08, 0.03, 0.04); r0, f not prime to pi (5e-3, 7e-3, 1e-3); zeta coefficients (0.07, 0.05, 0.06) |
| s6-T3 | exponent bookkeeping of OAI78 Part I, claimed unchanged for F_CM (Section B.3, route of OAI78; prop:s6-shared) | exact fractions: C_I(s) = s - 2/3 and C_II(s) = s - 11/16 from eq:stage-mellin-exponent; R(delta); the three pieces of E_I(1/2) on a grid of 601 values of delta; E_I(1/2) - Delta_1 <= -1/24 - Delta_1/2 under delta <= 5/6 + 2 Delta_1 (200 values of Delta_1); slope range [33/50, 62/75]; m_c = 1021/25000; small rows -43/300; large rows < -1 at v = 4000; m_w, m_z; 46 eps_* < 1/4800; E_I(d) = E_{11/12}(d; R, 0); Part II bootstrap Delta <= 1/24, kappa <= 5/6 | all exact identities hold | without the bin ceiling (delta = 1) the comparison fails |
| s6-T4 | "u = 1 is the only row with a principal numerator" (OAI78 L5502-5514; prop:s6-shared (d)) | for each unit u != 1, a primary prime p with chi_p(u) != 1 among the first 40; the six unit characters are distinct | pass | u = 1 is principal |
| s6-T5 | eq. (B.1): alpha(pi) nu(pi) = conj(psi((pi))); psi = conj(alpha) eta has infinity type z-bar/|z| and is not of finite order | 300 primes; 125 elements a = 1 mod 3 pi; arguments of psi at 50 primes are not multiples of 2 pi/M for M <= 1000. This is a sanity check of the definitions: the first identity holds by construction | pass | infinity type alpha (2.0) |

## Not tested numerically

- **The analytic inequalities.** Proposition prop:s6-ms (the twisted mean square), the completed estimate, cube
  removal and the transfer estimate over G_S are bounds and cannot be checked on finite data. Their algebra is
  covered by T1 and T2: the exact Poisson reduction, and the paired identities in both directions.
- **Lemma lem:s6-reflection end to end at a general modulus.** That would mean T(X; Psi) against the transformation
  formula with phi periodic modulo L = 18 pi. Its cusp denominators reach norm 468, so with the precision of
  `tests_cm/proofs_2026-10-12/b3/t_b3_prop_endtoend.py` the dual tables would need norms near 10^7 (V#(Kx) is
  still 1e-4 at Kx = 1e4 and 8e-10 at Kx = 1e5). That is out of reach for the brute-force Gauss-sum tables.
  - What is covered is the cusp-level identity for denominators with a prime factor outside 6. The single-cusp
    test `tests_cm/part_B_lines_0550_0997/t5_reflection_j2.py` includes C = -2 - 3 omega (norm 7) and C = 4 - omega
    (norm 21), with relative errors 1e-10.
  - Steps 1-2 at a general L are finite Fourier analysis plus the cusp bookkeeping of OAI1112 App. A.2, which is
    stated for general L.
- **Proposition prop:s6-shared and Lemma lem:s6-sextic.** They are structural statements about OAI78's proofs; the
  checks are the line-by-line reading of OAI78 recorded in the s6 agent's final report.
