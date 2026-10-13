# Tests for agent s2 (statement 2: Lemma cm-growth and the CM logarithmic control)

All scripts are self-contained (Python 3, mpmath, numpy not needed) and use `cmhecke.py`.
Run from this directory: `python3 t1_theta_fe.py` etc.; outputs are in `*.out`. All four report ALL PASS.
Each script exits with status 1 if a check fails or a falsifier does not fail (added 2026-10-12; the falsifier
verdicts are summarized in one final "falsifiers:" line). With `CM_TESTS_FORCE_FAIL=1` the first check counts as
failed, as a self-test of this wiring.

Caveats (reviewer, 2026-10-12):
- The values of L(s, psi) in t1 (part (b)) and t2 are computed from the split-integral formula, which builds in the
  functional equation with the root number W. The left-line identity of t2 therefore only checks consistency. In
  t1 the split-point independence is equivalent to the functional equation with this W, and the comparison at
  s = 3 ties the formula to the Dirichlet series. These checks reject a wrong W (falsifier W -> -W) but still use
  the functional equation as input. The PARI `lfuncheckfeq` checks in `../s3` (t3_fe.py, S3-T3) start from the
  Euler coefficients, so they give stronger, independent evidence. The theta transformation (t1, part (a)) and the
  reflection identity (t3) compare two independently computed sides.
- The conductors of t4 are computed by enumerating all residues modulo the given modulus, so they are exact for
  the 18 cases tested. The 180 conductor cases of `../s3/t2_conductors.py` instead use randomized congruence
  sampling, so they are checks, not exhaustive exact conductor computations.

## Test characters (cmhecke.make_test_chars)
psi((x)) = phi(x) conj(x)/|x| with phi a primitive character mod f, phi(unit) = unit (Lemma s2-pres):

| name | phi | f | Q | root number W = -gamma(phi) sqrt(Q)/f |
|---|---|---|---|---|
| A | unit = x mod 3 (psi = alphabar) | 3 | 9 | +1 (= root number of 27a, y^2=x^3+16) |
| B | chi_pi, N(pi)=7 | pi | 7 | 0.99799+0.06333i |
| C | zeta(x)^5 chi_pi, N(pi)=13 | 3 pi | 117 | -0.96425-0.26501i |
| D | zeta(x)^3 chi_pi^4, N(pi)=7 | lambda pi | 21 | 0.96804+0.25079i |

Conductor exponents at lambda: 2 (A, C), 0 (B), 1 (D). All four are checked primitive and satisfy phi(zeta)=zeta.

## t1_theta_fe.py (Lemma s2-fe)
* |gamma(phi)| = 1 and |W| = 1 (errors 1e-16).
* Theta transformation Theta(v) = -gamma/(f sqrt Q) v^{-2} conj Theta(1/(Qv)) at 4 values of v: rel. error <= 3e-15.
  Falsifiers (B, C, D; for A, whose theta series is real, three of them coincide with the truth):
  no conjugation (rel. err >= 6.6e-4), conj(gamma) (>= 0.13), sign + (2.0), conj(f) (>= 0.38).
* Functional equation via the split-integral (incomplete gamma) formula, which is independent of the split point
  iff the FE holds with this W: independence at s = 0.5+3i, 0.2-7i, -0.1+11i, 0.8+0.3i to <= 1e-12 (60-digit
  arithmetic); and the formula reproduces the Dirichlet series at s=3 (err <= 3e-11).
  Falsifiers: W -> -W gives split dependence 0.57-1.43; dual coefficients psi instead of conj psi give 0.13-14.8
  (for the real character A this coincides with the truth).
* W(alphabar) = +1, and the coefficients of L(s, alphabar) equal a_p(E_16)/sqrt p for all 76 primes 5 <= p < 400
  (inert p: psi(-p) = -1, a_p = 0).

## t2_growth.py (Lemma s2-gamma, proof of Lemma cm-growth)
* G(sigma,t)=|Gamma(3/2-sigma-it)/Gamma(1/2+sigma+it)|: G/(3+|t|)^{1-2 sigma} stays in a bounded interval of (0,oo)
  for sigma in {-0.25,-0.1,-0.006,0,0.25,0.49,0.5,0.75,1.1,1.25}, t in [0,1e4]; G/|t|^{1-2 sigma} = 1.000000 at t=1e4;
  CM / finite-order quotient = 1.000000 at t=1e4. On sigma=-1/10: G/(3+|t|)^{6/5} in [0.107, 0.9996].
  Falsifier: exponents (1-2 sigma) +- 0.1 drift by a factor > 1.5 between t=1e2 and 1e4 (all sigma).
* Numerical constants: zeta_K(11/10) approximately 6.629847; the largest sampled value of
  G(-1/10,t)/(1.9^2+t^2) is approximately 0.315425. These give sampled boundary contributions
  C_left approximately 0.4455 and C_right approximately 0.6899. The sampling does not certify the supremum
  or prove a bound with constant 0.69; the guide uses only an unspecified absolute constant.
* Values of L(s,psi) (A-D) constructed from the split functional-equation formula on the grid
  sigma in {-0.1,0.25,0.5,0.75,1.1}, t in {0,2,-5,10,-20,30}:
  max |L|/(Q^{3/5}|s+2|^2) = 0.095 <= 0.69.
* Left-line FE identity |L(-1/10+it)| = (3Q)^{3/5}(2 pi)^{-6/5} G |L(11/10-it, conj psi)|: err 1e-55;
  falsifier (3Q)^{1/2}: err 0.26-0.44.
  This is an internal consistency check because the function is defined using its functional equation.
  Independent theta, Dirichlet-series and PARI checks are listed in t1 and the s3 suite.

## t3_reflection.py (Lemma s2-reflection)
* X^{-1/2} sum psi(a) W(Na/X) = W(psi) Y^{-1/2} sum conj psi(a) W#(Na/Y), Y = 3Q/((2pi)^2 X),
  M W#(u) = M W(1-u) Gamma(u+1/2)/Gamma(3/2-u), W(y) = exp(-(log y)^2): rel. err 8e-16, 6e-15, 2e-15 (B, C, D).
  Falsifiers: OAI78's finite-order kernel Gamma(u)/Gamma(1-u): rel. err 2.2-3.8; conj W: 0.13-0.53.
* W#(y)/y^{1/2} -> M W(3/2) with error ~ y (y = 1e-2, 1e-3, 1e-4: 1.4e-2, 1.4e-3, 1.4e-4): no constant term;
  the finite-order W#(y) -> M W(1) != 0.

## t4_conductors.py (Lemma s2-pres(iii), Lemma s2-rowcond)
* conductor of alphabar = (3) (norm 9).
* 18 finite-order psi_0 = zeta^k chi_pi^j (trivial on units; N pi = 7, 13, 19): the conductor of alphabar psi_0
  divides lcm((3), f0) and has the same exponent as f0 at pi; at lambda the pairs (v(f0), v(f*)) occurring are
  (0,2), (2,2), (2,1), (1,2), (2,0). Falsifier lcm((lambda), f0): fails for 10 of the 18.

## apply_and_check.py
Applies edits.md and inserts.tex to a copy of the guide: every OLD block and the insertion anchor occur exactly
once. The result (`../compile/cm-guide-s2.tex`) compiles with pdflatex (62 pages) with no undefined references
and no overfull boxes from the new text.
