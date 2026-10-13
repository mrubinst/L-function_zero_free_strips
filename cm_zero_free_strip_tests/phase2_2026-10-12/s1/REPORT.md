# S1 report (item 1: CM transformation formula at full generality), saved by the main agent from the agent's final message
Result: new appendix (inserts.tex) + edits E1-E11; patched guide 65 pp clean. 0 ERROR, 1 GAP (filled), 1 MISLEADING, 2 IMPRECISE.
GAP (item 1, L2586-2591): App. A covers modulus 72 and compact weights only; probe uses the Gaussian; marked lemmas need
 exact dependence of phases on moving primes -> appendix.
MISLEADING (Sec 12, L2305-2316): paragraph "The compensated probe needs more than this..." misplaced between
 prop:cm-completed and "[OAI78] proves this through..." -> E1, E4.
IMPRECISE (Sec 12, L2291-2304): B_m definition omits c, n primary with (cn,S)=1 (may share primes); V Gaussian (OAI78
 L3287-3300); m != 0, Q = Z^M, constant depends on a seminorm of V -> E1.
IMPRECISE (notation L2328, 2338, 2687, 2701-2702): dual index mu clashes with Moebius -> rename to ell (E2, E3, E10, E11).
Proved: Prop s1-reflection (CM completed reflection: general modulus supported on S, arbitrary local prime sets, all
 log-Schwartz tests; Mellin representation, entire continuation with polynomial growth, identity with alpha(l)^2,
 two-derivative kernel, scalar C = -conj(alpha(c))^4/81 (...) = -i conj(alpha(c))^2 x OAI78 scalar; line independence
 for Re t > -4/3; kernel bounds with OAI78 seminorms; App A Steps 3-6 reused). Cited identical (checked): OAI78
 L1889-2279 (never differentiates theta), derivative-free statements L2464-3311. Prop s1-marked (CM cor:marked-
 reflection); Lemma s1-unmarked (prop:cm-completed precise form); Lemmas s1-reflected-energy, s1-probe-row-norm;
 hybrid-energy unchanged. Properties (R1)-(R6) of the reflection used later in OAI78, each checked; no exponent changes.
Tests: end to end with active prime of norm 7, j = 0..5, both primes; non-compact tests (log-Gaussian 1.7e-12,
 e^{-y-1/y} 5.2e-8), ungrouped over 252 translates and grouped; phase bookkeeping 2.4e-12; marked formula 8e-13;
 6 falsifiers fail (0.43-4.4); local identities at 15 primes 1e-13; kernel constants, growth, line independence.
Open: with active prime only cusps (lambda,c)=1 tested; modulus with primes outside {2, lambda} not tested; energy
 lemmas audited not tested; slot normalization psi-bar(p) taken from Sec 16 table; assumed OAI78's three cusp
 coefficient functions = d_sigma of App A (tests consistent). E6 (Sec 16 lead) needs manual merge; E9 points item 4 to
 lem:s1-reflected-energy.
