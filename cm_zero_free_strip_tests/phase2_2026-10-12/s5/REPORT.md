# S5 report (item 5: plain fourth moment with alpha^e), saved by the main agent from the agent's final message
Result: proved for e = +-1, every kappa in [3/4,1], same generality/uniformity as OAI78 lem:plain; only zero-free
input OAI78 thm:main. New appendix (inserts.tex, before thebibliography), edits E1-E7. Test build 63 pp.
New (proved in full): Lemma s5-rowchar (primitive CM row character Xi_k inducing alpha^e psi_k^0; conductor agrees with
finite part's away from lambda, ratio in [1/9,9]; [alpha^e psi_k^0] = Xi_k 1_{(.,r_k)=1}; C_k^CM N(r_k) <= 27 C_k N(R_0k));
Lemma s5-reflection (direct Poisson over Z[omega] [QF Thm 19, Lemma 17] + Jacobi-Anger; planar transform of
(z/|z|)^e W(|z|^2) = pi e X (wbar/|w|)^e J_1-Hankel, zero at w = 0; X^{-1/2} sum Xi W = eps Y^{-1/2} sum Xi-bar W#,
Y = 3 N f/(4 pi^2 X), W#(y) = int W(x) J_1(2 sqrt(xy)) dx, |eps| = 1; eps = 1 for Xi = alpha (y^2=x^3+16));
Lemma s5-sharp (Mellin M W(1-s) Gamma(s+1/2)/Gamma(3/2-s); fits OAI78 kernel-seminorms, m holomorphic Re s > -1/2;
growth |t|^{2sigma-1}; W# = O(y^{1/2}) at 0); Lemma s5-lattice (zero frequency int alpha^e f = 0 exactly; masked
lattice sums O(Z^eps (1+|t|)^J) uniformly in X; OAI78 (2.18d-f) with c = 0); Lemma s5-coefficients (alpha(l1 l2)^e
survives every operation in OAI78 L12602-14778 as unimodular frozen scalars); Step 5: kappa = 3/4 via OAI78 thm:main.
Cited as identical: OAI78 L12581-14984 in 13 steps (mask deletion, induction order, comparisons, both row transforms,
amplifier, second transform and children, exceptional-row count, completion).
Findings: 0 ERROR; GAP L2525-2529, L2595 (statement 5 only sketched -> appendix, E6, E7); IMPRECISE L2499, L2518 ("Theta_0"
should be OAI78's Theta, L12511-12521, incl. reciprocity characters; E1, E4); L2510 (uses Lemma cm-growth, only sketched;
now proved directly; E2); L2513 (conductor change at lambda quantified [1/9,9]; E3); L2521 (uniform error bounds now
supplied by s5-lattice; E5).
Tests: 6 scripts, 106 checks pass, 98 falsifiers fail (Bessel-Mellin 2.5e-31; planar transform 8e-14; reflection for
12 primitive characters, e = +-1, 3 profiles, 3.3e-11; row reflection with deleted/restored Euler factors 1.3e-15;
lattice sums bounded and -> 0 up to X = 3e5 vs growth without alpha; kappa inequalities).
Open: nothing in statement 5. Merge: update "Six statements must still be written out" and caveat; appendix order;
move tests into tests_cm and list in Section 18.
