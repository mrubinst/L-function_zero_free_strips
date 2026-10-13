# S6 report (item 6: beta_* <= 11/12 for F_CM), saved by the main agent from the agent's final message

Result: Theorem 2 of a new Appendix B (inserts.tex, before \begin{thebibliography}); edits E1-E5. Test build 61 pp.
Main finding: beta_* <= 11/12 is a per-member statement (no uniformity in eta or conductor is needed in the conclusion),
so it is proved one character at a time by running Part I with the twist nu = conj(eta) in the columns, exactly as
OAI1112 proves Thm 1.1 from Prop 3.1 with a finite-order twist.
- B.2 Theorem 2: S enlarged to primes over 6 and the modulus of eta; class group G_S (classes mod L_S, = 1 mod 3);
  paired identity (7.9) with xi_1 = conj(nu); Prop 8 (twisted mean square) re-checks every written-out CM proof that
  mentions 6, 72 or G (Sections 7.3, 8, 9.1, 10, App. A); Lemma 13 = App. A at general modulus L supported on S
  (Steps 1-2 from OAI1112 App. A.2, stated for general L, paper2.tex L2928-2931; Steps 3-6 independent of L);
  extraction with rows pi^6 (pi not in S), Mellin inversion, Hecke's theorem for L(s, conj psi).
- B.3 Prop 9: CM forms of the (U2) shared lemmas (OAI78 L15495-16460), assuming items 2 and 3. Target enters in five
  ways: correction H_psi (P_0 chosen before the target); global reciprocal bound (CM log control, conductor << q_u);
  bin ceiling a <= beta_* (OAI78's theorem for finite-order labels); principal clauses vacuous; numerators finite order.
- B.3 Lemma 14: sextic row count with CM labels (witness coefficients mu alpha^{+-1} theta chi^{+-1}); sieve takes
  arbitrary coefficients, R(delta) unchanged.
Findings: 0 ERROR, 1 GAP, 1 MISLEADING, 3 IMPRECISE.
- MISLEADING L2557-2558: "statement about the whole family, so not supplied by Part I" -> per-member; E1.
- GAP Section 16: CM forms of the (U2) lemmas not listed among the six; now Prop 9; E2.
- IMPRECISE: table row "4687-6807: none beyond the above" and Section 14 item 1: CM witnesses carry alpha^{+-1}; E3, E4.
- IMPRECISE: Section 16 "each is obtained from the corresponding statement of OAI78": item 6 cheaper via OAI1112; E2.
- IMPRECISE (flagged for s1): Prop cm-completed rests on item 1 (general modulus, log-Schwartz V_G).
Tests (run_all.sh, ~6 s, all pass, falsifiers fail): T2a (7.9) over G_S, 144 characters, 1.2e-14; T2b (7.10) 1.1e-14;
T1 twisted first Poisson reduction end to end <= 1.5e-14 at (120,170),(120,40),(200,300); T3 exponent bookkeeping;
T4 unit rows nonprincipal; T5 definitions.
Open: Prop 9 and Lemma 14 depend on items 2 and 3; Lemma 13 overlaps item 1's general-modulus part (merge with s1);
end-to-end test of Lemma 13 at L = 18 pi_13 not run (cost); abstract/Section 16 "conditional" wording to update.
