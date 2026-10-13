# S2 report (item 2: Lemma cm-growth + CM logarithmic control), saved by the main agent from the agent's final message
Result: statement 2 fully proved (new appendix app:s2-hecke, inserts.tex before thebibliography; edits E1-E6). Nothing open.
Counts: 0 ERROR, 3 GAP, 1 MISLEADING, 3 IMPRECISE.
1. GAP L2366-2367: missing reflected-line bound |L(s,psi)| << (3Q)^{1/2-sigma}(3+|t|)^{1-2sigma}|L(1-s-bar..)| needed by
   OAI78's detector (L4354-4367, L4461-4497) -> E1.
2. GAP L2408: conductor bound Q_psi << N(u) and factorization of zero-extended L for CM members of X_u asserted, not
   proved -> Lemma s2-rowcond: conductor | lcm((3), conductor of finite part), Q <= 9 Q_finite; same exponent away from
   lambda; proves "conductor scale changes only at 3" (L2513); factor 9 needed (exponent at lambda can rise or fall).
3. GAP L2386-2388: log control needed also for conjugates alpha theta chi_.(u) in X_u and, globally, for finite-order
   members (input OAI78's main theorem b >= 7/8); correct range b >= max(beta_*, 7/8) (= OAI78's under its hypotheses).
4. MISLEADING L2388-2389: "only fixed polynomial orders in the height change ... depend on the target": nothing
   changes (strip bound Q^{3/5}(3+|t|)^2; Gamma(3/2-s)/Gamma(1/2+s) has exponent 1-2sigma like Gamma(1-s)/Gamma(s)).
5. IMPRECISE L2366: constant must be absolute (|L| <= 0.69 Q^{3/5}|s+2|^2).
6. IMPRECISE L2510: plain-moment reflection needs its own short proof (Lemma s2-reflection; dual profile has no
   constant term at 0).
7. IMPRECISE L2374-2377 (also L2709, Sec 17): "xi-bar" in the Poisson step only for the Euclidean pairing; with the
   guide's e(z) the transform is a multiple of xi (theta series of psi-bar). Conclusion unaffected; E1 removes sketch.
Zeros in boxes: OAI78 never counts zeros; only discreteness in compact rectangles (L4325-4327); Jensen remark added.
Proved: Lemma s2-pres; Lemma s2-fe (theta series with harmonic polynomial, Lambda entire, bounded in strips, functional
equation with explicit root number W = -gamma(phi) sqrt(Q)/f, via [QF] Thm 19, Ex 20, Lemma 17; W(conj alpha) = +1,
matches y^2 = x^3+16, conductor 27); Lemma s2-gamma (Stirling); proof of Lemma cm-growth; Lemma s2-logcontrol (OAI78
L1554-1600 uses exactly four properties, each checked; matches OAI Lean logarithmic_control hypotheses); s2-rowcond;
s2-reflection (replaces OAI78 L12700-12731); Sec B.7 map of every use of the growth lemmas in OAI78.
Cited: OAI78 L1503-1519 (damped maximum principle), L1554-1600 (complex analysis of log control), deleted Euler
factors lemma, finite-order conductor bound (L4244-4254), OAI78 main theorem for finite-order rows, Stirling.
Tests (s2/tests): t1 theta transformation/FE for 4 characters (1e-12..1e-15; falsifiers O(1)); L(s, conj alpha) vs
point counts of y^2=x^3+16 at 76 primes; t2 Stirling exponents and explicit constants (worst ratio 0.095 < 0.69);
t3 reflection identity 1e-15 (OAI78 finite-order kernel fails); t4 conductor claims (weaker version fails).
Merge notes: E5 marks item 2 done; Section 16 count and closing sentence need updating; statements 3 and 5 can cite
s2-logcontrol, s2-rowcond, s2-reflection.
