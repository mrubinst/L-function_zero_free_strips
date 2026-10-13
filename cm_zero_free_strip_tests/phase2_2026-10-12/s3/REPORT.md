# S3 report (item 3: detector with label set X_u), saved by the main agent from the agent's final message
Result: proved assuming Lemma cm-growth (item 2, now proved by s2). New Appendix "The zero detector for a CM target"
(inserts.tex), edits E1-E7, tests t1-t5 (18 checks pass). Patched guide 60 pp.
Written: label set Lambda = {(nu, s): nu in Theta_0 cup alpha^s Theta_0} (4|Theta_0| labels, independent of u) = X_u;
Lemma s3-labels (a)-(h): closure under conjugation; numerator/denominator/conjugates/slot characters are members;
complete multiplicativity; each CM member induced by a primitive Hecke character in F_CM or its conjugate (never
principal, L entire); conductor bound; ramification at every good p | u (deleted factors only over S); zeros
Re rho <= beta# = max(beta_*, beta_*^(0)); reflected-line bound with Gamma(s+1/2). Lemma s3-logcontrol. Prop s3-detector
(CM forms of OAI78 Lemma buffered-bins, bin ceiling, Lemma detector-dyads, Prop detector-witness; every
character-dependent step of OAI78 L4201-4686 cited by line and replaced). Downstream: sextic row count verbatim for CM
labels; CM witnesses co-oriented, normalised to s = 1.
Findings: 0 ERROR, 3 GAP, 5 IMPRECISE, 1 TYPO.
R1 GAP L2407-2409: detector input list incomplete (log control; CM presentations as Hecke characters of type
 (z/|z|)^{+-1}; conjugates alpha theta chi(u) not covered by cm-growth; conductor bound and deletion identity unproved)
 -> appendix, E3.
R2 GAP OAI78 L4373-4379: "every character inducing a member of X_u has finite order" false for CM; table implies
 verbatim -> Prop (ii), E5.
R3 IMPRECISE L2410-2411: bin ceiling needs beta_* >= 7/8 and conjugates -> E3, beta#.
R4 GAP IN OAI78 ITSELF, L4246-4248: Q << q_u "by reciprocity" but OAI78's lemma L660-661 states no bound on conductor
 exponents -> Lemma (e) by local class field theory: exponent exactly 1 at good p | u; <= 3 at 2; <= 4 at lambda;
 <= 1 at other primes of S; 72 u' is a modulus for the Kummer part; maxima 3, 4 attained; Q/N(u') reaches 5184 = N(72).
R5 IMPRECISE L2171-2172: principal clauses remain for finite-order numerators (u = 1, rows supported on S) -> E1.
R6 IMPRECISE L2279-2280: primitive character inducing the denominator lies in F_CM -> E2.
R7 IMPRECISE L2438: OAI78's bound U^{eps_1} P_i^{delta/2+O(e)}; guide drops power of U -> E4.
R8 IMPRECISE L2386-2388: log control also for conjugates -> E7. R9 TYPO L2408 Q_psi -> E3. Status: E6 marks item 3.
Answers: CM members never principal; finite-order members nonprincipal when u has a good prime factor (180/180);
S-supported rows removed. Conductor linear in N(u): Q <= N(m_A) N(u); Q_CM <= 9 Q_0; agree away from lambda (7/180
unramified at lambda). Bounded fibres not used. Co-oriented CM members suffice.
Tests (all pass): conductors and infinity types, 180 pairs; PARI lfuncheckfeq for 8 CM members (conductor norms 28..
117936, <= 2^-49; finite-order gamma, Q, 12Q fail); deletion identity; label-set properties (Theta_0 of order 216);
truncated-inverse coefficient identity; gamma quotient for |t| <= 1e12; 68 OAI78 line citations.
Merge: overlaps s2 (s3-logcontrol duplicates s2-logcontrol; s3-labels (d),(e),(h) overlap s2-pres, s2-rowcond,
s2-gamma); keep s3 (e) (s2 relies on OAI78's unproved bound R4). E1-E3, E5-E7 overlap s2 E2, r3 E17 etc.
