# S4 report (inverse moment eq:cm-inverse), saved by the main agent from the agent's final message

Result: statement 4 written out in full (inserts.tex, appendix sec:s4-inverse), conditional on Hypothesis (R) = the CM
form of OAI78 lem:reflected-energy (L7747-7993, with its width charge L7994-8030), which is part of statement 1.
12 new identity checks pass; 12 falsifiers fail. Compiles (62 pp with this insert).

Key points:
- CM classes = OAI78 classes times conj(alpha) on each column: conj(alpha) gamma_2 nu -> conj(alpha)^2 gamma_2 nu and
  mu nu -> mu conj(alpha) nu, for both signs epsilon (OAI78 conjugates the polynomial for epsilon=+1, L11712-11716).
- The descent carries conj(alpha) through: each of the three Poisson steps fixes (n1,n2) and is pointwise in the columns,
  so the CM identity is OAI78's times conj(alpha(n1)) alpha(n2) (lemma in the appendix); the factor splits into a factor
  on new columns and a unimodular outer factor (conj(alpha(A1)) alpha(A2) in the first transformation, 1 in the others),
  bounded only in absolute value by OAI78 (L10563-10568).
- Class switches are eq:signal times conj(alpha): closure conj(alpha)^2 gamma_2 -> mu conj(alpha) -> conj(alpha)^2 gamma_2.
- Masks, punctures, scales, margins, counts and seminorm propagation (L9348-9712) do not involve the coefficients.
- Exact identities written out: cube inversion (conj(alpha)^3 -> ^6); cube coefficient beta, polynomial P_a, outer
  product e_sigma of the first transformation; a_0 = conj(alpha)^2 gamma_2 and the child polynomial of the second; in
  the initialization the overlap identity (alpha(j)^eps; slot coefficients mu(p) conj(alpha(p))^eps a_i(p)), residual
  polynomial R_0, A_* = conj(alpha)^2 gamma_2 nu_*; the sixth-power column identity (alpha(d)^eps).
- Dual-coefficient bound in the tail exponent B_ref (old-eq:4.3) holds for the CM theta function by Appendix A.

Findings (0 ERROR, 1 GAP, 6 IMPRECISE; guide line numbers of the version before integration):
1. GAP L2594, L2453-2485: eq:cm-inverse was only sketched; now written out conditional on (R).
2. IMPRECISE L2473: "for eps=+1 they are the complex conjugates" -> after OAI78's conjugation the classes are
   mu conj(alpha) <-> conj(alpha)^2 gamma_2 for both signs.
3. IMPRECISE L2467-2471: CM switches are eq:signal-cm and companion multiplied by conj(alpha(n)); "checked at 60 primes"
   was primes only (now 366 squarefree n, 157 composite).
4. IMPRECISE L2479-2480: "never see alpha": the descent's identities carry conj(alpha) as a column factor; harmless.
5. IMPRECISE L2484-2485: anti-oriented case: descent would run; only the terminal step (needs theta itself) is not
   covered; nothing shows failure.
6. IMPRECISE L2488: statement 3's column identity picks up alpha(d)^eps.
7. IMPRECISE L2455-2466: statement of eq:cm-inverse omits OAI78's parameter-range and uniformity clauses (L9251-9268).

Open: Hypothesis (R) (part of statement 1). Anti-oriented case not covered and not needed.
