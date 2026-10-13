# b3 tests: new identities and inequalities in b3.tex and b3_edits.md

Scripts written for this block are in `drafts_cm_fill/compile_b3/tests/`. All were run, and all pass.
- `b3_symbolic_checks.py` (about 1 minute; mpmath at 40 digits, sympy) covers the constants of Steps 3-6, the kernel
  growth and the chain rule.
- `t_b3_prop_endtoend.py` (45 s; output in `t_b3_prop_endtoend.out`) is an end-to-end check of Proposition
  `prop:b3-reflection` against the definition eq:T. It uses copies of `fastgauss.py` and `zw.py` from
  `tests_cm/part_B_lines_0550_0997`.

Existing scripts that already cover a claim are named by path. "B47" etc. are the claim ids of
`tests_cm/part_B_lines_0550_0997/claims.py`. I reran `t5_reflection_j2.py` (PASS, 41 cusp configurations, worst
relative error 2.3e-8, all five falsifiers fail as required) and an earlier independent script (OK for j = 1..4).

| id | claim (location) | kind | test | status |
|---|---|---|---|---|
| N1 | $d_\sigma(\ell)=0$ unless $\ell\in\lambda^{-4}\OO$; only $\overline\theta$ at the identity cusp has a constant term $\frac{3^{5/2}}2v^{2/3}$ (INSERT 1) | cited [QF Thm 23] | part B `t4_theta.py` B40, B52 | covered |
| N2 | $\widehat V(s)=\int_0^\infty V(y)y^{s-1}dy$ (E2) | definition | used with this normalization in `t5_reflection_j2.py` and `t_b3_prop_endtoend.py` (Gaussian $\widehat V$) | covered |
| N3 | Proposition `prop:b3-reflection`: $T(X;\Psi)$ = sum of terms (eq:reflection), with $C=-\overline{\alpha(c)}^{\,4}C'/81$ (E7, INSERT 2) | identity | **new** `t_b3_prop_endtoend.py`. For $\mathcal P=\emptyset$ and three characters $\xi_0$, it compares $T(X;\Psi)$ computed from eq:T with $-\frac1{81}\sum_{h_0}\widehat\phi(h_0)\overline{\kappa(g_1)}\,\overline{\alpha(c)}^{\,4}\sum_\ell d_\sigma(\ell)\alpha(\ell)^2\N(\ell)^{-1/2}\ce(-\delta\ell/c)V^\sharp(\N(\ell)X/\N(c)^2)$, i.e. Steps 1-6 before grouping. 36 translates per twist; all three cusp types and $v_\lambda(c)\in\{0,1,2\}$. Relative error at most $1.1\cdot10^{-12}$. Controls fail: $\overline{\alpha(c)}^{\,2}$ (error 1.0, at twist B), zeta prefactor $-i/81$ (1.4), one-derivative kernel (2.3). Also `t5_reflection_j2.py` (B53, B55: single cusps) and part D `t08_reflection_j.py` ($j=2$) | PASS |
| N3a | grouped form with $\mathcal P\ne\emptyset$: the sum over $h_\pi\ne0$ gives $\chi_\pi(\sigma_\pi)^{-2}\varpi_{\pi,j}B_{\pi,j}$ (Step 2(e), Step 6) | identity, cited [OAI (A.12)] | part B `t4_theta.py` B57 (local factors as finite sums). Suggested: rerun `t_b3_prop_endtoend.py` with $\mathcal P=\{\pi\}$, $\N(\pi)=7$, $j=0,\dots,5$, using $X\approx2\cdot10^4$ and tables to $B\approx10^6$ (cusps up to $\N(c)=252$) | not run (cost) |
| N4 | Lemma `lem:b3-uniform` | cited [OAI Lemma 6.3] | not testable directly ($\N(M^2)$ is about $10^{41}$) | cited |
| N5 | $R(t)$ holomorphic in $\Re t>-\frac43$; $\lvert R(\eta+iu)\rvert\ll(1+\lvert u\rvert)^{4\eta}$ on $-\frac14\le\eta\le A$ (eq:b3-stirling) | inequality | `b3_symbolic_checks.py` (6): fitted exponents $-1.0,0,2.0,8.0,20.0$ at $\eta=-\frac14,0,\frac12,2,5$; part B B61; part D `t04_kernels.py` | PASS |
| N6 | Lemma `lem:b3-decay`: $(x\partial_x)^jV^\sharp(x)\ll\min(x^{1/4},x^{-A})$ | asymptotic | size check only. Suggested: compute $V^\sharp$ for a $C_c^\infty$ bump $W$ by quadrature on $\Re t=0$ and fit slopes of $\log\lvert V^\sharp\rvert$ at $x\to0$ ($\ge1/4$) and $x\to\infty$ (steeper than any fixed power) | partial (part B B61) |
| N7 | remark after Lemma `lem:b3-decay`: $\N(\ell)\ge3^{-4}$ on $\lambda^{-4}\OO\setminus\{0\}$; $\prod_{\mathcal A}\lvert B_{\pi,j}\rvert\le\N(c)^{1/2}$; tail $O(D^{-B})$ | elementary / asymptotic | $\N(\lambda^{-4}m)=\N(m)/81$; $\lvert B_{\pi,4}\rvert\le\N(\pi)^{1/2}$ by the definition. Suggested quick check of $\max_x\lvert B_{\pi,j}(x)\rvert$ for $\N(\pi)\le200$ | elementary |
| N8 | Bessel-Mellin (eq:b3-bessel), cited [QF (5.21)] | identity | `b3_symbolic_checks.py` (1); part B B45 | PASS |
| N9 | $\phi$ periodic mod 72, $(\lambda/x)_3=\omega^{-b/3}$ for primary $x=a+b\omega$ (Step 1) | identity | `t_b3_prop_endtoend.py`: supplement checked at 300 primary primes (0 failures); part B B42 | PASS |
| N10 | eq:b3-coeff: coefficient of $\Theta_\Psi$ at $\lambda^{-3}nb^3$ is $3^{5/2}\lvert b\rvert\gamma_2(n)\Psi(n)\Psi(b)^3$ | identity | part B B43, B47; implicitly N3 | covered |
| N11 | eq:b3-translates, $\sum_{\mathbf h}c_{\rm F}=0$, $C_{\pi,j}(0)=0$ ($j\ne0$), $C_{\pi,0}(0)=1-\N(\pi)^{-1}$, at most $5184\cdot2^{\lvert\mathcal P\rvert}$ groups | identity, cited [OAI (A.2)-(A.5)] | N3 tests the Fourier expansion (period 18). $C_{\pi,j}(0)$ by orthogonality. $5184=\N(72)$ | covered |
| N12 | Step 2(a): $\lambda^2h_0/72=h_0/(8\lambda^2)$, so $c_0\mid8\lambda^2$; $a\equiv\sigma_\pi h_\pi\pmod\pi$ | identity | `b3_symbolic_checks.py` (9). The congruence is elementary (Step 2(a) in the map); suggested quick random check | PASS / elementary |
| N13 | eq:b3-inversion (Mellin inversion of eq:T) | identity | implicit in N3. Suggested direct check: integrate on $\Re w=1.5$ with Gaussian $\widehat V$ | covered implicitly |
| N14 | eq:b3-mode: one mode gives $-\frac\pi2(2\pi)^{-2s}\Gamma(s+\frac13)\Gamma(s+\frac23)\overline{\alpha(\ell)}^{\,2}\N(\ell)^{1/2-s}$ | identity | `b3_symbolic_checks.py` (1): 9 values of $(\ell,s)$, relative error $10^{-40}$; falsifier one derivative: part B B45 | PASS |
| N15 | $\alpha(\lambda^{-3})=i$; $\overline{\alpha(\ell)}^{\,2}=-\overline{\alpha(n)}^{\,2}\overline{\alpha(b)}^{\,6}$; $\N(\ell)^{1/2-s}$ factorization; $3^{5/2}27^{-1/2}=3$; $\lvert b\rvert\N(b)^{3/2-3s}=\N(b)^{2-3s}$ | identity | `b3_symbolic_checks.py` (7); part B B46 | PASS |
| N16 | eq:b3-J: $\mathcal J(s)=\frac{3\pi}2(\frac{27}{(2\pi)^2})^s\Gamma(s+\frac13)\Gamma(s+\frac23)\mathcal T(s-\frac12)$ for $\Re s>\frac32$ | identity | `b3_symbolic_checks.py` (2) (constant, sympy); implicit in N3. Suggested: numerical $\int f(v)v^{2s-1}dv$ at $s=2$ from the theta series against the Dirichlet series | PASS (constant) |
| N17 | termwise integration for $\Re s>\frac32$; dual series for $\Re s'>\frac32$; the bound $\sum\lvert d_\sigma\rvert\N(\ell)^{1/2-\Re s'}\le6\cdot27\sum_m3^{m(2/3-\Re s')}\cdots$ | inequality | exponent bookkeeping: $3^{m/6}\cdot3^{m(1/2-\sigma')}=3^{m(2/3-\sigma')}$, $\lvert b\rvert\N(b)^{3(1/2-\sigma')}=\N(b)^{2-3\sigma'}$. Suggested: partial sums from the part B tables | bookkeeping |
| N18 | Step 4: the six derivatives of $z',\bar z',v'$; only $\partial_{\bar z}z'=-(cv)^{-2}$ survives at $z=0$ | identity | `b3_symbolic_checks.py` part 2 (all True); an earlier independent script parts (1)-(2); part B B48-B49 | PASS |
| N19 | eq:chain for generic $F$ (proof in Step 4) | identity | `b3_symbolic_checks.py` part 2: generic $F$, $j=1,2,3$, True; an earlier independent script ($j\le4$) OK; part B B50-B51. **Note:** an earlier independent script part (3) prints False because of a test bug (it sets $z=0$ but not $\bar z=0$) | PASS |
| N20 | eq:b3-cusp with actual $\theta$ | identity | an earlier independent script ($j=2$); part D `t02_chain_rule.py` | covered |
| N21 | eq:b3-G: $\partial_{z'}^2$ kills $v'^{2/3}$ and gives $(2\pi i\ell)^2=-4\pi^2\N(\ell)\alpha(\ell)^2$ | identity | part B B44, B52 | covered |
| N22 | $\lvert\ell\rvert\ge3^{-3/2}$ on $\lambda^{-3}\OO\setminus\{0\}$, $\lvert\ell\rvert\ge1/9$ on $\lambda^{-4}\OO\setminus\{0\}$ | elementary | $\lvert\lambda\rvert=\sqrt3$ | elementary |
| N23 | eq:b3-dual (dual Mellin, $\ell$ in place of $\bar\ell$) | identity | as N14 with $\bar\ell\to\ell$ (same code path) | PASS |
| N24 | substitution $v=1/(\N(c)v')$ and $c^{-4}\N(c)^2=\overline{\alpha(c)}^{\,4}$ | identity | `b3_symbolic_checks.py` (4): two values of $s$, relative error $10^{-41}$ | PASS |
| N25 | eq:b3-fe (functional equation) | identity | end to end in N3 and `t5_reflection_j2.py` | covered |
| N26 | $\mathcal T(w)=\mathcal J(w+\frac12)/(\frac{3\pi}2(\cdots)^{w+1/2}\Gamma(w+\frac56)\Gamma(w+\frac76))$; growth $(1+\lvert\Im w\rvert)^{2-4\sigma_2}$ on $\Re w=\sigma_2$ | identity / inequality | `b3_symbolic_checks.py` (5) (shifts) and (8): fitted exponents $6.0,4.0,2.4$ at $\sigma_2=-1,-0.5,-0.1$ | PASS |
| N27 | assembly constant $\frac{-\frac\pi2(2\pi)^{-2-2t}}{\frac{3\pi}2(27/(2\pi)^2)^{1-t}}=-\frac1{81}((2\pi)^4/27)^{-t}$; $\frac32-w=1+t$, $w+\frac56=\frac43-t$, $w+\frac76=\frac53-t$ | identity | `b3_symbolic_checks.py` (3), (5); part B B54 | PASS |
| N28 | eq:b3-C and $\lvert C\rvert\le1/81$ | identity / inequality | the CM part $-\overline{\alpha(c)}^{\,4}/81$ is checked in N3 and B55. The remaining factors are those of [OAI App. A.2], and $\lvert\widehat\phi\rvert\le1$ | covered |
| N29 | E4: $\Psi=\xi_0\prod\chi_\pi^{j_\pi}$ with $j_\pi\equiv v_\pi(k)+4v_\pi(f)\pmod6$, $\xi_0$ a character of $\mathfrak G$ | identity, cited [OAI (6.17)] | part B B56, B59 | covered |

Suggested additions, in order of value:
1. N3a, the grouped proposition with one active prime and $j=0,\dots,5$. It would exercise the local transform and
   the constant $C$ of eq:b3-C, including $\kappa_0$, $\varpi_{\pi,j}$ and $\chi_\pi(\sigma_\pi)^{-2}$. To build the
   cusp data, follow OAI's $\delta$ congruences. Use larger tables.
2. N16, a direct check of the Mellin identity from the theta series at $s=2$.
3. N6, slope fits for $V^\sharp$ with a compactly supported $W$.
