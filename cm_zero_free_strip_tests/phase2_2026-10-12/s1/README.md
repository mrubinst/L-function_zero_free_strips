# s1 tests: Appendix app:s1-general (item 1 of Section 16)

Scripts in this folder; outputs in `../out/`. Helpers `fastgauss.py`, `zw.py` are unchanged copies from
`tests_cm/proofs_2026-10-12/b3/` (originally `tests_cm/part_B_lines_0550_0997/`). Python 3, numpy, scipy, mpmath,
sympy. All scripts were run; all pass.

| script | what it tests | run time | result |
|---|---|---|---|
| `t_s1_endtoend_active.py [gauss\|bessel] [0\|1]` | Proposition prop:s1-reflection end to end with an **active local prime** and **non-compact log-Schwartz tests**; grouped form; cross-prime bookkeeping; Proposition prop:s1-marked | ~2 min each | PASS (4 runs) |
| `t_s1_local_phases.py` | the finite local identities of OAI78 L2012-2279 and L7443-7451 used by the proposition and the marked form | 25 s | PASS |
| `t_s1_kernel.py` | constants, the kernel $R_2$, the reflected test $V^\sharp$ for two log-Schwartz tests, growth of $\mathcal T_2$ | ~3 min | PASS |
| (removed) `apply_edits_check.py` | checks that every OLD block of edits.md and the INSERT anchor occur exactly once, and applies them | 1 s | OK; the result compiles cleanly |

## t_s1_endtoend_active.py

Data: base twist $\Psi_0=\overline{(\lambda/\cdot)_3}$ (twist 'A' of `t_b3_prop_endtoend.py`, so
$\phi=\ind_{x\equiv1(3),\,2\nmid x}$, period 6, used with $L=18$); one local prime $p$ of norm 7 (both primes of norm
7: `pi_index` 0 and 1) with exponent $j_p=0,\dots,5$; $X=3000$; tables to norm $1.5\cdot10^5$.
Tests: `gauss` $V(y)=\exp(-(\log y)^2/(2s^2))$, $s=\frac12$; `bessel` $V(y)=e^{-y-1/y}$ with
$\widehat V(t)=2K_t(2)$ (Mellin transform decays only like $e^{-\pi|u|/2}$).

For each $j$ it computes:
1. $T_V(X;\Psi)$ from its definition (squarefree $n$, all primary $b$ prime to 6, $\Psi(nb^3)$ with all zeros,
   $\overline{\alpha(n)}^{\,2}\overline{\alpha(b)}^{\,6}$);
2. the **ungrouped** right side: all 252 translates $z_{\mathbf h}=\lambda^2(h_0/18+h_p/p)$ with
   $\widehat\phi(h_0)\ne0$, each with its own reduced cusp $a/c$, brute-force matrix $g=GH$, Kubota multiplier
   $\kappa(G)$, $-\frac1{81}\overline{\alpha(c)}^{\,4}\sum_\ell d_\sigma(\ell)\alpha(\ell)^2\N(\ell)^{-1/2}
   \ce(-\delta'\ell/c)V^\sharp(\N(\ell)X/\N(c)^2)$;
3. the **grouped** right side (the statement): for each $h_0$ and the active branch, the scalar
   $C=-\frac1{81}\overline{\alpha(c)}^{\,4}\widehat\phi(h_0)\overline{\kappa_F}\chi_p(\sigma_p)^{-2}\varpi_{p,j}$ with
   OAI78's $\sigma_p=\lambda^2c/p$, $\epsilon_p=-((\lambda^3c/p)\sigma_p)^{-1}$, $\varpi_{p,j}$ (L1832-1840), local factor
   $B_p(\lambda^4\ell)$ (L1786-1792), and the inactive branch with $1-\N(p)^{-1}$ for $j=0$. The combination
   $\overline{\kappa_F}\vartheta(\lambda^4\ell)=\overline{\kappa(G_{\mathbf h})}\ce(-\delta'\ell/c)\chi_p(\sigma_ph_p)^2
   e(-\epsilon_ph_p^{-1}\lambda^4\ell/p)$ is computed from every active translate and must not depend on $h_p$
   (this is OAI78 eq. reflection-multiplier + eq. reflection-additive-phase, the cross-prime bookkeeping);
4. the marked sum $T_{V,\{p\}}$ (insert $\ind_{p\mid nb^3}$, base twist) against Proposition prop:s1-marked.

Results (`../out/endtoend_{gauss,bessel}_pi{0,1}.txt`):

| run | worst rel. error (ungrouped / grouped, j=0..5) | spread of $\overline{\kappa_F}\vartheta$ over $h_p$ | marked | smallest control error |
|---|---|---|---|---|
| gauss, pi0 = (-2,-3) | 1.7e-12 | 2.4e-12 | 8.1e-13 | 0.47 |
| gauss, pi1 = (1,3) | 1.7e-12 | 2.2e-12 | 8.2e-13 | 0.47 |
| bessel, pi0 | 5.2e-8 | 2.4e-12 | 1.8e-8 | 0.43 |
| bessel, pi1 | 5.2e-8 | 2.2e-12 | 1.8e-8 | 0.43 |

Cusp types covered: $(\sigma,\N(c))\in\{0,+,-\}\times\{1,4,7,28\}$ (the case $(\lambda,c)=1$ of OAI78 L2160-2183;
the cases $\lambda\mid c$ are covered with $\mathcal P=\emptyset$ by `tests_cm/proofs_2026-10-12/b3/t_b3_prop_endtoend.py`).
The bessel runs are limited by the dual truncation (slower decay of $V^\sharp$), not by the identity.

Controls (each must fail, and does, for every $j$ and both tests; relative errors over all runs):
$\overline{\alpha(c)}^{\,2}$ instead of $\overline{\alpha(c)}^{\,4}$ (0.49-2.0); the zeta prefactor $-i/81$ (1.41);
the one-derivative kernel $R_1$ (0.72-4.3); $B_p$ with $\chi_p^{-j}$ instead of $\chi_p^{-j-2}$, i.e. the multiplier
$\chi_p(a)^2$ forgotten (0.47-4.4); $\varpi_{p,j}$ with $\epsilon_p\mapsto3\epsilon_p$ (0.43-2.0; not applicable for
$j=4$, where $\varpi_{p,4}$ does not involve $\epsilon_p$); the marked formula without its active branch (0.69-0.95).

## t_s1_local_phases.py

At the 15 good primes of norm at most 61 (split, and the inert $-5$), for all $j$ and all units $u=c/p$ (first 12):
(a) $C_{p,j}(h)$ (OAI78 eq. reflection-local-fourier) and $|\mathfrak g^\pm_{p,m}|=1$; (b) the local sum
$\sum_{h\ne0}C_{p,j}(h)\chi_p(\sigma_ph)^{-2}e(\epsilon_ph^{-1}x/p)=\chi_p(\sigma_p)^{-2}\varpi_{p,j}B_p(x)$ for all $x$;
(c) the cross-prime identity: $\chi_p(\sigma_p)^{-2}\varpi_{p,j}/\chi_p(u)^{2j+2}$ independent of $u$; (d) the pair
identity $\chi_p(q)^{2j_p+2}\chi_q(p)^{2j_q+2}=(q/p)_3^{j_p+j_q+2}$ over 600 ordered pairs and all $j_p,j_q$; (e) the
active zero-mask transform of $\ind_{q\mid A}$ gives $\chi_q(\sigma_q)^{-2}(-\varpi_{q,0})\N(q)^{-1/2}\chi_q(x)^{-2}$.
Max errors: 1.5e-15, 6.3e-14, 5.5e-15, 3.5e-15, 6.6e-15. Falsifiers (exponent $-j$ in (b), $2j$ in (c), $j_p+j_q$ in
(d)) deviate by at least 1.73. PASS (`../out/local_phases.txt`).

## t_s1_kernel.py

Output: `../out/kernel.txt`. PASS.
(1) one mode with two $\bar z$-derivatives (and with two $z$-derivatives, giving $\alpha(\ell)^2$): relative error
4e-31. (2) the constant of eq. s1-direct (sympy: ratio 1) and the assembly ratio $-\frac1{81}((2\pi)^4/27)^{-t}$: 4e-31;
falsifier $-\frac i{81}(\cdots)$ deviates by 1.41. (3) the Stirling exponent of $R_2$ equals $4\Re t$ at
$\Re t\in\{-\frac14,0,\frac12,2,5\}$ (to $2\cdot10^{-6}$ relative, an $O(1/u)$ effect); pole at $-\frac43$, finite on
$\Re t=-1.3$. (4) $V^\sharp$ is computed by the trapezoid rule (step 0.05 on $|u|\le80$) from the Mellin factor
$\widehat V(-t)R(t)$ tabulated on each line: (4a) the value is the same on the lines $\Re t=0,\frac12,2,-\frac14,-1$
(relative 4e-10 for the Gaussian, 3.5e-8 for $e^{-y-1/y}$; double-precision cancellation on the line $\Re t=2$);
(4b) small $x$: $V^\sharp(x)$ equals the residues at $t=-\frac43,-\frac53$ up to $O(x)$ (relative 4.3e-4 and 2.3e-5 at
$x=10^{-6}$); the local slope tends to $\frac43$, while the one-derivative kernel $R_1$ gives $\frac56$ (falsifier);
(4c) $|(x\partial_x)^jV^\sharp(x)|/\min(x^{1/4},(1+x)^{-A})$ is bounded on $[10^{-8},10^6]$ for $j\le2$: maxima 0.65
($A=2$) and 3.5e6 ($A=4$) for the Gaussian, 57 and 3.3e8 for $e^{-y-1/y}$ (the constants of (e) grow with $A$, as
they should). (5) the growth exponent of the gamma quotient on $\Re w=\sigma_2$ is $2-4\sigma_2$ for
$\sigma_2=-1,-\frac12,-0.1$.

A first version used `mp.quad` with infinite endpoints; it returned a wrong value for $e^{-y-1/y}$ on the line
$\Re t=2$ at $x=10^{-3}$ (0.00401 instead of 0.00410). A finite partition gives the same value on all five lines
(`dbg_bes.py`); the current script uses the tabulated trapezoid rule, which is also much faster.
