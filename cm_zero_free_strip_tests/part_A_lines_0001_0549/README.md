# Part A audit: lines 1–549 of `cm-eleven-twelfths-guide.tex`

This folder tests every equality, inequality, bound, count and numerical value in lines 1–549 of
`../../cm-eleven-twelfths-guide.tex`. That range covers the abstract, "The method in six sentences", Higher weight,
Section 1 (Why start with 11/12) and Sections 2–5: point counts, the triangular lattice, the average, and sixth powers.
The paper was not modified.

```sh
python3 run_all.py      # about 1 minute; one line per claim; outputs also go to outputs/*.txt
```
Exit status (2026-10-12): `cmcore.report` records failures and an atexit hook makes each script exit 1 if a check
failed or a falsifier was not detected; `run_all.py` exits 1 if any script exits nonzero or a claim fails.
`CM_TESTS_FORCE_FAIL=1` makes the first check of each script fail, as a self-test.

The scripts are `t0_intro.py`, `t2_integer.py`, `t3_lattice.py`, `t4_average.py` and `t5_sixth.py`. `cmcore.py` holds the
shared helpers. `zw.py` and `fastgauss.py` are copied unchanged from `tests_guide/`. The tests need Python 3 with numpy and
sympy, and PARI/GP (`gp`) for `ellap`, `ellan`, conductors, minimal models and `lfuncheckfeq`.

**Conventions (the paper's):**
- $n=A+B\omega$ is stored as the pair `(A,B)`.
- "Primary" means $\equiv1\pmod 3$.
- $\chi_\pi(u)$ is the sixth root of unity $\equiv u^{(N\pi-1)/6}\pmod\pi$, or $0$ when $\pi\mid u$.
- $\alpha(n)=n/|n|$, and $u_d=(4d)^5$.

The tests compute symbols in two ways: with exact Euler-criterion arithmetic in $\mathbb Z[\omega]$ (`zw.sym6_prime`),
and with vectorised lookup tables. The two methods were cross-checked on 5700 random cases.

**Exactness.** Every identity involving $\alpha$ is tested exactly in $\mathbb Z[\omega]$ after multiplying by
$\sqrt m$. For example, $\sqrt m\,q_d(m)=\sum_{N n=m}\mu_K(n)\,n\,\chi_n(u_d)$ is compared, as an element of
$\mathbb Z[\omega]$, with the integer $\sqrt m\,q_d(m)$, which is built from PARI's or brute-force $a_p$.

**Status words.**
- **PASS** (E/I): the claim held on the whole range, *and* a deliberately wrong variant (the falsifier, listed in the table)
  failed.
- **PASS** (A/C): an explicit finite instance or a spot check held.
- **SANITY**: sizes are printed for an asymptotic claim, with no pass/fail criterion.
- **UNTESTED**: the claim is not finitely testable; the reason is given in the table.

## Result

| class | claims | PASS | FAIL | SANITY | UNTESTED |
|---|---|---|---|---|---|
| E (exact identity) | 60 | 60 | 0 | – | – |
| I (explicit inequality / exponent arithmetic) | 13 | 13 | 0 | – | – |
| A (asymptotic / analytic) | 16 | 1 (toy instance of eq:extraction) | 0 | 10 | 5 |
| C (cited) | 4 | 3 (spot checks) | 0 | – | 1 |
| **total** | **93** | **77** | **0** | **10** | **6** |

**No genuine errors were found in lines 1–549.** The notes after the table record one small gap in the presentation and three
points of convention.

## Claim → test table

The range column gives the tested range. "fals." names the falsifier: the wrong variant that the test confirms fails.

| ID | line | label | claim | class | test / falsifier | range | result |
|---|---|---|---|---|---|---|---|
| A0.01 | 33–35 | | $E_d$ are exactly the $j=0$ curves; CM by $\mathbb Z[\omega]$ | E | PARI `E.j`=0; sympy: $j=0\iff a=0$; $(x,y)\mapsto(\omega x,y)$ is an automorphism. fals.: $y^2=x^3+x+1$ has $j\ne0$ | $0<\lvert d\rvert\le200$ | PASS |
| A0.02 | 38–39, 55–57 | | $L(E_d,s)\ne0$ for $\Re s>17/12$, for $\Re s>11/8$, and $L(s,\bar\alpha\eta)\ne0$ for $\Re s>7/8$ | A | main (conditional) theorems; not finitely testable | – | UNTESTED |
| A0.03 | 41–42, 172–174 | | critical strip $1/2<\Re s<3/2$; critical line 1 (arithmetic) or 1/2 (shifted) | E | PARI: weight 2, `lfuncheckfeq` ≈ −140 bits. fals.: weight-1 functional equation gives −3 bits | $d=1,7$ | PASS |
| A0.04 | 43–45, 104, 149 | | $11/12+1/2=17/12$, $7/8+1/2=11/8$ | I | Fraction. fals.: $11/12+1$ | – | PASS |
| A0.05 | 45–46 | | the Hecke character has infinite order | E | $\pi\not\sim\bar\pi$, so $\alpha(\pi)$ is not a root of unity. fals.: inert $-p$ does have $\alpha=-1$ | 1124 split $p<2\cdot10^4$ | PASS |
| A0.06 | 64–65, 120–130 | | $\eta(3z)^8$ (weight 4) and $\eta(2z)^3\eta(6z)^3$ (weight 3) are CM by $\mathbb Q(\sqrt{-3})$, character $\bar\alpha^k\eta$ | E | q-expansions: $a_p=0$ at inert $p$; $a_p=\mathrm{Tr}(\zeta^j\bar\pi^k)$; Ramanujan bound. fals.: $\eta(z)^2\eta(11z)^2$ | $q^{3000}$ | PASS |
| A0.07 | 72–74 | | $1/L_p(E_d,s+\frac12)$ has coefficients $1,-\lambda,1,0,0,\dots$ | E | sympy series. fals.: $+X^2\to-X^2$ | – | PASS |
| A0.08 | 76–77, 225 | | $c(m)$ are the coefficients of $1/(\zeta L(\chi_{-3}))$ | E | Dirichlet convolution equals $\delta$. fals.: drop $\chi_{-3}$ | $m\le2\cdot10^4$ | PASS |
| A0.09 | 82, 95 | | zeta sum: $\sum_{Nn=m}\mu_K(n)=c(m)$ | E | exact. fals.: $\mu(m)$ | $m\le2\cdot10^4$, $(m,6)=1$ | PASS |
| A0.10 | 87–88 | | $D^{2+\vartheta+\varepsilon}=D^{1+\varepsilon}H$ | I | sympy. fals.: $H^2$ | – | PASS |
| A0.11 | 86–88, 378 | | about $H$ rows, about $D$ columns | A | see A4.01 | – | SANITY |
| A0.12 | 91–94 | | $H^{1/6}\lvert S\rvert^2\lesssim DH$ gives $S\lesssim D^{11/12}$ at $H=D$ | I | Fraction. fals.: $H^{1/3}$ | – | PASS |
| A0.13 | 91–93 | | $\approx H^{1/6}/\log H$ rows; error $O(D/N\pi)$ | A | see A5.03, A5.05–A5.07 | – | SANITY |
| A0.14 | 96–104 | | Poisson/theta/descent summary | A | forward references to §§6–8 (other chunks). Text check only: matches $a_\xi=\bar\alpha^2\gamma_2\xi$ (l.661–674) | – | UNTESTED |
| A0.15 | 57, 122–124 | | weight $k+1$ ↔ $\bar\alpha^k\eta$, while Part I uses $\psi_d=\alpha\chi$ | E | $\sum n\chi_n(u_d)=\sum\bar n\,\overline{\chi_n(u_d)}=a_m$. fals.: $\alpha\bar\chi$ | $d=7$, $m\le3000$ | PASS (convention) |
| A0.16 | 66–67, 129–130 | | Sym$^k$ factors into Hecke $L$-functions | C | sympy local factors at split and inert $p$. fals.: wrong inert factor | $k\le6$ | PASS (spot); zero-free consequence UNTESTED |
| A2.01 | 158–162 | | $a_p=p+1-\#E_d(\mathbb F_p)=-\sum_x\left(\frac{x^3+d}{p}\right)$ | E | brute-force count = Legendre sum = PARI `ellap`. fals.: sign | 1065 pairs $(d,p)$, $p<600$ | PASS |
| A2.02 | 165–166 | | $\#\{y\}=1+\left(\frac{c}{p}\right)$ | E | exhaustive. fals.: $1-(\cdot)$ | $p<400$ | PASS |
| A2.03 | 168–171 | | $L_p(E_d,s+\frac12)^{-1}=1-\lambda p^{-s}+p^{-2s}$ | E | sympy. fals.: unshifted $p^{1-2s}$ | – | PASS |
| A2.04 | 172–174 | | critical lines $\Re s=1/2$ and $\Re s=1$ | E | see A0.03 | – | PASS |
| A2.05 | 177–183 | eq:qd | definition of $q_d$; $q_d\ne\mu$ | E | $q_1(25)=1\ne\mu(25)=0$ (the full identity is A2.07) | – | PASS |
| A2.06 | 185–193 | eq:integer-convolution | $q_d=\sum_{rt^2=m,(r,t)=1}\mu(r)\mu(t)^2\Lambda_d(r)$ | E | 32382 values, max err $3\cdot10^{-15}$. fals.: drop $(r,t)=1$; drop $\mu(t)^2$ | $m\le2\cdot10^4$, 5 values of $d$ | PASS |
| A2.07 | 198–208 | eq:reciprocal | $\sum q_d m^{-s}=\prod_{p\nmid6d}(\cdots)=B_d/L(E_d,s+\frac12)$ | E | exact integer convolution with PARI `ellan`, including non-minimal $d$ such as $5\cdot2^6$ and $7\cdot5^6$. fals.: $Q(p^2)=-p$ | $m\le6000$, 10 values of $d$ | PASS |
| A2.08 | 209–211 | | $B_d$ holomorphic and nonzero for $\Re s>0$, also for non-minimal models | E | every pole has $\lvert p^{-s}\rvert\ge1$. fals.: unshifted factor | $\lvert d\rvert\le300$ plus non-minimal $d$ | PASS |
| A2.09 | 218–220 | eq:target | $S_d(D)\ll D^{11/12+\varepsilon}$ | A | $\lvert S\rvert/D^{1/2}\le0.46$ | $D\le10^5$ | SANITY |
| A2.10 | 226–228 | | $\chi_{-3}$ is the Kronecker symbol $(-3/\cdot)$ | E | PARI `kronecker`. fals.: swapped signs | $t\le3000$ | PASS |
| A2.11 | 230–237 | table | $c(p),c(p^2),q_d(p),q_d(p^2)$, including $a_p=0$ at inert $p$ | E | exact. fals.: $q_d(p^2)=-1$ | $p<3000$, 6 values of $d$ | PASS |
| A2.12 | 238–240 | | $a_5(1)=0$, $a_7(1)=-4$, $q_1(5)=0$, $q_1(25)=1$, $q_1(7)=4/\sqrt7$ | E | PARI. fals.: $-4/\sqrt7$ | – | PASS |
| A3.01 | 248–250 | | $N(n)=\lvert n\rvert^2=A^2-AB+B^2$ | E | fals.: $A^2+AB+B^2$ | $\lvert A\rvert,\lvert B\rvert\le60$ | PASS |
| A3.02 | 252–253 | | six units $\pm1,\pm\omega,\pm\omega^2$; UFD | E (UFD cited) | norm-1 enumeration; 20000 Euclidean divisions. fals.: $A^2+B^2=1$ count | – | PASS |
| A3.03 | 253–258 | | unique primary associate; primary $\iff A\equiv1,B\equiv0\pmod3$ | E | fals.: "$\equiv\pm1$" | $\lvert A\rvert,\lvert B\rvert\le40$ | PASS |
| A3.04 | 260–261 | | coprime to 2 $\iff$ $A,B$ not both even | E | fals.: "$A$ even" | $\lvert A\rvert,\lvert B\rvert\le50$ | PASS |
| A3.05 | 264–271 | | $\mu_K(-5)=-1$; $\mu_K(7)=1$; $7=(1+3\omega)(-2-3\omega)$; $\mu_K\ne\mu\circ N$ | E | exact. fals.: $(1+3\omega)(2+3\omega)\ne7$ | – | PASS |
| A3.06 | 273–278 | | Euler criterion gives exactly one root of unity; $\chi$ multiplicative in $u$; $\chi_1=1$ | E | 4000 random cases. fals.: additivity | – | PASS |
| A3.07 | 287–290, 300–303 | eq:ud | $\overline{\chi_n(4d)}=\chi_n(4d)^5=\chi_n((4d)^5)$, zeros kept | E | index arithmetic. fals.: $(4d)^2$ | 12000 pairs $(d,n)$ | PASS |
| A3.08 | 291, 339–343 | eq:normgroup | $q_d(m)=\sum_{Nn=m}\mu_K\alpha\chi_n(u_d)$; $S_d=A_{u_d}$ | E | exact in $\mathbb Z[\omega]$; sums agree to $10^{-14}$. fals.: $u=4d$; $u=d^5$ | $m\le2\cdot10^4$, 8 values of $d$ | PASS |
| A3.09 | 304–306 | | $\chi_n(uv^6)=\chi_n(u)\mathbf1_{(n,v)=1}$ | E | fals.: $v^3$ | 3000 random cases | PASS |
| A3.10 | 306–308 | | factors at 2 and 3 do not matter; $(4d)^5=16d^5\cdot2^6$ | E | fals.: row $4d^5$ | 3000 columns | PASS |
| A3.11 | 309–310 | | $d=16$: $2^{24}=(2^4)^6$, conductor 27, row $u=1$; $E_1$ is row 16 | E | PARI; exact. fals.: row 1 vs $E_1$ | $m\le2\cdot10^4$ | PASS |
| A3.12 | 311–312 | | every rational row $u$ is $E_{4^5u^5}$ | E | PARI `ellan`. fals.: $d=4u^5$ | $0<\lvert u\rvert\le30$, $m\le3000$ | PASS |
| A3.13 | 312–315 | | rows that are not a rational times a sixth power have non-real coefficients in general | E | examples found at $m=7,13$; rows $r v^6$ are real away from primes of $v$. fals.: detector stays silent on a real row | $m\le3000$ | PASS |
| A3.14 | 295–297 | | disc $=1-4a$; $J(\kappa,\rho)=\kappa(4)J(\kappa,\kappa)$ | E | 414 cases. fals.: $\kappa(2)$, $\bar\kappa(4)$ | $p<3000$ | PASS |
| A3.15 | 318–319 | | $p=\pi\bar\pi$, both primary | E | fals.: some element has norm 5 | $p<10^5$ | PASS |
| A3.16 | 320–322 | eq:point-identity | $a_p=\pi\overline{\chi_\pi(4d)}+\bar\pi\chi_\pi(4d)$ | E | exact for both $\pi$. fals.: $\chi$ unconjugated; $\chi_\pi(d)$; $-\pi$ | 27528 pairs $(d,\pi)$, $\lvert d\rvert\le60$, $p<1500$ | PASS |
| A3.17 | 324–325 | IR Ch.18 | the IR convention ($\pi\equiv2$) changes the sign of $\pi$ | C | cited form holds with $\pi'=-\pi$. fals.: same form with $\pi$ | 5 values of $d$, 11 primes | PASS |
| A3.18 | 326 | | "checked for 480 pairs" | C | statement about a past check; superseded by A3.16 | – | UNTESTED |
| A3.19 | 327 | | $\chi_{\bar\pi}(x)=\overline{\chi_\pi(x)}$ for rational $x$ | E | fals.: non-rational $x$ | $p<600$ | PASS |
| A3.20 | 328–332 | | $(1-\psi\pi X)(1-\psi\bar\pi X)=1-\lambda X+X^2$ | E | fals.: unconjugated $\chi$ | $p<3000$, 5 values of $d$ | PASS |
| A3.21 | 333–337 | | inert $p$: $-p$ is primary of norm $p^2$; $\chi_{-p}(x)=1$; $\alpha(-p)=-1$; factor $1+p^{-2s}$ | E | fals.: non-rational $x$ | $p<400$ | PASS |
| A3.22 | 336–337 | | $p\mid d$ makes the symbols vanish | E | fals.: $\chi_{(1+3\omega)}(5)\ne0$ | – | PASS |
| A3.23 | 344–345 | | the sum is empty when $(m,6)>1$ | E | fals.: drop "not both even" | $N\le10^5$ | PASS |
| A3.24 | 348–361 | | $a_m(E_d)/\sqrt m=\sum_{Nn=m}\alpha(n)(u_d/n)_6$ (the paper checked $d\in\{1,2,-3,7\}$, $m\le400$) | E | exact against PARI. fals.: $u=4d$ | 51264 coefficients, $m\le2\cdot10^4$, 8 values of $d$ | PASS |
| A3.25 | 358–359 | | Frobenius eigenvalues are $\pi(u_d/\pi)_6$ and its conjugate | E | trace $a_p$, norm $p$. fals.: $u=4d$ | $p<5000$ | PASS |
| A3.26 | 362–369 | | Kummer/Artin interpretation | C | spot check: $(u/\pi)_6=1\iff x^6=u$ has 6 roots mod $p$ | $p<800$ | PASS (spot) |
| A4.01 | 378 | | about $H$ rows and about $D$ columns | A | rows/$H\to2\pi/\sqrt3$; columns/$X\to0.3023$ | $H\le10^6$ | SANITY |
| A4.02 | 379–402 | eq:coordinates, eq:double | coordinate form of $M$ = double sum | E | literal per-element evaluation vs vectorised $S_H$, relative error $10^{-15}$. fals.: $\alpha(n_1)$ unconjugated; $\chi_{n_1}$ unconjugated | 3 pairs $(D,H)$ | PASS |
| A4.03 | 403–405 | | $\bar\alpha_1\alpha_2=\bar n_1n_2/(\lvert n_1\rvert\lvert n_2\rvert)$ | E | fals.: conjugate the other factor | 2000 pairs | PASS |
| A4.04 | 408–411 | | $A_u(D)=\sum_m q_u(m)W(m/D)$ | E | fals.: drop $\mu$ | 5 rows, $D\le3000$ | PASS |
| A4.05 | 412–425 | eq:allrows | $q_u$ multiplicative; zero at $2,3$; local table; $q_u(p^{j\ge3})=0$; reduces to eq:qd | E | exact in $\mathbb Z[\omega]$. fals.: sign at split $p$; sign at inert $p^2$; non-coprime product | 30 rows, $m\le12000$ | PASS |
| A4.06 | 426–427 | | the symbols distinguish $\pi$ from $\bar\pi$ | E | examples. fals.: detector stays silent on rational rows | $p=7$ | PASS |
| A4.07 | 430–437 | eq:MS | $M(D,D^{1+\vartheta})\ll D^{2+\vartheta+\varepsilon}$ | A | $M/D^{2.1}\approx0.27$–$0.29$ | $D\le1600$, $\vartheta=0.1$ | SANITY |
| A4.07b | 437–439 | | constant depends on derivatives of $W$; independent of $d$ | A | structural remark | – | UNTESTED |
| A4.08 | 441–443 | | $S_H(n,n)=\#\{u:(u,n)=1\}\ll H$ | E | fals.: $\sum\chi$ instead of $\sum\lvert\chi\rvert^2$ | 150 columns, $H=2\cdot10^4$ | PASS |
| A4.09 | 445–451 | | $S_H(n,n)=\frac{2\pi}{\sqrt3}H\prod(1-N\pi^{-1})+O(\tau(n)(\sqrt H+1))$ | A | error ratio $\le0.35$ (93 without the Euler product) | $H\le4\cdot10^5$ | SANITY |
| A4.10 | 451–452 | | $2\pi/\sqrt3$ = disc area / covolume | I | sympy (covolume $\sqrt3/2$). fals.: covolume 1 | – | PASS |
| A4.11 | 452–453 | | diagonal is $O_W(HD)$ | A | diag/$(HD)\approx0.27$ | $D\le1600$ | SANITY |
| A4.12 | 453–455 | | what (MS) bounds | A | remark | – | UNTESTED |
| A5.01 | 459 | | $N_0=\lvert4d\rvert^{10}$ | E | fals.: $\lvert4d\rvert^5$ | $\lvert d\rvert\le200$ | PASS |
| A5.02 | 460–464 | | $\chi_n(u_0v^6)=\chi_n(u_0)\mathbf1_{(n,v)=1}$ | E | fals.: $v^2$ | 3000 random cases | PASS |
| A5.03 | 472–473 | | $J\asymp Y/\log Y$ | A | $J/(Y/\log Y)\approx0.50$–$0.52$; $J=0$ for $Y<7$ | $Y\le10^6$ | SANITY |
| A5.04 | 473 | | rows $u_0\pi^6$ distinct, of norm $\le H$ | E | fals.: associates repeat rows | $Y=3000$ | PASS |
| A5.05 | 474–478 | | $A_{u_0}-A_{u_0\pi^6}=\sum_{\pi\mid n}\cdots$ | E | fals.: sum over $\pi\nmid n$ | $D=3000$ | PASS |
| A5.06 | 479–481 | | difference $\le\max\lvert W\rvert\cdot\#$multiples; coefficients $\le1$ | I | fals.: count multiples of $\pi^2$ | $D\le8000$ | PASS |
| A5.07 | 479–484 | eq:approx | $\#$multiples $=O(D/N\pi)$; $\lvert A_{u_0}-A_{u_0\pi^6}\rvert\ll D/Y$ | A | ratios 0.60 and 0.36 | $D\le8000$, $Y\le200$ | SANITY |
| A5.08 | 485 | | $\lvert z\rvert^2\le2\lvert w\rvert^2+2\lvert z-w\rvert^2$ | I | sympy: the difference is $\lvert z-2w\rvert^2$; $10^6$ random pairs. fals.: constant 1 | – | PASS |
| A5.09 | 485–492 | eq:extraction | $J\lvert S\rvert^2\le2M+O(JD^2/Y^2)$ | A | toy cases with the actual $M$ ($H$ up to $10^6$ rows); explicit chain holds. fals.: dropping the error term fails | 3 toy cases | PASS |
| A5.10 | 494–502 | eq:exponents | $DH^{5/6}=D^{11/6+5\vartheta/6}$; $D^2H^{-1/3}=D^{5/3-\vartheta/3}$ | I | sympy, also re-derived from $2M/J$ and $D^2/Y^2$. fals.: $\vartheta/6$ | – | PASS |
| A5.11 | 503–507 | eq:Sbound | $S_d\ll D^{11/12+\eta}$ | I | exact Fractions over 399 values of $\eta$. fals.: $5/6+\eta$ | – | PASS |
| A5.12 | 508–513 | | $N_0^{1/12}D^{11/12+5\vartheta/12+\varepsilon}+N_0^{1/6}D^{5/6-\vartheta/6}$; $N_0^{1/12}=\lvert4d\rvert^{5/6}$ | I | sympy. fals.: $N_0^{1/6}$ in the first term | – | PASS |
| A5.12b | 514 | remark | "valid once $Y\ge2$" | I | $J=0$ for $Y<7$, but the bound then holds trivially (see note 1). fals.: "$J\ge1$ for $Y\ge2$" | – | PASS (note) |
| A5.13 | 514 | | $D^{1+\vartheta}\ge2^6N_0\iff Y\ge2$ | I | fals.: $2^3$ | integer grid | PASS |
| A5.14 | 514–515 | | nontrivial only for $D\gtrsim\lvert d\rvert^{10}$ | I | requires $D^{1-5\vartheta}>\lvert4d\rvert^{10}$. fals.: $1+5\vartheta$ | – | PASS |
| A5.15 | 515–516 | | $E_{dt^6}\cong E_d$ | E | substitution; PARI minimal models. fals.: $E_{dt^2}$ | 4 values of $d$, 3 of $t$ | PASS |
| A5.16 | 516–519 | | $f_p=2$ for $p\ge5$, $p\mid d$; $\lvert d\rvert\ll N^{5/2}$; $\lvert d\rvert^{10}\ll N^{25}$ | E | 8708 primes; $\lvert d\rvert\le6^5N^{5/2}$ (max ratio 0.11). fals.: $f_p=1$ | 5979 values of $d$ | PASS |
| A5.16b | 518 | | sharpness of the exponent $5/2$ | A | $\lvert d\rvert/N^2$ grows along $d=p^5$ | – | SANITY |
| A5.18 | 522–524 | | $D^{2-1/6}=D^{11/6}$ | I | fals.: $2-1/3$ | – | PASS |
| A5.19 | 528–531 | | $\chi_n(v^3)=(v/n)_2$, $\chi_n(v^2)=(v/n)_3$ | E | Euler criterion. fals.: $(v/n)_3$ for $v^3$ | 2000 random cases | PASS |
| A5.20 | 532–533 | | $(v/n)_2=(v/N n)$ (Jacobi) for rational $v$ | E | fals.: $-v$ | $\lvert v\rvert\le60$, $N\le5000$ | PASS |
| A5.21 | 534–537 | eq:quadratic-twists | $A_{u_dv^3}=\sum q_d(m)(v/m)W$ | E | exact coefficients. fals.: $v^2$ | 135000 coefficients | PASS |
| A5.22 | 538–540 | | $\mathbb F_p^*$ consists of squares in $\mathbb F_{p^2}$ | E | fals.: $x+\omega$ | $p<500$ | PASS |
| A5.23 | 542–543 | | rational slice has $2\lfloor\sqrt H\rfloor\asymp H^{1/2}$ rows | E | fals.: $\lfloor\sqrt H\rfloor$ | – | PASS |
| A5.24 | 544–545 | | sixth powers number $\asymp H^{1/6}$ (lattice) and $\asymp H^{1/12}$ (rational) | E | exact counts $=\#v/6$. fals.: $\#v/3$ | $H^{1/6}\le800$ | PASS |
| A5.25 | 546–548 | | remarks on lattice primes and rational numerators | A | not numerical | – | UNTESTED |

## Notes: no errors, one gap in presentation, three points of convention

1. **Lines 513–514, "valid once $D^{1+\vartheta}\ge2^6N_0$, so that $Y\ge2$"** (gap in presentation, not an error).
   There is no primary prime coprime to 6 of norm below 7, so $J=|\mathcal P_Y|=0$ for $2\le Y<7$. Also, $J\gg Y/\log Y$
   holds uniformly only for $Y$ large. In that range the "steps above", which divide by $J$, give nothing.

   The stated bound is still true there. When $Y<7$, $N_0>H/7^6$, so the first term is at least
   $D^{1+\vartheta/2+\varepsilon}/\sqrt7$. This exceeds the trivial bound $|S_d(D)|\ll_W D\log D$, which follows from
   $|q_d(m)|\le2^{\omega(m)}$ (Hasse; checked for $m\le2\cdot10^4$). The same holds when the primes of $d$ deplete
   $\mathcal P_Y$, since in the nontrivial range $D^{1-5\vartheta}>|4d|^{10}$ they number $O(1/\vartheta)$.

   A one-line remark ("for smaller $Y$ the bound is trivial") would close this.
2. **Lines 57 and 124 vs line 288** (convention). Part I writes $\psi_d=\alpha\,\chi_\cdot(u_d)$, while Part II and the
   Higher weight paragraph write $\bar\alpha^k\eta$. For $E_d$ both give the same $L$-function, because the coefficients
   are real (A0.15). In general the two conventions differ by complex conjugation. The test for the eta products (A0.06)
   cannot distinguish them.
3. **Line 82** (convention). The zeta sum $\sum_n\mu_K(n)W(N n/D)$ over primary $n$ coprime to 6 equals
   $\sum_m c(m)W(m/D)$ only after the Euler factors at 2 and 3 are removed, e.g. $c(3)=-1$ and $c(4)=-1$. This agrees
   with "up to finitely many Euler factors" (l.72).
4. **Lines 312–313** (convention). Rows $r v^6$ with $r$ rational and $v\notin\mathbb Z$ have non-real coefficients at the
   primes of $v$ (391 instances). This agrees with l.304–306: only finitely many Euler factors change.

## Not tested, and why

- The main theorems (17/12, 11/8, 7/8) and the mean-square estimate (MS) are conditional or analytic. Only sanity sizes
  were computed for (MS): $M/D^{2.1}\approx0.28$.
- The Poisson/theta/descent summary (l.96–104) refers forward to Sections 6–8, which other chunks cover.
- The 480-pair check claimed on l.326 was not re-run as such; A3.16 covers about 13700 pairs.
- The Artin/Kummer statement and unique factorisation are cited; only spot checks were made.
- The remarks on l.437–439, 453–455 and 546–548 are not numerical claims.
