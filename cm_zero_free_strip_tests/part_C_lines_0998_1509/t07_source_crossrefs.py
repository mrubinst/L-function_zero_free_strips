"""t07: restatements of OAI78 (build/paper.tex) and G78 (seven-eighths-guide.tex) in guide lines 998-1509, checked against
the source text.  Each check searches for the exact TeX of the source formula in a stated line window.  A restatement
whose source text is found and agrees (up to the CM change the guide declares) is reported CITED; otherwise CITED-X.
Also checks every line number quoted by the guide (L2331, L3871, L4056, L2003-2279, the table of ranges L1470-1480).
The falsifier for this file is a deliberately wrong needle per block, which must not be found.
Run time: < 1 s.
"""
import os, sys, re
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from cmreport import claim, cited, info, untested, src_has, read_lines, OAI, G78

OL = read_lines(OAI)
def has(needle, lo=1, hi=None, path=OAI):
    return src_has(path, needle, lo, hi)

def chk(cid, gline, label, desc, needles, extra_ok=True, note=""):
    found = []; ok = extra_ok
    for nd, lo, hi in needles:
        f, at = has(nd, lo, hi)
        ok &= f
        found.append(f"L{at}" if f else f"MISSING[{nd[:40]}]")
    return cited(cid, gline, label, desc, ok, "OAI78 " + ", ".join(found) + (("; " + note) if note else ""))

# falsifier: a wrong needle must not be found
wrong = [("\\zeta=-\\frac1{81}\\bar\\alpha(c)^4", 1, None), ("r+2z\\le m-c_2", 1, None), ("\\frac{49}{440641}", 1, None),
         ("C(s)=s-\\frac{11}{17}", 1, None), ("Q^{3/5}(3+|t|)^3", 1, None)]
fw = [has(n, lo, hi)[0] for n, lo, hi in wrong]
claim("C62", 1002, "source-check falsifier", "C", "deliberately wrong needles (CM scalar in OAI78, wrong constants) are absent",
      not any(fw), f"found: {fw}", falsifier_failed=not any(fw))

n = len(OL) - (1 if OL[-1] == "" else 0)
cited("C63", 1040, "16,700 lines", "OAI78 source has about 16,700 lines", abs(n - 16700) < 100, f"{n} lines")
chk("C64", 1075, "thm:main", "finite-order Hecke L-functions over Q(sqrt-3) have no zero in Re s>7/8",
    [("\\begin{theorem}\\label{thm:main}", 100, 110), ("Every finite-order Hecke $L$-function over", 100, 112),
     ("has no zero in $\\Re s>7/8$", 100, 112)])
chk("C65", 1082, "continuation criterion", "Prop. lem:continuation-criterion; reciprocal Euler product bounded on Re s>=2; "
    "sigma_0=7/8, C(s)=s-11/16", [("\\begin{proposition}[Continuation from a common signal]", 395, 402),
    ("The reciprocal Euler product is absolutely and uniformly", 440, 460),
    ("\\sigma_0=\\frac78,\\qquad C(s)=s-\\frac{11}{16},", 510, 516)])
chk("C66", 1107, "probe weight", "eta enters the probe only through gamma_2(c) abar(A) Psi(A), Psi(a) = eta(a) x fixed phases",
    [("\\gamma_2(c)\\overline{\\alpha(cn^3)}\\Psi_{m,s,\\eta}(cn^3)", 3380, 3390),
     ("\\eta(a)\\Xi(a)^{-1}\\mathcal R(a,s)\\chi_a(m)", 3355, 3366), ("\\label{eq:general-probe}", 3390, 3400)])
chk("C67", 1129, "L2331", "OAI78 L2331: the zbar-derivative removes the constant term C_H(v)",
    [("the derivative removes $C_H(v)$ from every cusp expansion", 2331, 2331)])
untested("C68", 1126, "cusp constant terms", "C",
         "theta has a nonzero constant term at the identity cusp and none at the other two cusps",
         "needs the cusp expansions of Dunn-Radziwill App. A (F_1, F_19, F_10); OAI78 L1936-1940 only defines C_H(v)")
chk("C69", 1140, "rows", "row u = zeta^S(6z) L^S(w,chi(u)) / L^S(x, eta conj chi(u)) H_{eta,u} (eta -> psi in the guide)",
    [("\\frac{\\zeta_F^S(6z)L^S(w,\\chi_\\bullet(u))}", 4030, 4055), ("{L^S(x,\\eta\\overline{\\chi_\\bullet(u)})}", 4030, 4055),
     ("\\label{eq:probe-poisson-identity}", 3831, 3831)])
chk("C70", 1148, "L4056, H = 1 + O(P_0^-4/5)", "second region (L4056); H << q_u^eps; H_{eta,1} = 1+O(P_0^{-c_H}), c_H = 4/5",
    [("x_r\\ge7/8,\\quad z_r\\ge33/200,\\quad w_r\\ge19/20.", 4056, 4056), ("\\mathcal H_{\\eta,u}\\ll_\\epsilon q_u^\\epsilon", 4058, 4060),
     ("$c_H=4/5$", 4060, 4064)])
chk("C71", 1154, "L3871", "OAI78 L3871: no ray-class restriction on p", [("No ray-class restriction on $p$ is imposed.", 3871, 3871)])
chk("C72", 1160, "eq:signal", "gamma1 gamma2 = mu alpha G for squarefree c (rearranged in eq:signal-cm)",
    [("\\gamma_1(c)\\gamma_2(c)=\\mu(c)\\alpha(c)G(c)", 965, 972), ("For squarefree $c$ the complete Gauss identities are", 960, 966)])
chk("C73", 1168, "a_p, b_p", "a_p = abar^3 eta^3, b_p = abar^2 eta^2 conj chi_p(4), b_p^3 = a_p^2; J_0 = -eta(p) rho^-1 Q^-x + ...",
    [("a_p=\\overline{\\alpha(p)}^3\\eta(p)^3", 3860, 3870), ("b_p=\\overline{\\alpha(p)}^2\\eta(p)^2\\overline{\\chi_p(4)}", 3860, 3870),
     ("b_p^3=a_p^2", 3860, 3870), ("0&-\\eta(p)\\rho^{-1}Q^{-x}+W_{\\rm loc}R", 4018, 4026)])
chk("C74", 1182, "c_S", "principal data: H_eta = H_{eta,1}(s,1,1/6), c_S = W1^(1) M(1/6) (Res zeta^S)^2/6 (positive iff M(1/6)>0)",
    [("H_\\eta(s)=\\mathcal H_{\\eta,1}(s,1,1/6)", 5518, 5522), ("c_S=\\widehat W_1(1)M(1/6)", 5518, 5522)])
chk("C75", 1192, "B_m = C_V with abar -> abar^2", "OAI78's completed row C_V(X;m,nu)",
    [("\\frac{\\gamma_2(c)\\overline{\\alpha(cn^3)}\\nu(cn^3)\\chi_{cn^3}(m)}", 2818, 2826)])
chk("C76", 1202, "lem:unmarked-completed-moment", "bound Z^{O/2+max(M-O,2M-O-N)+eps} (=> Z^eps(Q+Q^2/Z), see C48)",
    [("\\label{lem:unmarked-completed-moment}", 3139, 3139), ("\\ll Z^{O/2+\\max\\{M-O,\\,2M-O-N\\}+\\epsilon}", 3155, 3162)])
chk("C77", 1207, "marked lemmas", "lem:hybrid-energy, lem:reflected-energy, lem:probe-row-norm, Cor. cor:marked-reflection exist",
    [("\\label{lem:hybrid-energy}", 7510, 7510), ("\\label{lem:reflected-energy}", 7747, 7747),
     ("\\label{lem:probe-row-norm}", 8111, 8111), ("\\begin{corollary}[Completed reflection with whole-index marks]", 7360, 7362)])
chk("C78", 1221, "table, OAI78 column", "dbar theta; Gamma(s+1/3)Gamma(s+2/3); R(t)=prod Gamma(1+t+-1/6)/Gamma(1-t+-1/6); "
    "-i/81 abar(c)^2; abar(n)abar(b)^3; alpha(mu)",
    [("Since $\\partial_{\\bar z}\\check e(\\mu z)=2\\pi i\\bar\\mu\\check e(\\mu z)$", 2290, 2300),
     ("J(s)=\\frac i4(2\\pi)^{-2s}\\Gamma(s+1/3)\\Gamma(s+2/3)", 2295, 2302),
     ("R(t)=\\prod_{\\pm}\\frac{\\Gamma(1+t\\pm1/6)}{\\Gamma(1-t\\pm1/6)}", 1738, 1744),
     ("\\zeta=-\\frac i{81}\\bar\\alpha(c)^2", 1830, 1836), ("\\frac{\\gamma_2(n)\\overline{\\alpha(nb^3)}\\Psi(nb^3)}", 1720, 1730),
     ("\\frac{d(\\mu)\\alpha(\\mu)\\vartheta(\\lambda^4\\mu)}{\\sqrt{q_\\mu}}", 1764, 1772)],
    note="CM column verified independently in t03 (C28-C32)")
chk("C79", 1232, "kernel bounds", "|R(sigma+iu)| << (1+|u|)^{4 sigma} on -1/4<=sigma<=A; first pole -5/6; min(x^{1/4},(1+x)^{-A})",
    [("|R(\\sigma+iu)|\\ll_A(1+|u|)^{4\\sigma}.", 2412, 2416), ("On $-1/4\\le\\sigma\\le A$, Stirling's formula gives", 2410, 2414),
     ("The first negative pole of the numerator of $R(t)$ is $-5/6$.", 2405, 2410),
     ("\\ll_{A,j}\\min\\{x^{1/4},(1+x)^{-A}\\}", 1855, 1862)],
    note="guide writes x^{-A} for (1+x)^{-A} (equivalent for the min)")
chk("C80", 1242, "marks", "chi_p(x)^0 = 1_{p not | x} (zero convention)", [("$\\chi_p(x)^0=1_{p\\nmid x}$", 1670, 1675)])
chk("C81", 1247, "Z^{3/16}", "compensated low estimate Z^{l_x/2+b/12+eps} = Z^{3/16+eps}",
    [("\\ll Z^{l_x/2+b/12+\\epsilon}=Z^{3/16+\\epsilon}.", 8565, 8572), ("C(s)=s-\\frac23,", 505, 511)])
chk("C82", 1266, "lem:hecke-strip-growth", "|L(sigma+it,psi)| << Q^{3/5}(3+|t|)^2 on [-1/10,11/10]; FE with Gamma(s); "
    "|L(-1/10+it)| << Q^{3/5}(3+|t|)^{6/5}",
    [("|L(\\sigma+it,\\psi)|\\ll Q^{3/5}(3+|t|)^2,", 1428, 1432), ("-1/10\\le\\sigma\\le11/10.", 1428, 1432),
     ("\\Lambda(s,\\psi)=(3Q)^{s/2}(2\\pi)^{-s}\\Gamma(s)L(s,\\psi),", 1438, 1444),
     ("|L(-1/10+it,\\psi)|\\ll Q^{3/5}(3+|t|)^{6/5}.", 1444, 1450)],
    note="CM: Gamma(s) -> Gamma(s+1/2), verified in t04 (C37-C40)")
chk("C83", 1286, "lem:logarithmic-control", "logarithmic control lemma, with principal clauses",
    [("\\begin{lemma}[Logarithmic control]\\label{lem:logarithmic-control}", 1531, 1531), ("For the principal", 1531, 1540)])
chk("C84", 1292, "X_u and Theta", "Theta = <eta, T^>, X_u = {nu chi(u)^sigma}, Q_psi << q_u; slot characters lie in X_u",
    [("\\Theta=\\langle\\eta,\\widehat T\\rangle", 4220, 4226), ("\\bigl\\{\\psi_{u,\\nu,\\varsigma}(n)=\\nu(n)\\chi_n(u)^\\varsigma:", 4230, 4236),
     ("Q_\\psi\\ll_{\\mathcal A}q_u,", 4246, 4250), ("Every resulting prime character belongs to \\(\\mathcal X_u\\).", 15029, 15029)])
chk("C85", 1338, "lem:bin-prime-bound", "|Q_i(u;z)| << U^eps P_i^{a-1/2+O(e)}",
    [("\\label{lem:bin-prime-bound}", 15015, 15015), ("U^{\\epsilon_1}P_i^{a-1/2+O(e)}.", 15018, 15024)])
chk("C86", 1340, "old-eq:3.5", "uniform slot bound Z^{kappa z+eps} under beta_* <= (1+kappa)/2, used only for slot characters "
    "(only zero-free use in lem:plain: L12833-12898)",
    [("\\le C Z^{\\kappa z+\\epsilon_1}p_j(\\boldsymbol W)^b", 12870, 12878), ("For \\(\\kappa<1\\), put \\(s_\\kappa=(1+\\kappa)/2\\).", 12833, 12833),
     ("For \\(\\kappa<1\\) the hypothesis remains", 12895, 12898)])
chk("C87", 1350, "lem:sextic-large-sieve", "arbitrary coefficients c_n independent of the row; bound (KD)^eps{K+D+(KD)^{2/3}}",
    [("any sequence of complex numbers indexed by the squarefree ideals outside $S$.", 4708, 4712),
     ("\\bigl\\{K+D+(KD)^{2/3}\\bigr\\}", 4712, 4720)])
chk("C88", 1354, "lem:marked", "conditions r+2z<=m-c_1, 2r+8z<=3m-c_2; M_u, Q_u with psi_u = nu chi(u)^eps on n and p; "
    "bound Z^{m+eps}; either sign",
    [("r+2z\\le m-c_1,\\qquad 2r+8z\\le3m-c_2.", 9254, 9256), ("M_u(Z^r;W)=Z^{-r/2}\\sum_n\\mu(n)\\psi_u(n)W(q_n/Z^r).", 9228, 9234),
     ("\\psi_u(n)=\\nu(n)\\chi_n(u)^{\\varepsilon_\\chi}.", 9224, 9228), ("a_i(p)\\psi_u(p)W_i(q_p/Z^{z_i}),", 9236, 9242),
     ("|M_u(Z^r;W)Q_u|^2\\ll Z^{m+\\epsilon}.", 9257, 9261), ("for either common sign", 9262, 9266)])
chk("C89", 1387, "lem:inverse-amplification", "input is lem:marked with no prime slots",
    [("Lemma~\\ref{lem:marked}, with no prime slots.", 12363, 12366)])
chk("C90", 1393, "lem:plain", "3/4<=kappa<=1; z=0 no restriction; n1+n2+6 kappa z<=M; beta_*<=(1+kappa)/2; M=m+q; uniform in moving "
    "moduli, masks, frozen labels",
    [("Let \\(3/4\\le\\kappa\\le1\\).", 12532, 12533), ("n_1+n_2+6\\kappa z=A+(6\\kappa-1)z\\le M.", 12560, 12563),
     ("\\(\\beta_*\\le(1+\\kappa)/2\\).", 12564, 12566), ("effective width \\(M=m+q\\)", 12592, 12592),
     ("moving moduli, admissible masks, and frozen outer labels", 12574, 12579)])
chk("C91", 1431, "eq:common-high-exponent", "E(d) = C_0 + 2delta/3 + q/6 - h(1-R) + (d-h)(R+delta/2-17/50), C_0 = -1/48",
    [("&=C_0+\\frac23\\delta+\\frac q6-h(1-R)", 15728, 15733), ("C_0=-\\frac1{48}.", 15729, 15734)])
chk("C92", 1444, "lem:balanced-endpoint", "-E_* >= 49/440640 for 0<=delta<=5/6, 0<=x<=1/2, q = x delta, R = R_*",
    [("\\label{lem:balanced-endpoint}", 15994, 15994), ("For \\(0\\le\\delta\\le5/6\\) and \\(0\\le x\\le1/2\\)", 15995, 15996),
     ("-E_*\\ge\\frac{49}{440640}>\\frac1{10000}.", 15998, 16002)])
chk("C93", 1450, "eq:adaptive-endpoint", "E_actual(h) - Delta <= -49/440640 - (1-h/4)Delta + eps_total = ... - (51/64)Delta",
    [("E_{\\rm actual}(h)-\\Delta", 16095, 16100), ("\\le-\\frac{49}{440640}-\\left(1-\\frac h4\\right)\\Delta", 16096, 16100),
     ("=-\\frac{49}{440640}-\\frac{51}{64}\\Delta", 16098, 16102), ("\\label{eq:adaptive-endpoint}", 16104, 16104)])
chk("C94", 1452, "h Delta/4 from kappa", "kappa = 3/4+2Delta; capacity comparison raises the row exponent by at most Delta/4; "
    "coefficient of R is d", [("\\kappa=\\frac34+2\\Delta,\\qquad \\alpha=\\frac56.", 15924, 15928),
    ("raises the row exponent by at most \\(\\Delta/4\\)", 15958, 15962), ("is \\(d\\), so at \\(d=h\\) the capacity replacement contributes", 16085, 16092)])

# ---------------------------------------------------------------- table of line ranges (L1470-1480) and L2003-2279
bounds = [1648, 3313, 4687, 6807, 9209, 12477, 14985]
starts_ok = all(OL[b - 1].startswith("\\section") or OL[b - 1].startswith("\\part") for b in bounds)
contents = [(401, 1, 1648), (965, 1, 1648), (1425, 1, 1648), (1677, 1648, 3313), (3139, 1648, 3313), (3397, 3313, 4687),
            (3831, 3313, 4687), (4201, 3313, 4687), (4707, 4687, 6807), (6493, 4687, 6807), (6810, 6807, 9209),
            (8061, 6807, 9209), (9250, 9209, 12477), (12362, 9209, 12477), (12531, 12477, 14985), (15015, 14985, 16677),
            (15994, 14985, 16677)]
in_range = all(lo <= l <= hi for l, lo, hi in contents)
last = n == 16677
cited("C95", 1470, "table of OAI78 line ranges", "eight consecutive parts 1-1648,...,14985-16677 start at section/part "
      "boundaries and contain the cited lemmas", starts_ok and in_range and last and len(bounds) + 1 == 8,
      f"boundaries {bounds} are \\section/\\part lines: {starts_ok}; cited labels inside the stated parts: {in_range}; last line {n}")
f, at = has("eq:reflection-cross-prime-phase", 2003, 2279)
f2, at2 = has("cubic supplementary laws", 2003, 2006)
cited("C96", 1488, "L2003-2279", "cross-prime phase bookkeeping lies in OAI78 L2003-2279", f and f2,
      f"supplementary-law periodicity at L{at2}, eq:reflection-cross-prime-phase at L{at}")

# ---------------------------------------------------------------- G78 section numbers
GL = read_lines(G78)
sec = 0; sub = 0; where = {}
for i, l in enumerate(GL, 1):
    if l.startswith("\\section{"):
        sec += 1; sub = 0; where[(sec, 0)] = (i, l)
    elif l.startswith("\\subsection{"):
        sub += 1; where[(sec, sub)] = (i, l)
props = [i for i, l in enumerate(GL, 1) if "\\begin{proposition}" in l]
def g(s_, ss, needle):
    i, l = where.get((s_, ss), (0, ""))
    return needle in l, i
checks = [("G78 Prop. 1", props[:1] == [137] and "prop:criterion" in GL[136], 137),
          ("G78 3.1 probe formula", *g(3, 1, "The formula")), ("G78 3.2 first evaluation", *g(3, 2, "First evaluation")),
          ("G78 3.3 second evaluation", *g(3, 3, "Second evaluation")), ("G78 4.2 compensated probe", *g(4, 2, "compensated probe")),
          ("G78 4.4 high side", *g(4, 4, "High side")), ("G78 5.1 three averages", *g(5, 1, "Three averages")),
          ("G78 6 final inequality", *g(6, 0, "The final inequality"))]
ok = all(c[1] for c in checks)
f1, a1 = src_has(G78, "Q_i(u)=P_i^{-1/2}\\sum_{p\\in\\text{slot }i}\\Omega_i\\Bigl(\\frac{\\N(p)}{P_i}\\Bigr)", 470, 480)
f2, a2 = src_has(G78, "to each row a bin $a\\le\\beta_*$", 320, 326)
f3, a3 = src_has(G78, "\\label{eq:witness}", 326, 330)
f4, a4 = src_has(G78, "A large sieve for sextic characters", 336, 342)
cited("C97", 1083, "G78 references", "G78 Prop. 1, Sections 3.1, 3.2, 3.3 (bin ceiling, witnesses, sextic sieve), 4.2, 4.4 "
      "(Q_i(u)), 5.1, 6 are the cited items", ok and f1 and f2 and f3 and f4,
      "; ".join(f"{c[0]}: L{c[2]}" for c in checks) + f"; Q_i(u) at G78 L{a1}; bin a<=beta_* L{a2}; witnesses L{a3}; sieve L{a4}")

# ---------------------------------------------------------------- C105 counts in the status section (L1463, L1484)
from cmreport import GUIDE
GLs = read_lines(GUIDE)
i0 = next(i for i, l in enumerate(GLs, 1) if l.startswith("Six statements must still be written out"))
i1 = next(i for i in range(i0, len(GLs)) if GLs[i].startswith("\\end{enumerate}"))
items = sum(1 for l in GLs[i0:i1] if l.lstrip().startswith("\\item"))
j0 = next(i for i, l in enumerate(GLs, 1) if "lines & content & CM change" in l)
rows_ = sum(1 for l in GLs[j0:j0 + 12] if re.match(r"^\d+--\d+ &", l))
claim("C105", 1484, "six statements; eight parts", "E", "the status list has 6 items and the table 8 consecutive parts",
      items == 6 and rows_ == 8, f"items {items} (L{i0}-{i1 + 1}), table rows {rows_}; falsifier: 7 items",
      falsifier_failed=items != 7)

f1, a1 = src_has(G78, "Part~I uses only $M_r$ with", 330, 345)
f2, a2 = src_has(G78, "the size\n$U^{\\delta/2}$", 570, 595)
f2, a2 = src_has(G78, "$U^{\\delta/2}$ of the numerator", 570, 595) if not f2 else (f2, a2)
cited("C106", 1351, "only row count of Part I; numerator U^{delta/2}", "G78: Part I uses only M_r (sextic large sieve); "
      "G78 6: the size U^{delta/2} of the numerator on the reflected line", f1 and f2, f"G78 L{a1}, L{a2}")
