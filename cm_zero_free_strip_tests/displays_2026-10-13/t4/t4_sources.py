"""t4 / S: source comparisons for the definitions (D) and quoted statements (C) of t4's inventory, and a consistency
check that every display of the inventory is still verbatim in the current guide.

S0  every 'source' string of sources/inventory/inventory_t4.json occurs verbatim in the current guide, and the
    printed numbers in sources/cm_zero_free_strip.aux agree with the inventory.
S1  eq:d081 (13.5) is character-for-character G78 eq:amplitude (whitespace ignored).
S2  eq:d080 (13.4) is G78 eq:witness with psi(n) -> conj(alpha(n)) theta(n) conj(chi_n(u)) (the CM row character)
    (whitespace and thin spaces ignored).
S3  needle checks against OAI78 (paper.tex of 2026-09-30) at the stated lines:
      eq:d074 beta_*        <-> old-eq:1.1b (L378-384), family 'primitive finite-order Hecke character' -> F_CM;
      eq:d075 row u         <-> eq:probe-high (L4043-4048) with eta -> psi;
      eq:d076 B_m           <-> eq:unmarked-central-completion C_V (L2819-2826) with abar -> abar^2;
      eq:d199 M_u, Q_u      <-> lem:marked set-up (L9225-9244), psi_u = nu chi_n(u)^eps;
      eq:cm-inverse         <-> old-eq:4.1 (L9258-9262);
      eq:d082 companion     <-> L10774-10776;
      eq:signal-cm          <-> eq:signal (L965-970);
      eq:d083 S_k(n_j)      <-> old-eq:1.1 (L704-707) and psi_k = tau chi_n(k), M = m+q (L12493-12497);
      eq:d084               <-> old-eq:2.1 (L12549-12554);
      eq:d085               <-> lem:balanced-endpoint (L15994-16000);
      eq:d086               <-> eq:adaptive-endpoint (L16096-16104).
    Falsifier: altered needles (e.g. 49/440641, q_n^2, chi_n(u)^{2}) must NOT be found.
Exit status 1 on any failure.  Run time: < 1 s.
"""
import os, re, sys, json

HERE = os.path.dirname(os.path.abspath(__file__))
_d = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(_d, "cmsources.py")) and os.path.dirname(_d) != _d:
    _d = os.path.dirname(_d)
sys.path.insert(0, _d)
import cmsources            # locations of the guide and of the cited sources (tests_cm/sources/)
GUIDE = cmsources.GUIDE
AUX = cmsources.GUIDE_AUX
INV = cmsources.inventory("inventory_t4.json")
OAI = cmsources.OAI78
G78 = cmsources.G78
NFAIL = [0]
def report(tid, ok, msg, word="CITED"):
    if not ok:
        NFAIL[0] += 1
    print(f"[{word if ok else 'FAIL'}] {tid} {msg}", flush=True)

guide = open(GUIDE).read(); oai = open(OAI).read().split("\n"); g78 = open(G78).read()
inv = json.load(open(INV)); aux = open(AUX).read()

# ------------------------------------------------------------------ S0
miss = [e["labels"][0] for e in inv if e["source"] not in guide]
badnum = []
for e in inv:
    m = re.search(r"\\newlabel\{" + re.escape(e["labels"][0]) + r"\}\{\{([^}]*)\}", aux)
    if not m or m.group(1) != e["numbers"][0]:
        badnum.append(e["labels"][0])
report("S0", not miss and not badnum and len(inv) == 29,
       f"{len(inv)} inventory displays: verbatim in the current guide (missing: {miss or 'none'}); aux numbers agree "
       f"(mismatch: {badnum or 'none'})", word="PASS")

def body(text, label):
    m = re.search(r"\\begin\{equation\}\\label\{" + re.escape(label) + r"\}(.*?)\\end\{equation\}", text, re.S)
    return m.group(1)
norm = lambda s: re.sub(r"\s+|\\,", "", s)

# ------------------------------------------------------------------ S1, S2
a = norm(body(guide, "eq:d081")); b = norm(body(g78, "eq:amplitude"))
report("S1", a == b and norm(body(g78, "eq:amplitude").replace("z-1", "z")) != a,
       "eq:d081 = G78 eq:amplitude character for character (falsifier: exponent z for z-1 differs)")
w = body(g78, "eq:witness")
sub = w.replace(r"\psi(n)", r"\overline{\alpha(n)}\theta(n)\overline{\chi_n(u)}").replace(
      r"\psi(l)", r"\overline{\alpha(l)}\theta(l)\overline{\chi_l(u)}")
a = norm(body(guide, "eq:d080")); b = norm(sub)
fal = norm(w.replace(r"\psi(n)", r"\alpha(n)\theta(n)\overline{\chi_n(u)}").replace(
      r"\psi(l)", r"\alpha(l)\theta(l)\overline{\chi_l(u)}"))
report("S2", a == b and a != fal, "eq:d080 = G78 eq:witness with psi -> conj(alpha) theta conj(chi.(u)) "
       "(falsifier: the anti-oriented substitution alpha theta conj(chi) does not match)")

# ------------------------------------------------------------------ S3
def has(lo, hi, needles):
    txt = "\n".join(oai[lo - 1:hi])
    return all(n in txt for n in needles)
CHECKS = [
    ("eq:d074", 378, 384, [r"\beta_* = \sup\left(\{1/2\}\cup", r"1/2\le\Re\rho\le1,\quad L_F(\rho,\eta)=0",
                           "primitive finite-order Hecke character"]),
    ("eq:d075", 4043, 4048, [r"\frac{\zeta_F^S(6z)L^S(w,\chi_\bullet(u))}", r"{L^S(x,\eta\overline{\chi_\bullet(u)})}",
                             r"\mathcal H_{\eta,u}(x,w,z)"]),
    ("eq:d076", 2819, 2826, [r"\gamma_2(c)\overline{\alpha(cn^3)}\nu(cn^3)\chi_{cn^3}(m)", r"{\sqrt{q_c}\,q_n}",
                             r"V\!\left(\frac{q_cq_n^3}{X}\right)", r"c\ {\rm sf},\ n"]),
    ("eq:d199", 9225, 9244, [r"\psi_u(n)=\nu(n)\chi_n(u)^{\varepsilon_\chi}", r"M_u(Z^r;W)=Z^{-r/2}\sum_n\mu(n)\psi_u(n)W(q_n/Z^r)",
                             r"Q_u=\prod_{i\in I}Z^{-z_i/2}", r"\sum_{p\in\mathcal P_i}a_i(p)\psi_u(p)W_i(q_p/Z^{z_i})"]),
    ("eq:cm-inverse", 9258, 9262, [r"0<q_u\ll Z^m", r"|M_u(Z^r;W)Q_u|^2\ll Z^{m+\epsilon}"]),
    ("eq:d082", 10774, 10776, [r"\mu(z)\gamma_{-1}(z)", r"=\bar\alpha(z)\gamma_2(z)\chi_z(-1)\overline{G(z)}"]),
    ("eq:signal-cm", 965, 970, [r"\gamma_1(c)\gamma_2(c)=\mu(c)\alpha(c)G(c)", r"|G(c)|=1"]),
    ("eq:d083", 704, 707, [r"S_\psi(n;W)=Z^{-n/2}\sum_l\psi(l)W(q_l/Z^n)"]),
    ("eq:d083 (rows)", 12493, 12497, [r"\psi_k(n)=\tau(n)\chi_n(k),\qquad M=m+q"]),
    ("eq:d084", 12549, 12554, [r"\sum_{k\in\mathcal R_z}", r"\left|S_{\psi_k}(n_1;W_1)S_{\psi_k}(n_2;W_2)Q\right|^2",
                               r"\ll Z^{M+\epsilon}"]),
    ("eq:d085", 15994, 16000, [r"\(0\le\delta\le5/6\) and \(0\le x\le1/2\)", r"\(q=x\delta\) and \(R=R_*\)",
                               r"-E_*\ge\frac{49}{440640}"]),
    ("eq:d086", 16096, 16104, [r"E_{\rm actual}(h)-\Delta", r"\le-\frac{49}{440640}-\left(1-\frac h4\right)\Delta",
                               r"+\epsilon_{\rm total}", r"=-\frac{49}{440640}-\frac{51}{64}\Delta"]),
]
FALS = [("eq:d085", 15994, 16000, [r"-E_*\ge\frac{49}{440641}"]), ("eq:d076", 2819, 2826, [r"{\sqrt{q_c}\,q_n^2}"]),
        ("eq:d199", 9225, 9244, [r"\psi_u(n)=\nu(n)\chi_n(u)^{2}"])]
ok = all(has(lo, hi, nd) for _, lo, hi, nd in CHECKS)
bad = [lab for lab, lo, hi, nd in CHECKS if not has(lo, hi, nd)]
fal_ok = not any(has(lo, hi, nd) for _, lo, hi, nd in FALS)
report("S3", ok and fal_ok, f"{len(CHECKS)} OAI78 source passages found at the stated lines "
       f"({', '.join(f'{l} L{a}-{b}' for l, a, b, _ in CHECKS)}); not found: {bad or 'none'}; "
       f"falsifier needles absent: {fal_ok}")
print(f"# t4_sources: {NFAIL[0]} failure(s)")
sys.exit(1 if NFAIL[0] else 0)
