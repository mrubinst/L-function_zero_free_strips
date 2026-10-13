"""t4 / F: L-functions of members of F_CM and of higher-weight CM characters (PARI/GP Grossencharacters).

Guide displays (cm-eleven-twelfths-guide.tex, 2026-10-13 numbering):
  eq:d077 (13.1)  Lambda(s,psi) = (3Q)^{s/2} (2 pi)^{-s} Gamma(s+1/2) L(s,psi) satisfies Lambda(s,psi) = eps Lambda(1-s,conj psi),
                  |eps| = 1, Q = N(conductor of psi).
  eq:d078 (13.2)  |L(sigma+it,psi)| << (3Q)^{1/2-sigma} (3+|t|)^{1-2 sigma} |L(1-sigma-it, conj psi)| on compact I in (-1/2,3/2).
  eq:d071 (11.1)  L(s,psi) != 0 for Re s > 7/8 (finite-height sanity check only).
  eq:d087 (17.1)  expected: L(f,s) != 0 for Re s > k/2 + 7/8, weight k+1 (finite-height sanity check only).

Characters: PARI gcharinit over K = Q(sqrt -3), moduli 3 and 84; those with local component [k,0] at the complex place
(infinity type of conj(alpha)^k).  Six members of F_CM (k = 1, conductor norms Q = 9, 7, 28, 112, 441, 1008, most of
them not self-dual), four of higher weight (k = 2, 3), and the eta products eta(3z)^8, eta(2z)^3 eta(6z)^3.

Checks
  F0  parameters: gchar [1,0] mod 3 equals conj(alpha) on primary primes (order of psi*alpha^k is 1) and
      L(5/2+3i, conj alpha) agrees with an independent Dirichlet sum over primary n (numpy); every psi * alpha^k
      has finite order on primary primes (so psi = conj(alpha)^k eta, eta of finite order); PARI's conductor is 3Q
      and its gamma factor is Gamma_R(s+k/2)Gamma_R(s+k/2+1) (= const (2 pi)^{-s} Gamma(s+k/2)); lfuncheckfeq <= -100.
  F1  eq:d077: Lambda(s)/conj(Lambda(1-conj s)) is the same unimodular eps at five points s (incl. Re s = -0.6 and
      Re s = 1.7), for the six members.  Falsifiers: Q for 3Q; Gamma(s) for Gamma(s+1/2); pi^{-s} for (2 pi)^{-s}.
  F2  eq:d078 (size check): the ratio |L(s)| / ((3Q)^{1/2-sigma}(3+|t|)^{1-2 sigma}|L(1-sigma-it, conj psi)|) equals
      (2pi)^{2sigma-1}|Gamma(3/2-s)/Gamma(1/2+s)|/(3+|t|)^{1-2sigma} (independent of psi and Q) at sigma in
      {-2/5,-1/10,1/5,4/5,11/10,7/5}, |t| <= 100; that function is bounded on I = [-2/5,7/5] x R (sup on a grid
      to |t| = 1e6) and tends to (2pi)^{2sigma-1}.  Falsifier: exponent (1-2sigma)/2 makes the ratio grow in t.
  F3  eq:d071 (sanity): argument principle on 0.88 <= Re s <= 2, |Im s| <= 30: no zeros, for the six members of F_CM.
      Falsifier: the same winding count on the box 0.4 <= Re s <= 0.6, 0.5 <= Im s <= 15 equals the number of
      zeros found by lfunzeros.
  F4  eq:d087 (sanity): the same in the unitary normalization for conj(alpha)^2, conj(alpha)^3 (= eta(3z)^8 shifted)
      and two twisted characters of weights 3 and 4, and for the eta products on w/2+0.38 <= Re s <= w/2+1.5.
Exit status 1 on any failure.  Run time: about 20 s.
"""
import os, sys, math, cmath, subprocess, time
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()
NFAIL = [0]
def report(tid, ok, msg):
    if not ok:
        NFAIL[0] += 1
    print(f"[{'PASS' if ok else 'FAIL'}] {tid} {msg}", flush=True)

out = subprocess.run(["gp", "-q", "-f"], stdin=open(os.path.join(HERE, "t4_lfunctions.gp")), capture_output=True,
                     text=True, timeout=1200)
if out.returncode != 0 or "***" in out.stdout + out.stderr:
    print(out.stdout, out.stderr)
    report("GP", False, "gp script failed"); sys.exit(1)
rows = {}
for l in out.stdout.splitlines():
    f = l.split("|")
    rows.setdefault(f[0], []).append(f[1:])
num = lambda s: float(s.replace(" E", "e").replace(" ", ""))

# ------------------------------------------------------------------ F0
ok = True; det = []
for name, k, Q, NP, Vga, feq, order, dev in rows["ID"]:
    k, Q, NP, feq, order = int(k), int(Q), int(NP), int(feq), int(order)
    vg = [float(eval(x)) for x in Vga.strip("[]").split(",")]
    good = (NP == 3*Q and abs(vg[0] - k/2) < 1e-12 and abs(vg[1] - k/2 - 1) < 1e-12 and feq <= -100
            and 1 <= order <= 60 and num(dev) < 1e-20)
    if name == "abar_mod3":
        good &= order == 1
    ok &= good
    det.append(f"{name}: k={k} Q={Q} N_PARI={NP} ord={order}")
# independent Dirichlet sum for L(s, conj alpha), s = 5/2+3i: primary n = a + b omega, a = 1, b = 0 mod 3
s = 2.5 + 3j; X = 3.0e5; R = int(math.sqrt(4*X/3)) + 3
A = np.arange(-R, R + 1); a, b = np.meshgrid(A, A, indexing="ij")
msk = ((a - 1) % 3 == 0) & (b % 3 == 0)
a = a[msk].astype(float); b = b[msk].astype(float)
Nn = a*a - a*b + b*b; keep = (Nn > 0) & (Nn <= X); a, b, Nn = a[keep], b[keep], Nn[keep]
z = a + b*complex(-0.5, math.sqrt(3)/2)
Lpy = np.sum(np.conj(z)/np.abs(z)*Nn**(-s))
lv = rows["LVAL"][0]; Lgp = complex(num(lv[1]), num(lv[2]))
tail = 2*math.pi/(3*math.sqrt(3))*X**(-1.5)/1.5      # crude size of the omitted tail
Lfal = np.sum((np.conj(z)/np.abs(z))**2*Nn**(-s))     # falsifier: conj(alpha)^2 (alpha itself gives the same sum,
                                                     # since the primary elements are closed under conjugation)
ok &= abs(Lpy - Lgp) < 10*tail and abs(Lfal - Lgp) > 1e-3
report("F0", ok, "parameters: " + "; ".join(det) + f". gchar [1,0] mod 3 = conj(alpha) on primary primes; "
       f"L(5/2+3i, conj alpha): PARI {Lgp:.10f}, Dirichlet sum (N(n)<=3e5) {Lpy:.10f}, diff {abs(Lpy-Lgp):.1e} "
       f"(tail bound {tail:.0e}); falsifier conj(alpha)^2: diff {abs(Lfal-Lgp):.2e}")

# ------------------------------------------------------------------ F1
ok = True; spread = 0; unim = 0; fmin = [1e9]*3
for r in rows["FE"]:
    v = [num(x) for x in r[2:]]
    spread = max(spread, v[0]); unim = max(unim, v[1])
    for i in range(3):
        fmin[i] = min(fmin[i], v[2 + 2*i])
ok = spread < 1e-25 and unim < 1e-25 and min(fmin) > 0.1
report("F1", ok, f"eq:d077: Lambda(s)/conj(Lambda(1-conj s)) constant over 5 points s, for {len(rows['FE'])} members of "
       f"F_CM (Q = {', '.join(r[1] for r in rows['FE'])}): max spread {spread:.1e}, ||eps|-1| <= {unim:.1e}; "
       f"falsifiers (min spread, must be large): Q for 3Q {fmin[0]:.2f}, Gamma(s) for Gamma(s+1/2) {fmin[1]:.2f}, "
       f"pi^-s for (2pi)^-s {fmin[2]:.2f}")

# ------------------------------------------------------------------ F2
mp.mp.dps = 20
def pred(sg, t):
    s = mp.mpc(sg, t)
    return float((2*mp.pi)**(2*sg - 1)*abs(mp.gamma(1.5 - s)/mp.gamma(0.5 + s))/(3 + abs(t))**(1 - 2*sg))
grid_sup = 0.0; lim_dev = 0.0
for sg in np.linspace(-0.4, 1.4, 19):
    vals = [pred(sg, t) for t in [0] + list(np.logspace(-2, 6, 60))]
    grid_sup = max(grid_sup, max(vals))
    lim_dev = max(lim_dev, abs(vals[-1]/float((2*mp.pi)**(2*sg - 1)) - 1))
relmax = max(num(r[3]) for r in rows["GR"]); rmax = max(num(r[2]) for r in rows["GR"])
fgrow = min(num(r[4]) for r in rows["GR"])
ok = relmax < 1e-20 and grid_sup < 1e4 and lim_dev < 1e-4 and fgrow > 10
report("F2", ok, f"eq:d078 size check: measured ratio = (2pi)^(2sigma-1)|Gamma(3/2-s)/Gamma(1/2+s)|/(3+|t|)^(1-2sigma) "
       f"to rel. err {relmax:.1e} for {len(rows['GR'])} characters (Q = 7..1008), sigma in [-2/5,7/5], |t|<=100; "
       f"max ratio {rmax:.0f} (at sigma=7/5, t=0), the same for every Q; sup of that function on [-2/5,7/5] x [0,1e6] "
       f"{grid_sup:.0f}, limit (2pi)^(2sigma-1) reached to {lim_dev:.0e}; falsifier exponent (1-2sigma)/2: ratio grows "
       f"by x{fgrow:.1f} from |t|=3 to 100 at sigma=-2/5")

# ------------------------------------------------------------------ F3, F4
def judge(rs):
    good = True; d = []
    for name, kind, w1, mn, step, w2, nz in rs:
        w1, mn, step, w2, nz = num(w1), num(mn), num(step), num(w2), int(nz)
        g = abs(w1) < 0.01 and mn > 1e-3 and step < 0.5
        f = abs(w2 - nz) < 0.01 and nz > 0
        good &= g and f
        d.append(f"{name}: winding {w1:+.0e}, min|L| {mn:.2f}; crit. box {w2:.2f} vs {nz}")
    return good, d
k1 = [r for r in rows["AP"] if r[1] == "k1"]; hw = [r for r in rows["AP"] if r[1] != "k1"]
ok, d = judge(k1)
report("F3", ok, "eq:d071 sanity (argument principle, 0.88<=Re s<=2, |Im s|<=30): no zeros; " + "; ".join(d) +
       " (falsifier: critical-box winding = lfunzeros count, so the count detects zeros)")
ok, d = judge(hw)
report("F4", ok, "eq:d087 sanity (unitary normalization 0.88<=Re s<=2, |Im s|<=30; eta products w/2+0.38<=Re s<=w/2+1.5): "
       "no zeros; " + "; ".join(d))
print(f"# t4_lfunctions: {time.time() - T0:.1f} s, {NFAIL[0]} failure(s)")
sys.exit(1 if NFAIL[0] else 0)
