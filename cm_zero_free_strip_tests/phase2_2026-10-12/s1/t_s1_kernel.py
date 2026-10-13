"""s1 checks of the analytic constants and of the CM reflected test V# for log-Schwartz tests (no compact support).

 (1) one mode, two zbar-derivatives: int_0^oo (2 pi i conj l)^2 v K_{1/3}(4 pi|l|v) v^{2s-1} dv
       = -(pi/2)(2 pi)^{-2s} G(s+1/3)G(s+2/3) conj(alpha(l))^2 N(l)^{1/2-s};  the same with (2 pi i l)^2 gives alpha(l)^2.
 (2) J(s) = (3 pi/2)(27/(2pi)^2)^s G(s+1/3)G(s+2/3) T(s-1/2)  [constant: 3^{5/2} 27^{-1/2} (-1)(-pi/2) ...];
     assembly ratio  [-(pi/2)(2pi)^{-2(3/2-w)}] / [(3pi/2)(27/(2pi)^2)^{w+1/2}] = -(1/81) K^{-t},  t = 1/2-w, K=(2pi)^4/27.
     Falsifier: OAI78's one-derivative ratio is -(i/81)K^{-t}, not -(1/81)K^{-t}.
 (3) R_2(t) = G(4/3+t)G(5/3+t)/(G(4/3-t)G(5/3-t)): Stirling exponent 4 Re t; holomorphic for Re t > -4/3 (pole at -4/3).
 (4) V# for V_gauss (Vhat(t)=e^{t^2}, OAI78's Gaussian) and V_bes(y)=exp(-y-1/y) (Vhat(t)=2K_t(2)):
     (a) independent of the line Re t = s0 for s0 in {0, 1/2, 2, -1/4, -1} (CM: any s0 > -4/3);
     (b) small x: V#(x) = Res_{t=-4/3} + Res_{t=-5/3} + O(x^{7/3}) = Vhat(4/3) G(1/3)/(G(8/3)G(3)) (Kx)^{4/3}\n         + Vhat(5/3) G(-1/3)/(G(3)G(10/3)) (Kx)^{5/3} + ...;
         falsifier: the one-derivative kernel R_1 gives exponent 5/6 instead;
     (c) |(x d/dx)^j V#(x)| / min(x^{1/4}, (1+x)^{-A}) stays bounded on 1e-8 <= x <= 1e6 (j = 0,1,2; A = 2, 4).
 (5) growth of T^CM on Re w = s2 < 0: |G(11/6-w)G(13/6-w)/(G(w+5/6)G(w+7/6))| ~ |Im w|^{2-4 s2}.
Usage: python3 t_s1_kernel.py   (a few minutes; the Bessel tables dominate)
"""
import os
import mpmath as mp, sympy as sp
mp.mp.dps = 30
third = mp.mpf(1)/3
K = (2*mp.pi)**4/27
out = {}
# (1)
worst = 0
for l in [mp.mpc(0.3, 0.7), mp.mpc(-1.1, 0.2)]:
    for s in [mp.mpf(1.7), mp.mpc(2.1, 0.8)]:
        I = mp.quad(lambda v: v*mp.besselk(third, 4*mp.pi*abs(l)*v)*v**(2*s-1), [0, 1, mp.inf])
        al = l/abs(l); Nl = abs(l)**2
        pred = -(mp.pi/2)*(2*mp.pi)**(-2*s)*mp.gamma(s+third)*mp.gamma(s+2*third)*Nl**(mp.mpf(1)/2-s)
        worst = max(worst, abs((2j*mp.pi*mp.conj(l))**2*I-pred*mp.conj(al)**2)/abs(pred),
                    abs((2j*mp.pi*l)**2*I-pred*al**2)/abs(pred))
out['(1) mode Mellin, two derivatives'] = worst
# (2)
s, t, w = sp.symbols('s t w')
c1 = sp.simplify((sp.pi/2)*3**sp.Rational(5, 2)*27**(s-sp.Rational(1, 2))*(2*sp.pi)**(-2*s)/((3*sp.pi/2)*(27/(2*sp.pi)**2)**s))
print("(2) J-constant ratio (must be 1):", c1)
ratio = lambda wv: (-(mp.pi/2)*(2*mp.pi)**(-2*(mp.mpf(3)/2-wv)))/((3*mp.pi/2)*(27/(2*mp.pi)**2)**(wv+mp.mpf(1)/2))
worst = 0; fals = 1e9
for tv in [mp.mpf('0.3'), mp.mpc(0.6, 2.2), mp.mpc(1.4, -5)]:
    wv = mp.mpf(1)/2-tv
    worst = max(worst, abs(ratio(wv)/(-K**(-tv)/81)-1))
    fals = min(fals, abs(ratio(wv)/(-1j*K**(-tv)/81)-1))
out['(2) assembly ratio vs -(1/81)K^-t'] = worst
print("(2) falsifier -(i/81)K^-t deviation:", mp.nstr(fals, 3))
# (3)
R2 = lambda t: mp.gamma(4*third+t)*mp.gamma(5*third+t)/(mp.gamma(4*third-t)*mp.gamma(5*third-t))
R1 = lambda t: mp.gamma(mp.mpf(7)/6+t)*mp.gamma(mp.mpf(5)/6+t)/(mp.gamma(mp.mpf(7)/6-t)*mp.gamma(mp.mpf(5)/6-t))
for eta in [-0.25, 0, 0.5, 2, 5]:
    u1, u2 = mp.mpf(1e3), mp.mpf(1e5)
    ex = (mp.log(abs(R2(eta+1j*u2)))-mp.log(abs(R2(eta+1j*u1))))/(mp.log(u2)-mp.log(u1))
    print(f"(3) Re t={eta}: Stirling exponent of R_2 = {mp.nstr(ex, 8)} (4 Re t = {4*eta})")
    out[f'(3) exponent at {eta}'] = abs(ex-4*eta)*(u1/10)    # O(1/u) correction scaled
print("(3) |R_2(-4/3+1e-8)| =", mp.nstr(abs(R2(-4*third+mp.mpf('1e-8'))), 5), " (pole);  R_2 on Re t=-1.3:",
      mp.nstr(max(abs(R2(mp.mpc(-1.3, u))) for u in [0, 1, 5, 20]), 5), " finite")
# (4)  V#(x) by the trapezoid rule on a fixed grid of each line (spectrally accurate: analytic, rapidly decaying
#      integrand); the Mellin factor Vhat(-t)R(t) is tabulated once per line with mpmath.
import numpy as np
VH = {'gauss': lambda z: mp.e**(z*z), 'bes': lambda z: 2*mp.besselk(z, 2)}
H = 0.05; U = np.arange(-80, 80+1e-9, H)
TAB = {}
def table(kind, s0, R):
    key = (kind, s0, R)
    if key not in TAB:
        TAB[key] = np.array([complex(VH[kind](-mp.mpc(s0, u))*R(mp.mpc(s0, u))) for u in U])
    return TAB[key]
Kf = float(K)
def Vs(x, kind, s0=0, R=R2, j=0):
    t = s0+1j*U
    return complex(np.sum((-t)**j*table(kind, s0, R)*np.exp(-t*np.log(Kf*float(x))))*H/(2*np.pi))
for kind in ('gauss', 'bes'):
    worst = 0
    for x in [1e-3, 0.7, 30.0]:
        ref = Vs(x, kind, 0)
        for s0 in [0.5, 2, -0.25, -1]:
            worst = max(worst, abs(Vs(x, kind, s0)-ref)/abs(ref))
    out[f'(4a) line independence {kind}'] = worst
    # (4b) small-x asymptotics: residues at -4/3 and -5/3
    res = complex(VH[kind](4*third)*mp.gamma(third)/(mp.gamma(8*third)*mp.gamma(3)))
    res2 = complex(VH[kind](5*third)*mp.gamma(-third)/(mp.gamma(3)*mp.gamma(10*third)))
    xs = [1e-4, 1e-6]
    errs = [abs(Vs(x, kind)/(res*(Kf*x)**(4/3)+res2*(Kf*x)**(5/3))-1) for x in xs]
    slope = (np.log(abs(Vs(xs[1], kind)))-np.log(abs(Vs(xs[0], kind))))/(np.log(xs[1])-np.log(xs[0]))
    slope1 = (np.log(abs(Vs(xs[1], kind, R=R1)))-np.log(abs(Vs(xs[0], kind, R=R1))))/(np.log(xs[1])-np.log(xs[0]))
    print(f"(4b) {kind}: V#(x)/[residue terms] - 1 at x=1e-4,1e-6: {errs[0]:.3g}, {errs[1]:.3g}"
          f" (two residues -4/3,-5/3; O(x) expected); slope {slope:.6g} (-> 4/3); falsifier R_1 slope {slope1:.6g} (-> 5/6)")
    out[f'(4b) small-x {kind}'] = errs[1]*1e2
    # (4c) bounds: line 0 for x <= 10, line 6 for x > 10 (same value by (4a))
    for A in (2, 4):
        mx = 0
        for j in (0, 1, 2):
            for e in range(-8, 7):
                x = 10.0**e
                v = Vs(x, kind, s0=(6 if e > 1 else 0), j=j)
                mx = max(mx, abs(v)/min(x**0.25, (1+x)**(-A)))
        print(f"(4c) {kind}: max over x in [1e-8,1e6], j<=2 of |(x d/dx)^j V#| / min(x^1/4,(1+x)^-{A}) = {mx:.5g}")
    print(f"(4a) {kind}: line independence (Re t = 0.5, 2, -1/4, -1 vs 0): {worst:.2e}")
# (5)
for s2 in [-1, -0.5, -0.1]:
    f = lambda y: abs(mp.gamma(mp.mpf(11)/6-(s2+1j*y))*mp.gamma(mp.mpf(13)/6-(s2+1j*y))/(mp.gamma(s2+1j*y+mp.mpf(5)/6)*mp.gamma(s2+1j*y+mp.mpf(7)/6)))
    ex = (mp.log(f(1e5))-mp.log(f(1e3)))/(mp.log(1e5)-mp.log(1e3))
    print(f"(5) sigma_2={s2}: exponent {mp.nstr(ex, 7)} vs 2-4sigma_2={2-4*s2}")
    out[f'(5) at {s2}'] = abs(ex-(2-4*s2))*100
for k, v in out.items(): print(f"   {k}: {mp.nstr(v, 3)}")
ok = all(v < 1e-6 for k, v in out.items() if not k.startswith('(3)') and not k.startswith('(5)') and not k.startswith('(4b)')) \
     and all(v < 1 for k, v in out.items()) and fals > 0.1
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':   # self-test of the exit status (added 2026-10-12)
    ok = False; print("CM_TESTS_FORCE_FAIL=1: the check below is counted as failed")
print("PASS" if ok else "FAIL")

# A printed failure must also fail the suite.
raise SystemExit(0 if (ok) else 1)
