"""b3 symbolic/numeric checks of the CM-specific constants in Steps 3-6 and Lemma b3-decay, and of the chain rule (Step 4).
Part 1: constants (mpmath/sympy). Part 2: generic-F chain rule at z = zbar = 0 and the explicit derivatives of Step 4.
Usage: python3 b3_symbolic_checks.py   (about 1 minute)"""
import mpmath as mp, sympy as sp
import os, sys
mp.mp.dps=40
CHECKS=[]   # (name, ok); exit status added 2026-10-12, see the end of the file
# (1) single-mode two-derivative Mellin: int_0^inf (2 pi i conj l)^2 v K_{1/3}(4 pi |l| v) v^{2s-1} dv
#     = -(pi/2)(2pi)^{-2s} Gamma(s+1/3)Gamma(s+2/3) conj(alpha(l))^2 N(l)^{1/2-s}
worst=0
for l in [mp.mpc(0.3,0.7), mp.mpc(-1.1,0.2), mp.mpc(0.05,-0.4)]:
  for s in [mp.mpf(1.7), mp.mpc(2.1,0.8), mp.mpc(1.6,-3)]:
    lhs=(2j*mp.pi*mp.conj(l))**2*mp.quad(lambda v: v*mp.besselk(mp.mpf(1)/3,4*mp.pi*abs(l)*v)*v**(2*s-1),[0,1,mp.inf])
    al=l/abs(l); Nl=abs(l)**2
    rhs=-(mp.pi/2)*(2*mp.pi)**(-2*s)*mp.gamma(s+mp.mpf(1)/3)*mp.gamma(s+mp.mpf(2)/3)*mp.conj(al)**2*Nl**(mp.mpf(1)/2-s)
    worst=max(worst,abs(lhs-rhs)/abs(rhs))
print("(1) mode Mellin, two derivatives: worst rel err", mp.nstr(worst,3))
CHECKS.append(("(1) worst rel err <= 1e-20", worst <= 1e-20))
# (2) constant in J(s)=(3pi/2)(27/(2pi)^2)^s G G T(s-1/2):
s=sp.symbols('s')
expr=-(sp.pi/2)*(2*sp.pi)**(-2*s)*(-1)*27**(s-sp.Rational(1,2))*3**sp.Rational(5,2)
r2=sp.simplify(expr/((3*sp.pi/2)*(27/(2*sp.pi)**2)**s))
print("(2) J constant ratio simplifies to 1:", r2)
CHECKS.append(("(2) ratio == 1", r2 == 1))
# (3) assembly constant: [-(pi/2)(2pi)^{-2-2t}] / [(3pi/2)(27/(2pi)^2)^{1-t}] = -(1/81)((2pi)^4/27)^{-t}
t=sp.symbols('t')
lhs=(-(sp.pi/2)*(2*sp.pi)**(-2-2*t))/((3*sp.pi/2)*(27/(2*sp.pi)**2)**(1-t))
rhs=-sp.Rational(1,81)*((2*sp.pi)**4/27)**(-t)
r3=sp.simplify(sp.expand_power_base(sp.powsimp(lhs/rhs,force=True),force=True))
print("(3) assembly constant ratio:", r3)
CHECKS.append(("(3) ratio == 1", r3 == 1))
for tv in [0.3, 1.7+2j, -0.4+5j]:
    r3n=complex(sp.N((lhs/rhs).subs(t,tv)))
    print("    numeric t=",tv, r3n)
    CHECKS.append((f"(3) numeric ratio at t={tv} within 1e-12 of 1", abs(r3n-1) <= 1e-12))
# (4) substitution v=1/(N(c)v'): (cv)^{-4} v^{2s-1} dv -> c^{-4} N(c)^{4-2s} v'^{3-2s} dv' ; c^{-4}N(c)^2=conj(alpha(c))^4
c=mp.mpc(2.3,-1.4); N=abs(c)**2
e4=abs(c**-4*N**2-mp.conj(c/abs(c))**4)
print("(4) c^-4 N(c)^2 vs conj(alpha(c))^4:", mp.nstr(e4,3))
CHECKS.append(("(4) c^-4 N(c)^2 = conj(alpha(c))^4 within 1e-30", e4 <= 1e-30))
G=lambda vp: mp.exp(-vp)*vp**2   # test function
for sv in [mp.mpc(0.2,1.0), mp.mpc(-0.7,0.3)]:
    L=mp.quad(lambda v: (c*v)**-4*G(1/(N*v))*v**(2*sv-1),[0,mp.mpf(1)/N,1,mp.inf])
    R=c**-4*N**(4-2*sv)*mp.quad(lambda vp: G(vp)*vp**(3-2*sv),[0,1,mp.inf])
    e4s=abs(L-R)/abs(R)
    print("    substitution check rel err", mp.nstr(e4s,3))
    CHECKS.append((f"(4) substitution rel err <= 1e-20 at s={sv}", e4s <= 1e-20))
# (5) variables: t=1/2-w gives 3/2-w=1+t, w+5/6=4/3-t, w+7/6=5/3-t ; Gamma(s+1/3) at s=w+1/2 is Gamma(w+5/6)
w=sp.symbols('w'); tt=sp.Rational(1,2)-w
r5=(sp.simplify(sp.Rational(3,2)-w-(1+tt)), sp.simplify(w+sp.Rational(5,6)-(sp.Rational(4,3)-tt)), sp.simplify(w+sp.Rational(7,6)-(sp.Rational(5,3)-tt)))
print("(5)", *r5)
CHECKS.append(("(5) three shifts == 0", all(x == 0 for x in r5)))
# (6) Stirling exponent: |Q(eta+iu)| ~ |u|^{4 eta}
Q=lambda t: mp.gamma(mp.mpf(4)/3+t)*mp.gamma(mp.mpf(5)/3+t)/(mp.gamma(mp.mpf(4)/3-t)*mp.gamma(mp.mpf(5)/3-t))
for eta in [-0.25,0,0.5,2,5]:
    u1,u2=mp.mpf(1e3),mp.mpf(1e5)
    ex=(mp.log(abs(Q(eta+1j*u2)))-mp.log(abs(Q(eta+1j*u1))))/(mp.log(u2)-mp.log(u1))
    print(f"(6) eta={eta}: growth exponent {mp.nstr(ex,8)} vs 4eta={4*eta}")
    CHECKS.append((f"(6) eta={eta}: fitted exponent within 1e-3 of 4eta", abs(ex-4*eta) <= 1e-3))
# (7) alpha(lambda^-3)=i, conj(alpha(lambda^-3))^2=-1; 3^{5/2} 27^{-1/2}=3; |b| N(b)^{3/2-3s}=N(b)^{2-3s}
lam=mp.mpc(0,mp.sqrt(3)); a=lam**-3/abs(lam**-3)
print("(7) alpha(lambda^-3)=",mp.nstr(a,5)," conj^2=",mp.nstr(mp.conj(a)**2,5), " 3^{5/2}/27^{1/2}=",mp.nstr(mp.mpf(3)**2.5/mp.sqrt(27),10))
CHECKS.append(("(7) alpha(lambda^-3) = i, conj^2 = -1, 3^{5/2}/27^{1/2} = 3 (within 1e-12)",
               abs(a-1j) <= 1e-12 and abs(mp.conj(a)**2+1) <= 1e-12 and abs(mp.mpf(3)**2.5/mp.sqrt(27)-3) <= 1e-12))

# ---------------- Part 2 ----------------
# Corrected version of part (3) of an earlier independent chain-rule check: evaluate at z = zbar = 0 (both), generic F.
z,zb,v,c,cb,d,db=sp.symbols('z zb v c cb d db')
claimZ=-d/c - zb/(c**2*(v**2+z*zb)); claimZb=-db/cb - z/(cb**2*(v**2+z*zb)); claimV=v/(c*cb*(v**2+z*zb))
F=sp.Function('F'); Zs,Zbs,Vs=sp.symbols('Zs Zbs Vs')
for j in (1,2,3):
    expr=sp.diff(F(claimZ,claimZb,claimV),zb,j).subs({z:0,zb:0}).doit()
    target=((-1/(c**2*v**2))**j*sp.diff(F(Zs,Zbs,Vs),Zs,j)).subs({Zs:-d/c,Zbs:-db/cb,Vs:1/(c*cb*v)}).doit()
    okj=sp.simplify(expr-target)==0
    print(f"j={j}: generic-F chain rule at z=zbar=0:", okj)
    CHECKS.append((f"chain rule j={j}", okj))
# explicit first and second zbar-derivatives used in the guide insert
D=lambda e,k: sp.simplify(sp.diff(e,zb,k))
dzp=(sp.simplify(D(claimZ,1)-(-v**2/(c**2*(v**2+z*zb)**2)))==0, sp.simplify(D(claimZ,2)-(2*z*v**2/(c**2*(v**2+z*zb)**3)))==0)
print("dz'   :", *dzp)
dzbp=(sp.simplify(D(claimZb,1)-(z**2/(cb**2*(v**2+z*zb)**2)))==0, sp.simplify(D(claimZb,2)-(-2*z**3/(cb**2*(v**2+z*zb)**3)))==0)
print("dzbar':", *dzbp)
dvp=(sp.simplify(D(claimV,1)-(-z*v/(c*cb*(v**2+z*zb)**2)))==0, sp.simplify(D(claimV,2)-(2*z**2*v/(c*cb*(v**2+z*zb)**3)))==0)
print("dv'   :", *dvp)
CHECKS += [("dz' first and second", all(dzp)), ("dzbar' first and second", all(dzbp)), ("dv' first and second", all(dvp))]
# (8) growth of the gamma ratio on Re w = sigma_2 (Step 5): |Gamma(11/6-w)Gamma(13/6-w)/(Gamma(w+5/6)Gamma(w+7/6))| ~ |Im w|^{2-4 sigma_2}
for s2 in [-1,-0.5,-0.1]:
    f=lambda y: abs(mp.gamma(mp.mpf(11)/6-(s2+1j*y))*mp.gamma(mp.mpf(13)/6-(s2+1j*y))/(mp.gamma(s2+1j*y+mp.mpf(5)/6)*mp.gamma(s2+1j*y+mp.mpf(7)/6)))
    ex=(mp.log(f(1e5))-mp.log(f(1e3)))/(mp.log(1e5)-mp.log(1e3))
    print(f"(8) sigma_2={s2}: exponent {mp.nstr(ex,7)} vs 2-4sigma_2={2-4*s2}")
    CHECKS.append((f"(8) sigma_2={s2}: fitted exponent within 1e-3 of 2-4sigma_2", abs(ex-(2-4*s2)) <= 1e-3))
# (9) lambda^2/72 = 1/(8 lambda^2): c_0 | 8 lambda^2
lam=sp.I*sp.sqrt(3); r9=sp.simplify(lam**2/72*8*lam**2); print("(9) lambda^2/72*8*lambda^2 =", r9)
CHECKS.append(("(9) == 1", r9 == 1))

# Exit status (added 2026-10-12): 1 unless every check recorded in CHECKS holds (exact results equal their targets;
# mpmath errors at 40 digits <= 1e-20 or 1e-30; fitted growth exponents within 1e-3).  CM_TESTS_FORCE_FAIL=1 counts
# the first check as failed (self-test of the wiring).
if os.environ.get('CM_TESTS_FORCE_FAIL') == '1':
    CHECKS[0] = (CHECKS[0][0] + ' [CM_TESTS_FORCE_FAIL=1]', False)
bad = [n for n, ok in CHECKS if not ok]
print('EXIT STATUS %d: ' % (1 if bad else 0) + ('%d of %d checks failed: %s' % (len(bad), len(CHECKS), '; '.join(bad))
                                                if bad else 'all %d checks passed' % len(CHECKS)))
sys.exit(1 if bad else 0)
