\\ t4_lfunctions.gp -- PARI/GP part of t4_lfunctions.py (guide displays eq:d071, eq:d077, eq:d078, eq:d087).
\\ Prints lines "TAG|name|v1|v2|..." that the Python driver parses and judges.  Read-only; no files written.
\\ Characters: PARI Grossencharacters (gcharinit) of K = Q(sqrt(-3)).  A character whose local component at the
\\ complex place is [k,0] has infinity type that of conj(alpha)^k (checked for k = 1 below: gchar [1,0] mod 3 is
\\ conj(alpha) on primary elements), so it is conj(alpha)^k * eta with eta of finite order.

default(realprecision, 38);
bnf = bnfinit(y^2+y+1, 1);
emb(h) = subst(lift(h), y, (-1+sqrt(-3))/2);
prim(pr) = {my(g, h); g = nfbasistoalg(bnf, bnfisprincipal(bnf, pr)[2]);
  for(k = 0, 5, h = lift(g*Mod(-y, y^2+y+1)^k); if((polcoeff(h,0)-1)%3 == 0 && polcoeff(h,1)%3 == 0, return(h)));
  error("no primary generator")};

gc3 = gcharinit(bnf, 3);
gc84 = gcharinit(bnf, 84);
\\ first character (in forvec order) of infinity type [k,0] mod 84 with conductor norm q
pick(k, q) = {my(r = 0); forvec(v = [[0,5],[0,5],[0,5],[0,1]], my(c = [v[1],v[2],v[3],v[4],k,0]~);
  if(idealnorm(bnf, gcharconductor(gc84, c)[1]) == q, r = c; break)); if(r == 0, error("none")); r};

\\ ---------------------------------------------------------------- the characters
CH = List();
listput(CH, ["abar_mod3", gc3, [1,0]~, 1]);
foreach([7, 28, 112, 441, 1008], q, listput(CH, [Str("k1_Q", q), gc84, pick(1, q), 1]));
HW = List();
listput(HW, ["abar2_mod3", gc3, [2,0]~, 2]);
listput(HW, ["abar3_mod3", gc3, [3,0]~, 3]);
listput(HW, ["k2_Q112", gc84, pick(2, 112), 2]);
listput(HW, ["k3_Q28", gc84, pick(3, 28), 3]);

\\ ---------------------------------------------------------------- 1. identification and parameters
\\ ID|name|k|Q|N_PARI|Vga|feqbits|order(psi*alpha^k on primary primes)|max|psi(pi)|-1|
{
for(i = 1, #CH + #HW, my(c = if(i <= #CH, CH[i], HW[i-#CH]), gc = c[2], chi = c[3], k = c[4], L, Q, P, ordr, dev = 0.,
    vals = List());
  L = lfuncreate([gc, chi]); Q = idealnorm(bnf, gcharconductor(gc, chi)[1]);
  if(gcharlocal(gc, chi, 1) != [k, 0], error("infinity type"));
  forprime(p = 5, 300, P = idealprimedec(bnf, p);
    for(j = 1, #P, my(pr = P[j], h = prim(pr), a = emb(h), v);
      if(idealval(bnf, 84, pr) > 0, next);
      v = gchareval(gc, chi, pr)*(a/abs(a))^k; dev = max(dev, abs(abs(v) - 1)); listput(vals, v)));
  ordr = 0; for(n = 1, 60, if(vecmax(apply(v -> abs(v^n - 1), Vec(vals))) < 1e-25, ordr = n; break));
  print("ID|", c[1], "|", k, "|", Q, "|", lfunparams(L)[1], "|", lfunparams(L)[3], "|", lfuncheckfeq(L), "|", ordr, "|",
        dev));
}
\\ L(s, conj(alpha)) at s = 5/2+3i, for comparison with an independent Dirichlet sum in Python
print("LVAL|abar_mod3|", real(lfun(lfuncreate([gc3,[1,0]~]), 5/2+3*I)), "|", imag(lfun(lfuncreate([gc3,[1,0]~]), 5/2+3*I)));

\\ ---------------------------------------------------------------- 2. eq:d077: Lambda(s) = eps conj(Lambda(1 - conj(s)))
\\ Lambda_v(s) = (A Q)^{s/2} (B pi)^{-s} Gamma(s + g) L(s); variant v = [A, B, g]; guide: [3, 2, 1/2]
Lam(L, Q, v, s) = (v[1]*Q)^(s/2) * (v[2]*Pi)^(-s) * gamma(s + v[3]) * lfun(L, s);
SPTS = [0.3+2*I, 0.8-5*I, 1.7+11*I, -0.6+3*I, 0.5+0.25*I];
{
for(i = 1, #CH, my(c = CH[i], L = lfuncreate([c[2], c[3]]), Q = idealnorm(bnf, gcharconductor(c[2], c[3])[1]), out = []);
  foreach([[3,2,1/2], [1,2,1/2], [3,2,0], [3,1,1/2]], v,
    my(e = vector(#SPTS, j, my(s = SPTS[j]); Lam(L, Q, v, s)/conj(Lam(L, Q, v, 1 - conj(s)))));
    out = concat(out, [vecmax(apply(z -> abs(z - e[1]), e)), abs(abs(e[1]) - 1)]));
  print("FE|", c[1], "|", Q, "|", strjoin(apply(x -> Str(x), out), "|")));
}

\\ ---------------------------------------------------------------- 3. eq:d078: growth ratio
\\ r(sigma,t) = |L(sigma+it)| / ((3Q)^{1/2-sigma} (3+|t|)^{1-2sigma} |L(1-sigma-it, conj psi)|),
\\ |L(1-sigma-it, conj psi)| = |L(1-sigma+it, psi)|; pred = (2pi)^{2sigma-1}|Gamma(3/2-s)/Gamma(1/2+s)|/(3+|t|)^{1-2sigma}
default(realprecision, 28);
SIG = [-2/5, -1/10, 1/5, 4/5, 11/10, 7/5];
TT = [0, 3, -3, 10, -10, 30, -30, 100, -100];
{
for(i = 1, #CH, my(c = CH[i], L, Q = idealnorm(bnf, gcharconductor(c[2], c[3])[1]), rmax = 0., relmax = 0., f3 = 0., f100 = 0.);
  L = lfuninit(lfuncreate([c[2], c[3]]), [1/2, 1, 101]);
  foreach(SIG, sg, foreach(TT, t, my(s = sg + I*t, num, den, r, pr, rf);
    num = abs(lfun(L, s)); den = (3*Q)^(1/2 - sg)*(3 + abs(t))^(1 - 2*sg)*abs(lfun(L, 1 - sg + I*t));
    r = num/den; pr = (2*Pi)^(2*sg - 1)*abs(gamma(3/2 - s)/gamma(1/2 + s))/(3 + abs(t))^(1 - 2*sg);
    rmax = max(rmax, r); relmax = max(relmax, abs(r/pr - 1));
    if(sg == -2/5 && t == 3, f3 = num/((3*Q)^(1/2 - sg)*(3 + abs(t))^((1 - 2*sg)/2)*abs(lfun(L, 1 - sg + I*t))));
    if(sg == -2/5 && t == 100, f100 = num/((3*Q)^(1/2 - sg)*(3 + abs(t))^((1 - 2*sg)/2)*abs(lfun(L, 1 - sg + I*t))))));
  print("GR|", c[1], "|", Q, "|", rmax, "|", relmax, "|", f100/f3));
}

\\ ---------------------------------------------------------------- 4. argument principle
default(realprecision, 19);
wind(L, a, b, c, e, n) = {my(P = [], v, tot = 0., prev, cur, d, mn, dmax = 0.);
  P = concat(P, vector(n, k, a + c*I + (b - a)*(k - 1)/n));
  P = concat(P, vector(4*n, k, b + c*I + (e - c)*I*(k - 1)/(4*n)));
  P = concat(P, vector(n, k, b + e*I - (b - a)*(k - 1)/n));
  P = concat(P, vector(4*n, k, a + e*I - (e - c)*I*(k - 1)/(4*n)));
  P = concat(P, [P[1]]);
  v = lfun(L, P[1]); prev = arg(v); mn = abs(v);
  for(k = 2, #P, v = lfun(L, P[k]); mn = min(mn, abs(v)); cur = arg(v); d = cur - prev;
    while(d > Pi, d -= 2*Pi); while(d < -Pi, d += 2*Pi); tot += d; dmax = max(dmax, abs(d)); prev = cur);
  [tot/(2*Pi), mn, dmax]};
\\ AP|name|kind|box winding|min|L| on boundary|max arg step|critical-box winding|#lfunzeros in the critical box
{
for(i = 1, #CH + #HW, my(c = if(i <= #CH, CH[i], HW[i-#CH]), L, r1, r2, nz);
  L = lfuninit(lfuncreate([c[2], c[3]]), [1.2, 0.85, 31]);
  r1 = wind(L, 0.88, 2.0, -30., 30., 1500);
  r2 = wind(L, 0.4, 0.6, 0.5, 15., 600); nz = #lfunzeros(L, [0.5, 15]);
  print("AP|", c[1], "|", if(i <= #CH, "k1", "hw"), "|", r1[1], "|", r1[2], "|", r1[3], "|", r2[1], "|", nz));
}
\\ the two eta products of Section 17 (arithmetic normalization, centre w/2, w = weight = k+1)
{
foreach([["eta(3z)^8", Mat([3,8])], ["eta(2z)^3eta(6z)^3", [2,3;6,3]]], c,
  my(E = lfunetaquo(c[2]), w = lfunparams(E)[2], L, r1, r2, nz);
  L = lfuninit(E, [w/2 + 0.7, 0.85, 31]);
  r1 = wind(L, w/2 + 0.38, w/2 + 1.5, -30., 30., 1500);
  r2 = wind(L, w/2 - 0.1, w/2 + 0.1, 0.5, 15., 600); nz = #lfunzeros(L, [0.5, 15]);
  print("AP|", c[1], "|eta w=", w, "|", r1[1], "|", r1[2], "|", r1[3], "|", r2[1], "|", nz));
}
