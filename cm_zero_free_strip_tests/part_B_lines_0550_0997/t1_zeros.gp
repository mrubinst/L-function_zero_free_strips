\\ PARI/GP part of t1_zeros.py (claims B02, B04, B07).  Prints machine-readable lines "KEY|...".
\\ For each d: L(E_d,s+1/2) * sum_{m<=M} q_d(m) m^{-s}  versus  B_d(s) = prod_{p | 6d} L_p(E_d, s+1/2),
\\ with L_p the local factor of a MINIMAL model (good: 1/(1-a_p p^-x+p^(1-2x)), bad: 1/(1-a_p p^-x)).
default(realprecision, 38);
M = 60000;
qtab(E, d, M) =
{
  my(v = vector(M), f, r, p, e, ap = Map());
  forprime(p = 2, M, if((6*d) % p, mapput(ap, p, ellap(E, p))));
  for (m = 1, M,
    if (gcd(m, 6*d) > 1, v[m] = 0; next);
    f = factor(m); r = 1.;
    for (i = 1, #f~,
      p = f[i,1]; e = f[i,2];
      if (e == 1, r *= -mapget(ap, p)/sqrt(p), if (e == 2, r *= 1, r = 0; break)));
    v[m] = r);
  v
};
Bd(Emin, d, s) =
{
  my(N = ellglobalred(Emin)[1], x = s + 1/2, r = 1., ap);
  fordiv(6*d, p, if (isprime(p),
     ap = ellap(Emin, p);
     if (N % p, r /= (1 - ap*p^(-x) + p^(1-2*x)), r /= (1 - ap*p^(-x)))));
  r
};
{
foreach([1, 2, -3, 7, 16, -432, 5], d,
  E = ellinit([0,0,0,0,d]); Emin = ellminimalmodel(E); N = ellglobalred(Emin)[1];
  q = qtab(E, d, M);
  Lf = lfuncreate(E);
  foreach([3, 2.5 + 4*I, 1.7 - 2*I], s,
     S = sum(m = 1, M, q[m]*m^(-s));
     lhs = lfun(Lf, s + 1/2)*S;
     rhs = Bd(Emin, d, s);
     f1 = abs(lhs - 1);                      \\ falsifier: B_d omitted
     f2 = abs(lfun(Lf, s)*S - rhs);          \\ falsifier: no shift by 1/2
     print("MZ|", d, "|", N, "|", s, "|", abs(lhs - rhs), "|", abs(rhs), "|", f1, "|", f2));
  \\ B07: local polynomials at p | 6d for the minimal model: roots in X = p^{-s-1/2}
  fordiv(6*d, p, if (isprime(p),
     ap = ellap(Emin, p);
     if (N % p, rr = polroots(p*'X^2 - ap*'X + 1); print("LR|", d, "|", p, "|good|", ap, "|", abs(rr[1]), "|", abs(rr[2]), "|", p^(-1/2)),
                rr = if (ap, abs(1/ap), 0); print("LR|", d, "|", p, "|bad|", ap, "|", rr, "|", rr, "|", 0))));
  \\ B02 sanity: S_d(D) = sum_m q_d(m) W(m/D), W = bump on (1,2)
  foreach([1000, 10000, 30000], D,
     S = sum(m = D+1, min(2*D-1, M), q[m]*exp(-1/((m/D-1)*(2-m/D))));
     print("SD|", d, "|", D, "|", abs(S)))
);
}
quit;
