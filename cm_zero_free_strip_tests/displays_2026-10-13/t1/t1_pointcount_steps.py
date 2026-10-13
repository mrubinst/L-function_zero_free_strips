"""Agent t1, new test T1-P: the character facts used in the proof of Proposition prop:b1-points that the earlier
test b1/t_points.py does not check separately.

  eq:d013 (3.11)  kappa = chi^2 = (./pi)_3,  rho = chi^3 = (./pi)_2  (the Legendre symbol mod p), kappa(-1) = 1
  eq:d014 (3.12)  #{y : y^2 = U} = 1 + rho(U),   #{x : x^3 = V} = 1 + kappa(V) + conj(kappa(V))
  eq:d016 (3.15)  J(kappa,kappa) = sum_a kappa(a)(1 + rho(1-4a)) = sum_a kappa(a) rho(1-4a) = conj(kappa(4)) J(kappa,rho),
                  including #{t : t^2 - t + a = 0} = 1 + rho(1-4a)
  eq:b1-S (3.16)  S = sum_V rho(V+d) kappa(V) = -rho(d) kappa(4d) pi  and  a_p(d) = Tr(rho(d) kappa(4d) pi), checked
                  directly (b1/t_points.py checks it only through the chain S = rho(d)kappa(d)J(kappa,rho), J(kappa,kappa) = -pi)
                  for 0 < |d| <= 30, p not dividing d, a_p(d) = -sum_x ((x^3+d)/p) by brute force

chi = chi_pi is the sextic symbol of a primary prime pi of norm p = 1 mod 3, computed by Euler's criterion
(tests_cm/part_A_lines_0001_0549/cmcore.sym_table).  All sums are exact in Z[omega]: a sum of sixth roots of unity is
stored as the integer pair sum_j c_j * zeta^j with zeta^j = UNITS[j].
Range: every prime 7 <= p < 1200 with p = 1 mod 3, both primary primes above p, every residue class mod p.
"""
import math
import numpy as np
from t1lib import check
from cmcore import sym_table, primes_above, primes_upto, UNITS

PMAX = 1200
U = np.array(UNITS, dtype=np.int64)            # zeta^j as (A, B)


def zsum(idx, coef=None):
    """Exact sum of zeta^idx[k] (times integer coef[k]) over k with idx[k] >= 0, as an (A, B) pair."""
    idx = np.asarray(idx)
    c = np.ones(len(idx), dtype=np.int64) if coef is None else np.asarray(coef, dtype=np.int64)
    m = idx >= 0
    cnt = np.bincount(idx[m] % 6, weights=c[m], minlength=6).astype(np.int64)
    return (int(cnt @ U[:, 0]), int(cnt @ U[:, 1]))


def zmul(x, y):
    a, b = x; c, d = y
    return (a*c - b*d, a*d + b*c - b*d)


def zconj(x):
    a, b = x
    return (a - b, -b)


cases = 0
bad = dict(kap=0, rho=0, kneg=0, sq=0, cube=0, quad=0, J1=0, J2=0, J3=0, S=0, ap=0)
nS = 0
falS = dict(plus=0, kd=0)
fal = dict(rho_is_chi=0, kap_is_chi4=0, cube_is_1rho=0, cube_with_chi=0, disc_plus=0, disc_2a=0, kap4_unconj=0,
           kap2=0)
for p in primes_upto(PMAX):
    if p < 7 or p % 3 != 1:
        continue
    t = np.arange(p, dtype=np.int64)
    sqcount = np.bincount((t*t) % p, minlength=p)              # #{y : y^2 = U}
    cubecount = np.bincount((t*t % p)*t % p, minlength=p)      # #{x : x^3 = V}
    # #{t : t^2 - t + a = 0}, i.e. a = t - t^2 = t(1-t)
    tt = (t*(1 - t)) % p
    quadcount = np.bincount(tt, minlength=p)
    leg = np.array([0] + [1 if pow(int(x), (p - 1)//2, p) == 1 else -1 for x in range(1, p)])
    for pi in primes_above(p):
        cases += 1
        mod, w0, T, kind = sym_table(pi)
        assert kind == 'split' and mod == p
        chi = T[t]                                   # index j of chi(t), -1 at t = 0
        nz = chi >= 0
        kap = np.where(nz, (2*chi) % 6, -1)          # index of kappa = chi^2
        rho_idx = np.where(nz, (3*chi) % 6, -1)      # index of rho = chi^3, in {0, 3}
        rho = np.where(nz, np.where(rho_idx == 0, 1, -1), 0)
        img = (U[:, 0] + U[:, 1]*w0) % p             # image of zeta^j in O/pi = F_p
        # eq:d013: kappa is the cubic residue symbol (value in {1, omega, omega^2} congruent to t^((p-1)/3) mod pi),
        # rho is the Legendre symbol, kappa(-1) = 1
        cub = np.array([pow(int(x), (p - 1)//3, p) for x in t])
        bad['kap'] += int(np.sum(nz & ((kap % 2 != 0) | (img[kap % 6] != cub))))
        bad['rho'] += int(np.sum(rho != leg))
        bad['kneg'] += int(kap[p - 1] != 0)
        fal['rho_is_chi'] += int(np.any(np.where(nz, np.where(chi == 0, 1, np.where(chi == 3, -1, 99)), 0) != leg))
        kap4 = np.where(nz, (4*chi) % 6, -1)
        fal['kap_is_chi4'] += int(np.any(nz & (img[kap4 % 6] != cub)))
        # eq:d014: 1 + kappa + conj(kappa) is 3 (nonzero cube), 0 (non-cube) or 1 (V = 0), exactly in Z[omega]
        for V in range(p):
            k = kap[V]                                       # k in {0, 2, 4}, or -1 at V = 0
            kv = (0, 0) if k < 0 else (int(U[k, 0]), int(U[k, 1]))
            tot = (1 + kv[0] + zconj(kv)[0], kv[1] + zconj(kv)[1])
            if tot != (int(cubecount[V]), 0):
                bad['cube'] += 1
            if int(cubecount[V]) != 1 + rho[V]:
                fal['cube_is_1rho'] += 1
            c = chi[V]
            cv = (0, 0) if c < 0 else (int(U[c, 0]), int(U[c, 1]))
            if (1 + cv[0] + zconj(cv)[0], cv[1] + zconj(cv)[1]) != (int(cubecount[V]), 0):
                fal['cube_with_chi'] += 1
        bad['sq'] += int(np.sum(sqcount != 1 + rho))
        # eq:d016
        bad['quad'] += int(np.sum(quadcount != 1 + rho[(1 - 4*t) % p]))
        Jkk = zsum(np.where((kap >= 0) & (kap[(1 - t) % p] >= 0), (kap + kap[(1 - t) % p]) % 6, -1))
        mid1 = zsum(kap, 1 + rho[(1 - 4*t) % p])
        mid2 = zsum(kap, rho[(1 - 4*t) % p])
        Jkr = zsum(kap, rho[(1 - t) % p])
        k4 = (int(U[kap[4], 0]), int(U[kap[4], 1]))
        right = zmul(zconj(k4), Jkr)
        bad['J1'] += int(Jkk != mid1)
        bad['J2'] += int(mid1 != mid2)
        bad['J3'] += int(mid2 != right)
        # (1+4a in place of 1-4a is NOT a valid falsifier of the sum: a -> -a and kappa(-1) = 1 give the same value;
        #  it is a falsifier of the root count.)  Sum-level falsifier: discriminant 1-2a.
        fal['disc_plus'] += int(np.any(quadcount != 1 + rho[(1 + 4*t) % p]))
        fal['disc_2a'] += int(Jkk != zsum(kap, 1 + rho[(1 - 2*t) % p]))
        fal['kap4_unconj'] += int(Jkk != zmul(k4, Jkr))
        k2 = (int(U[kap[2], 0]), int(U[kap[2], 1]))
        fal['kap2'] += int(Jkk != zmul(zconj(k2), Jkr))
        # eq:b1-S
        for d in range(-30, 31):
            if d == 0 or d % p == 0:
                continue
            nS += 1
            S = zsum(kap, rho[(t + d) % p])
            c = zmul((int(rho[d % p]), 0), (int(U[kap[(4*d) % p], 0]), int(U[kap[(4*d) % p], 1])))
            target = zmul(c, pi)
            ap = -int(np.sum(leg[(t*t % p*t + d) % p]))
            bad['S'] += int(S != (-target[0], -target[1]))
            bad['ap'] += int(ap != 2*target[0] - target[1])          # Tr(A + B omega) = 2A - B
            falS['plus'] += int(S == target)
            kd = zmul((int(rho[d % p]), 0), (int(U[kap[d % p], 0]), int(U[kap[d % p], 1])))
            falS['kd'] += int(S == tuple(-x for x in zmul(kd, pi)))

check('T1-P1', 'eq:d013', '3.11', bad['kap'] == 0 and bad['rho'] == 0 and bad['kneg'] == 0,
      f"{cases} pairs (p, pi), 7<=p<{PMAX}, all residues: chi^2 = cubic symbol t^((p-1)/3) mod pi (values in mu_3), "
      f"chi^3 = Legendre (t/p), kappa(-1)=1; failures {bad['kap']}, {bad['rho']}, {bad['kneg']}; "
      f"falsifiers: chi itself is not Legendre in {fal['rho_is_chi']}/{cases}, chi^4 is not the cubic symbol in "
      f"{fal['kap_is_chi4']}/{cases} cases",
      falsifier=fal['rho_is_chi'] == cases and fal['kap_is_chi4'] == cases)
check('T1-P2', 'eq:d014', '3.12', bad['sq'] == 0 and bad['cube'] == 0,
      f"{cases} pairs (p, pi): #{{y:y^2=U}} = 1+rho(U) and #{{x:x^3=V}} = 1+kappa(V)+conj kappa(V) exactly for every "
      f"U, V in F_p; failures {bad['sq']}, {bad['cube']}; falsifiers: 1+rho(V) misses {fal['cube_is_1rho']} cube "
      f"counts, chi in place of kappa misses {fal['cube_with_chi']}",
      falsifier=fal['cube_is_1rho'] > 0 and fal['cube_with_chi'] > 0)
check('T1-P3', 'eq:d016', '3.15', all(bad[k] == 0 for k in ('quad', 'J1', 'J2', 'J3')),
      f"{cases} pairs (p, pi), exact in Z[omega]: #{{t:t^2-t+a=0}} = 1+rho(1-4a) (failures {bad['quad']}); "
      f"J(k,k) = sum k(a)(1+rho(1-4a)) ({bad['J1']}) = sum k(a)rho(1-4a) ({bad['J2']}) = conj k(4) J(k,rho) "
      f"({bad['J3']}); falsifiers fail in: root count 1+rho(1+4a) {fal['disc_plus']}, sum with 1-2a {fal['disc_2a']}, "
      f"k(4) unconjugated {fal['kap4_unconj']}, conj k(2) {fal['kap2']} of {cases}",
      falsifier=fal['disc_plus'] > 0 and fal['disc_2a'] > 0 and fal['kap4_unconj'] > 0 and fal['kap2'] > 0)
check('T1-P4', 'eq:b1-S', '3.16', bad['S'] == 0 and bad['ap'] == 0,
      f"{nS} triples (p, pi, d), 7<=p<{PMAX}, 0<|d|<=30, p not | d, exact in Z[omega]: S = -rho(d)kappa(4d)pi "
      f"(failures {bad['S']}) and brute-force a_p(d) = Tr(rho(d)kappa(4d)pi) (failures {bad['ap']}); falsifiers hold in: "
      f"+rho(d)kappa(4d)pi {falS['plus']}, -rho(d)kappa(d)pi {falS['kd']} of {nS}",
      falsifier=falS['plus'] == 0 and falS['kd'] < nS)
