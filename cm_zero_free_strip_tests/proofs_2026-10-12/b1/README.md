# Block b1 — testable claims introduced in b1.tex / b1_edits.md

Scripts (copies kept in `compile_b1/tests/`; run with `python3 -I <script>` and `gp -q t_bd.gp`):
- `t_points.py` — every step of the proof of Proposition `prop:b1-points`, exactly in Z[ω] (needs `zw.py` from
  `tests_cm/part_A_lines_0001_0549/`, imported by absolute path).
- `t_bd.gp` — local factors at p | 6d (Lemma `lem:b1-euler`(c)), PARI/GP.
- `t_extract.py` — counts and inequalities of Lemmas `lem:b1-extract`, `lem:b1-uniform` (numpy, sympy).
- `apply_b1.py` — applies b1_edits.md and b1.tex to `compile_b1/guide_orig.tex`, producing `guide_b1.tex`.

Status: **all tests below were run on 2026-10-12 and passed; every falsifier failed as required.**

| ID | Where | Claim | Class | Test (falsifier) | Range | Result | Existing test |
|---|---|---|---|---|---|---|---|
| b1.01 | prop:b1-points (b), eq:b1-count | #E_d(F_p) − 1 = Σ_V (1+ρ(V+d))(1+κ(V)+κ̄(V)) = p + S + S̄ | E | exact in Z[ω] vs brute-force count | 4412 pairs (d,π): 0<\|d\|≤30, 7≤p<400, p≡1 (3), both π | PASS | — |
| b1.02 | prop:b1-points (b) | S = ρ(d)κ(d)J(κ,ρ) (substitution V = −dt, κ(−1) = 1) | E | exact | same | PASS | — |
| b1.03 | prop:b1-points (b) | J(κ,κ) = conj(κ(4)) J(κ,ρ) (discriminant 1−4a) | E | exact | same | PASS | A3.14 (414 cases, form J(κ,ρ)=κ(4)J(κ,κ)) |
| b1.04 | prop:b1-points (b) | J(κ,κ) = −π for primary π (input [QF, Thm 15]) | E/C | exact (fals.: J = +π never holds: 0 hits) | same | PASS | an earlier independent check; QF numerics |
| b1.05 | prop:b1-points (b), eq:b1-S | ρ(d)κ(4d) = conj χ_π(4d) | E | exact (fals.: ρ(d)κ(4d) = χ_π(4d) for non-real values: 0 hits) | same | PASS | A3.07 (related) |
| b1.06 | prop:b1-points (b) | a_p(d) = Tr(π conj χ_π(4d)), both π | E | exact (fals.: unconjugated symbol agrees in only 1536/4412 cases; χ_π(d) in 1320/4412) | same | PASS | A3.16 (27 528 pairs), an earlier independent check (24 347 pairs) |
| b1.07 | prop:b1-points | \|a_p\| ≤ 2√p, i.e. \|λ_d(p)\| ≤ 2 | I | exact integer check a_p² ≤ 4p | same | PASS | A3.25 |
| b1.08 | prop:b1-points (a) | a_p = 0 for p ≡ 2 (3), p∤d | E | — | — | — | A2.11, A3.21 |
| b1.09 | INSERT 1, first paragraph | j(y²=x³+Ax+B) = 1728·4A³/(4A³+27B²); j = 0 ⇔ A = 0 | E | suggested: sympy from c4, Δ (fals.: 1728·4A³/(4A³−27B²)) | symbolic | not run (standard) | A0.01 (j = 0 for E_d) |
| b1.10 | lem:b1-euler (a) | split/inert Euler factors equal 1 − λX + X²; ψ_d(−p) = −1; ψ_d = 0 above p \| d | E | — | — | — | A3.20, A3.21, A3.22 |
| b1.11 | lem:b1-euler (b), eq:b1-euler | the five-term chain for Re s>1; (3.4); a_m(E_d)/√m = Σ_{N n=m} ψ_d(n) for (m,6d)=1 | E | — | — | — | A2.07, A3.08, A3.24, B04; an earlier independent check |
| b1.12 | lem:b1-euler (c) | at p \| 6d: bad ⇒ a_p ∈ {0,±1} (here a_p = 0); good ⇒ \|a_p\| ≤ 2√p; roots of L_p(E_d,s+½)^{-1} in X = p^{-s} have \|X\| ≥ 1 | E/I | PARI on minimal models (fals.: unshifted factor 1 − a_pX + pX² has a root inside the unit disc at every good p) | 20 918 pairs (d, p\|6d): \|d\| ≤ 3000, d = 5⁶k, 7⁶k, 2⁶3⁶5⁶k; 96 good p∈{2,3}, 58 good p≥5 with p\|d | PASS | A2.08, B07 |
| b1.13 | lem:b1-euler (c), remark | p ≥ 5, p \| d: good ⇔ 6 \| v_p(d), and then a_p(E_d) = a_p(E_{d₁}), d₁ = d p^{−v_p(d)}; p ∈ {2,3} good ⇒ \|a_p\| ≤ p | E/I | PARI (same run) | same | PASS | A5.15 (E_{dt⁶} ≅ E_d) |
| b1.14 | lem:b1-extract (i) | #{m ∈ Z[ω]: N(m) ≤ x} ≤ 1 + 8x for all x ≥ 0 | I | exact counts at every norm value ≤ 1.2·10⁵ (fals.: 1 + 3x is violated) | x ≤ 0.75·400² | PASS | A5.06–A5.07 (sizes only) |
| b1.15 | lem:b1-extract (iii) | #{n primary, (n,6)=1, π \| n, N(n) ≤ X} ≤ 1 + 8X/N(π) | I | exact, X = 3000, all primary primes of norm ≤ 200 | — | PASS | A5.06 |
| b1.16 | lem:b1-extract (ii) | rows u₀π⁶ distinct; N ≤ H | E | — | — | — | A5.04 |
| b1.17 | lem:b1-extract (iii) | Y ≤ H^{1/6} = D^{(1+ϑ)/6} ≤ D for D ≥ 1, ϑ ≤ 1/10, N₀ ≥ 1 | I | exact fraction | — | PASS | — |
| b1.18 | lem:b1-extract (iv) | J\|S\|² ≤ 2M + 2(1+16b)²‖W‖²JD²/Y² | I | suggested: rerun the toy extraction of A5.09 with the explicit constant | — | not run | A5.08, A5.09 (with O(·)) |
| b1.19 | proof of (5.4) | J ≥ Y/(5 log Y) for Y ≥ Y₀(d) | A | sanity: P(Y) ≥ 0.4·Y/log Y for all integers 97 ≤ Y ≤ 10⁶ (P = primary primes ∤6 with Y/2 < N ≤ Y); min of P/(Y/log Y) on [10³,10⁶] is 0.497 | Y ≤ 10⁶ | SANITY | A5.03 |
| b1.20 | lem:b1-uniform | no primary prime coprime to 6 has norm < 7, so P(Y) = J = 0 for Y < 7; P(7) = 2 | E | exact | — | PASS | A5.03 note |
| b1.21 | lem:b1-uniform | #{π \| d primary prime: N(π) > Y/2} < 2 log\|d\| / log(Y/2) for Y > 2 | I | 20 000 random (d,Y), \|d\| ≤ 10⁷ (fals.: without the factor 2, violated at d = p split, Y = 2p − 1) | — | PASS | — |
| b1.22 | lem:b1-uniform | log\|d\| < (1/10) log N₀ ≤ (1/10) log H when Y ≥ 1 | I | immediate from N₀ = \|4d\|¹⁰ and N₀Y⁶ = H | — | (by inspection) | A5.01 |
| b1.23 | lem:b1-uniform, Case 2 | Y ≥ Y₁ ≥ 7 and c₁Y/log Y < log H/(5 log(Y/2)) ⇒ Y < 2 log H/(5c₁); log(Y/2) ≥ ½ log Y for Y ≥ 4 | I | 2·10⁵ random (c₁, Y, log H); grid for Y ∈ [4,10⁶] | — | PASS | — |
| b1.24 | lem:b1-uniform | with Y⁶ = D^{1+ϑ}/N₀: D^{1+ϑ/2}Y^{−1/2} = N₀^{1/12}D^{11/12+5ϑ/12}, D/Y = N₀^{1/6}D^{5/6−ϑ/6} | I | sympy (fals.: N₀^{1/6} in the first term) | — | PASS | A5.12 |
| b1.25 | lem:b1-uniform | the uniform bound holds for all D ≥ 1 | A | suggested: for d ∈ {1, 2, 7, 5⁶} and D ≤ 10⁵, compare \|S_d(D)\| with N₀^{1/12}D^{11/12+5ϑ/12} (sizes only) | — | not run | A2.09 (sizes of S_d) |
| b1.26 | lem:b1-mellin (i) | \|q_d(m)\| ≤ r(m) = #{n primary, (n,6)=1, N n = m}; Σ r(m)m^{−σ} ≤ (1−3^{−σ})ζ_K(σ) | I | by inspection: triangle inequality in (3.4); both sides multiplicative with \|q_d(p)\| ≤ 2 = r(p) at split p (fals.: \|q_d(m)\| ≤ 1 fails at m = 7 for d = 1, \|q_1(7)\| = 4/√7) | — | PASS (by inspection) | A3.08 |
| b1.27 | lem:b1-mellin (i) | F_W(s) = Ŵ(s)Σ q_d(m)m^{−s} for Re s > 1 | E | — | — | — | B04 (7 values of d, 3 values of s) |
| b1.28 | lem:b1-mellin (ii) | F_W Ŵ′ = F_{W′} Ŵ on Re s > 1 | E | suggested: two weights, s = 2, 1.5 + 3i, via the series of b1.27 | — | not run (algebraic consequence of b1.27) | B04 |
| b1.29 | E9 | \|d\| ≤ 6⁵ N(E_d)^{5/2} for sixth-power-free d; f_p = 2 for p ≥ 5, p \| d | E/I | — | — | — | A5.16 (5979 values of d, max ratio 0.11) |

Not testable numerically: the existence of the absolute constants c₁, Y₁ (prime number theorem; only the sanity check
b1.19), the cited theorems (Hasse; Deuring; Hecke; [QF, Thm 15]), and the analytic statements of Lemma b1-mellin
beyond b1.27–b1.28.
