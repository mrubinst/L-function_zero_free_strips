"""t8_tex_compare.py -- symbol-by-symbol comparison of displays of Appendices E and G of cm-eleven-twelfths-guide.tex
with their sources (OAI78 = sources/openai/qrh-september-30-2026-paper.tex, OAI1112 = sources/openai/qrh-october-5-2026-
paper2.tex) and with the guide's own displays that they restate.

Each comparison normalizes both TeX strings (spacing commands, \\Bigl/\\Bigr, alignment marks, \\frac{A}{B} -> (A)/(B) with the parentheses
dropped around a single atom such as q_n or Z^r,
\\N(x) and \\NK(x) -> q_x, \\bar\\alpha(x) -> \\overline{\\alpha(x)}, \\sum_{x\\ {\\rm sf}} -> \\ssum_x, 1_{ -> \\ind_{,
\\varepsilon -> \\epsilon, \\mathcal O and \\OO -> O, \\muk -> \\mu, final punctuation), then applies to the source
an explicit list of CM edits (the changes the guide says it makes) and requires the two strings to be equal.
A falsifier for each comparison: the source WITHOUT the CM edits must differ from the guide (where there are edits).

Usage: python3 t8_tex_compare.py [guide.tex] [OAI78 paper.tex] [OAI1112 paper2.tex]   (defaults: the project paths)
Exit status 1 if a comparison fails.
"""
import os, sys, re
HERE = os.path.dirname(os.path.abspath(__file__))
_d = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(_d, "cmsources.py")) and os.path.dirname(_d) != _d:
    _d = os.path.dirname(_d)
sys.path.insert(0, _d)
import cmsources            # locations of the guide and of the cited sources (tests_cm/sources/)
GUIDE = sys.argv[1] if len(sys.argv) > 1 else cmsources.GUIDE
O78 = sys.argv[2] if len(sys.argv) > 2 else cmsources.OAI78
O1112 = sys.argv[3] if len(sys.argv) > 3 else cmsources.OAI1112
TXT = {k: open(p).read() for k, p in (("G", GUIDE), ("O78", O78), ("O1112", O1112))}
LINES = {k: v.split('\n') for k, v in TXT.items()}

def env_of_label(src, label):
    """The display environment containing \\label{label} (equation/align), without begin/end and labels."""
    t = TXT[src]; i = t.index('\\label{' + label + '}')
    b = max(t.rfind('\\begin{equation}', 0, i), t.rfind('\\begin{align}', 0, i))
    e = min(x for x in (t.find('\\end{equation}', i), t.find('\\end{align}', i)) if x >= 0)
    return t[b:e]

def lines(src, a, b):
    return '\n'.join(LINES[src][a - 1:b])

def frac_to_slash(s):
    out = ''; i = 0
    while True:
        j = s.find('\\frac', i)
        if j < 0: return out + s[i:]
        out += s[i:j]; k = j + 5
        def arg(k):
            if s[k] == '{':
                d = 0
                for m in range(k, len(s)):
                    d += (s[m] == '{') - (s[m] == '}')
                    if d == 0: return s[k + 1:m], m + 1
            if s[k] == '\\':
                m = k + 1
                while m < len(s) and s[m].isalpha(): m += 1
                return s[k:m], m
            return s[k], k + 1
        A, k = arg(k); B, k = arg(k)
        out += '(' + frac_to_slash(A) + ')/(' + frac_to_slash(B) + ')'; i = k

def norm(s):
    s = re.sub(r'%.*', '', s)
    s = re.sub(r'\\label\{[^}]*\}|\\tag\{[^}]*\}|\\begin\{[a-z*]*\}|\\end\{[a-z*]*\}|\\hspace\{[^}]*\}', '', s)
    s = frac_to_slash(s)
    for a, b in ((r'\Bigl', ''), (r'\Bigr', ''), (r'\Biggl', ''), (r'\Biggr', ''), (r'\biggl', ''), (r'\biggr', ''),
                 (r'\left', ''), (r'\right', ''), ('{}&', ''), ('&', ''), (r'\\', ' '), (r'\qquad', ' '),
                 (r'\quad', ' '), (r'\,', ''), (r'\;', ''), (r'\!', ''), (r'\times', ''), (r'\varepsilon', r'\epsilon'),
                 (r'\mathcal O', 'O'), (r'\OO', 'O'), (r'\muk', r'\mu'), (r'\NK(', r'\N('), (r'\mathcal ', r'\mathcal')):
        s = s.replace(a, b)
    s = re.sub(r'\\N\((\w)\)', r'q_\1', s)
    s = re.sub(r'\\bar\\alpha\((\w)\)', r'\\overline{\\alpha(\1)}', s)
    s = re.sub(r'\\sum_\{(\w)\\ \{\\rm sf\}\}', r'\\ssum_\1', s)
    s = re.sub(r'(?<![\\\w])1_\{', r'\\ind_{', s)
    s = re.sub(r'\s+', '', s)
    # drop the parentheses that \frac{A}{B} -> (A)/(B) put around a single atom (letter, optional subscript and power),
    # so that \frac{\N(n)}{Z^r} and OAI78's q_n/Z^r agree; products keep their parentheses
    br = r'\{(?:[^{}()]|\{[^{}()]*\})*\}'                       # a brace group, one nesting level
    atom = r'[A-Za-z](?:_\w|_' + br + r')?(?:\^\w|\^' + br + r')?'
    s = re.sub(r'\((' + atom + r')\)(?=/)', r'\1', s)
    s = re.sub(r'(?<=/)\((' + atom + r')\)', r'\1', s)
    return s.rstrip('.,;')

RESULTS = []
def compare(name, guide_tex, src_tex, edits, extra_src_edits=()):
    """edits: CM edits (old, new) applied to the normalized source; extra_src_edits: notational edits (not CM)."""
    g = norm(guide_tex); s0 = norm(src_tex)
    for a, b in extra_src_edits:
        assert a in s0, (name, 'notation', a, s0)
        s0 = s0.replace(a, b)
    s = s0
    for a, b in edits:
        assert a in s, (name, a, s)
        s = s.replace(a, b)
    ok = (g == s)
    fals = (g != s0) if edits else None
    RESULTS.append((name, ok, fals))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {len(edits)} CM edit(s)"
          + ("" if fals is None else f"; source without the edits differs from the guide: {'yes (falsifier fails as required)' if fals else 'NO'}"))
    if not ok:
        print("   guide :", g); print("   source:", s)

A = r'\overline{\alpha(%s)}'
# ---------------------------------------------------------------- Appendix E
compare("eq:s4-M vs OAI78 L9233 (M_u)", env_of_label("G", "eq:s4-M"), lines("O78", 9233, 9233),
        [(r'\mu(n)\psi_u(n)', r'\mu(n)\alpha(n)^{\epsilon}\psi_u(n)')])
# eq:s4-M and eq:d171 restate the two lines of eq:d199 (Section sec:78moments), with psi_u = nu chi_.(u)^eps
d199 = norm(env_of_label("G", "eq:d199"))
gM = norm(env_of_label("G", "eq:s4-M")).replace(r'M_u(Z^r;W)', 'M_u').replace(r'\psi_u(n)', r'\nu(n)\chi_n(u)^{\epsilon}')
gQ = norm(env_of_label("G", "eq:d171")).replace(r'\psi_u(p)', r'\nu(p)\chi_p(u)^{\epsilon}').replace(r'\prod_{i\inI}', r'\prod_{i}') \
    .split(r'z=\sum')[0].rstrip(',')
ok = d199.startswith(gM + ',') and (gQ in d199)
RESULTS.append(("eq:s4-M, eq:d171 = eq:d199", ok, None))
print(f"[{'PASS' if ok else 'FAIL'}] eq:s4-M and eq:d171 are the two lines of eq:d199 (psi_u = nu chi_.(u)^eps written out)")
if not ok: print("   d199:", d199); print("   M   :", gM); print("   Q   :", gQ)
compare("eq:d171 vs OAI78 L9240-9243 (Q_u)", env_of_label("G", "eq:d171"), lines("O78", 9241, 9243), [])
compare("eq:d172 vs OAI78 L9254 (hypotheses of lem:marked)", env_of_label("G", "eq:d172"), lines("O78", 9254, 9254), [])
compare("eq:s4-marked vs OAI78 old-eq:4.1", env_of_label("G", "eq:s4-marked"), env_of_label("O78", "old-eq:4.1"), [])
compare("eq:s4-a (first form) vs OAI78 L7293 (a(n))", env_of_label("G", "eq:s4-a").split(r'=\varphi(n)')[0],
        lines("O78", 7293, 7293).split(r',\qquad')[0],
        [(r'a(n)=' + A % 'n', r'a(n)=' + A % 'n' + r'^{2}')],
        )
# the second form of eq:s4-a: phi(n) times the OAI78 coefficient
g_a = norm(env_of_label("G", "eq:s4-a")).split('=')
s_a = norm(lines("O78", 7293, 7293).split(r',\qquad')[0]).split('=')[1]
ok = g_a[2] == r'\varphi(n)\cdot' + s_a
RESULTS.append(("eq:s4-a second form = phi(n) x OAI78 coefficient", ok, None))
print(f"[{'PASS' if ok else 'FAIL'}] eq:s4-a second form = \\varphi(n)\\cdot (OAI78 L7293 coefficient)")
compare("eq:s4-E vs OAI78 eq:canonical-energy", env_of_label("G", "eq:s4-E"), env_of_label("O78", "eq:canonical-energy"),
        [(A % 'n' + r'\gamma_2(n)\nu(n)\rho(n)', 'a(n)')])
compare("eq:s4-T vs eq:s1-marked-completion (guide)", env_of_label("G", "eq:s4-T"),
        env_of_label("G", "eq:s1-marked-completion"), [(r'\mathcalT^{(2)}_', r'\mathcalT_')])
compare("eq:s4-T vs OAI78 eq:marked-completion", env_of_label("G", "eq:s4-T"), env_of_label("O78", "eq:marked-completion"),
        [(A % 'n' + r'\gamma_2', A % 'n' + r'^{2}\gamma_2'), (A % 'b' + '^3', A % 'b' + '^{6}')],
        extra_src_edits=[('V_*(q_nq_b^3/X)', 'V_*((q_nq_b^3)/X)')])      # notation only: \frac vs slash
compare("eq:s4-inversion vs OAI78 eq:marked-completion-inversion", env_of_label("G", "eq:s4-inversion"),
        env_of_label("O78", "eq:marked-completion-inversion"), [(A % 'h' + '^3', A % 'h' + '^{6}')])
compare("eq:s4-Pa vs OAI78 eq:inverse-first-polynomial", env_of_label("G", "eq:s4-Pa"),
        env_of_label("O78", "eq:inverse-first-polynomial"), [(r'\mu(x)\nu_a(x)', r'\mu(x)' + A % 'x' + r'\nu_a(x)')])
compare("eq:s4-a0 vs OAI78 L10807 (a_0)", env_of_label("G", "eq:s4-a0"), lines("O78", 10807, 10807).split(r',\qquad')[0],
        [(A % 'z', A % 'z' + '^{2}')])
compare("eq:s4-child vs OAI78 eq:canonical-child-polynomial", env_of_label("G", "eq:s4-child"),
        env_of_label("O78", "eq:canonical-child-polynomial"), [(A % 'n', A % 'n' + '^{2}')])
compare("eq:s4-R0 vs OAI78 eq:initial-residual-polynomial", env_of_label("G", "eq:s4-R0"),
        env_of_label("O78", "eq:initial-residual-polynomial"), [(r'\mu(c)\nu_0(c)', r'\mu(c)' + A % 'c' + r'\nu_0(c)')])
compare("eq:s4-Astar vs OAI78 L11812-11813 (A_*, nu_*)", env_of_label("G", "eq:s4-Astar"),
        lines("O78", 11813, 11813).split(r',\qquad')[0] + r',\qquad' + lines("O78", 11812, 11812).split(r',\qquad')[0],
        [(A % 'z', A % 'z' + '^{2}')])
compare("eq:s4-sixth vs OAI78 eq:sixth-power-column-identity", env_of_label("G", "eq:s4-sixth"),
        env_of_label("O78", "eq:sixth-power-column-identity"), [(r'\mu(d)\psi_u(d)', r'\mu(d)\alpha(d)^{\epsilon}\psi_u(d)')])
# ---------------------------------------------------------------- Appendix G
g190 = norm(env_of_label("G", "eq:d190"))
s_high = norm(lines("O78", 6587, 6588)); s_low = norm(lines("O78", 6692, 6693).split('=')[0])
ok = (s_low in g190) and (s_high in g190) and (r'\omega=\frac{\Delta_1}{2}' in lines("O78", 6606, 6606)) \
     and (r'C_{\mathrm I}(s)=s-\frac23' in lines("O78", 5547, 5547))
RESULTS.append(("eq:d190 vs OAI78 L6587-6588, L6692, L6606, L5547", ok, None))
print(f"[{'PASS' if ok else 'FAIL'}] eq:d190: both estimates occur verbatim in OAI78 (L6587-6588, L6692); omega = Delta_1/2 (L6606);"
      " C_I(s) = s - 2/3 (L5547)")
compare("eq:d193 vs OAI1112 L673 (A_u)", env_of_label("G", "eq:d193"), lines("O1112", 673, 673),
        [(r'A_u(D):=\sum_{(n,S)=1}\mu(n)\nu(n)', r'A_u^\nu(D)=\sum_n\mu(n)\alpha(n)\nu(n)')])
g194 = norm(env_of_label("G", "eq:d194")); s194 = norm(env_of_label("O1112", "eq:ms"))
g194c = g194.replace(r'A_u^\nu', 'A_u').replace(r'\ll_{\psi,W,\vartheta,\epsilon}', r'\ll')
s194c = s194.replace(r'\ll_{\nu,S,I,\vartheta,\epsilon}(\max_{0\lej\lek}\|W^{(j)}\|_\infty)^2', r'\ll')
ok = g194c == s194c
RESULTS.append(("eq:d194 vs OAI1112 eq:ms", ok, None))
print(f"[{'PASS' if ok else 'FAIL'}] eq:d194 vs OAI1112 eq:ms (Prop. 3.1): same sum and exponent D^(2+theta+eps); the constant's "
      "dependence (finite C^k norm of W) is stated in words after eq:d194")
if not ok: print("   guide :", g194c); print("   source:", s194c)
# eq:d196 vs OAI1112 L742-746: same Mellin identity with alpha nu and the Euler product inserted
s196 = norm(lines("O1112", 742, 747)); g196 = norm(env_of_label("G", "eq:d196"))
ok = (r'\widehatW(s)\sum_{(n,S)=1}(\mu(n)\nu(n))/q_n^s=(\widehatW(s))/(L_K^S(s,\nu))' in s196.replace(r'\mathcalM_W(s)=', '')) \
     and g196.startswith(r'F_W(s)=\widehatW(s)\sum_n(\mu(n)\alpha(n)\nu(n))/q_n^s=') \
     and g196.endswith(r'=(\widehatW(s))/(L^S(s,\overline\psi))')
RESULTS.append(("eq:d196 vs OAI1112 L742-747", ok, None))
print(f"[{'PASS' if ok else 'FAIL'}] eq:d196 vs OAI1112 L742-747: same Mellin identity, with alpha(n)nu(n) and L^S(s, conj psi) "
      "for mu(n)nu(n) and L_K^S(s, nu); the Euler product is the middle member")
if not ok: print("   guide :", g196); print("   source:", s196)
# eq:d197 vs eq:s4-M of the guide with nu = theta, Z^r = D (and psi_u(n) = theta(n) chi_n(u)^eps)
g197 = norm(env_of_label("G", "eq:d197")).split(r'\epsilon=\pm1')[0].rstrip('.,;')
g4M = norm(env_of_label("G", "eq:s4-M")).replace(r'M_u(Z^r;W)=Z^{-r/2}', r'M_u(D;W)=D^{-1/2}') \
    .replace(r'\psi_u(n)', r'\theta(n)\chi_n(u)^{\epsilon}').replace('q_n/Z^r', 'q_n/D')
ok = g197 == g4M
RESULTS.append(("eq:d197 = eq:s4-M with nu = theta, Z^r = D", ok, None))
print(f"[{'PASS' if ok else 'FAIL'}] eq:d197 is eq:s4-M with nu = theta in Theta_0, Z^r = D, psi_u = theta chi_.(u)^eps")
if not ok: print("   d197:", g197); print("   s4-M:", g4M)

print()
nok = sum(1 for _, ok, f in RESULTS if ok and f is not False)
print(f"SUMMARY: {nok}/{len(RESULTS)} comparisons pass" + ("" if nok == len(RESULTS) else "  <-- PROBLEM"))
sys.exit(0 if nok == len(RESULTS) else 1)
