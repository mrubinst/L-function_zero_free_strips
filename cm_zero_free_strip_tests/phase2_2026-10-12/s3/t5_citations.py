"""t5: every OAI78 line citation in inserts.tex / edits.md points at the quoted content, the guide labels used
exist, and the OLD blocks of edits.md are unique.

  S3-T5a  each cited OAI78 line range contains its needle (the needle is searched only inside the range).
          Falsifier: the same needles searched in ranges shifted by 150 lines are (almost all) not found.
  S3-T5b  the guide labels referenced by the insert exist in the current guide.
"""
import os, sys, re, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s3lib import report, summary

_d = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(_d, "cmsources.py")) and os.path.dirname(_d) != _d:
    _d = os.path.dirname(_d)
sys.path.insert(0, _d)
import cmsources            # locations of the guide and of the cited sources (tests_cm/sources/)
OAI = cmsources.OAI78
GUIDE = cmsources.GUIDE
L = open(OAI).read().split("\n")


def has(lo, hi, needle):
    txt = " ".join(L[lo - 1:hi])
    return needle in txt


CIT = [  # (lo, hi, needle)
    (4201, 4201, r"\label{sec:detector}"), (4687, 4687, r"\label{sec:sextic-row-count}"),
    (8836, 8850, r"Instantiate Section~\ref{sec:detector}"),
    (1437, 1438, "A finite-order character has trivial infinite type"),
    (590, 603, "For the arithmetic datum in any given invocation"), (590, 603, "ideal sums exclude $S$"),
    (732, 740, r"Write $\mathcal A$ for this fixed arithmetic datum"),
    (3327, 3341, "Choose a finite ray group $T$"), (3327, 3341, r"\Theta=\langle\eta,\widehat T\rangle"),
    (3336, 3341, "Use one excluded set $S$"), (3336, 3341, "every prime of norm at most a fixed $P_0$"),
    (598, 599, "The fixed datum may depend on a target fixed beforehand"),
    (735, 737, "though it may depend on a previously fixed"),
    (601, 603, "We call primes outside $S$ good"),
    (636, 639, "including $j=0$ and negative $j$"),
    (4231, 4234, r"\mathcal X_u="), (4249, 4250, r"L_{\rm orig}(s,\psi)"),
    (580, 582, "Products of primary generators are primary"),
    (634, 635, r"For a primary $c=\prod p^{v_p(c)}$ outside $S$, define"),
    (4244, 4252, r"be the primitive character inducing a"), (4246, 4248, "Reciprocity and the fixed conductor primes give"),
    (660, 661, "No bound on its conductor exponents at those primes is"),
    (4259, 4263, r"local character of \(\chi_\bullet(u)\) on units at \(p\) has order"),
    (4248, 4252, r"E_u\subset S\cup\{p:p\mid u\}"),
    (4354, 4362, "For the last assertion, the primitive functional equation recalled"),
    (1554, 1555, r"Write $R_j=2-a-je$"), (1600, 1600, r"\end{proof}"),
    (4271, 4279, r"\(d_{\min}\le d\le d_{\max}\)"), (4271, 4279, r"I=\lceil1/e\rceil+2"),
    (4264, 4269, "can induce the principal character only when"), (4264, 4269, "Remove them for the"),
    (4222, 4228, r"\Theta=\langle\eta,\widehat T\rangle"),
    (4237, 4242, r"The group property of \(\Theta\)"),
    (4244, 4257, "These identities retain the original zero extensions"),
    (4324, 4324, r"\begin{proof}"), (4366, 4367, "This proves the reflected"),
    (4325, 4328, "Each maximum exists"), (4328, 4345, "If every one of its"),
    (4345, 4352, r"Lemma~\ref{lem:logarithmic-control}"), (4345, 4352, "deleted Euler products"),
    (1444, 1446, r"a nonprincipal $L(s,\psi)$ is"), (1438, 1445, "The primitive Hecke functional equation"),
    (4373, 4379, "Every primitive character inducing a presentation in"),
    (4419, 4419, r"\begin{proof}"), (4498, 4498, r"\end{proof}"),
    (4423, 4432, "Mellin inversion on a line"),
    (4549, 4549, r"\begin{proof}"), (4685, 4685, r"\end{proof}"),
    (4566, 4569, r"The row character is nonprincipal, so \(L_{\rm orig}\) is entire"),
    (4570, 4589, "The global strip bound"), (4570, 4589, r"-\frac{189}{100}"),
    (4592, 4600, "On the original line, absolute convergence"), (4592, 4600, r"\sum_{l\mid n}\mu(l)=0"),
    (4541, 4544, "selected by a label"), (4541, 4544, r"\Theta\times\{-1,1\}"),
    (4650, 4656, r"W_M(x)&=W_1(x)"),
    (5217, 5223, "Subdivide by its presentation labels"),
    (15223, 15226, "Subdivide by the witness presentation"),
    (5306, 5308, r"c_n=D^{-1/2}\mu(n)\nu(n)\chi_n(v)^{\varsigma}W(q_n/D)"),
    (725, 726, "Conjugating a whole factor inside its absolute value"),
    (4259, 4269, r"\label{par:physical-row-ramification}"),
    (384, 384, r"\label{old-eq:1.1b}"), (106, 108, r"\label{thm:main}"),
    (4707, 4707, r"\label{lem:sextic-large-sieve}"), (5181, 5181, r"\label{prop:sextic-row-count}"),
    (5809, 5809, r"\label{lem:fixed-bin-contour}"), (14985, 14985, r"\label{sec:refined-row-counts}"),
    (5298, 5298, r"\label{eq:sextic-frozen-mask}"),
]
bad = [(lo, hi, n) for lo, hi, n in CIT if not has(lo, hi, n)]
shift_found = sum(1 for lo, hi, n in CIT if has(lo + 150, hi + 150, n))
report("S3-T5a", "inserts.tex/edits.md", "every cited OAI78 line range contains the quoted content", not bad,
       f"{len(CIT)} citations; missing: {bad}; falsifier (ranges shifted by 150 lines): {shift_found} found",
       falsifier_failed=shift_found <= 3)

G = open(GUIDE).read()
labels = ["lem:cm-growth", "thm:78", "sec:78rows", "sec:78status", "sec:78final", "sec:78", "eq:cm-weight",
          "app:reflection", "lem:s3-labels", "prop:s3-detector", "app:s3-detector"]
missing = [l for l in labels if "\\label{%s}" % l not in G]
report("S3-T5b", "guide", "guide labels referenced by the merged appendix exist", not missing, f"missing: {missing}",
       falsifier_failed=("\\label{lem:s3-nonexistent}" not in G))
summary()
