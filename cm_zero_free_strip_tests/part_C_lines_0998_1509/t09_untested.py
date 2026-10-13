"""t09: claims of guide lines 998-1509 that are not finitely testable here, with the reason and the tests that
cover their finitely checkable ingredients.  No computation."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
from cmreport import untested

U = [
 ("U01", 1011, "thm:78", "A", "L(s, abar*eta) != 0 for Re s > 7/8, every finite-order eta (hence L(E_d,s) != 0 for Re s > 11/8)",
  "the theorem is labelled 'claimed'; it is conditional on the six statements of L1484-1496. Only finite-height "
  "numerics (C98) and the exponent bookkeeping (C51-C61) are tested."),
 ("U02", 1046, "section-by-section check", "C", "each proof of OAI78 uses only properties CM characters share; no new idea needed",
  "requires re-reading 16,677 lines of OAI78; not redone here (prior audits: notes78/audit_*.md)."),
 ("U03", 1071, "no poles", "A", "every L(s,psi), psi in F_CM, is entire",
  "Hecke's theorem; the mechanism (no constant term in the theta transform) is tested in C41 and C101."),
 ("U04", 1082, "continuation criterion", "C", "OAI78 Prop. lem:continuation-criterion holds verbatim with eta -> psi",
  "proof inspection; its three inputs are tested: |1/L^S| <= zeta_K(2) (C42), nonvanishing local factors (C05), beta_* (def.)."),
 ("U05", 1146, "Prop. cm-rows properties", "C", "H_{psi,u} converges in OAI78's regions, H << N(u)^eps, H_{psi,1} = 1+O(P_0^{-4/5})",
  "needs the full local Euler factor of OAI78 lem:local-euler; restatement checked (C70), formal substitution checked (C15, C16, C27)."),
 ("U06", 1151, "local table formal in psi(p)", "C", "every entry of OAI78's local table holds with eta(p) -> psi(p)",
  "only the entries quoted by the guide are tested (C16, C27, C73); the full six-row table J_0..J_5 is not re-derived."),
 ("U07", 1202, "Prop. cm-completed", "A", "sum_{N(m)<=Q}|B_m|^2 << Z^eps(Q+Q^2/Z) for the CM weight abar^2",
  "an analytic estimate (needs the CM transformation formula for d_zbar^2 theta at general moduli: item 1 of L1486); "
  "exponent algebra C48, sanity numbers C49, reflection table C28-C36."),
 ("U08", 1205, "marked lemmas (CM)", "A", "CM forms of lem:hybrid-energy, lem:reflected-energy, lem:probe-row-norm",
  "stated by the guide as not written out (item 1)."),
 ("U09", 1236, "dual lengths, sieve", "A", "dual lengths unchanged; alpha(mu)^2 absorbed by the quadratic large sieve",
  "argument only; unimodularity of alpha(mu)^2 is trivial; no finite test of the sieve."),
 ("U10", 1246, "low estimates", "A", "|I_psi(Z^1/2,Z^1/2,Z)| << Z^{1/4+eps} and |I_mod(Z)| << Z^{3/16+eps} for the CM target",
  "analytic estimates; exponents checked in C51, C52, C81."),
 ("U11", 1286, "logarithmic control for F_CM", "A", "OAI78 lem:logarithmic-control holds for F_CM with the same exponents",
  "item 2 of L1491; its growth input is tested (C37-C40)."),
 ("U12", 1306, "detector with X_u", "A", "X_u, Lemma cm-growth, Q_psi << N(u) and non-principality are all the detector uses",
  "item 3 of L1492; X_u properties tested (C45, C46), conductor sizes (C102)."),
 ("U13", 1362, "eq:cm-inverse", "A", "sum_{N(u)<<Z^m}|M_u Q_u|^2 << Z^{m+eps} for co-oriented coefficients",
  "item 4 of L1493; class identities tested (C23-C25), sanity numbers without slots (C103)."),
 ("U14", 1379, "descent never sees alpha", "C", "the descent, masks, punctures and seminorm propagation never see alpha",
  "proof inspection of OAI78 L9209-12477; not redone."),
 ("U15", 1402, "CM plain fourth moment", "A", "sum_k |S_k(n1)S_k(n2)Q|^2 << Z^{M+eps} with alpha(l)^e on the plain variables",
  "item 5 of L1494; ingredients tested: dual profile (C43), conductor change only at lambda (C44), vanishing lattice main "
  "terms (C47), kappa = 3/4 (C55, C86)."),
 ("U16", 1439, "no CM gamma factor in E(d)", "A", "the four sources of E(d) do not involve the gamma factor of a CM L-function",
  "argument only (the gamma factor enters through fixed powers of the height); E(d) itself checked (C56, C91)."),
 ("U17", 1455, "beta_* <= 11/12 for F_CM", "A", "the balanced probe with the sextic row count gives beta_* <= 11/12 for F_CM",
  "item 6 of L1495; not proved in the notes (Part I treats one curve at a time)."),
]
for cid, line, label, cls, desc, reason in U:
    untested(cid, line, label, cls, desc, reason)
