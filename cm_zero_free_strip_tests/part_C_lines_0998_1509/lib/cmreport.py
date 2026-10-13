"""Reporting helpers for the Part C audit (guide lines 998-1509).

Every claim is reported on one line:

    [STATUS] <claim id> L<line> <label> {<class>} <description> :: <details>

STATUS is one of
    PASS      the claim was tested and holds on the tested range, AND its falsifier failed;
    FAIL      the claim failed, or the falsifier unexpectedly passed (test has no power);
    CITED     a restatement of OAI78/G78 checked against the source text (OK);
    CITED-X   a restatement that does NOT match the source text;
    INFO      sanity numbers for an asymptotic bound with an unspecified constant (no criterion);
    UNTESTED  not finitely testable here (reason given).

Paths: the OAI78 source and the G78 guide are read from tests_cm/sources/ (see cmsources.py).
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
_d = os.path.dirname(os.path.abspath(__file__))
while not os.path.exists(os.path.join(_d, "cmsources.py")) and os.path.dirname(_d) != _d:
    _d = os.path.dirname(_d)
sys.path.insert(0, _d)
import cmsources            # locations of the guide and of the cited sources (tests_cm/sources/)
OAI = cmsources.OAI78
G78 = cmsources.G78
# Part C refers to the guide of 2026-10-09 (line numbers L998-L1509); its frozen copy is read, not the live guide, whose
# status section no longer has the "Six statements must still be written out" list that C105 counts (changed 2026-10-13).
GUIDE = os.path.join(ROOT, "..", "cm-eleven-twelfths-guide_2026-10-09.tex")

_results = []

# Exit status (added 2026-10-12): a script that reports a failed check (or a falsifier that did not fail) through
# the helpers below exits with status 1 (atexit hook).  Self-test: with CM_TESTS_FORCE_FAIL=1 in the environment the
# first check reported here counts as failed.
import atexit
_FORCE_FAIL = [os.environ.get("CM_TESTS_FORCE_FAIL") == "1"]
_NFAILED = [0]

def _forced():
    """True exactly once (for the first check) when CM_TESTS_FORCE_FAIL=1."""
    if _FORCE_FAIL[0]:
        _FORCE_FAIL[0] = False
        return True
    return False

def _exit_status_hook():
    if _NFAILED[0]:
        try:
            print(f"EXIT STATUS 1: {_NFAILED[0]} check(s) failed", file=sys.stderr)
            sys.stdout.flush(); sys.stderr.flush()
        finally:
            os._exit(1)

atexit.register(_exit_status_hook)


def claim(cid, line, label, cls, desc, ok, details="", falsifier_failed=None, status=None):
    """Report a tested claim.  ok: claim held.  falsifier_failed: True if the deliberately wrong variant failed
    (None means no falsifier applies, e.g. CITED/INFO/UNTESTED lines, which must pass status explicitly)."""
    if status is None:
        if falsifier_failed is None:
            raise ValueError("tested claims need a falsifier")
        status = "PASS" if (ok and falsifier_failed) else "FAIL"
        if ok and not falsifier_failed:
            details += " [falsifier did NOT fail: test has no power]"
        elif falsifier_failed:
            details += " [falsifier failed as required]"
    if _forced():   # self-test: the first reported line of any kind becomes a FAIL
        status = "FAIL"; details = "CM_TESTS_FORCE_FAIL=1 (forced failure) " + details
    if status in ("FAIL", "CITED-X"):
        _NFAILED[0] += 1
    s = f"[{status}] {cid} L{line} {label} {{{cls}}} {desc} :: {details}"
    print(s, flush=True)
    _results.append((status, cid))
    return status


def cited(cid, line, label, desc, ok, details=""):
    return claim(cid, line, label, "C", desc, ok, details, status=("CITED" if ok else "CITED-X"))


def info(cid, line, label, cls, desc, details=""):
    return claim(cid, line, label, cls, desc, True, details, status="INFO")


def untested(cid, line, label, cls, desc, reason):
    return claim(cid, line, label, cls, desc, True, "reason: " + reason, status="UNTESTED")


def read_lines(path):
    with open(path) as f:
        return f.read().split("\n")


def src_has(path, needle, lo=None, hi=None):
    """True if the string needle occurs in lines lo..hi (1-based, inclusive) of path; returns (bool, first line)."""
    L = read_lines(path)
    lo = 1 if lo is None else lo
    hi = len(L) if hi is None else hi
    for i in range(lo - 1, min(hi, len(L))):
        if needle in L[i]:
            return True, i + 1
    # allow needles that span a line break: search the joined text
    txt = "\n".join(L[lo - 1:hi])
    j = txt.replace("\n", " ").find(needle)
    if j >= 0:
        return True, lo + txt[:j].count("\n")
    return False, None
