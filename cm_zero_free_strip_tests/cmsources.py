"""Locations of the files, other than the test scripts, that some checks read.

The tests compare displays with the text of the CM notes and with their sources, so a few checks read TeX files:
  GUIDE      the CM notes, cm_zero_free_strip.tex (the version tested), and GUIDE_AUX, its .aux file (equation numbers);
  OAI78      OpenAI, "The Quasi-Riemann Hypothesis" (September 30, 2026), build/paper.tex;
  OAI1112    OpenAI, "The Quasi-Riemann Hypothesis" (October 5, 2026), build/paper2.tex;
  ZG, G78    our companion guides to the 11/12 and 7/8 arguments for zeta(s);
  INVENTORY  the list of the numbered displays of GUIDE (one JSON file per test group, and inventory.json for all).
All of them are in the folder sources/ next to this file (see sources/README.md).  Any one can be replaced by setting
the environment variable named in brackets below to a path.

Import from any test script with
    sys.path.insert(0, <path of the tests folder>); from cmsources import GUIDE, OAI78, ...
or call require(path, ...) to stop with a clear message when a file is missing.
"""
import os
import sys

TESTS = os.path.dirname(os.path.abspath(__file__))
SOURCES = os.environ.get("CM_SOURCES") or os.path.join(TESTS, "sources")


def _path(env, *rel):
    return os.environ.get(env) or os.path.join(SOURCES, *rel)


GUIDE = _path("CM_GUIDE_TEX", "cm_zero_free_strip.tex")                              # [CM_GUIDE_TEX]
GUIDE_AUX = _path("CM_GUIDE_AUX", "cm_zero_free_strip.aux")                          # [CM_GUIDE_AUX]
GUIDE_OCT9 = os.path.join(TESTS, "cm-eleven-twelfths-guide_2026-10-09.tex")          # frozen version of October 9
OAI78 = _path("CM_OAI78_TEX", "openai", "qrh-september-30-2026-paper.tex")           # [CM_OAI78_TEX]
OAI1112 = _path("CM_OAI1112_TEX", "openai", "qrh-october-5-2026-paper2.tex")         # [CM_OAI1112_TEX]
ZG = _path("CM_ZG_TEX", "zeta_11_12_guide.tex")                                      # [CM_ZG_TEX]
G78 = _path("CM_G78_TEX", "zeta_7_8_guide.tex")                                      # [CM_G78_TEX]
INVENTORY = _path("CM_INVENTORY_DIR", "inventory")                                   # [CM_INVENTORY_DIR]


def inventory(name="inventory.json"):
    """Path of a display inventory, e.g. inventory('inventory_t4.json')."""
    return os.path.join(INVENTORY, name)


def require(*paths):
    """Exit with status 2 and a clear message if any of the given files is missing."""
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        print("[FAIL] missing source file(s): " + ", ".join(missing))
        print("       See sources/README.md in the tests folder, or set the environment variable named in cmsources.py.")
        sys.exit(2)
    return paths[0] if len(paths) == 1 else paths
