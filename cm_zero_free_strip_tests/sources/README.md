# Source files read by the tests

Some checks compare displays of the CM notes with the text of the notes and of the papers they cite, so they read
the TeX files below. `../cmsources.py` gives their locations; any one can be replaced through the environment
variable named there.

| File | What it is |
|---|---|
| `cm_zero_free_strip.tex` | The CM notes, in the version the tests were run on. It is identical to `../../cm_zero_free_strip.tex` in the repository at the time of posting. |
| `cm_zero_free_strip.aux` | Its LaTeX `.aux` file, from which the equation numbers are read. |
| `openai/qrh-september-30-2026-paper.tex` | OpenAI, *The Quasi-Riemann Hypothesis* (September 30, 2026), the 7/8 paper (cited as OAI78). |
| `openai/qrh-october-5-2026-paper2.tex` | OpenAI, *The Quasi-Riemann Hypothesis* (October 5, 2026), the 11/12 paper (cited as OAI1112). |
| `openai/LICENSE` | The Apache License 2.0, under which the two OpenAI files are distributed. |
| `zeta_11_12_guide.tex` | Source of our guide to the 11/12 argument for ζ(s) (cited as ZG). |
| `zeta_7_8_guide.tex` | Source of our guide to the 7/8 argument for ζ(s) (cited as G78). |
| `inventory/inventory*.json` | The numbered displays of the notes: label, number, source text and context, one file per test group t1–t9, plus `inventory.json` for all of them. |

The two OpenAI files are unmodified copies of `preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex`
and `preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex` from
<https://github.com/openai/math>, commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a`, copyright OpenAI, distributed
under the Apache License 2.0 (`openai/LICENSE`). Line numbers cited in the tests (L1234) refer to these files.
