# Quadratic number fields, with Q(√−3) worked out in detail

Background notes for reading the 11/12 argument. The compiled PDF is
[`../quadratic_fields_notes.pdf`](../quadratic_fields_notes.pdf) (33 pages). Contents:

1. **Quadratic fields in general:** rings of integers, units, ideals and class groups, splitting of
   primes, ζ_K = ζ·L(s, χ_D).
2. **ℚ(√−3) and the Eisenstein integers ℤ[ω]:** complete elementary proofs, primary elements in the
   convention α ≡ 1 (mod 3), and ideals as lattices, with figures.
3. **Residue symbols:** quadratic, cubic and sextic symbols, cubic reciprocity, and cubic Gauss and
   Jacobi sums.
4. **Poisson summation over ℤ[ω]:** additive characters and duality, Gauss-sum identities, and Poisson
   summation with a character, as used in the 11/12 argument.
5. **The cubic theta function of Kubota and Patterson:** the properties that the 11/12 argument
   imports, and a sketch of how it uses them.
6. **Tables.**

The notes use λ = 1 − ω. OpenAI's paper and Dunn–Radziwiłł write λ for √−3 = 1 + 2ω = ωλ;
Sections 4 and 5 write √−3 explicitly.

## Files

| File | Contents |
|---|---|
| `quadratic-fields.tex` | LaTeX source. It inputs the generated tables and figures from `tables/`. |
| `tables/` | Generated tables and TikZ figures. `lattice.tex` (Figure 1) was drawn once and is not regenerated. |
| `scripts/zomega.py` | Generates the tables of Sections 3 and 6 and checks the numerical claims of Sections 2–3. |
| `scripts/ideal_figures.py` | Generates Figures 2–4 and checks the claims of the subsection on ideals. |
| `scripts/poisson_theta_checks.py` | Numerical checks for Sections 4 and 5 (about 100 seconds). |
| `scripts/revision_checks.py` | Checks the cusp construction of Section 5 and the index of Γ₂ exactly, and the derivative and Mellin substitution formulas symbolically; needs sympy. |

## Building and checking

```
latexmk -pdf quadratic-fields.tex             # build the PDF
python3 scripts/zomega.py tables              # regenerate the tables and run their checks
python3 scripts/ideal_figures.py tables       # regenerate Figures 2-4 and run their checks
python3 scripts/poisson_theta_checks.py       # checks for Sections 4-5; needs mpmath
python3 scripts/revision_checks.py            # cusp construction, index of Gamma_2, Mellin substitution; needs sympy
```

Each script stops with an `AssertionError` if a check fails. `poisson_theta_checks.py` and `revision_checks.py` print one
`[PASS]` line per check.

These checks test identities on finite samples. They do not prove the imported theorems, among them
Kubota's theorem that the Kubota character is a homomorphism, the theorem of Kubota and Patterson on
the cubic theta function, Hasse's formula for the cube of a cubic Gauss sum, and the Hasse–Davenport
product formula.
