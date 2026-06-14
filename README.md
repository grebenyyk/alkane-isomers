# Alkane isomer counts (OEIS A000602): asymptotics and finite exclusion tests

Computational study of the sequence **A000602** — the number of structural
isomers of the alkane CₙH₂ₙ₊₂, equivalently the number of unlabeled, unrooted
trees on *n* vertices with maximum degree 4.

**Provenance: the mathematics, code, and manuscript in this repository were
assisted by Claude Fable 5 model.** See [PROVENANCE.md](PROVENANCE.md) before citing. The numerical claims and finite certificate searches are reproducible by running the scripts in this.

## Results

**1. Asymptotic expansion with computed constants**:

    a(n) ~ C · α^n · n^(−5/2) · (1 + c₁/n + c₂/n² + c₃/n³ + ...)

    α  = 2.81546003317615074652661677824269954253651...   (growth constant; = 1/ρ)
    ρ  = 0.355181742314377392882244473647632636708747...   (= OEIS A261340)
    C  = 0.65631869584183475...
    c₁ = 0.0300686950399...
    c₂ = 4.7882749807...
    c₃ ≈ 22.18

The growth constant and ρ are published (Kotesovec 2015 on A000598; A261340).
The unrooted amplitude C extends Finch's ~10 digits (§5.6.1). I found no OEIS
entry/formula line or searchable published source for the correction
coefficients c₁–c₃; corrections and prior references are welcome.

Empirical error bounds on the checked range 25 ≤ n ≤ 600: the 1-, 2-, 3-term
truncations have relative error ≤ 0.276/n, 6.16/n², 35.5/n³ respectively.

**2. Finite non-existence certificates for standard formula classes.** Exhaustive,
exact-arithmetic (modular) searches prove the following statements within the
reported finite search regions:

- **no polynomial-coefficient recurrence** (homogeneous or inhomogeneous) in a
  staircase region up to order 2/degree 150, order 10/degree 46, order 304 with
  constant coefficients — excluding Wilf–Zeilberger-class
  (binomial/Γ/hypergeometric finite-sum) formulas whose canonical annihilator
  lies in that region;
- **no algebraic equation** P(x, A(x)) = 0 with deg_y ≤ 28 (staircase), and no
  first-order algebraic ODE of degree ≤ 7 per variable;
- **no small tested digit-automatic structure**: the 2-kernel of a(n) mod 2
  grows maximally (≥ 31 states through depth 4), and no eventual periodicity
  mod 2..7 with period ≤ 200 appears in the tested window.

**3. Structural comparison with the labeled problem.** The *labeled* analogue
has an exact two-index finite-sum closed form (derived and verified herein), so
the degree-4 constraint itself is not the main difficulty. The unlabeled quotient
introduces a Burnside sum over all partitions of n, and no analogue of the
Rademacher modular-circle mechanism is evident for the mixed (2,3)-Mahler
system arising here.

Full statements, proofs, tables, and references: `noA_certificates.pdf`.

## Repository contents

| File | Purpose |
|---|---|
| `asym.py` | Exact terms a(0..620); Newton solve for (ρ, T(ρ)) to 55 digits; Richardson extrapolation of C, c₁, c₂, c₃; error table |
| `bounds.py` | Certified error-bound scan over 25 ≤ n ≤ 600 (depends on `asym.py`) |
| `solveA.py` | The certificate battery: P-recurrence / algebraic / ADE / automaticity searches with decoy validation and second-prime re-checks; labeled-tree closed form verified against brute force |
| `crosscheck.py` | Pipeline validation: reproduces Kotesovec's published A000598 constant to 22 digits |
| `verify_bfile.py` | Term-by-term comparison against the independent OEIS b-file (all 621 values agree) |
| `b000602.txt` | OEIS b-file for A000602 (P. von Brömssen, terms 0..1000) |
| `noA_certificates.tex` / `.pdf` | The full write-up (~13 pp.) |

## Reproducing

Requires Python ≥ 3.9 with `numpy` and `mpmath`:

```sh
python asym.py          # constants + error table            (~1 min)
python bounds.py        # empirical bounds scan              (~1 min)
python solveA.py        # finite non-existence certificates  (~2 min)
python verify_bfile.py  # data certification vs OEIS
python crosscheck.py    # pipeline validation vs published constant
```

Each certificate in the paper corresponds to a "full rank" line in
`solveA.py` output; the searches are validated by planted decoys (Catalan,
Thue–Morse) that the same code *does* find.

## Status of the results

- The leading asymptotic form is classical (Pólya 1937; constants:
  Robinson–Harary–Balaban 1976, Finch 2003). Computed here: additional digits
  for the unrooted amplitude, correction coefficients c₁–c₃, empirical error
  bounds on 25 ≤ n ≤ 600, and finite non-existence certificates in the reported
  search regions.
- The finite rank certificates are rigorous modulo correctness of the exact
  integer/modular arithmetic; full rank modulo a single prime proves
  non-existence over ℂ for the corresponding finite-dimensional ansatz (see
  paper, §3).
- The natural-boundary and global non-holonomicity statements are
  conditional/structural and clearly flagged as such in the paper (§7).
- The manuscript has not been peer reviewed.

## License

Code and text: MIT (see [LICENSE](LICENSE)). The manuscript PDF/TeX may be
reused under the same terms with provenance retained.
