# Provenance statement

This repository's mathematical content, code, and manuscript were produced in an
AI-assisted (Claude Fable 5) interactive session on June 12, 2026, directed by the repository owner.

## Division of labor

**AI-assisted:** problem analysis; choice of methods; Python scripts; the
derivations (Wilf–Zeilberger reduction, Otter cancellation lemma, labeled-tree
closed form and its proof); the LaTeX manuscript including prose; the numerical
constants and certificate tables as outputs of the scripts.

**Human (repository owner):** posing the problem; directing the investigation
across its stages; running the computations on local hardware; checking that the
scripts execute and reproduce the reported outputs; curation and publication of the repo.

## What this means for citation and trust

- The repository owner does **not** claim sole mathematical authorship of the
  results.
- No claim in this repository should be taken on authority. The data is checked
  against the independent OEIS b-file; the finite certificate searches are
  validated on decoy sequences with known answers; the analytic pipeline
  reproduces an independently published constant (Kotesovec, A000598) to 22
  digits. Please check the scripts.
- The manuscript (`noA_certificates.pdf`) has **not** been peer-reviewed.
- Prior art: the leading-order asymptotics are classical (Pólya 1937;
  Robinson–Harary–Balaban 1976; Finch, *Mathematical Constants*, §5.6.1;
  Flajolet–Sedgewick, *Analytic Combinatorics*, p. 477–478). A literature
  check on June 12--13, 2026 found no OEIS formula line or searchable published
  source for the correction terms c₁, c₂, c₃ for A000602, and no published
  finite non-existence certificates of the kinds presented here; this check was
  not exhaustive and corrections are welcome.

## Errors

Please report errors or prior references as issues. The main numerical claims
and finite certificate searches are tied to reproducible scripts, so many
questions about the computations can be checked directly by rerunning them.
