# FINAL_MODEL_DISCOVERY_REPORT.md

Autonomous discovery investigation, 2026-11-05, branch
`arena/1aa803b8-negative-reflection`. Artifacts under `research_novelty/`.
G4 remains OPEN; no production code touched.

## A. Research coverage

~14 literature clusters searched (plate ZGV/NR, piezo backward waves, Stoneley
and generalized interfaces, gratings, shunts/metasurfaces, Willis, spatiotemporal
modulation, dynamic/impedance edges, topological, nonlocal interfaces), access
mostly abstract/bibliographic (stated per source in PRIOR_ART_EQUATION_MATRIX.md).
12 candidate model families examined (C1–C3, E, D, M1, M2, N1–N3, N5).

## B. Rejection findings

- **C2 (Arena's selected model) REJECTED and numerically falsified**: the
  Al6061|PZT-5H Stoneley branch is exactly scale-free — column-normalized
  interface determinant hits ~1e-45 at (k,ω) and at (2k,2ω)
  (`EXPLORATORY_EXPERIMENTS/exp1_stoneley_vs_plate.py`), so ω/k is constant:
  nondispersive, no backward branch, no ZGV; moreover a uniform planar interface
  provides no momentum-transfer scatterer, so no reflection event exists.
  Closest prior: Murdoch 1977.
- **E/D/M1 edge networks REJECTED as governing-law novelty**: eliminating the
  internal variables from the edge Lagrangian yields K_eff(ω,kx) = rational in
  ω² plus one A kx² term (`DERIVATIONS/edge_elimination.md`) — the established
  dynamic-impedance boundary-condition class (Acta Acustica 2025; QJMAM 2024;
  Vinh–Giang 2011 for nonlocal variants). Also double-counts plate piezo coupling.
- **C1, N1**: reduce to Milton–Willis effective media. **C3**: inadmissible.
  **N2**: published active modulation. **N3**: topological edge ≠ NR. **M2**:
  published plate mechanism (Germano 2012; Smirnov 2023).

## C. Strongest surviving model

None survives as a *new governing law*. The strongest surviving *phenomenon*
route is M2 (piezoelectric plate edge/step conversion into a backward branch) —
but its governing equations are published; it is an application/validation
baseline, not a novel model.

## D. Mathematical evidence

- E1 (numerical observation, provenance exp1 script): Stoneley branch exactly
  nondispersive for Al|PZT → C2 impossible.
- Edge elimination (derived, sympy-verified): K_eff(ω,kx) rational + A kx² →
  standard dynamic BC → M1/D/E not a new law.

## E. Novelty gates

See NOVELTY_GATE_RESULTS.md. No candidate passes N1+N4+N6; C2 fails N4/N5,
edge family fails N1/N6.

## F. Confidence

**High** that no defensible novel governing model exists in the investigated
space: every candidate reduced to a published class under elimination,
scaling, or limiting-case analysis, and the one Arena-selected candidate was
numerically falsified.

## G. Decision

**NO-GO** — no defensible novel governing model found within the investigated
space. (Equivalently: "NO DEFENSIBLE NOVEL GOVERNING MODEL FOUND WITHIN THE
INVESTIGATED SPACE.")

## H. Next necessary proof

If the programme continues, the single highest-value bounded investigation is a
*proof-level* classification: show for a fully nonlocal segmented-electrode edge
(memory kernel in x, not just kx²) whether the eliminated boundary operator is
equivalent to published nonlocal interface laws (Vinh–Giang class) or genuinely
new — with an existence/nonexistence theorem for a backward edge branch. Any
candidate must also pass the G4 energy identity first (currently OPEN).

## I. Repository evidence

Research-only artifacts committed under `research_novelty/` on
`arena/1aa803b8-negative-reflection`; production code, tests and blocker
status unchanged.
