# CANDIDATE_REGISTER.md

Status classes: REJECTED / UNRESOLVED / PROMISING / CONDITIONALLY ACCEPTABLE.

| id | candidate | family | decisive check | status |
|---|---|---|---|---|
| C1 | continuum shunt-resonator field ψ | bulk effective | reduces to Milton–Willis effective medium | REJECTED |
| C2 | backward Stoneley surface branch, half-space interface | interface | E1: det(λk,λω)~1e-45 → exactly nondispersive; no scatterer → no reflection event | REJECTED (numerically falsified, exp1) |
| C3 | relaxed electric displacement D+τ_D Ḋ | constitutive | charge conservation/entropy inadmissible; Cattaneo relabelling | REJECTED |
| E | distributed electromechanical edge (edge Lagrangian) | dynamic BC | symbolic elimination → K_eff(ω,kx) rational + A kx² = standard dynamic impedance edge | REJECTED (reduces to established BC class) |
| D | segmented-electrode edge | dynamic BC | same elimination; segmentation beyond kx² unresolved but nonlocal impedance BCs published | REJECTED as stated; UNRESOLVED only for fully nonlocal segmentation |
| M1 | electrically dressed plate edge | dynamic BC | same as E/D | UNRESOLVED→ not acceptable without nonlocal proof |
| M2 | piezo plate step/edge scattering (benchmark) | plate | published (Smirnov 2023; Gérardin 2016) | REJECTED as novelty; kept as oracle |
| N1 | piezoelectric Willis material | bulk effective | Willis + piezo = composition of published laws | REJECTED |
| N2 | spatiotemporally modulated shunt edge | active modulation | published 2021–2026; active, nonreciprocal — not passive NR | REJECTED |
| N3 | topological phononic domain-wall edge | topological | published topological acoustics; chiral edge is not negative reflection | REJECTED |
| N4 | dynamic impedance edge with resonator (elastic "impedance patch") | dynamic BC | Acta Acustica 2025 / QJMAM 2024 — established | REJECTED |
| N5 | nonlocal (memory-kernel) interface spring | nonlocal interface | Vinh–Giang 2011 class — published | REJECTED |

## Exploratory evidence

- `EXPLORATORY_EXPERIMENTS/exp1_stoneley_vs_plate.py`: Al6061|PZT-5H interface,
  column-normalized 8×8 interface determinant, dps 30. Roots found at every k
  with rel|det| ≈ 7e-45; scale test det(2k, 2ω) ≈ 1.25e-45 ⇒ ω/k exactly
  constant ⇒ Stoneley branch nondispersive ⇒ no backward branch, no ZGV ⇒ C2
  impossible for this pair. (Numerical observation, provenance: script above.)
- `DERIVATIONS/edge_elimination.md`: elimination of Q and r from the Arena edge
  Lagrangian gives K_eff(ω,kx)/Kc = −2A kx²(CLω²−1)/(ACL kx²ω² − A kx² − C Kc L ω² + CLmω⁴ + Kc + Lα²ω² − mω²),
  a rational dynamic impedance plus an ordinary A kx² gradient term — the
  established dynamic-boundary-condition class.
