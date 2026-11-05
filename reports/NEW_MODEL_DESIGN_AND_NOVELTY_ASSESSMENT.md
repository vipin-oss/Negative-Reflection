# NEW GOVERNING MODEL DESIGN AND NOVELTY ASSESSMENT — NEGATIVE REFLECTION

Branch `arena/1aa803b8-negative-reflection`, checkpoint 6ef33e9.
Status of existing work: the 6mm Lord–Shulman (LS) piezothermoelastic
half-space model, the repaired scattering solver, the 2U-W normal-incidence
completeness certificate are the **established foundation**, not the novelty.
G4 (energy identity) remains **OPEN** and is not hidden here; G5, SCT02
unresolved. No production code modified by this report.

---

## PHASE A — Competing published models (governing-equation level)

| # | Model family | Actual governing structure | Relation to NR |
|---|---|---|---|
| A1 | Backward Lamb waves / ZGV mode conversion in plates (Mindlin 1955; Tolstoy–Usdin 1957; Prada–Clorennec 2012; Gérardin 2018/2019 Sci. Rep.) | Rayleigh–Lamb dispersion of a **plate of finite thickness**; S1↔S2b conversion at a free edge or thickness step; ρ = −1 at the ZGV point where k(S1)=k(S2b) | Direct NR, but requires a plate branch with negative slope (descending onset from a thickness resonance) |
| A2 | Grating/phononic interface negative reflection (APL 2023, "strange case") | Helmholtz for out-of-plane shear + periodic hole grating; NR from isofrequency-contour folding, bisector rotated ~45° | NR via periodic interface, not a continuum model |
| A3 | Locally resonant / shunted piezoelectric metamaterials (Milton–Willis effective mass; negative-capacitance shunts) | Effective ρ_eff(ω), K_eff(ω) from unit-cell resonances; double negativity ⇒ negative index | Negative **refraction** mostly; effective-medium, not a new continuum PDE |
| A4 | Nonlocal / gradient / flexoelectric / fractional / porous / rotating models | Add Eringen kernel, ∇²u, flexo coupling, ∂^α, Biot porosity, Ω×u terms to classical elasticity | Excluded by the task as mere add-ons |
| A5 | Thermodynamically consistent coupled theories: LS, Green–Lindsay, Green–Naghdi II/III, DPL, MGT, Guyer–Krumhansl | Differ in the heat-flux law q + τq q̇ = −κ∇θ (LS), thermal displacement α̇_T=θ (GN), dual lags, etc. | Set the **thermodynamic frame** the new model must respect; LS is our frozen foundation |

Key gap identified: every published NR mechanism is either (i) a **plate**
(geometric dispersion, thickness resonances, ZGV) or (ii) a **periodic
effective medium**. No published continuum model produces negative
reflection at a **planar interface between two homogeneous half-spaces**
through a **backward surface-wave branch** of the half-space itself. That
is the design target.

## PHASE B — Three differentiated candidates

### Candidate C1 — "Continuum shunt-resonator electromechanical field"
- Architecture: 6mm piezothermoelastic half-space in which the piezoelectric
  layer is bonded to a **distributed shunt lattice**, homogenized not as
  effective parameters but as a bulk internal field ψ(z,t) (shunt charge)
  with its own kinetic and potential energy.
- Fields: u, φ, θ (foundation) + internal variable ψ.
- Potential: Helmholtz free energy Ψ(ε, E, θ, ψ, ∇ψ) with an added
  quadratic ψ-energy and a coupling term linear in ψ·div D.
- Governing PDEs: foundation equations + one new second-order-in-time
  equation for ψ (a damped oscillator field) + modified Gauss law.
- Backward branch: the ψ-resonance gives a backward electromechanical
  branch; NR at a material step via matched backward branches.
- Closest prior: A3 (shunt metamaterial effective models). **Exact
  difference**: prior art replaces the lattice by frequency-dependent
  effective constants; C1 keeps ψ as a genuine continuum field with a
  variational energy, so interface conditions and energy balance are
  derived, not fitted.
- Objections: homogenization validity (subwavelength) must be proved;
  effective-medium referee will say "this is Milton–Willis with extra
  steps"; thermodynamic consistency of the ψ-energy must be shown (G4!).

### Candidate C2 — "Backward piezothermoelastic surface branch at a half-space interface" (SELECTED)
- Architecture: two bonded homogeneous 6mm LS piezothermoelastic
  half-spaces (Al6061 dielectric surrogate | PZT-5H class), planar
  interface; a surface (Stoneley-type) wave localized at the interface is
  the carrier. No plate, no thickness resonance, no periodicity.
- Fields: u, φ, θ per half-space (exactly the frozen foundation).
- Potential: the foundation free energy; LS relaxation via the standard
  Coleman–Noll procedure (heat-flux relaxation only — the frozen frame).
- Governing PDEs: unchanged foundation PDEs in each half-space; new content
  is the **interface eigenproblem** and the surface-branch dispersion.
- Backward branch: the pyroelectric coupling p and the thermal relaxation
  modify the surface-wave determinant; near a degeneracy of the surface
  branch with a bulk branch, the surface dispersion ω(k) can bend
  downward (negative dω/dk) — a backward **surface** branch.
- Negative reflection: an incident forward surface wave at the interface
  converts into the backward surface branch of the second half-space with
  matched k (analogue of the ZGV matching, but in a half-space surface
  eigenproblem, not the Rayleigh–Lamb plate equation).
- Closest prior: A1 (backward Lamb/ZGV in plates). **Exact difference**:
  plate (two stress-free faces, thickness resonances, Rayleigh–Lamb
  determinant) vs **half-space interface surface wave** (one interface,
  Stoneley-type determinant with piezoelectric + pyroelectric + LS terms).
- Objections: backward surface branches in half-spaces are not guaranteed —
  must be found numerically; if none exists the model fails (falsifiable);
  the mechanism is still "mode conversion with matched k", so a referee may
  say it is the plate mechanism transplanted. Mitigation: the eigenproblem
  is genuinely different (surface, not plate) and includes thermal/
  pyroelectric terms absent from A1.

### Candidate C3 — "Charge-relaxation (Cattaneo-type) electric displacement"
- Architecture: foundation model + a relaxation law for the electric
  displacement: D + τ_D Ḋ = e ε + ε E + p θ (electromagnetic inertia /
  charge-carrier relaxation), instead of relaxing only the heat flux.
- Fields: u, φ, θ; D becomes an independent dynamic field.
- Potential: extended free energy with a ½τ_D|Ḋ|² kinetic-like term
  (must be justified thermodynamically — the hard part).
- Backward branch: the τ_D relaxation adds a second timescale to the
  electromechanical dispersion; near it, backward branches can appear.
- Closest prior: Cattaneo–Maxwell/LS applied to heat flux (A5); applying
  relaxation to **D** rather than q is not standard in piezoelectricity.
- Objections: thermodynamic admissibility of a relaxed D is doubtful
  (D is not a flux; relaxing it may violate the entropy inequality —
  likely referee kill); risk of being "Cattaneo with different letters".

## PHASE C — Ranking and selection

| criterion | C1 | C2 | C3 |
|---|---|---|---|
| Mathematical originality | high | medium | medium-high |
| Physical credibility | medium | high | low (thermo risk) |
| Likelihood of actual NR | medium | high | low-medium |
| Realistic parameters | medium (shunt values) | high (foundation tables) | low |
| Independent validation feasibility | medium | high (existing solver + oracle) | low |
| Referee strength | medium | high | low |

**Selected: C2** — backward piezothermoelastic surface branch at a
half-space interface, because it keeps the frozen, validated foundation
(LS 6mm, existing characteristic matrix, existing energy oracle), targets
negative **reflection** (not refraction) directly, uses only measured
material tables, and is falsifiable: if no backward surface branch exists
in the interface eigenproblem, the model is dead.

## PHASE D — Derivation of C2 (sketch at governing-equation level)

Free energy per half-space (foundation, 6mm + LS):

  Ψ = ½Cijkl εij εkl − e_kij E_k εij − β_i θ ε_i − ½ε_ij E_i E_j − ½ρcE θ²/T0 + p θ E_3,

with σ = ∂Ψ/∂ε, D = −∂Ψ/∂E, entropy s = −∂Ψ/∂θ; LS law
(1+τq∂t)q = −κ∇θ; momentum ρü = div σ; Gauss div D = 0.

Surface-wave ansatz in each half-space: fields ∝ exp(i(kx + ξ_α z − ωt))
with Re(ξ_α) > 0 decay into that half-space — the existing 8-root
characteristic polynomial of the foundation provides ξ_α(k, ω).

Interface conditions (bonded, no free charge, thermal contact):
[u]=0, [σ·n]=0, [φ]=0, [D·n]=0, [θ]=0, [q·n]=0.

Dispersion eigenproblem: det M_surface(k, ω) = 0 with M_surface the 8×8
interface matrix built from the 4 admissible roots per side (exactly the
machinery of the repaired production solver, now applied to a
**surface/localization** selection instead of the k=0 scattering basis).

Backward-branch mechanism: surface dispersion ω_s(k) near an avoided
crossing with a bulk branch; negative dω_s/dk segment ⇒ phase and group
velocity antiparallel ⇒ backward surface wave. Negative reflection:
incident forward surface branch in medium 1 at (ω, k) matching the
backward branch of medium 2 at the same k ⇒ near-total conversion with
ρ ≈ −1, energy returning on the incident side (group velocity of the
converted wave points back) — the NR signature.

Boundary conditions at an edge/interface terminating the half-space:
the standard 8 conditions above plus, at a free edge, σ·n = 0, D·n = 0.

## PHASE E — Decisive proof of novelty / falsification tests

Minimum results to distinguish C2 from A1 (plate Lamb/ZGV) and A2 (grating):

1. Compute the surface-branch dispersion ω_s(k) for Al|PZT and PZT|BaTiO3
   half-space interfaces; exhibit a segment with dω_s/dk < 0. (If none:
   model falsified.)
2. Reflection eigenproblem: incident forward surface branch, compute the
   amplitude into the backward branch; show |ρ| → 1 with negative real
   part at the matched k, and energy flux (group velocity) reversed.
3. Classical limits: τq→0, p→0, e→0 must reduce to the classical elastic
   Stoneley problem (no backward branch — confirming the branch is due to
   the new coupling terms, not numerics); plate limit d→∞ must NOT
   reproduce Rayleigh–Lamb (different determinant).
4. Thermodynamic check: the surface-branch energy flux must satisfy the
   G4 identity once G4 is accepted (currently OPEN — flagged, not hidden).
5. Falsification tests: swap ε11/ε33, flip sign of p, remove thermal
   relaxation — each must destroy the backward segment (the G4-1 mutation
   machinery is reusable).

## Critical risks

- No backward surface branch may exist in the half-space eigenproblem
  (most likely failure mode; falsifiable by test E1).
- Referee objection "plate mechanism transplanted": answer is the different
  determinant and the thermal/pyroelectric terms absent from A1.
- G4 energy identity unresolved — the model's energy claims must wait for
  G4 acceptance; this report makes no energy-conservation claim.
- No experimental parameters invented; all parameters from the frozen
  material tables.
