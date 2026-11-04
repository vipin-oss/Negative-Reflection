# 07 — Five Research Problems

**Document status:** Phase 5 deliverable
**What this document is:** five genuinely distinct research problems, each specified against the
18 required elements. Distinctness is achieved by varying the **mechanism**, not the material.
**What it is not:** a proposal. Ranking is in `09`; only the top one or two are developed fully, in
`10`.

---

## The 18 required elements (identical for all five)

| # | Element | Why it is required |
|---|---------|--------------------|
| 1 | Title | — |
| 2 | Research question | Must be answerable |
| 3 | Central hypothesis | Must be falsifiable |
| 4 | Physical mechanism claimed | Must be a mechanism, not a material |
| 5 | Differentiation from ordinary specular reflection | The brief's core demand |
| 6 | Governing equations and boundary conditions | Complete, incl. higher-order conditions |
| 7 | Governing dimensionless group(s) | What actually controls the answer |
| 8 | Material system and parameter values | Sourced, realistic |
| 9 | Analytical programme | What is done by hand |
| 10 | Numerical programme | What is done by machine |
| 11 | Validation strategy | Against what independent reference |
| 12 | Falsification / go-no-go test | Cheap test that can kill it |
| 13 | Expected results — positive and null | Both must be stated in advance |
| 14 | Novelty claim and prior art | With the exact residual claim |
| 15 | Application / impact | Concrete |
| 16 | Feasibility risk register | Ranked |
| 17 | Integrity statement | What is explicitly *not* claimed |
| 18 | Outputs and staged timeline | 18 months |

---

# PROBLEM P1 — Bias-tuned negative reflection in a polarized piezoelectric plate: does forward↔backward near-degeneracy survive electromechanical coupling?

### 1. Title
*Electromechanical coupling and the near-degeneracy of forward/backward Lamb modes: a modal
scattering theory for the DC-biased piezoelectric plate edge.*

### 2. Research question
In a piezoelectric plate poled through the thickness, the S₁/S₂ᵦ zero-group-velocity fold moves by
−35 % relative to the elastically identical uncoupled plate (computed, `05` §3.1). **Does the
near-degeneracy of the two mode shapes — the condition that makes the free-edge reflection
coefficient approach −1 — survive this shift, and does an applied DC bias provide continuous control
of the negative-reflection window?**

### 3. Central hypothesis
**H1.** The mode pair carrying the fold in a coupled piezoelectric plate is *not* the same pair as in
the elastic plate; the pair is selected by the electrical boundary condition and by the
piezoelectric stiffening, not by the elastic branch repulsion alone.
**H2.** A DC bias, which changes the effective elastic constants through the converse piezoelectric
effect, tunes f_ZGV continuously and monotonically while preserving the near-degeneracy, because it
acts on the *dispersion* rather than on the *mode-shape symmetry*.
**H1 and H2 are separately falsifiable.** H1 fails if the pair is unchanged and the mode shapes are
unchanged; H2 fails if bias destroys the near-degeneracy before it usefully moves the fold.

### 4. Physical mechanism claimed
Forward↔backward mode conversion at a traction-free edge (mechanism **M1** of `03`), modulated by
**converse-piezoelectric stiffening**. The bias is a *stiffness* control, not a *boundary-condition*
control — this is the essential distinction from the open/short switch, which was computed to be
inert for this mode (−0.0006 %).

### 5. Differentiation from ordinary specular reflection
Ordinary reflection: the reflected beam is on the opposite side of the normal, and the tangential
Poynting flux reverses with the wave vector. Here the reflected mode is a **backward** mode at the
same wavenumber, so the tangential flux *keeps its sign* while the beam returns along the incident
path. The operational test is the flux test stated in `05` §1, applied to the electromechanical flux,
not to Re(ρ).

### 6. Governing equations and boundary conditions

Quasi-static piezoelectricity, time-harmonic e^{−iωt}:

$$\rho\,\ddot u_i = \sigma_{ij,j},\qquad D_{i,i}=0,$$
$$\sigma_{ij}=c_{ijkl}\,u_{k,l}+e_{kij}\,\varphi_{,k},\qquad D_i = e_{ijk}\,u_{j,k}-\varepsilon_{ij}\varphi_{,j},$$

with c, e, ε the *biased* effective constants (linearised about the bias point; the bias enters
through the third-order constants if these are available for the material, and otherwise through the
measured bias dependence of c^E — the proposal must state which, because this is the main modelling
assumption).

Plate: |x₂| ≤ h, free faces x₂ = ±h: σ₂ⱼ = 0; electrical condition: φ = 0 (short) or D₂ = 0 (open),
**and** a bias circuit condition that must be specified (a bias field with a blocked charge is not
the same as an open circuit).

Edge at x₁ = 0: σ₁ⱼ = 0, plus the electrical edge condition, whose correct form for a *cut* edge
(electrically open across the cut face) is one of the things this problem must settle — it is not
obvious and no source was found that states it for the modal-scattering problem.

### 7. Governing dimensionless group
$$\Pi_1=\frac{e_{31}^2}{\varepsilon_{33}\,c_{11}}\ \ (\text{electromechanical coupling, }\approx 0.026\text{ for PZT-5H}),\qquad
\Pi_2=\frac{\omega h}{c_T}\ \ (\text{normalised frequency}),\qquad
\Pi_3=\frac{E_{\text{bias}}\,e_{33}}{c_{11}}\ \ (\text{bias strength}).$$

**Note the structure:** Π₁ contains no length. This is why pure piezoelectricity cannot produce a
size effect, and it is why Π₃ (also dimensionless, also length-free) can tune but cannot create a
size-dependent phenomenon.

### 8. Material system and parameters
PZT-5H (c₁₁ = 126, c₁₂ = 79.5, c₁₃ = 84.1, c₃₃ = 117, c₄₄ = 23.0 GPa; e₃₁ = −6.5, e₃₃ = 23.3,
e₁₅ = 17.0 C/m²; ε₁₁ = 1700ε₀, ε₃₃ = 1470ε₀; ρ = 7500 kg/m³). These are the standard
manufacturer/published values used in `code/plate_spectral.py` and already validated. Plate
thickness 0.5 mm and 1.5 mm. Frequency range 0.5–4 MHz.
**Bias range** must be quoted from measured bias dependence of the elastic constants for PZT — this
is a data-sourcing task in month 1; if no reliable data exist for PZT-5H under bias, the problem
falls back to LiNbO₃, for which bias-dependent behaviour is better characterised.

### 9. Analytical programme
(i) Derive the modal orthogonality and the electromechanical flux for the *biased* plate.
(ii) Derive the completeness relation for the mixed propagating/evanescent spectrum at the edge.
(iii) Obtain the near-degeneracy condition (C2) as an explicit overlap integral of the two mode
shapes in the traction norm.
(iv) Perturbation theory in Π₃ for the shift of the fold.
(v) **Explicitly prove or disprove** the claim that the pair is unchanged by coupling — this is H1
and it is an analytical question, not a numerical one, if the perturbation is small.

### 10. Numerical programme
Spectral Chebyshev solver (exists, `code/plate_spectral.py`), extended with:
(i) bias-dependent constants and a bias circuit condition;
(ii) bias continuation;
(iii) the edge scattering matrix by the modal projection method with pseudo-inverse of an
ill-conditioned evanescent basis;
(iv) energy-budget audit of every computed reflection matrix.

### 11. Validation strategy
Four independent checks, in increasing strength:
- cutoffs against exact n·c/(4h) analytic values (already 6.3 × 10⁻⁷ here);
- dispersion against an independent Rayleigh–Lamb root finder (already 5.3 × 10⁻⁷);
- **bias → 0 limit must reproduce the unbiased result exactly** (a strong consistency test);
- **reciprocity**: ρ_ij = ρ_ji in the flux-normalised basis;
- **energy**: Σ|ρ|²(P/P_inc) = 1 to machine precision for the lossless case.

### 12. Falsification / go-no-go (month 2)
Compute the mode shapes of the pair at the coupled fold and their overlap in the traction norm.
**If the overlap is not near-degenerate (i.e. if the pair cannot approximately satisfy the free-edge
condition alone), the problem is terminated** and the negative result reported. This is cheap — one
eigenvalue run across the fold — and it is scheduled before any scattering code is written.

### 13. Expected results
*Positive:* a modal scattering matrix for a biased piezoelectric plate edge; a demonstrated,
continuously tunable negative-reflection window; a criterion for when coupling destroys the effect.
*Null:* the pair is unchanged by coupling and bias moves the fold without affecting conversion
efficiency — in which case the honest conclusion is that electromechanical coupling is a
*dispersion* modifier, not a *reflection* modifier, which is itself a clean and useful statement.

### 14. Novelty claim and prior art
Prior art: M1 in isotropic plates (Germano 2012; Gérardin 2016; Veres 2016; Gérardin 2019);
backward waves in piezoelectric plates with ZGV (Smirnov et al. 2023). **Residual claim:** the first
modal edge-scattering calculation for a *coupled* plate, and the first test of whether
electromechanical coupling preserves the near-degeneracy that M1 depends on. Not claimed: the
existence of backward waves in piezoelectric plates (published), or a size effect (computed to be
absent).

### 15. Application
An electrically tunable, passive acoustic mirror for guided waves — a "varifocal" reflector with no
moving parts, relevant to guided-wave structural health monitoring (beam steering without
repositioning transducers) and to acoustic delay lines.

### 16. Feasibility risk register (ranked)
1. **Bias-dependent elastic constants may not exist for the chosen material** (high) — fallback to
   LiNbO₃ or to a purely theoretical bias treatment with stated linearisation.
2. **Near-degeneracy destroyed by coupling** (medium) — then the problem converts into a negative
   result, still publishable but smaller.
3. **Edge electrical condition for a cut face is not uniquely defined** (medium) — must be resolved
   in month 1; will be stated as a modelling choice with a sensitivity study.
4. **Evanescent-basis ill-conditioning** (medium, standard) — mitigated by the pseudo-inverse with
   truncation and by a convergence study in the number of retained modes.
5. **Effect too small to be interesting** (low–medium).

### 17. Integrity statement
This problem does **not** claim that negative reflection has been observed in a piezoelectric plate.
It does not claim a size effect. It does not transplant any elastic identity to the coupled case
without re-derivation. The open/short switch has been computed to be inert for this mode and that
negative finding will be reported, not hidden.

### 18. Outputs and timeline
M1–M3: bias-dependent constitutive data, edge-condition resolution, near-degeneracy gate.
M4–M8: modal scattering code + validation. M9–M12: bias sweeps, energy audit. M13–M15: manuscript.
M16–M18: response to review, extension to a second material. Outputs: one journal paper, one open
solver release, one dataset.

---

# PROBLEM P2 — Does material grading preserve the ρ ≈ −1 condition? Thickness-graded plates as a controllable negative-reflection window

### 1. Title
*Functional grading as a mode-shape control: preserving and widening forward/backward conversion at
a plate edge.*

### 2. Research question
Grading changes the mode shapes continuously through the thickness. Since the ρ ≈ −1 condition is a
statement about **mode-shape near-degeneracy**, not about the dispersion relation alone: **can a
through-thickness property gradient be designed to widen the frequency/angle window over which
negative reflection is efficient, and does a sufficiently strong gradient destroy it?**

### 3. Central hypothesis
**H.** The negative-reflection bandwidth is controlled by the *derivative* of the mode-shape
difference across the fold, and grading can flatten that derivative over a finite band, producing a
wider and flatter ρ ≈ −1 window than a homogeneous plate — up to a gradient strength beyond which
the pair ceases to be near-degenerate.

### 4. Physical mechanism claimed
M1, with **mode-shape engineering by grading**. The mechanism is not "a new material" but "a designed
spatial variation of the constitutive law that targets a specific algebraic condition".

### 5. Differentiation from ordinary specular reflection
Same flux criterion as P1. Additionally, grading means that a naive plane-wave/Snell analysis is
invalid: the local wavenumber varies through the thickness, and only a full modal treatment can
distinguish backward conversion from a mere angle change.

### 6. Governing equations and boundary conditions
Elastic (or piezomagnetic) plate with x₂-dependent c(x₂), ρ(x₂):
$$\rho(x_2)\ddot u_i = \big(c_{ijkl}(x_2)u_{k,l}\big)_{,j},$$
free faces σ₂ⱼ = 0 at x₂ = ±h, free edge σ₁ⱼ = 0 at x₁ = 0. No higher-order conditions — this is the
crucial simplification relative to P3, and it is why P2 is the low-risk member of the family.

Grading law: power-law or exponential, c(x₂) = c₀·f(x₂/h; n), with n the grading index. Must state
whether the mass density is graded identically (it must be, physically, for a real FG material, and
grading ρ and c independently would be a modelling error).

### 7. Governing dimensionless group
$$\Pi_g=\frac{h}{L_c},\quad L_c=\Big|\frac{\mathrm{d}\ln c}{\mathrm{d}x_2}\Big|^{-1}\ \text{(gradient length)},\qquad
\Pi_2=\frac{\omega h}{c_T}.$$
The claim is that the bandwidth scales with a function of Π_g alone at fixed Π₂.

### 8. Material system and parameters
Two-material FG systems with tabulated properties, e.g. Al/Al₂O₃ or Ti/TiB, and — for the coupled
variant — a graded BaTiO₃/CoFe₂O₄ MEE plate whose constituent constants are well tabulated. Real
FG plates achieve Π_g ~ 0.1–10. Grading index n ∈ [0.2, 5].

### 9. Analytical programme
(i) Derive the modal orthogonality for a graded plate (non-self-adjoint in general? — must be
checked; the operator remains self-adjoint for real c(x₂), which is a useful check).
(ii) Asymptotic analysis in the thin-gradient limit Π_g → 0 to recover the homogeneous result.
(iii) Derive the bandwidth–curvature relation from the flux identity, then show how grading modifies
the curvature of the fold.
(iv) State a design criterion: the gradient profile that minimises d(ΔΦ)/dω across the band, where
ΔΦ is the mode-shape difference in the traction norm.

### 10. Numerical programme
Extend the spectral solver to x₂-dependent coefficients (Chebyshev differentiation with variable
coefficients — straightforward). Compute: dispersion and folds vs n; mode shapes and overlaps;
edge scattering matrix; bandwidth map over (n, Π_g).

### 11. Validation strategy
- Π_g → 0 must reproduce the homogeneous solution to 10⁻⁶.
- Independent validation of the graded dispersion against a Riccati/compound-matrix (Thomson–Haskell
  with many sublayers) solution — the standard independent method for graded media.
- Energy audit and reciprocity, as in P1.

### 12. Falsification / go-no-go (month 3)
Compute the bandwidth vs Π_g for three grading profiles. **If bandwidth is monotone decreasing in
Π_g with no interior optimum**, the hypothesis is false and the problem reduces to "grading degrades
NR", which is a smaller but still reportable result.

### 13. Expected results
*Positive:* a design rule for graded negative-reflection mirrors, with a demonstrated bandwidth gain.
*Null:* grading monotonically degrades conversion — a clean statement that the M1 condition is
fragile to inhomogeneity, which constrains all real devices.

### 14. Novelty claim and prior art
Prior art: graded plates with ZGV points (`10.1121/1.4983296`); M1 in homogeneous plates; Dr. Gupta's
FG piezomagnetic wave work. **Residual claim:** the first study to treat grading as a *mode-shape
control* for the negative-reflection condition, and the first bandwidth criterion derived from
dispersion curvature.

### 15. Application
Wideband passive acoustic retroreflectors; graded buffers that preserve or suppress back-reflection
by design.

### 16. Feasibility risk register
1. **Gain may be marginal** (medium–high) — the homogeneous case is already fairly broadband.
2. **Manufacturability of the designed profile** (medium) — real FG plates have limited profile
   control; the design rule must be robust to profile error (this becomes a sensitivity study).
3. **Numerical cost of variable-coefficient spectral solves** (low).
4. **Novelty may be judged incremental** (medium) — mitigated by the bandwidth *criterion*, which is
   analytical and new.

### 17. Integrity statement
No claim that grading creates negative reflection (it exists already at Π_g = 0); the claim is only
about its *window*. No claim of experimental validation.

### 18. Outputs and timeline
M1–M4: variable-coefficient solver + independent Thomson–Haskell validation. M5–M8: bandwidth map
and go/no-go. M9–M13: design rule and optimisation. M14–M16: manuscript. M17–M18: robustness to
profile error.

---

# PROBLEM P3 — Higher-order boundary conditions and the fate of forward↔backward conversion in a strain-gradient plate

**Lead candidate. Highest novelty, highest risk.**

### 1. Title
*Do higher-order boundary conditions destroy negative reflection? A modal scattering theory for the
strain-gradient elastic plate edge.*

### 2. Research question
Negative reflection at a free edge works because **two** modes suffice to satisfy the **two** classical
traction conditions. A strain-gradient (or flexoelectric) theory requires additional natural
boundary conditions — a higher-order traction, and in some formulations an essential condition on the
normal derivative of the displacement. **Can those additional conditions still be satisfied, and at
what cost to the conversion efficiency?**

### 3. Central hypothesis
**H.** The additional boundary conditions cannot be met by the forward/backward pair alone; they are
satisfied by an **edge boundary layer of evanescent modes whose thickness scales with the internal
length ℓ**. Consequently the reflection coefficient's modulus is reduced from ≈ 1 by a factor that
depends on ℓ/λ, and negative reflection survives only for ℓ/h below a computable threshold.

**This hypothesis is quantitatively falsifiable, and its negation is equally interesting:** if the
higher-order conditions collapse onto the classical ones in the limit (a degenerate formulation), the
model is inadmissible and that must be reported.

### 4. Physical mechanism claimed
M1 **modified by a structural change in the boundary-value problem**. This is the only one of the
five problems in which the coupling changes the *number and type* of conditions, not merely their
coefficients — which is why it is the strongest scientific claim available.

### 5. Differentiation from ordinary specular reflection
Two layers of differentiation:
(a) the flux/backward test of `05` §1, which distinguishes NR from any sign of a coefficient;
(b) **the boundary-layer structure itself**: in a gradient theory, reflection is no longer described
by a coefficient between plane waves but by a coefficient between *modal fields with a boundary
layer*. Any analysis that treats the gradient term as a perturbation of the reflection coefficient
(rather than of the boundary-value problem) is wrong in principle, and demonstrating that is part of
the contribution.

### 6. Governing equations and boundary conditions
Mindlin's strain-gradient elasticity (chosen over Eringen integral nonlocality because it has a
well-defined variational structure and therefore well-defined natural boundary conditions):

$$\sigma_{ij}=c_{ijkl}\varepsilon_{kl},\qquad
\tau_{ijk}=\ell^2\,\partial_i\sigma_{jk}\ \ (\text{simplified form}),$$
$$\big(\sigma_{ij}-\tau_{ijk,k}\big)_{,j}+\rho\,\omega^2 u_i=0,$$
with natural conditions on the boundary:
$$n_j\big(\sigma_{ij}-\tau_{ijk,k}\big)-D_j\big(n_k\tau_{ijk}\big)+\dots=0
\quad\text{(classical traction)},$$
$$n_i n_j \tau_{ijk}=0\quad\text{(double traction)}.$$

The exact form depends on the strain-gradient variant, and **the first analytical task is to fix the
variant and derive the conditions from the variational principle rather than quoting them** — the
brief forbids transplanting equations between theories without verification, and this is exactly the
point where that rule bites.

Eringen-type nonlocal elasticity is to be treated in a **separate, clearly-labelled section** (if at
all), because its consistent boundary conditions are an open modelling question and mixing the two
would compromise the rigour of the main result.

### 7. Governing dimensionless group
$$\Pi_\ell=\frac{\ell}{h}\quad\text{(internal length to half-thickness)},\qquad
\Pi_2=\frac{\omega h}{c_T}.$$
The central prediction is a **collapse of the conversion efficiency onto a function of Π_ℓ alone**
(in the long-wavelength limit), and a threshold Π_ℓ* beyond which ρ < 1 measurably.

### 8. Material system and parameters
Two tiers:
(i) **Generic**: ℓ scanned over 0–0.5 h, no material claim. This is the honest core of the work and
requires no material data.
(ii) **Material-specific**: aluminium and a spinodal alloy, ℓ ~ 1–10 µm from the strain-gradient
literature. **These values must be sourced per material, and the range is wide** — the proposal must
state the uncertainty explicitly rather than quoting a single number.
Physical accessibility: ℓ/h ~ 0.01 requires h ~ 100 µm to 1 mm for ℓ ~ 1–10 µm, i.e. **the effect is
in the micrometre-to-millimetre regime at MHz frequencies** — numerically accessible, experimentally
demanding.

### 9. Analytical programme
(i) Variational derivation of the complete boundary conditions for the chosen gradient variant.
(ii) Proof (or disproof) that the classical limit ℓ → 0 is non-singular for the edge problem — this
is a known hazard and must be settled analytically before any numerics.
(iii) Asymptotic boundary-layer analysis: show that the extra conditions are met by a layer of
thickness O(ℓ) with an explicitly computable exponential decay rate.
(iv) Derive the correction to ρ as a function of Π_ℓ, at leading order, using the layer solution.
(v) State the threshold Π_ℓ* explicitly.

### 10. Numerical programme
Extend the solver to the fourth-order system (the Chebyshev method handles this naturally by
increasing the state dimension; the difficulty is the extra boundary rows, which must be built from
the derived conditions, not from the classical ones). Compute:
dispersion vs Π_ℓ; folds; mode shapes and overlaps including the layer; edge scattering matrix;
the limit Π_ℓ → 0 as a convergence test.

### 11. Validation strategy
- **Reduction test:** Π_ℓ → 0 must reproduce the classical M1 result (Germano/Gérardin geometry) to
  within the asymptotic error.
- **Independent method:** a collocation or finite-element solution of the same boundary-value problem
  for one parameter point.
- **Energy:** the full gradient energy (including the double-stress contribution) must be conserved
  for the lossless case — this is a much stricter test than the classical one and is a strong check
  on the derived boundary conditions.
- **Mesh convergence** in both the spectral order and the retained evanescent modes.

### 12. Falsification / go-no-go (months 2–3, two gates)
**Gate 1 (analytical):** is the ℓ → 0 limit regular? If singular, the formulation must be changed
before proceeding.
**Gate 2 (numerical):** compute the conversion efficiency at Π_ℓ = 0.05, 0.1, 0.2. **If the change is
below 1 % at Π_ℓ = 0.2**, the effect is not physically meaningful and the problem is reframed as a
*negative result* — which is still publishable and is explicitly pre-registered here as an acceptable
outcome.

### 13. Expected results
*Positive:* a quantified degradation law ρ(Π_ℓ), an explicit threshold, and the demonstration that
higher-order boundary conditions are met by an edge layer — the first modal-scattering treatment of a
gradient plate edge.
*Null:* the gradient conditions are satisfied with negligible effect (or the formulation degenerates),
reported as a caution to the size-dependent-waves community, whose dispersion-only treatments
implicitly assume this question away.

### 14. Novelty claim and prior art
Prior art: gradient/size-dependent *dispersion* for plates and Lamb waves (extensive, including
Dr. Gupta's own); gradient electro-mechanics via Stroh formalism (`10.1063/5.0196873`); M1 in
classical plates. **No paper was found that solves a modal edge-scattering problem in a gradient
plate.** Residual claim: the first such treatment, and the first statement of whether higher-order
boundary conditions permit negative reflection at all.

### 15. Application
Determines whether negative-reflection devices can be scaled into the micro/nano regime where
gradient effects matter — i.e. whether the phenomenon is a *macro* curiosity or an *integrated-device*
possibility. It also gives the size-dependent-waves community a boundary-condition result it currently
lacks.

### 16. Feasibility risk register (ranked)
1. **The gradient formulation's boundary conditions may be contested** (high) — mitigated by
   deriving them from the variational principle, choosing one variant explicitly, and treating
   alternatives in a labelled section.
2. **The ℓ → 0 limit may be singular** (medium–high) — Gate 1 exists for exactly this.
3. **The effect may be numerically small** (medium) — Gate 2, and pre-registered as an acceptable
   negative result.
4. **Material data for ℓ are uncertain by an order of magnitude** (medium) — addressed by making the
   dimensionless scan the core and the material values an illustration.
5. **Fourth-order numerics with extra boundary rows are error-prone** (medium) — mitigated by the
   reduction test and the full gradient-energy check.

### 17. Integrity statement
This problem does not assert that nonlocal/strain-gradient theory is "more correct". It does not
transplant the classical modal-expansion method into the gradient setting without re-deriving
completeness for the fourth-order operator. It does not quote boundary conditions from a secondary
source. It commits in advance to reporting a null result. It makes **no experimental claim**.

### 18. Outputs and timeline
M1–M3: formulation + Gate 1. M4–M7: fourth-order solver + reduction validation. M8–M10: Gate 2 and
the Π_ℓ scan. M11–M14: asymptotic layer theory and the ρ(Π_ℓ) law. M15–M17: manuscript. M18: response
and the nonlocal-variant comparison. Outputs: one substantial paper (target: *Wave Motion*, *JSV*, or
*Proceedings of the Royal Society A*), one solver, one dataset.

---

# PROBLEM P4 — Stiffness-controlled negative reflection at an imperfectly bonded joint

### 1. Title
*Tunable negative reflection at a spring-type joint between two elastic waveguides: the reflection
matrix as a function of interfacial stiffness.*

### 2. Research question
At a perfectly free edge the conversion is fixed. At a *bonded joint* with finite stiffness the
boundary condition is a continuous family between free and clamped. **Does negative reflection
survive over a finite stiffness range, and is there a threshold stiffness below which it abruptly
fails?**

### 3. Central hypothesis
**H.** Because M1 depends on a delicate near-degeneracy, the conversion efficiency is **not** a
smooth monotone function of interfacial stiffness: it is approximately unity up to a threshold
(where the joint is soft enough to behave as free) and then collapses. The threshold defines a
practical tolerance for real bonded structures.

### 4. Physical mechanism claimed
M1 with a **mixed boundary condition**. The mechanism claim is about *robustness*: how much
imperfection the phenomenon tolerates.

### 5. Differentiation from ordinary specular reflection
Same flux/backward criterion. Additionally: a naive treatment would use a scalar reflection
coefficient at the joint; the correct object is a **matrix** in the modal basis, and the "negative"
character can hold for one mode pair while failing for another at the same stiffness. Demonstrating
that the matrix (not a scalar) is required is part of the contribution.

### 6. Governing equations and boundary conditions
Two semi-infinite plates joined at x₁ = 0, with the joint modelled by distributed springs:
$$\sigma^{(1)}_{1j}=K_j\big(u^{(2)}_j-u^{(1)}_j\big),\qquad
\sigma^{(2)}_{1j}=\sigma^{(1)}_{1j},$$
with K_j → 0 recovering the free edge and K_j → ∞ the rigid bond. Normal and tangential stiffness
treated separately (they have very different physical values and very different effects).

### 7. Governing dimensionless group
$$\kappa=\frac{K_\parallel}{\rho c_T^2/h}\quad\text{(normalised interfacial stiffness)},$$
the standard dimensionless joint parameter. The prediction is a threshold κ*.

### 8. Material system and parameters
Aluminium or Duralumin plates (to connect directly to the M1 literature), 1–3 mm, with κ spanning
10⁻²–10². Real joints: adhesive bonds give κ ~ 10¹–10³; kissing bonds (disbonded but contacting)
give κ ~ 10⁻¹–10¹. **These ranges must be sourced**; the disbond literature has measured values.

### 9. Analytical programme
(i) Derive the modal scattering matrix for the spring joint, including the evanescent basis on both
sides.
(ii) Show that the free-edge result is recovered as κ → 0 and the rigid-bond result as κ → ∞
(two consistency limits, both analytical).
(iii) Derive the threshold κ* from the near-degeneracy condition.
(iv) Analyse the special case of a purely tangential spring (which allows slip but prevents
separation) — physically the most relevant imperfection and the least studied.

### 10. Numerical programme
Modal projection with pseudo-inverse, swept over κ, both spring components, and frequency.

### 11. Validation strategy
- Recovery of both limits (free and rigid) — this is a very strong check because both are
  independently known.
- Reciprocity of the joint scattering matrix.
- Energy conservation for the lossless spring.

### 12. Falsification / go-no-go (month 4)
Sweep κ over four decades. **If the efficiency varies smoothly and monotonically** with no threshold,
the hypothesis is false; the result becomes "NR degrades gracefully", which weakens the application
but is still a legitimate joint-characterisation result.

### 13. Expected results
*Positive:* a measurable tolerance criterion for negative-reflection-based monitoring of bonded
joints.
*Negative:* graceful degradation, no threshold — the practical message becomes "NR is robust".

### 14. Novelty claim and prior art
Prior art: spring interfaces in dispersion problems (extensive); imperfect joints in edge-resonance
studies; M1 at free edges. Residual claim: the joint reflection *matrix* for forward↔backward
conversion, and the threshold criterion.

### 15. Application
Guided-wave structural health monitoring: since disbonds change κ, a negative-reflection measurement
becomes a *quantitative* bond-quality probe — a concrete, industrially relevant use of the
phenomenon, and one that is much easier to deploy than free-edge retroreflection.

### 16. Feasibility risk register
1. **Result may be a monotone parameter study** (high) — this is the main risk; the threshold
   hypothesis is what makes it a science paper rather than an engineering report.
2. **Spring model validity at high frequency** (medium) — the quasi-static spring assumption fails
   when the wavelength approaches the joint dimension; must be stated as a validity limit.
3. **Crowded literature** (medium) — mitigated by the NR-specific framing.

### 17. Integrity statement
No claim that a specific joint has been measured. No claim that κ values are known for a specific
adhesive without sourcing them. The spring model's validity range is stated as a limitation.

### 18. Outputs and timeline
M1–M4: matrix formulation + limit checks. M5–M9: κ sweeps, threshold analysis. M10–M13: tangential-only
case. M14–M16: manuscript. M17–M18: sensitivity to the spring-model validity limit.

---

# PROBLEM P5 — Directional negative reflection in an anisotropic coupled plate

### 1. Title
*Anisotropy, group-velocity misalignment and the negative-reflection locus in a piezoelectric plate of
LiNbO₃.*

### 2. Research question
In an anisotropic plate the ZGV condition defines a **locus** in (k, θ) rather than a point, and the
group velocity is **not** parallel to the wave vector. **What does "negative reflection" mean when
the energy direction and the wave-vector direction differ, and can the resulting reflection be made
directional — i.e. retroreflective for one incidence direction and ordinary for another?**

### 3. Central hypothesis
**H.** The correct criterion in an anisotropic plate is on the **tangential component of the energy
velocity**, not on the wave vector, and the two disagree over a finite angular range around the
ZGV locus. Consequently there is an angular window in which a beam is *apparently* negatively
reflected (by wave-vector analysis) but is not (by energy analysis), and a correct experiment would
detect this.

### 4. Physical mechanism claimed
M1 in an **anisotropic** medium, where the geometric condition must be reformulated. This is a
*foundational clarification* problem rather than a new-effect problem, and its value is that a
correct criterion is a prerequisite for any anisotropic device claim.

### 5. Differentiation from ordinary specular reflection
Strongly differentiated: in the anisotropic case even the *definition* of the reflection angle is
ambiguous (wave-vector angle vs energy angle). The contribution is to state which definition is
physical and to show where they disagree.

### 6. Governing equations and boundary conditions
Full 3-D anisotropic elasticity reduced to the plate problem, using the existing
`Plane(k)` 6mm machinery (`code/plate_spectral.py`, validated via exact cutoffs). Free faces, free
edge, and — for LiNbO₃ — the piezoelectric coupling with short/open electrical conditions.

### 7. Governing dimensionless group
$$\Pi_2=\frac{\omega h}{c_T},\qquad \Psi=\angle(\mathbf{v}_e,\mathbf{k})\ \text{(misalignment angle)},$$
with the claim that the usable negative-reflection window is bounded by the locus Ψ = 0.

### 8. Material system and parameters
YX-LiNbO₃, for which backward waves and ZGV points are **experimentally established**
(`10.3390/s23020648`) — this is the decisive advantage of this problem: the physical existence of the
required branches is not in doubt. Standard constants; h = 0.3–1 mm.

### 9. Analytical programme
(i) Formal definition of the negative-reflection criterion in terms of the energy-velocity angle.
(ii) Show the equivalence with the isotropic criterion when isotropy holds (consistency).
(iii) Perturbative treatment of weak anisotropy around the isotropic fold.
(iv) Map the locus where the wave-vector and energy-angle criteria disagree.

### 10. Numerical programme
Anisotropic dispersion over (k, θ) with the existing rotation machinery; ZGV locus extraction;
edge scattering matrix as a function of incidence angle; energy-angle computation from the flux.

### 11. Validation strategy
- Reproduce the isotropic criterion as a limit (rotate a cubic crystal to isotropy).
- Compare the computed ZGV locus with the experimentally reported backward-wave frequency for LiNbO₃.
- Reciprocity and energy audit.

### 12. Falsification / go-no-go (month 5)
Compute Ψ over the ZGV neighbourhood. **If Ψ is negligible (< 0.1°) everywhere**, the ambiguity is
vacuous and the problem collapses; the criterion reduces to the isotropic one.

### 13. Expected results
*Positive:* a corrected criterion, an angular window of ambiguity, and a directional reflector design.
*Null:* the misalignment is negligible and the isotropic criterion transfers unchanged — a useful
clarification but a short paper.

### 14. Novelty claim and prior art
Prior art: anisotropic ZGV (Prada et al. 2009 region); backward waves in LiNbO₃ (Smirnov et al. 2023);
M1 in isotropic plates. Residual claim: the energy-angle criterion for anisotropic NR and the
directionality it enables.

### 15. Application
Direction-selective acoustic mirrors for anisotropic substrates — relevant to SAW/BAW device design
and to directional acoustic sensing.

### 16. Feasibility risk register
1. **Effect may be vacuous** (medium–high) — Gate at month 5.
2. **Computational cost of the angular sweep with full anisotropy** (medium).
3. **Algebraic complexity of the analytical part** (medium) — mitigated by the perturbative route.
4. **Overlap with existing anisotropic-ZGV work** (medium) — the NR framing is the differentiator.

### 17. Integrity statement
No claim that anisotropic negative reflection has been observed. The LiNbO₃ material choice is
justified by a published experimental result, and the citation is given, not assumed.

### 18. Outputs and timeline
M1–M5: anisotropic dispersion + ZGV locus + Gate. M6–M10: energy-angle criterion and ambiguity map.
M11–M14: edge scattering vs angle. M15–M17: manuscript. M18: device-design implications.

---

## Distinctness audit

The brief requires the five problems to be genuinely distinct. They differ in **mechanism**, not
merely in material:

| Problem | What varies | Type of claim | Risk | Novelty |
|---|---|---|---|---|
| P1 | Electrical bias (stiffness) | Does coupling preserve the condition? | Medium | High |
| P2 | Spatial grading (mode shape) | Can the condition be *engineered*? | Low | Medium |
| P3 | Higher-order BCs (structure of the BVP) | Does the condition *survive at all*? | High | **Very high** |
| P4 | Interfacial stiffness (mixed BC) | How *robust* is the condition? | Low–medium | Medium |
| P5 | Anisotropy (geometry of the criterion) | Is the criterion even *correct*? | Medium | Medium–high |

No two share a mechanism, a dimensionless group, or a mathematical difficulty. Each could be
executed without the others. Each has a pre-registered cheap falsification test. Each commits in
advance to reporting a null result.

## Problems considered and **not** included, with reasons

- **Thermoelastic / MGT control (candidate E):** physics is a perturbation plus damping; low novelty.
  Retained as a robustness section of P1.
- **Porosity / viscoelasticity (candidate I):** an ablation, not a problem. Retained as the loss check
  in P2.
- **Flexoelectric/flexomagnetic (candidate H):** structurally the same question as P3 with worse
  material data. Subsumed by P3.
- **Frequency-selective NR (candidate G):** a section of P2.
- **Non-Hermitian/EP (candidate N):** a different mechanism; already investigated as an NR problem.
- **Layered composites (candidate J):** already investigated in substance.
