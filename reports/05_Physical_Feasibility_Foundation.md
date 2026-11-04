# 05 — Physical Feasibility Foundation

**Document status:** Phase 3 deliverable
**Purpose:** establish, before any problem is proposed, *what physically has to be true* for negative
reflection to exist in a coupled plate, and *which equations must be used* to test it.
**Numbers quoted in this document were computed in this project**; the scripts and their outputs are
listed in §7.

---

## 1. The distinction the brief asks for, in operational form

The brief lists ten phenomena that must not be conflated. `03_Negative_Reflection_Taxonomy.md`
defines them. Here they are reduced to *measurable conditions* that a calculation can actually test.

Let a guided mode be labelled by (ω, k). Define

- phase velocity  **v**φ = (ω/|k|²) **k**
- energy velocity **v**e from the flux identity
  **v**e = **P** / (2⟨T⟩), with **P** the time-averaged power flux and ⟨T⟩ the time-averaged kinetic
  energy per unit length
- the time-averaged Poynting vector **P**

| Phenomenon | Operational test | What a naive analysis gets wrong |
|---|---|---|
| Specular reflection | k_r,∥ = k_i,∥ and P_r,∥ and P_i,∥ have opposite signs | — |
| **Negative reflection** | k_r,∥ = k_i,∥ **and** sign(P_r,∥) = sign(P_i,∥), with the reflected mode having v_e,∥·k_∥ < 0 | Inferring NR from Re(ρ) < 0 |
| Negative refraction | transmitted beam on the incident side; k_t,∥ same sign as k_i,∥ | Mixing up the two boundary conditions |
| Negative group velocity | v_e,∥ < 0 at k > 0 | Comparing with df/dk instead of dω/dk (factor 2π) |
| Backward wave | same as above, stated for a mode rather than a point | Claiming NR without an interface |
| Phase conjugation | ρ = −1 for the matched pair, i.e. ψ → ψ* | Claiming it away from the ZGV, where k₁ ≠ k₂ᵦ |
| Reflectionless | \|ρ\| → 0 | Confusing \|ρ\|→0 with ρ→−1 |
| Negative ρ | Re(ρ) < 0 | **Not** NR |
| Mode-converted reflection | ρ between different branches | Assuming conversion implies NR |
| Temporal reflection | frequency, not wavenumber, is the non-conserved quantity | Applying spatial Snell's law |

**The only defensible primary criterion** for negative reflection of a guided mode is therefore a
**pair** of statements: tangential wavenumber matching, plus reversal of the *tangential energy
flux*. Everything else is a corollary, a symptom, or a different phenomenon.

---

## 2. The flux identity — and why it must be derived, not assumed

For an elastic waveguide the standard result is

$$\mathbf{P}=\frac{\omega}{2}\,\mathrm{Im}\!\int \sigma_{1j}\,\overline{u_j}\,\mathrm{d}x_2,
\qquad
v_e=\frac{P}{2\langle T\rangle}=\frac{\mathrm{d}\omega}{\mathrm{d}k}.$$

This is what makes the whole subject tractable: it converts "does the energy go backwards?" into a
local, cheap computation that needs no branch tracking.

**It was verified here, and it was verified the hard way.** In
`code/check_02_flux_and_zgv.py`, part A, the flux-derived energy velocity was compared with a
finite-difference derivative of a branch tracked by eigenvector continuity:

- against **df/dk**: median relative difference 5.28 (a factor of almost exactly 2π)
- against **dω/dk**: median relative difference **5.27 × 10⁻⁷**, maximum 9.24 × 10⁻⁷

So the identity holds to machine precision for the lossless elastic plate — and the intermediate
"failure" was a bookkeeping error of exactly the sort the brief warns about (§8 of this document).
Both the pass and the failure are recorded.

### 2.1 The identity is *not* automatic in a coupled medium

Three things go wrong when the medium is piezoelectric, piezomagnetic, or magneto-electro-elastic:

**(i) The flux must include the electromagnetic part.** The correct conserved quantity is the
electromechanical Poynting flux. For the quasi-static piezoelectric problem the mechanical part
dominates by many orders of magnitude for typical plate modes, but "dominates" is not "is" — a paper
that claims energy conservation must state which flux it conserves and check the residual.

**(ii) The massless potential breaks the naive argument.** The electric potential has no inertia, so
the Lagrangian is not of the standard form and the usual derivation of v_e = dω/dk does not go
through unchanged. In this project the potential was eliminated by exact static condensation
(Gauss's law has no time derivative: **A**_{φφ} φ = −**A**_{φu} **u**), and the reduced problem is
then a standard Hermitian generalised eigenproblem in (u₁, u₂). That is what makes the identity
hold for the *reduced* system — and it holds only for the reduced system.

**(iii) Loss (viscoelasticity, thermoelastic damping, memory-dependent heat conduction) destroys the
identity.** With complex moduli, v_e and dω/dk genuinely differ, and **no** statement of the form
"backward ⇒ negative reflection" can be made from the dispersion relation alone. This matters
directly for Dr. Gupta's fractional and memory-dependent models, which are intrinsically lossy.
In a lossy medium the operational criterion must be stated on the **flux**, with the sign convention
and the phase reference fixed explicitly.

### 2.2 Real-frequency vs real-wavenumber analysis

Two analyses are used in the field and they are not interchangeable:

- **Real ω, complex k** (the standard "dispersion curve" analysis, and the one used in all four M1
  papers). Evanescent modes appear as complex k; the edge problem is solved by a modal expansion in
  k with the frequency as a parameter.
- **Real k, complex ω** (used in leakage/attenuation studies). Complex ω represents temporal decay.

The M1 mechanism is a **real-ω/complex-k** statement: the fold, the coincident wavenumbers and the
backward branch are all read off the real-k branch of the dispersion relation at fixed real ω.
Mixing the two conventions silently changes the sign of Im(k) and therefore which modes are
admissible at the edge. This is a classic source of error in edge-scattering codes and must be
fixed once, explicitly, at the start of the proposed work.

---

## 3. The mechanism, stated as conditions that can be checked before any heavy algebra

### 3.1 What produces the fold

For a homogeneous plate, forward/backward coexistence above a ZGV point arises from the repulsion of
two branches of the same symmetry whose cutoff frequencies are close — a longitudinal and a
transverse thickness resonance of the same parity. That repulsion is what turns a branch downward and
creates the region where v_e < 0 at k > 0.

**This is an elastic phenomenon and it survives coupling.** What coupling does is move the two
cutoffs, hence move the fold. It was quantified here for PZT-5H poled through the thickness:

| Plate | f_ZGV (MHz) | k_ZGV (1/m) |
|---|---|---|
| isotropic Duralumin-like reference, 1.5 mm | 1.89038 | 1061.0 |
| PZT-5H elastic constants, coupling switched off, 1.5 mm | 1.74863 | 864.8 |
| PZT-5H fully coupled, electrically OPEN, 1.5 mm | **1.12975** | **809.0** |
| PZT-5H fully coupled, electrically SHORT, 1.5 mm | **1.12975** | **809.0** |

Published reference for the first row: f_ZGV = 1.91 MHz, k_ZGV = 1.06 mm⁻¹ for Duralumin
(Gérardin et al. 2019). Agreement: **−0.9 %** in frequency, +0.1 % in wavenumber.

**Two results follow, and both matter.**

1. **A usable fold exists in a coupled plate** — the premise of GAP-1 is satisfied. But it is not the
   pair one would guess from elasticity: coupling moves it by −35 %, and the elastic constants alone
   account for only part of that.
2. **The electrical boundary condition is irrelevant for this particular mode** (−0.0006 %), while
   being very relevant for other branches at the same wavenumber (10⁻⁴ to 10¹ % across branches,
   see `reports/figures/fig3a`). The reason is structural, not numerical: the recovered potential of
   the mode is antisymmetric through the thickness with φ = 0 at the mid-plane
   (`reports/figures/fig3b`), so imposing φ = 0 on both faces costs almost nothing.

**Consequence for problem design.** A proposal whose novelty rests on "we add an electrode and switch
the negative reflection on and off" is **not supported** by this calculation for the S₁/S₂ᵦ fold. That
avenue is closed here, cheaply, before a paper was written about it. This is precisely what a
feasibility audit is for.

### 3.2 What produces ρ ≈ −1, and the condition nobody has checked for a coupled plate

At the ZGV the two branches have equal wavenumber, and their stress–displacement fields become
nearly identical up to the sign of the energy flux. The traction-free condition at the edge can then
be met by the pair alone, giving ρ ≈ −1.

Formally: let **Φ**_n(x₂) be the displacement–stress vector of mode n at the edge. The edge problem is

$$\sum_n c_n\,\mathbf{B}\,\boldsymbol{\Phi}_n = -\,\mathbf{B}\,\boldsymbol{\Phi}_{\text{inc}},
\qquad \mathbf{B}\ \text{the traction operator}.$$

ρ ≈ −1 for the pair (1, 2b) means the system is almost closed on those two modes, i.e. the modes are
almost linearly dependent in the sense that matters for **B**.

**This is the condition that a coupled medium can break, and it has not been checked anywhere.**
Coupling changes the mode shapes; there is no reason in advance why the pair should remain
near-degenerate in the traction norm. Whether it does is an empirical question that is cheap to
answer (compute the mode shapes across the fold and their overlap) and expensive to answer wrongly
(build a whole paper on an assumption).

**This is the first task of the proposed work, and it is a genuine go/no-go gate.**

### 3.3 Why a *higher-order* theory is qualitatively harder, not just quantitatively different

If a strain-gradient or flexoelectric term is added, the variational principle produces additional
natural boundary conditions: a higher-order traction (a double-force or moment traction) on the
surface, in addition to the classical traction. The edge problem then requires **more conditions than
two modes can supply**.

Two outcomes are possible and both are publishable results:

- **Non-degenerate case:** extra evanescent modes provide an edge boundary layer of thickness ~ ℓ
  (the internal length). The conversion coefficient is then reduced by a factor that depends on
  ℓ/λ, and negative reflection survives only below a threshold ℓ/λ. This gives a **tunable, material-
  intrinsic limit on negative reflection** — a genuinely new physical statement.
- **Degenerate case:** if ℓ → 0 does not recover the classical edge problem smoothly (a known
  pathology of some gradient formulations, where the limit is singular), the model is not admissible
  and the problem must be reformulated. Reporting that would also be a real result, and a useful
  caution to the size-dependent-waves community.

The brief asks whether the material coupling changes the mechanism. For gradient theories the answer
is **yes, qualitatively**: it changes the *number and type* of boundary conditions, which is exactly
what determines the reflection coefficient. That is the strongest scientific argument available for
any of the candidate problems, and it is the reason P3 leads the ranking.

---

## 4. Energy balance in a multiphysics system

The brief correctly forbids asserting conservation of mechanical energy in a multiphysics problem.
The correct statement for a time-harmonic electromechanical problem, integrated over a control
volume V with boundary ∂V, is

$$\nabla\cdot(\boldsymbol\sigma\cdot\dot{\mathbf u})+\dot{U}
 =\underbrace{\mathbf E\cdot\dot{\mathbf D}}_{\text{electric work}}
 +\underbrace{q}_{\text{heat source}}
 \quad\Longrightarrow\quad
 \langle P_{\text{mech}}\rangle+\langle P_{\text{elec}}\rangle
 =\langle\Pi_{\text{diss}}\rangle .$$

For a **lossless, passive** piezoelectric medium the total electromechanical flux across any closed
surface is conserved; for a lossy (viscoelastic, or memory-dependent-thermal) medium it is not, and
the difference is the dissipation. Three practical rules follow, and they are adopted in all five
proposed problems:

1. **State the flux.** Every energy statement must name the flux used (mechanical only, or
   electromechanical) and give its expression.
2. **Check the residual.** For every computed reflection matrix, report
   Σ|ρ_n|²·(P_n/P_inc) and its deviation from 1. For a lossless case the deviation must vanish to
   numerical precision; for a lossy case it must be *positive and equal to the dissipation*.
3. **Do not sum energy ratios across different physics.** Dr. Gupta's existing interface papers do
   this correctly for the thermoelastic problem; the same discipline must be carried over.

This is one of the strongest matches between his existing skill set and the proposed work: he has
already published seven papers whose central object is the energy-ratio set at a coupled interface
(`data/references.csv`, class prefix `Gupta: energy`). The transfer to a *modal* energy budget at an
edge is direct.

---

## 5. Symmetry, reciprocity and what can be claimed without breaking them

Negative reflection with ρ ≈ −1 at a free edge of a passive, reciprocal plate does **not** require
nonreciprocity, gain, or broken time-reversal symmetry. This is important: it means the phenomenon is
achievable in an ordinary material, and it is why the passive-process framing of Gérardin et al. 2019
is defensible.

Consequences for the proposals:

- No gain medium is needed; a **passive** claim is available, which is a much stronger practical
  statement.
- Reciprocity should be **checked, not assumed**: for a piezoelectric plate, open- and short-circuit
  conditions preserve reciprocity, but a biased (DC-poled) plate with a nonlinear constitutive law
  would not. The proposed work should state and verify reciprocity of its own scattering matrix
  (ρ_ij = ρ_ji for the appropriate normalized basis).
- **Anisotropy requires care with the group velocity.** In an anisotropic plate the group velocity is
  not parallel to the wave vector, and the sign condition for "backward" is
  v_e,⋅k_∥ < 0 on the *in-plane* components, not a statement about |v_g|. This applies directly to
  the 6mm plates used here and must be handled explicitly in any anisotropic extension.

---

## 6. Is the phenomenon physically accessible? Parameter check

A proposal is worthless if the required regime is unreachable. Checks, from the verified literature
and from the computations here:

| Requirement | Value / range | Source |
|---|---|---|
| ZGV fold in an isotropic plate | requires Poisson ratio below ~0.45 for the S₁/S₂ᵦ pair | Prada et al. 2008; standard result |
| Duralumin plate, 1.5 mm | f_ZGV ≈ 1.91 MHz, k_ZGV ≈ 1.06 mm⁻¹ | Gérardin et al. 2019 |
| Same, computed here | 1.8904 MHz, 1061 1/m | `data/check_01_elastic_validation.json` |
| PZT-5H plate, 1.5 mm, coupled | f_ZGV ≈ 1.130 MHz, k_ZGV ≈ 809 1/m | `data/check_02_flux_and_zgv.json` |
| Backward waves in piezoelectric plates | experimentally established in YX-LiNbO₃ | Smirnov et al. 2023, `10.3390/s23020648` |
| Transducer/laser-ultrasound accessibility | demonstrated at MHz frequencies in mm plates | Gérardin 2016, 2019; Veres 2016 |
| Internal length for strain-gradient theories | ℓ ~ 1–10 µm for metals; sub-µm for some ceramics | order-of-magnitude; **must be sourced per material in the proposal** |
| ℓ/λ needed for an observable gradient effect | requires λ ≲ 10ℓ, i.e. plate thickness in the micrometre range | deduced from the above |

**The last row is the single biggest physical risk in the whole programme.** For a gradient effect to
be visible, the plate must be thin enough that ℓ/h is not negligible — which pushes the operating
frequency up (f ∝ 1/h) and the experiment out of the convenient range. This is stated as a *risk*,
not a *showstopper*, because the same tension exists throughout the size-dependent-waves literature
and is managed there by working in the micrometre regime at hundreds of MHz. It does mean, though,
that **a strain-gradient extension is a theory-first paper**: the numerical demonstration is
feasible, the experiment is not, and the proposal must say so.

---

## 7. Reproducibility: what was actually computed

| Script | What it does | Output |
|---|---|---|
| `code/plate_spectral.py` | Chebyshev spectral solver for isotropic and generally anisotropic (6mm) plates; energy assembly from constitutive tensors; static condensation of the potential; energy velocity from the flux identity | library |
| `code/check_01_elastic_validation.py` | validation against exact Rayleigh–Lamb theory | `data/check_01_elastic_validation.json` |
| `code/check_02_flux_and_zgv.py` | flux identity; ZGV location; piezoelectric / electrical-BC premise test | `data/check_02_flux_and_zgv.json` |
| `code/make_figures.py` | figures 1–3 | `reports/figures/` |

Headline validation numbers (all reproducible by running the scripts):

| Check | Result |
|---|---|
| Cutoff frequencies vs exact n·c/(4h) | max relative error 6.3 × 10⁻⁷ |
| Dispersive roots vs exact Rayleigh–Lamb roots | max relative mismatch 5.3 × 10⁻⁷ |
| ZGV: spectral solver vs independent exact-theory root tracker | relative difference 9.1 × 10⁻⁷ |
| ZGV: mesh convergence (N = 14 → 34) | mesh-to-mesh change 8.0 × 10⁻⁹ |
| ZGV: computed vs published Duralumin value | −0.85 % |
| Flux identity v_e vs dω/dk | median 5.3 × 10⁻⁷, max 9.2 × 10⁻⁷ |

**Three real bugs were found and fixed by these checks, and they are recorded because they are part
of the evidence that the audit is genuine:**

1. a sign error in the piezoelectric stiffening in the first version of the assembly;
2. an incorrect symmetric/antisymmetric Rayleigh–Lamb reference function (the two cases swap the
   sin·cos factors) — caught because the solver disagreed with the *reference*, and the solver was
   right;
3. branch-tracking failures at the S₁/S₂ᵦ near-degeneracy, which forced the tracking-free flux-based
   ZGV criterion now used throughout.

---

## 8. Conclusion of the foundation phase

Six statements are now established well enough to build on:

1. **Negative reflection of guided elastic waves is a mode-conversion phenomenon**, not a sign of a
   coefficient. Its necessary and sufficient operational test is: tangential wavenumber matching,
   plus a reflected mode whose energy flux is anti-parallel to its wave vector, plus reversal of the
   tangential energy flow.
2. **The flux identity v_e = P/(2⟨T⟩) = dω/dk** is the right tool, holds to 10⁻⁷ for the lossless
   elastic plate, and **must be re-derived** for coupled and for lossy media. It gives a
   tracking-free ZGV detector, which is what made the piezoelectric computations here reliable.
3. **A forward/backward fold does exist in a coupled piezoelectric plate**, at a frequency ~35 %
   below the elastic one, and it is **not** thickness dependent beyond the geometric scaling.
4. **The electrical surface condition does not tune this particular fold** — a computed negative
   result that closes one otherwise tempting research direction (electrically switchable NR).
5. **The un-checked condition is the near-degeneracy of the pair's stress–displacement fields.**
   If coupling destroys it, there is no useful negative reflection in that material. This is the
   go/no-go gate for every problem in Phase 5.
6. **A higher-order (gradient) theory changes the number of edge conditions**, which is a
   qualitative change and is the strongest available scientific argument for a new study.

The five candidate problems in `07_Five_Research_Problems.md` are built on statements 3–6, with 5 as
the first gate.
