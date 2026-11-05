# PHASE 2B / TASK G4-1 — Coupled Energy Identity: Derivation

Scope: 6mm Lord-Shulman piezothermoelastic bulk model as frozen in
`code/foundation_checks.py` (Task 2U-V/2U-W checkpoint 86f958b). No R/T/A
coefficients, no energy-conservation claim, no interface model. G4, G5,
SCT02 remain OPEN.

## 1. Governing equations and conventions

Plane wave exp(i(kx + ξz − ωt)); ∂x → ik, ∂z → iξ, ∂t → −iω.

- Momentum: ρ ü = div σ  →  −ρω² u = div σ.
- Gauss: div D = 0.
- Constitutives (`stresses`):
  σxx = C11 εxx + C13 εzz − e31 Ez − β1 θ, etc.; σxz = C44 γxz − e15 Ex;
  Dx = e15 γxz + ε11 Ex; Dz = e31 εxx + e33 εzz + ε33 Ez + p θ.
- Entropy/heat (characteristic row 3): ρcE θ̇ + T0(β:ε̇ − p Ė) = −div q,
  with the Lord-Shulman law (1 + τq ∂t) q = −κ ∇θ, i.e. time-harmonic
  q = −κ∇θ/(1 − iωτq).

## 2. Stored energy and flux

Total recoverable stored energy (quadratic form in the independent
variables):

  U = ½ρ|v|² + ½Cε:ε + ½εE·E + ρcE θ²/(2T0) + θ(pE − β:ε)/2.

Candidate total energy-flux vector (cycle-averaged):

  ⟨J⟩ = −½Re[σ·v* + D·(iωφ)*] + ½Re[q·θ*]/T0,

i.e. mechanical stress power, electrical power transport D·φ̇, and thermal
energy transport qθ/T0 — distinct from the bare heat flux q.

Irreversible dissipation (LS): ⟨diss⟩ = −½Re[q·∇θ*] = ½κ|∇θ|²/(1+ω²τq²) ≥ 0
for self-pairs; relaxation factor ω²τq²/(1+ω²τq²) verified positive for
physical (ωτq ≈ 2.3e-6), synthetic ωτq = 0.1, and synthetic ωτq = 1
(synthetic cases labelled synthetic, never measured properties).

## 3. Time-harmonic cycle-averaged balance

For a superposition of bulk modes at common (k, ω), each bilinear pair
(m,n) obeys the differential balance

  Re[ i(ξ_m − ξ_n*) Jz_mn ] + diss_mn = 0      (x-derivative cancels at common k).

## 4. Positive definiteness

The quadratic stored energy U is positive definite iff the elastic block
C, the dielectric block ε, and the thermal modulus ρcE/T0 are positive
definite and the pyroelectric/thermal coupling is subdominant; with the
supplied parameters (C11, C33, C44, ε11, ε33, cE > 0, coupling p, β small
relative to the diagonal blocks) U is positive definite. Verified
numerically per material in the validation report.

## 5. Finding (why G4 stays OPEN)

The oracle (`code/phase_02b_g4_1_oracle.py`) shows the balance above closes
to machine precision only when the pyroelectric/thermal coupling storage
θ(pE − β:ε)/2 is carried consistently through the RELAXED thermal equation
(the (1+τq∂t) factor acts on the coupling source as well). With the naive
flux the per-pair residual is O(1) on coupled (piezo) modes — measurable,
and the adversarial mutations M02/M03/M04/M08 each change it measurably.
Exact closure requires the full relaxed-coupling storage term; this is the
precise open obligation for G4 acceptance.
