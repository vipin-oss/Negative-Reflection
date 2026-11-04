# Phase 2A Task 2 — Two-Half-Space Scattering Assembly and Classical Validation

**Branch** `arena/a3ede4fe-negative-reflection`
**Start commit** `34df97733bcaad4e554a69b94d34ea51597663c6`
**Scope** classical elastic validation with an independent oracle; coupled
piezothermoelastic assembly as a **DIAGNOSTIC** result only.

---

## 1. Executive summary

| Question | Answer |
|---|---|
| Is the elastic reference trustworthy? | Yes — validated against closed-form Fresnel, Snell, acoustic and grazing limits; worst amplitude difference `3.99e-52` |
| Is `BLK-02A-04` closed on independent evidence? | Yes — 14/14 lossless branches, 0 failures on every check |
| Do classical P–SV amplitudes match an independent reference? | Yes — 88/88 benchmark checks, worst residual `5e-42` |
| Are incoming/outgoing modes selected correctly? | Yes for the elastic problem; **two defects found and fixed** in the coupled classifier (below) |
| Is the coupled piezothermoelastic R/T/A validated? | **No. G4 and G5 remain OPEN.** Every coupled amplitude is labelled DIAGNOSTIC |
| Did the fixes change published evidence? | They changed classification of 4 + 42 inventory rows. The Task 1R artifact is **unchanged**; the divergence is measured and published |

---

## 2. The mandatory prerequisite: an independent elastic oracle

Per the brief, the reference had to be derived independently of the
classification code. `code/scattering/elastic_reference.py` builds its own
modes, tractions and fluxes, and is validated against closed forms rather than
against the existing solver.

### 2.1 The legacy reference is defective — confirmed, not assumed

`code/elastic_psv.py` returns sign-flipped transmitted amplitudes at **every**
angle:

| θ | tP correct | tP existing | tS correct | tS existing |
|---|---|---|---|---|
| 0 | +0.5338019843 | −0.5338019843 | 0 | 0 |
| 15 | +0.5327599438 | −0.5327599438 | −0.0093422831 | +0.0093422831 |
| 30 | +0.5291290046 | −0.5291290046 | −0.0151649300 | +0.0151649300 |
| 45 | +0.5204900820 | −0.5204900820 | −0.0134289241 | +0.0134289241 |
| 60 | +0.4970348428 | −0.4970348428 | −0.0005768053 | +0.0005768053 |

The reflected amplitudes agree exactly. **Energy is blind to this defect**:
`|R + T − 1| < 1e-16` for every row, which is why `check_I1_elastic_psv` (which
checks only `R`, `T` and `|R+T−1|`) never caught it. The displacement jump at
the weld is **1.074** while the solver's own residual field reports `~1e-16`.

### 2.2 The corrected reference reproduces closed forms

| Quantity | Result |
|---|---|
| Normal-incidence Fresnel `r`, `t` | `0.466198015697`, `0.533801984303`, `|d| = 3.99e-52` |
| Physical-form boundary residuals | `ux = 0.0`, `uz = 2.5e-51`, `txz = 0.0`, `tzz = 0.0` |
| Oblique energy closure 0°–80° | `|R + T − 1| ≤ 2.7e-51` at every angle |
| Acoustic limit vs scalar Fresnel | worst `1.23e-9` over 0°–60° |
| Grazing | `rP = −0.99109` at 89.9°, `−0.99910` at 89.99° |
| S-incidence normal | `rS = 0.498728029328`, vs closed form `|d| = 3.07e-52` |
| Row-equilibrated condition number | 1.766 (0°), 2.005 (30°), 2.012 (60°), 2.535 (80°) |

Row equilibration is applied before the condition number is quoted: the raw
matrix mixes metres and pascals and reports `~2.5e14` for a well-behaved
problem, which is a statement about units, not about degeneracy.

### 2.3 `BLK-02A-04` — closed

`validation/scattering/lossless_oracle.py`, 14/14 branches:

| check | failures |
|---|---|
| `lossless_ok` | 0 |
| `flux_nonzero` | 0 |
| `role_ok` | 0 |
| `class_ok` | 0 |

Worst polarisation residual exactly `0.0`; worst dispersion residual `2.79e-51`.
All 14 are S-type, `Im ξ = 0` exactly, and the flux sign matches the recorded
`radiation_dir` **with no reference to `sign(Im ξ)`** — which is the
independence the brief demanded.

---

## 3. The clean physical state

`s = (u_x, u_z, φ, ϑ, σ_xz, σ_zz, D_z, q_z)^T`, derived from the governing
equations under `e^{-iωt}`, `e^{i(kx+ξz)}`, plane strain:

```
e_xx = i k u_x          e_zz = i ξ u_z        γ_xz = i(ξ u_x + k u_z)
E_x  = −i k φ           E_z  = −i ξ φ
σ_zz = C13 e_xx + C33 e_zz − e33 E_z − β3 ϑ
σ_xz = C44 γ_xz − e15 E_x
D_z  = e31 e_xx + e33 e_zz + ε33 E_z + p ϑ
q_z  = −κ33 (i ξ ϑ)/(1 − i ω τ_q)
⟨J_z⟩ = −½Re[σ_xz (iωu_x)* + σ_zz (iωu_z)* + iω D_z φ*] + ½Re[q_z ϑ*]/T0
```

All three limits are exact (`relative error 0.0`) for Al6061 and PZT-5H.

**One correction of substance.** PZT-5H initially failed the elastic reduction
with `tzz_rel = 0.0087`. This is **not a defect**: PZT-5H is transversely
isotropic (`C13 = 8.41e10 ≠ C12 = 7.95e10`, `C33 = 1.17e11 ≠ λ + 2μ = 1.261e11`),
so the isotropic Hooke identity does not apply to it. The check now gates on
isotropy and falls back to a hand-expanded 6 mm law
(`tzz_ref = C13 e_xx + C33 e_zz`), recording `isotropic_medium` and `basis` in
the output.

---

## 4. Classical elastic validation

`validation/scattering/classical_benchmarks.py` — **88 checks, 88 pass**.

| case | what it establishes | result |
|---|---|---|
| B1 | identical welded half-spaces, normal incidence | 3/3 |
| B2 | identical welded half-spaces, oblique | 15/15 |
| B3 | dissimilar media, **complex amplitudes** at normal incidence | 8/8 |
| B4 | oblique P–SV with mode conversion (amplitude + Snell + energy) | 12/12 |
| B5 | acoustic limit vs scalar Fresnel | 10/10 |
| B6 | grazing limit, per-angle stated bounds | 3/3 |
| B7 | critical angles on a pair that has them | 20/20 |
| B8 | invariance under acoustic-impedance scaling | 3/3 |
| B9 | reciprocity, stated correctly | 9/9 |
| B9b | scalar reciprocity control | 3/3 |
| B10 | legacy defect reproduced and pinned | 2/2 |

### 4.1 Four failures were defects in my own tests, not the solver

I record these because the brief requires honest diagnosis rather than
tolerance adjustment.

1. **`str()` round-trip noise.** `mp.nstr(x, 30)` then `mp.mpf(str)` re-introduces
   `2.7e-51`. Comparisons are now made on `mpf` objects directly.
2. **B4 Snell test compared two different quantities.** I computed
   `k/(ω/cL₁)` and compared it to `k/(ω/cL₂)` — the sines of the *incident* and
   *transmitted* angles, which differ by construction. The correct invariant is
   `k = (ω/cL₁) sin θ₁`, now tested directly.
3. **B7 tested a material pair with no critical angle.** Al→steel has
   `cL(Al) = 6312.8 > cL(steel) = 6002.1`, so no transmitted branch can go
   evanescent. The test now uses a soft solid → steel (`cL = 2000`), where the
   predicted critical angles are 19.464° and 38.086° and both fire exactly.
4. **Working precision was set after the media were constructed.** `AL`/`ST`
   were built at import time while mpmath was still at 15 digits, then `dps`
   was raised to 50 inside `main()`. Every "high-precision" agreement was
   really a 15-digit one, which is why the Fresnel check failed at `1.7e-17`.
   Setting `mp.mp.dps` before construction fixed all eight affected checks.

### 4.2 B9 — a refuted expectation, recorded rather than deleted

`r₁₂ + r₂₁ = 0` and `|R₁₂| = |R₂₁|` are **both false** for P–SV with mode
conversion. Measured `|R₁₂ − R₂₁| = 1.4e-3 / 6.0e-3 / 2.1e-3` at 10°/30°/45°,
with `R₁₂ = 0.2027923177` vs `R₂₁ = 0.2013830537`.

The cause is physical: the reflected field is a **coherent P+SV pair**, and the
way energy splits between the two branches inverts with the direction of
incidence. A scalar control settles it — in the acoustic limit, where there is
one reflected branch, `R₁₂ − R₂₁` falls to `1e-10` (case B9b). The residual is
retained as a stated, explained quantity, not widened away.

### 4.3 A real degeneracy: the acoustic limit

As `C44 → 0`, `|rS|` **does not converge** — it stays at `0.300508` while the
S energy share falls `3.18e-8 / 3.18e-10 / 3.18e-12` for
`C44 = 1e-2 / 1e-6 / 1e-10`. The S mode becomes a zero-traction null direction
and the `(u,t)` 4×4 loses rank. Only `rP` is meaningful in that limit, and the
record says so. Condition number is the wrong detector here (after equilibration
it stays near 2.2), so **energy share** is used.

---

## 5. The coupled assembly — degrees of freedom counted, not assumed

The brief forbids assuming an 8×8 system. The count is computed:

| configuration | θ | roots | gauges | admitted (1 / 2) | equations | unknowns | rank | square |
|---|---|---|---|---|---|---|---|---|
| Al6061 \| PZT-5H | 0° | 16 | 2 + 2 | 3 / 3 | 7 (of 8) | 6 | 6 | no |
| Al6061 \| PZT-5H | 30° | 16 | 0 + 0 | 4 / 4 | 8 | 8 | 8 | **yes** |
| PZT-5H \| PZT-5H | 0° | 16 | 2 + 2 | 3 / 3 | 7 (of 8) | 6 | 6 | no |
| PZT-5H \| PZT-5H | 30° | 16 | 0 + 0 | 4 / 4 | 8 | 8 | 8 | **yes** |
| Al6061 \| Al6061 | 0° | 16 | 2 + 2 | 3 / 3 | 4 (of 8) | 6 | 4 | no |

Where the physics supplies a square system, an 8×8 system is used. Where it does
not, the reduction is reported:

- **At normal incidence** the sagittal geometry sets `u_x = σ_xz = 0` for every
  mode, so those continuity equations read `0 = 0`. They are **identified and
  removed**, with the removal recorded (`rows_removed`), not kept with a
  tolerance and not repaired by adding unknowns.
- **Al6061 | Al6061 at normal incidence is genuinely under-determined**: only 4
  independent equations for 6 unknowns, condition number `1.1e21`. This is
  physical (an elastically identical interface is invisible), and it is
  reported as unresolved rather than regularised.

### 5.1 Two classification defects found and fixed

**(a) The decoupled Laplace electrostatic mode was being discarded as a gauge.**

`classify_mode` labelled any mechanically-decoupled mode `GAUGE`. But at
`k ≠ 0` the pair at `ξ = ∓ik` is the Laplace solution, with
`D_z = −ε₃₃ ξ φ = ∓ i ε₃₃ k φ ≠ 0` — measured exactly `−6.954e-9 = −ε₃₃k` for
Al6061 at 30°. That mode satisfies the `D_z` continuity equation like any other
mode. Calling it a gauge **deleted the electrostatic channel of an elastic
half-space while keeping the equivalent channel of a piezoelectric one** — an
exclusion of a physical channel, which the brief explicitly forbids.

The test is now `D_z = −ε₃₃ ξ φ`: zero ⇒ genuine gauge (only at `k = 0`,
verified); non-zero ⇒ `ELECTROSTATIC`, admissible.

**(b) `Im ξ ≠ 0` was treated as a cusp, which removed every incident channel
from PZT-5H.**

`classify_role` tested `Im ξ != 0` strictly. Thermomechanical coupling puts every
acoustic mode *just* off zero: PZT-5H's shear wave at 30° has
`Im ξ = 1.05e-9` against `Re ξ = 3.5e3` — an attenuation of **3e-13 per
wavelength** and a decay length of **9.6e8 m** for a 1.8 mm wave. The strict rule
declared it localised on the far side, so PZT-5H had **no incoming mode at all**
at `θ > 0` and could not be excited.

The rule is now a stated physical criterion, `ATTENUATION_PROPAGATING = 1e-4`
per wavelength: the largest affected case is `3e-13` and the weakest genuinely
damped mode is `6.7e-11`, so nothing is near the boundary. The threshold appears
in every reason string it produces.

### 5.2 A sign error that a converged solve did not reveal

The first version built `[B1 +B2]` while its own docstring said `[B1 −B2]`. The
solve then satisfied `s₁ + s₂ = −s_inc` instead of `s₁ − s₂ = −s_inc`, and the
linear-algebra residual was still `1e-37`, so it **looked converged**.
Recomputing the interface field independently exposed a residual of exactly
**2× the incident state** (measured `114.1` on `ϑ`). After the fix, interface
residuals are `2.7e-51`.

This is the concrete demonstration of why interface residuals are recomputed
from the modal expansion rather than read off the solver.

### 5.3 What the coupled solver now shows

| configuration | θ | eq×unk | cond (equilibrated) | interface residual | rescale verdict |
|---|---|---|---|---|---|
| Al6061 \| PZT-5H | 0° | 7×6 | 5.2e3 | 6.2e-45 | INVARIANT |
| Al6061 \| PZT-5H | 30° | 8×8 | 2.2e4 | 2.7e-51 | INVARIANT |
| Al6061 \| PZT-5H | 60° | 8×8 | 3.1e4 | 2.7e-51 | INVARIANT |
| PZT-5H \| PZT-5H | 0° | 7×6 | 4.1e3 | 3.3e-47 | INVARIANT |
| PZT-5H \| PZT-5H | 30° | 8×8 | 2.8e4 | 2.7e-51 | INVARIANT |
| PZT-5H \| PZT-5H | 60° | 8×8 | 3.1e4 | 2.7e-51 | INVARIANT |
| Al6061 \| Al6061 | 0° | 4×6 | 1.1e21 | — | DEFERRED_TO_CONDITIONING |
| Al6061 \| Al6061 | 30° | 8×8 | 1.3e4 | 2.7e-51 | INVARIANT |
| Al6061 \| Al6061 | 60° | 8×8 | 2.9e4 | 2.7e-51 | INVARIANT |

**Rescaling invariance** is the sharpest validity test available. A modal basis
is determined only up to a complex scalar per mode, so the physical field the
basis reconstructs must not move when every column is multiplied by an arbitrary
complex factor. Eight of nine configurations are **INVARIANT at 1e-38 or
better**. The ninth is deferred to its conditioning, explicitly, rather than
counted as a pass.

Two bugs in my *test* had to be fixed before this number was meaningful: the
comparison was originally made per-coefficient against a coefficient that is
exactly zero for identical media (reporting `54.2` for a quantity of order
`1e-37`), and an earlier version renormalised the scaled columns, which
cancelled the scaling and made the test a no-op.

---

## 6. Regression: restoring Task 1R

The two classification fixes changed labels for the frozen inventory. Rather
than rewrite a published artifact, the legacy rules were made reachable through
explicit keyword arguments (`legacy_gauge=True`, `strict_legacy=True`) and the
divergence was measured:

| change | rows |
|---|---|
| `spatial_class`: `GAUGE` → `ELECTROSTATIC` | 4 |
| side role: `EVANESCENT_NO_TRANSPORT` → `INCOMING` | 42 |
| side role: `GAUGE_EXCLUDED` → `OUTGOING` | 4 |
| side role: `GAUGE_EXCLUDED` → `EVANESCENT_NO_TRANSPORT` | 4 |

All 42 propagation-band rows have attenuation in `6.14e-14 … 4.57e-7` per
wavelength, i.e. decay lengths of `8.4e5 … 1.6e7` m. **The inventory file is
byte-identical to its published sha256** and `INV00` continues to enforce that.

`validation/scattering/convention_divergence.py` →
`validation/phase_02a_task2_convention_divergence.json`.

---

## 7. Gate status

| gate | status | basis |
|---|---|---|
| classical normal-incidence Fresnel amplitudes | **PASS** | `|d| = 3.99e-52` vs closed form |
| oblique P–SV complex amplitudes | **PASS** | 88/88, worst `5e-42` |
| incoming/outgoing selection | **PASS** | independent oracle, 0 role failures |
| gauge and degenerate handling | **PASS** | 4 + 42 fixes, both hazards as negative controls |
| boundary residuals and equation counts | **PASS** | 108/108 assembly checks |
| precision convergence | **PASS** | dps 30/50/80 agree to `1e-10` |
| every published field matches its source | **PASS** | INV01 120/120 under the CSV's own convention |
| **G4 — coupled energy conservation** | **OPEN** | no independent coupled energy gate passed |
| **G5 — coupled thermodynamic consistency** | **OPEN** | no independent coupled thermodynamic gate passed |

Coupled amplitudes remain **DIAGNOSTIC**. A converging boundary solve is not
evidence that `R`, `T` and `A` are validated piezothermoelastic quantities.

---

## 8. Files delivered

| file | role |
|---|---|
| `code/scattering/state.py` | physical state, dimensional scales, three limit checks |
| `code/scattering/elastic_reference.py` | independent elastic oracle, closed-form validated |
| `code/scattering/assembly.py` | degrees-of-freedom counting, mode selection, row reduction |
| `code/scattering/solver.py` | coupled solve, interface residuals, rescaling invariance |
| `code/scattering/modal.py` | GAUGE / ELECTROSTATIC separation (modified) |
| `code/scattering/radiation.py` | scale-aware propagation criterion (modified) |
| `validation/scattering/lossless_oracle.py` | `BLK-02A-04` closure evidence |
| `validation/scattering/classical_benchmarks.py` | 88 classical checks |
| `validation/scattering/assembly_validation.py` | 108 assembly and solver checks |
| `validation/scattering/convention_divergence.py` | measured divergence from Task 1R |
| `validation/phase_02a_task2_results.json` | machine-readable results and checksums |
| `validation/phase_02a_task2_benchmark_comparison.csv` | per-check benchmark records |
| `validation/phase_02a_task2_assembly_validation.csv` | per-check assembly records |
