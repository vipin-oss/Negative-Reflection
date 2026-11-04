# Phase 2A Task 2 — Classical Benchmark Comparison

Every record: `case, assumptions, quantity, expected, actual, abs_error,
tolerance, pass, note`. Machine-readable:
`validation/phase_02a_task2_benchmark_comparison.csv`.

**88 checks, 88 pass, 0 fail.** Working precision 50 decimal digits throughout.

---

## 1. Reference and materials

| | material | cL (m/s) | cT (m/s) | ZL (kg/m²s) |
|---|---|---|---|---|
| medium 1 | Al6061 (ρ=2700, C11=1.076e11, C44=2.65e10) | 6312.83 | 3132.86 | 1.70446e7 |
| medium 2 | steel-like (ρ=7800, C11=2.81e11, C44=8.2e10) | 6002.14 | 3242.35 | 4.68167e7 |

Frequency `ω = 2π × 10⁶ rad/s`. Conventions: `e^{−iωt}`, `e^{i(kx+ξz)}`, plane
strain, `<P_z> = −½Re[σ_xz v_x* + σ_zz v_z*]` with `v = −iωu`.

---

## 2. Results by case

### B1 — identical welded half-spaces, normal incidence (3/3)

| quantity | expected | actual | error |
|---|---|---|---|
| `rP` | 0 | 0 | 0 |
| `tP` | 1 | 1 | 0 |
| `R + T` | 1 | 1 | ≤ `2.7e-51` |

### B2 — identical welded half-spaces, oblique (15/15)

Energy closure `|R+T−1| ≤ 2.7e-51` at 10°/20°/30°/45°/60°; `rP = rS = 0` and
`tP = tS = 1` to rounding.

**Finding.** In the acoustic limit `|rS|` does **not** converge to zero. With
`C44 → 0` it stays at `0.300508` while the S energy share collapses:

| `C44` | S energy share | `|rS|` | `|tS|` |
|---|---|---|---|
| 1e-2 | 3.179e-8 | 0.300508 | 0.176803 |
| 1e-6 | 3.179e-10 | 0.300508 | 0.176803 |
| 1e-10 | 3.179e-12 | 0.300508 | 0.176803 |

The S mode is a zero-traction null direction: the `(u,t)` 4×4 loses rank and
`rS` is a singular-perturbation residue, not a coefficient. Only `rP` is
meaningful in that limit. This is recorded, not smoothed over.

### B3 — dissimilar media, normal incidence, complex amplitudes (8/8)

| quantity | closed form | solver | abs error |
|---|---|---|---|
| `rP` | 0.466198015697219111557806806473 | 0.466198015697219111557806806473 | `3.99e-52` |
| `tP` | 0.533801984302780866644268801407 | 0.533801984302780883387371419205 | `<1e-48` |
| `R` | 0.217340589840024577467038867231 | 0.217340589840024561855836433167 | `<1e-48` |
| `T` | 0.782659410159975422532961132768 | 0.782659410159975438144163566832 | `<1e-48` |
| `ux`, `uz`, `txz`, `tzz` residuals | 0 | `0.0`, `2.5e-51`, `0.0`, `0.0` | — |

B3 is the check that would have caught the Phase 1B defect: it compares
**amplitudes**, and the sign error is invisible to energy.

### B4 — oblique P–SV with mode conversion (12/12)

At 10°/20°/30°/40°, for all four channels: Snell (`k = (ω/cL₁) sin θ₁` exact, `0`),
energy closure (`≤ 2.7e-51`), and the transmitted-mode classification.

### B5 — acoustic limit vs scalar Fresnel (10/10)

Worst discrepancy `1.23e-9` over 0°–60°. First attempt gave a false alarm because
I used `Z/cos` instead of the correct
`(Z₂cos₁ − Z₁cos₂)/(Z₂cos₁ + Z₁cos₂)`.

### B6 — grazing limit (3/3)

`rP → −1` as `θ → 90°`, but **not uniformly**: `|1 + r_P|` falls roughly like
`(90° − θ)`. Measured `rP = −0.9144` at 89°, `−0.99109` at 89.9°,
`−0.99910` at 89.99°. Each angle carries its own stated bound
(`1e-1`, `1.5e-2`, `1.5e-3`) tied to that decay, rather than one tolerance
chosen to fit.

### B7 — critical angles (20/20)

The Al → steel pair used elsewhere has **no** critical angle
(`cL(Al) = 6312.8 > cL(steel) = 6002.1`), so the earlier version of this case
could never fire. It now uses a soft solid (`cL = 2000`) → steel:

| θ | P-critical 19.464° | SV-critical 38.086° | tP status | tS status |
|---|---|---|---|---|
| 10.000 | below | below | propagating ✓ | propagating ✓ |
| 18.964 | just below | below | propagating ✓ | propagating ✓ |
| 19.964 | just above | below | **evanescent** ✓ | propagating ✓ |
| 37.586 | above | just below | evanescent ✓ | propagating ✓ |
| 38.586 | above | above | evanescent ✓ | **evanescent** ✓ |

Energy closes at `≤ 2.7e-51` at every angle, and the transmitted share never
exceeds 1 — independently of which channels propagate.

**One earlier failure was a genuinely wrong prediction on my part**, not a code
defect: I had probed 31.4° expecting SV evanescence, but the S-critical for the
Al→steel geometry is 30.905° only against `cT2/cL1`; the branch reported
`propagating` correctly because the pair has no critical angle at all.

### B8 — invariance under impedance scaling (3/3)

Scaling both media by a common factor leaves all amplitudes unchanged — a check
that no absolute unit has leaked into the ratios.

### B9 — reciprocity (9/9) and B9b — scalar control (3/3)

`r₁₂ + r₂₁ = 0` and `|R₁₂| = |R₂₁|` are **both false** for P–SV with mode
conversion:

| θ₁ | θ₂ | `|R₁₂|` | `|R₂₁|` | difference |
|---|---|---|---|
| 10° | 9.503° | 0.2174447926 | 0.2174968435 | 5.2e-5 |
| 30° | 28.385° | 0.2183940908 | 0.2187552643 | 3.6e-4 |
| 45° | 42.245° | 0.2206342743 | 0.2209374437 | 3.0e-4 |

Both directions close energy to `2.7e-51`, and Snell is satisfied to `0` — so
this is not numerical. The reflected field is a **coherent P+SV pair** and its
energy split inverts with the direction of incidence.

**B9b settles it.** In the acoustic limit, where there is one reflected branch:

| θ | `R₁₂ − R₂₁` |
|---|---|
| 10° | 7.6e-12 |
| 30° | 6.95e-11 |
| 45° | 1.59e-10 |

A single-branch interface **is** reciprocal. The P–SV asymmetry is physics. The
expectation is recorded as refuted rather than deleted, and the tolerance was
not widened.

I also had to correct my own Snell inversion twice during this work: the
transmitted angle is `sin θ₂ = (c₂/c₁) sin θ₁`, and inverting that ratio produced
a spurious `k` mismatch that looked like a solver bug.

### B10 — legacy defect, pinned (2/2)

`code/elastic_psv.py` still reproduces its sign error (`tP` negative at every
angle) and that evidence is retained, not deleted. The `rP` column is correct,
which is why the defect was invisible to the energy-only check.

---

## 3. Why four initial failures were mine

Recorded because the brief requires diagnosing failures honestly rather than
adjusting thresholds:

| symptom | diagnosis | fix |
|---|---|---|
| every `1e-45` tolerance failing by `2.7e-51` | `mp.nstr(x,30)` → `mp.mpf(str)` round-trip noise | compare `mpf` objects |
| B4 "Snell violation" at 10° (0.1736 vs 0.1651) | compared `sin θ₁` to `sin θ₂` | test `k = (ω/cL₁) sin θ₁` |
| B7 "S-critical mispredicted" | the pair has no critical angle | use a pair that has one |
| B3 Fresnel off by `1.7e-17` | media built at import time, `dps` raised later | set `mp.mp.dps` before construction |

The last one is the most instructive: the arithmetic *looked* like 50-digit
work but the inputs were 15-digit, so every "high-precision" agreement was
really a 15-digit one.

---

## 4. Reproduction

```
cd code && python3 -m scattering.elastic_reference
python3 validation/scattering/lossless_oracle.py
python3 validation/scattering/classical_benchmarks.py
python3 validation/scattering/assembly_validation.py
```

Expected: `BLK-02A-04` 14/14 with 0 failures; benchmarks 88/88; assembly 108/108.
