# Phase 2A — Modal Classification Audit

**Experiment ID:** `PHASE-02A-MODAL-CORRECTED-1`
**Base commit:** `520fb5dfcf594898e94271db38bb9428821d209b`
**Scope:** Phase 2A task 1 only. Tasks 2–8 not started.

This report audits `validation/phase_01d_independent_comparison.csv` and the
classifier that produced it. It does **not** modify that file. The Phase 1D CSV
remains in the tree unchanged as historical evidence of what the Phase 1D
pipeline actually output, defects included; test `MAC00` pins its sha256
(`76c7174c…1eb56b`) so the record cannot drift.

---

## 1. The audited artifact

`validation/phase_01d_independent_comparison.csv` — 120 rows, 3 materials
(PZT-5H, Al6061, CdSe) × 5 angles (0/15/30/45/60°) × 8 branches.

| `mode_class` | count |
|---|---|
| `PROPAGATING` | 92 |
| `DEGENERATE_EXCLUDED` | 20 |
| `ELASTIC_DECOUPLED` | 8 |

A single `mode_class` field was the whole classification. It carried no decay
rate, no direction, no energy status, and no confidence — so a reader of the CSV
had no way to see that a "propagating" mode decays by a factor of *e* in a third
of a micrometre.

---

## 2. Findings

### D1 — `PROPAGATING` assigned to strongly decaying modes (severe)

Of the 92 `PROPAGATING` rows, **30 have |Im ξ| > 10³ m⁻¹**. The worst:

| material | branch | Im ξ (m⁻¹) | decay length | in transverse wavelengths |
|---|---|---|---|---|
| PZT-5H | 0, 7 (all angles) | ∓2.569e6 | 3.89e-7 m = 0.39 µm | **0.159 λ** |
| CdSe | 0, 7 | ∓8.659e5 | 1.15e-6 m | 0.159 λ |
| Al6061 | 0, 7 | ∓2.167e5 | 4.61e-6 m | 1.0 λ |
| PZT-5H | 3, 4 @ 45–60° | ∓1196 … ∓1466 | ~1e-3 m | — |

A mode whose amplitude falls to 1/e in 0.159 of a transverse wavelength is not
propagating by any convention. All thirty are **thermal-diffusive** branches:
their |ξ| matches the axial thermal wavenumber
`ξ_th = sqrt(ω ρ cE / κ33)` to within 1.6 %.

The label was produced by a rule effectively equivalent to "large |Re ξ| or a
non-zero root" rather than by the derived direction convention. It also caused a
downstream consequence: these thermal branches were available for selection as
"propagating" modes in a half-space, which would silently corrupt any
scattering matrix built on top of them.

### D2 — `PHYSICAL_GROWING` is self-contradictory in a half-space

The vocabulary emitted `PHYSICAL_GROWING`. In a half-space problem a mode either
decays away from the interface (admissible) or grows without bound (not
admissible). "Physical" and "growing" cannot both hold for a half-space
eigenmode; the growth has to be attributed to a direction change or excluded.

### D3 — one field conflated four independent questions

`mode_class` encoded spatial behaviour, radiation direction, energy flow, and
numerical confidence simultaneously. Any attempt to fix one — say, spatial
behaviour — unavoidably moved the others. This is a structural cause of D1.

### D4 — no regression test existed

Nothing in the Phase 1D suite would have failed if the classifier were replaced
by a constant. The `mode_class` field was generated but never checked.

### D5 — newly found in Phase 2A: the null-vector route is precision-dependent (severe)

Not visible from the CSV, found while rebuilding the classifier. `kernel_scaled`
intermittently returns the **wrong branch** as a function of working precision:

```
Al6061, 30°, branch 3 (a thermal branch)
  dps  60..85  ->  correct branch        (|θ|/max = 1.0)
  dps  90      ->  WRONG: wire gauge     (|φ|/max = 1.0, |θ| = 0)
  dps  95..105 ->  correct
  dps 110      ->  WRONG
  dps 120      ->  WRONG
  dps 125..130 ->  correct
  dps 135      ->  WRONG
```

**Mechanism.** `kernel_scaled` equilibrates matrix *columns* with
`cs[j] = max_i |A[i,j]|`, guarded by `if cs[j] > 0`. For a centrosymmetric
medium the piezoelectric coupling is exactly zero, so at a thermal branch the
electric-potential column of M is not small — it is **zero**. Whether the
computed column maximum is `8.04e-87` (roundoff) or `0.0` (guard fires) depends
on the last digits of the root and therefore on dps. If the guard fires, that
column is left at roundoff scale and a spurious 3×3 minor wins the free-variable
search; the routine returns a different null vector — the pure wire gauge —
without raising any error.

Note the consequence for Phase 1D: the BLK-U3 closure was a comparison against
an *independent construction*, and it agreed to ~1e-58. That agreement was real,
but both routes were evaluated at a single precision, so neither exposed this.
Agreement between two routes is not independence from a shared failure mode.

### D6 — the cross-check itself could produce a false positive

The first replacement cross-check compared the two null-vector routes after
normalising by the largest **mechanical** component. For the genuine wire-gauge
branch the mechanical content is ~1e-61, so that normalisation divided by
roundoff and amplified the two routes' noise by 1e61, reporting a disagreement
of `6.4e60`. The routes had actually agreed. A cross-check that reports failure
when it has succeeded is as dangerous as the reverse.

### D7 — the energy test was a tautology

The status test was written `|J_z| <= 1e-12 * |J_z|`, which for any non-zero
`J_z` is `1 <= 1e-12` — never true. It could only fire when `J_z` was exactly
`0.0`. As a result the pure electric wire gauges at Al6061 30° branches 3 and 4,
with `|J_z| = 2.787e-107`, were reported as `FLOWING`.

### D8 — half-space screening checked boundedness only

Admissibility was decided by `Im(ξ)` sign alone. A wire gauge is bounded, is a
genuine null vector of the operator, and carries no mechanical motion. Screening
on boundedness alone therefore admitted gauges as physical modes in both
half-spaces. Admissibility needs **both** boundedness and non-triviality.

### D9 — normal incidence had a hidden, unprincipled threshold

At normal incidence `k = 0`, so `|Im ξ|/|k|` is undefined. The fallback was:

```python
spatial = "WEAKLY_DAMPED" if abs(xi_i) < mp.mpf("1") else "EVANESCENT"
```

— an absolute threshold of `1 m⁻¹` with no physical justification, while the
comment directly above it claimed the comparison was against the thermal
wavenumber scale. It happened to give the right answer for these three materials
and would not generalise. Compounding it, the fallback set `confidence =
"MEDIUM"` which a later block then overwrote to `HIGH`, so 6 normal-incidence
rows were reported `HIGH` confidence with a blank attenuation field and no way
to tell how they had been classified.

Fixed: the fallback is `|Im xi| / |xi_th|`, the thermal wavenumber is the only
other length scale the branch set supplies, and the basis used is written into
every row as `attenuation_basis`. At `k = 0` the thermal branches land at ratio
`0.707 = 1/sqrt(2)`, which is the correct value for an axial diffusion branch in
that limit — a check on the fallback itself, not just on its bookkeeping.

---

## 3. The corrected convention (derived, not assumed)

Fields are `Re{ A exp(i(kx + ξz − ωt)) }`, i.e. `e^{−iωt}` with spatial
`e^{+iξz}`. Therefore:

- `|e^{iξz}| = e^{−Im(ξ) z}`
- **`Im(ξ) > 0` ⇒ decays toward `+z`**; `Im(ξ) < 0` ⇒ decays toward `−z`
- phase travels toward `+z` iff `Re(ξ) > 0`, but see below

This was re-derived numerically, not quoted: `verify_convention()` evaluates the
field at two z values and compares the measured decay rate against `−Im(ξ)`.
Result: **16/16 branches agree**, residual ~1e-50 at 60 digits.

### Why the sign of Re(ξ) is not a classifier

At PZT-5H 30°, branch 3 has `Re ξ = +1.36e-11` while `k = 785.4`. `Re ξ` is
1.4e-11 — effectively zero. A rule "`Re ξ > 0` ⇒ propagating toward +z" would
call this a forward-travelling wave with phase velocity 4.6e11 m/s, 1.5e6 times
*c*. The mode is actually evanescent: `|Im ξ|/|k| = 1.076`.

Separately, branch 3 has `J_z = −1.50e10` while branch 4 has `+1.50e10` — energy
flux and `Re ξ` do not have to agree in sign for a near-cutoff mode.

`sign(Re ξ)` is therefore retired as a classifier. Mutation test `M-C` restores
it and is caught by `MAC03` and `MAC07`.

---

## 4. The corrected classification

Four **independent** machine-readable fields replace the single label:

| field | values |
|---|---|
| `spatial_class` | `EVANESCENT` \| `WEAKLY_DAMPED` \| `THERMAL_DIFFUSIVE` \| `GAUGE` |
| `radiation_dir` | `DECAYS_TOWARD_PLUS_Z` \| `DECAYS_TOWARD_MINUS_Z` \| `NONE` |
| `energy_flux_status` | `FLOWING` \| `ZERO` \| `UNDEFINED` |
| `confidence` | `HIGH` \| `MEDIUM` \| `LOW` |

They are independent on purpose. A mode can be `EVANESCENT` and still carry a
non-zero flux (PZT-5H branch 3); it can be `THERMAL_DIFFUSIVE` and still have a
well-defined radiation direction (branch 0). `PROPAGATING` is retired and
`PHYSICAL_GROWING` is removed.

### Thresholds — all measured, none guessed

| constant | value | justification |
|---|---|---|
| `ATTENUATION_EVANESCENT` | 1.0 | `|Im ξ|/|k| = 1` is one neper per radian of transverse phase; reported alongside so borderline cases are visible |
| `ATTENUATION_WEAK` | 1e-3 | 0.6 % amplitude change per wavelength |
| `THERMAL_DISPERSION_TOL` | 0.1 | measured gap: thermal branches land at 1.000–1.016, acoustic at ≤6.1e-3, a factor >150 |
| `GAUGE_MECHANICAL` | 1e-40 | measured gap: true gauges are at **exactly 0** relative to scale; smallest physical mode is 4.01e-13 |

The attenuation ratio is reported together with the scale it was measured
against (`attenuation_basis`), because `|Im xi|/|k|` and `|Im xi|/|xi_th|` are
different quantities and a bare number invites the reader to compare across
them.

Two calibration errors were made and corrected while deriving these, both
recorded so they are not repeated:

1. **A unit-mixed vector norm is not a gauge detector.** `|u| ~ 1 m` against
   `|φ| ~ 1e9 V` makes every piezoelectric mode look gauge-dominated. The test
   must compare against the largest *component* (a pure normalisation scale).
2. **Amplitude composition does not identify a wave type.** Al6061 branch 2 is a
   genuine acoustic mode that carries a large temperature swing; judging by
   amplitude misclassified it. Only the **dispersion relation** separates
   thermal from acoustic branches.

---

## 5. Corrected results

PZT-5H @ 1 MHz, 30° (`k = 785.398 m⁻¹`):

| br | Re ξ | Im ξ | spatial class | radiation | flux |
|---|---|---|---|---|---|
| 0 | −2.5694e6 | −2.5694e6 | THERMAL_DIFFUSIVE | toward −z | FLOWING |
| 1 | −3497.7 | −1.05e−9 | WEAKLY_DAMPED | toward −z | FLOWING |
| 2 | −1139.0 | −7.62e−8 | WEAKLY_DAMPED | toward −z | FLOWING |
| 3 | +1.36e−11 | −845.2 | EVANESCENT | toward −z | FLOWING |
| 4 | −1.36e−11 | +845.2 | EVANESCENT | toward +z | FLOWING |
| 5–7 | conjugate roots | | WEAKLY_DAMPED ×2, THERMAL_DIFFUSIVE | toward +z | FLOWING |

Half-space screening: `z>0 → {4,5,6,7}`, `z<0 → {0,1,2,3}`, `non_decaying → []`.
Al6061: `z>0 → {5,7}`, `z<0 → {0,2}`, `non_decaying → {1,6}` (`Im ξ = 0`
exactly), `gauge_excluded → {3,4}`.

Corpus-wide (120 rows): `THERMAL_DIFFUSIVE` 30, `WEAKLY_DAMPED` 68,
`EVANESCENT` 8, `GAUGE` 14, `PROPAGATING` **0**. Confidence: 86 `HIGH`,
34 `MEDIUM`, 0 `LOW`. **78 of 120 labels changed** against Phase 1D.

The transitions are not all of one kind, and conflating them would overstate the
correction:

| transition | count | is this a physics correction? |
|---|---|---|
| `PROPAGATING → THERMAL_DIFFUSIVE` | 30 | **yes** — D1, the severe defect |
| `PROPAGATING → EVANESCENT` | 8 | **yes** — D1, genuinely evanescent |
| `DEGENERATE_EXCLUDED → GAUGE` | 14 | **yes** — D7/D8, gauges were miscoded as degenerate |
| `PROPAGATING → WEAKLY_DAMPED`, ratio ≥ 1e-6 | 12 | **yes** — ratio 0.47–0.96, wrongly called propagating |
| `DEGENERATE_EXCLUDED → WEAKLY_DAMPED` | 6 | no — already excluded from screening |
| `ELASTIC_DECOUPLED → WEAKLY_DAMPED` | 8 | no — decoupled media, correctly elastic |
| **total** | **78** | |

So of the 78 changed rows: **64 are physics corrections** (30 thermal, 8
evanescent, 14 gauge, 12 genuinely damped) and **14 are reclassifications of
rows that were already excluded from screening** in Phase 1D — a wording change,
not a change to any admissibility decision.

Separately, 48 rows moved `PROPAGATING → WEAKLY_DAMPED` with an attenuation
ratio below `1e-6` and are flagged `label_changed=False` on purpose: those modes
really were propagating, the old label was defensible, and only the vocabulary
moved. Counting them in would let this audit claim 126 corrections where 64
occurred, so the two sets are kept apart throughout.

---

## 6. Regression tests

`validation/tests/test_12_modal_classification.py` — 9 cases, all PASS.

| id | guards against |
|---|---|
| `MAC00` | Phase 1D evidence being edited (sha256 pin) |
| `MAC01` | D1 — spatial label contradicting its own attenuation ratio |
| `MAC02` | D1 — thermal branches labelled propagating (dispersion oracle) |
| `MAC03` | D3/D2 — fields not independent; `PHYSICAL_GROWING`; `sign(Re ξ)` rule |
| `MAC04` | convention silently flipping |
| `MAC05` | D5 — classification depending on working precision |
| `MAC06` | D7 — tautological energy test; gauge carrying flux |
| `MAC07` | D8 — admissibility checked on boundedness alone |
| `MAC08` | D9 — normal-incidence attenuation basis unstated |

### Mutation verification

A test that cannot fail is not evidence. Each defect was reintroduced and the
suite re-run:

| mutation | injected defect | caught by | result |
|---|---|---|---|
| M1 | `Im ξ ≠ 0 ⇒ PROPAGATING` (the Phase 1D rule) | MAC02 | **DETECTED** |
| M-A | restore `|J_z| <= 1e-12·|J_z|` | MAC06 | **DETECTED** |
| M-B | screening on boundedness only | MAC06, MAC07 | **DETECTED** |
| M-C | radiation by `sign(Re ξ)` | MAC03, MAC07 | **DETECTED** |
| M-D | remove structural-zero detection | MAC05 | **DETECTED** |

5 of 5 mutations detected.

---

## 7. Provenance

Artifacts written by this task:

| path | content |
|---|---|
| `code/scattering/modal.py` | classifier, four independent fields, screening |
| `code/scattering/nullspace.py` | precision-stable null vector, two cross-checked routes |
| `code/scattering/export_inventory.py` | inventory export under the new experiment ID |
| `code/scattering/__init__.py` | package exports |
| `validation/tests/test_12_modal_classification.py` | 9 regression cases |
| `validation/phase_02a_modal_inventory.csv` | 120 corrected rows |
| `validation/phase_02a_modal_inventory_manifest.json` | commit, precision, material set, sha256, screens |

Reproduce:

```bash
cd code && python3 -m scattering.modal              # corrected inventory table
cd code && python3 -m scattering.nullspace           # precision-stability check
cd code && python3 -m scattering.export_inventory    # regenerate the CSV
cd validation && python3 -m pytest tests/test_12_modal_classification.py -v
```

CSV sha256: `1f09b475d7133dc6a88e5b59283a75673d9f1ed36455f29c302faee27ee9499c`.

---

## 8. Consequence for later tasks

D5 is the important one for Task 2 onward. The coupled scattering system will
need null vectors and a mode-selection rule at the interface. If the underlying
null-vector route can silently return the wire gauge instead of a thermal branch
depending on the precision requested, then any scattering matrix built on it is
unreliable in a way that no residual check would reveal — the returned vector
*is* a valid null vector, just of the wrong branch. `code/scattering/nullspace.py`
is therefore the required entry point for Task 2, and it refuses to return a
vector when its two routes disagree or when the residual is not at the level the
arithmetic supports.

Two items remain open and are carried into the blockers report rather than
closed here:

- **The Phase 1D BLK-U3 closure did not test precision invariance.** It
  established agreement at one precision. That closure remains valid as stated,
  but it is weaker than it appeared, and Phase 2A does not treat it as covering
  the failure mode in D5.
- **`kernel_scaled` itself is unchanged.** It is still used elsewhere in the
  repository, including by the Phase 1D evidence. Retrofitting it would
  invalidate that evidence, so the fix lives in a new module and the old one is
  left as-is and documented.
