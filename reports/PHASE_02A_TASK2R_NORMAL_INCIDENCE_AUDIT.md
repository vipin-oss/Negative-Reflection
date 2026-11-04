# PHASE 2A — TASK 2R: NORMAL-INCIDENCE DEGENERACY AUDIT

**Status: COMPLETE — the Task 2 normal-incidence reduction was wrong, and it was
wrong in exactly the way the supervisor suspected.**

Branch `arena/a3ede4fe-negative-reflection`, based on `ca2254f`.

---

## 1. Headline result

The Task 2 statement

> "At normal incidence the sagittal geometry sets `u_x = σ_xz = 0` for every mode"

is **FALSE**, and it was load bearing. At normal incidence in an isotropic
elastic solid the shear (SV) mode carries

```
u_x = 1,  u_z = 0,  φ = 0,  ϑ = 0
σ_xz = i ξ C44 u_x = 5.3147726e13 j        (≠ 0)
```

The `u_x` and `txz` interface rows are therefore the **only** rows that
distinguish a shear mode from a longitudinal one. Task 2 removed them, which
deleted the physics that the problem is about.

The defect was in `assembly.py`'s trivial-row reduction: a row was dropped when
it was zero **across the admitted basis**, without requiring that it also be
zero on the incident right-hand side. `rhs_is_zero` was computed and recorded
but never used in the decision. Two lines of the return unit did not consult it.

That single omission produced upstream consequences that were then misdiagnosed:

* The aluminium normal-incidence system came out `rows < unknowns`.
* The solver reported `cond = 1.10e21` and the rescaling test returned
  `DEFERRED_TO_CONDITIONING`.
* This was recorded as **BLK-02A-14: "Al|Al normal incidence under-determined"**
  and treated as a candidate physical non-uniqueness.

**It was never non-uniqueness. It was a missing equation.**

---

## 2. What the exact normal-incidence states are

Measured at Al6061, `ω = 2π·1e6`, `k = 0` exactly, at 40 decimal digits.
All eight roots were verified against the characteristic matrix; every one
admits a null vector at the arithmetic floor.

| branch | role | ξ | u_x | u_z | φ | ϑ | t_xz | t_zz | D_z | q_z |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | OUTGOING | −216714.95 − 216712.02j | 0 | 1 | 0 | 6.26e9 | 0 | 3.48e11 | 0 | 3.21e17 |
| 1 | OUTGOING | −2005.5746 | **1** | 0 | 0 | 0 | **5.31e13** | 0 | 0 | 0 |
| 2 | OUTGOING | −979.77 − 1.6e−4j | 0 | 1 | 0 | 6.40e5 | 0 | 1.09e14 | 0 | 1.05e11 |
| 3 | GAUGE_EXCLUDED | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| 4 | GAUGE_EXCLUDED | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| 5 | INCOMING | 979.77 + 1.6e−4j | 0 | 1 | 0 | 6.40e5 | 0 | 1.09e14 | 0 | 1.05e11 |
| 6 | INCOMING | 2005.5746 | **1** | 0 | 0 | 0 | **5.31e13** | 0 | 0 | 0 |
| 7 | EVANESCENT | 216714.95 + 216712.02j | 0 | 1 | 0 | 6.26e9 | 0 | 3.48e11 | 0 | 3.21e17 |

Branches 1 and 6 are the shear modes. Their whole state is `u_x` and `σ_xz`.

### Which rows are genuinely zero

| component | vanishes for every admitted mode? | verdict |
|---|---|---|
| `u_x` | **NO** — 1.0 on br1, br6 | must be kept |
| `u_z` | no — 1.0 on br0, br2, br5, br7 | kept |
| `φ` | **yes** — 0 for all (Al6061 is non-piezoelectric) | removable **only because the RHS is also 0** |
| `ϑ` | no — 6.3e9 etc. | kept |
| `t_xz` | **NO** — 5.31e13 on br1, br6 | must be kept |
| `t_zz` | no — 1.09e14 | kept |
| `D_z` | **yes** — 0 for all (non-piezoelectric) | removable with zero RHS |
| `q_z` | no — 3.2e17 | kept |

Note the trap in the two removable rows: `φ` and `D_z` vanish because
**Al6061 has `e31 = e33 = e15 = 0`**, not because of normal incidence. The same
reduction run on PZT-5H keeps the `φ` row, and that is the correct behaviour.
Anything that attributes the `φ`/`D_z` removal to "sagittal decoupling" has the
causal story backwards.

---

## 3. The root-error hypothesis was wrong, and the way it was wrong matters

The investigation began from a measurement reported as:

```
xi_roots  = 2005.5745662580894209
exact     = 2005.5745662580893715       |Δ| = 4.93e-14
det M(xi_roots) = 0.009219 − 279.759j   (not zero)
```

This motivated a Newton polish against `det M`. **Both the measurement and the
polish are wrong.**

The probe built its "exact" reference as `cT = sqrt(m.C44 / m.rho)` with both
operands **float64**. That ratio is itself only accurate to `1.2e-13` as a
wavenumber. Measured properly, at the working precision:

| dps | min \|ξ_roots − w/cT\| |
|---|---|
| 30 | 2.0e-28 |
| 40 | 0.0 |
| 50 | 2.7e-48 |
| 60 | 0.0 |

**`xi_roots` is exact to working precision. There is no root error.**

And the polish built to repair that non-existent fault is not neutral — it is
destructive. From the exact root it moves ξ by `+2.0e-17j` and takes

```
M[0,0]:   1.65e-24   →   −2.13e-3
```

turning a structurally-zero row into a loud one and destroying the nullspace it
was supposed to protect. **The polish has been removed from `roots.py`** and the
reason is recorded in the module docstring so it is not reintroduced.

This is the third time in this phase that a false alarm came from the probe
rather than the code. It is recorded here rather than quietly corrected because
the "fix" would have introduced a real defect in place of an imaginary one.

---

## 4. The two defects that were real

### 4.1 The null-vector metric normalised by the global matrix max

`nullspace._rel_residual` computes

```
max|A v| / (max|A| · max|v|)
```

`M`'s rows are equations in different units and span `1e-5` to `1e19`. A
candidate violating the electric-potential row by `3.6e-5` scores
`3.6e-5 / 1.99e19 = 1.8e-24` and **passes**. That is how the pure gauge vector
`(0, 0, 1, 0)` — no mechanical content at all — was returned as "the mode" on a
branch whose physical content is `(1, 0, 0, 0)`.

### 4.2 Row equilibration amplified a roundoff row to unit scale

`_two_sided_equilibrate` rescales rows. At `k = 0` on the shear branch, row 0 of
`M` evaluates to `1.65e-24` — that is **roundoff**, the residue of subtracting
two `1.07e17` numbers. Equilibration lifts it to unit magnitude, after which it
dominates the eigenvector selection and excludes the true shear direction.

Measured: row-0 entry equals the sqrt-then-square roundoff and tracks `10^-dps`
exactly — `2.8e-14` at dps 30, `1.7e-24` at 40, `1.1e-44` at 60.

### 4.3 Confirmed independently of the root

At the **exact** root, the old path still returns the wrong vector:

```
existing null_vector  →  (0, 1, 0, −2.16e6)
gram_residual         →  3.004e-7
row_relative_residual →  4.954e-7
true answer           →  (1, 0, 0, 0)
```

Both metrics agree it is not a null vector. **The defect is in the algorithm,
not in the root** — which is why removing the polish alone was not enough.

---

## 5. The repair

### `code/scattering/roots.py` (new)

* `row_relative_residual` — normalises **each row by itself**, and skips rows
  whose largest entry is at the working-precision floor. That last clause is
  essential: without it, the shear root's row 0 (value `1.65e-24`, equal to its
  own row max) is scored as `1.0` and **rejects the correct physical answer**.
* `null_vector_verified` — finds the null direction by **enumeration** over
  (dropped row, free column), scoring every candidate with the row-relative
  metric. Eigenvector methods cannot do this: whichever way `A` is scaled, the
  smallest singular direction of `AᴴA` is dominated by the loudest row.
  Returns `ok=False` rather than a compressed direction.
* `matrix_rank_nullity` — on a row-equilibrated matrix, with the tolerance
  reported alongside the rank. `mp.svd`/`mp.cond` fail on complex input;
  `mp.eighe` on `AᴴA` is used instead.
* `verified_roots` — every root checked; a root with no demonstrated null vector
  is returned with `ok=False` and must not be used to build a mode.

Measured floors: every branch of both media verifies at `≤ 1.8e-24`.
`NULL_TOL = 1e-20` sits four orders above the floor.

### `code/scattering/assembly.py`

* `Side._build` now calls `null_vector_verified` and **honours `ok`**. The old
  code called `nullspace.null_vector` and discarded the `agreement` flag it
  returns — a guard no caller honours is not a guard.
* Trivial-row reduction now requires **both** conditions:

```python
if row_zero and rhs_zero:      # removable
elif row_zero and not rhs_zero: # kept, and flagged INCONSISTENT
```

* New `rows_inconsistent` field surfaces contradictory systems instead of
  silently dropping the equation.

### `code/scattering/modal.py`

Routed through the same verified path, so its diagnostic main agrees with the
assembly rather than diverging from it.

---

## 6. Effect on the identical-medium problem (BLK-02A-14)

Al6061 | Al6061 at `k = 0`, both polarizations, after the repair:

| incidence | n_eq | n_unk | rank | square | max\|reflected field\| | continuity error |
|---|---|---|---|---|---|---|
| P (br5) | 6 | 6 | 6 | yes | 1.1e-37 | 1.4e-41 |
| SV (br6) | 6 | 6 | 6 | yes | 1.0e-41 | 0.0 |

* Reflected coefficients zero in every channel: `0, −5.2e-42, 0`.
* Transmitted coefficient **exactly 1.0** on the shear channel.
* All eight physical components satisfy `incident + reflected = transmitted`.
* `cond = 1.10e21` no longer appears; the system is square and full rank.

**BLK-02A-14 is CLOSED as an artifact.** The identical-medium problem was never
under-determined — see the companion rank-and-gauge report for why
`rows < unknowns` was never evidence of non-uniqueness.

The independent check: for matched impedances the closed-form Fresnel amplitude
is exactly zero, and `elastic_reference.fresnel_normal_amplitudes` gives
`r = 0, t = 1` exactly.

---

## 7. Regression tests

`validation/tests/test_15_scattering_task2r.py` — 7 tests, all PASS.

| id | pins |
|---|---|
| SCTR01 | a shear boundary row is never removed when an admitted mode is nonzero on it |
| SCTR02 | a P-wave zero component is not generalised to a modal zero |
| SCTR03 | trivial-row removal does not change the rank |
| SCTR04 | no duplicate or proportional modal columns |
| SCTR05 | gauge freedom is not reported as physical non-uniqueness |
| SCTR06 | identical media reflect nothing |
| SCTR07 | reduction is verified on the physical basis and the incident data |

Negative-control mutations MR01–MR06 prove each detects its defect. See the
companion test-results JSON.

---

## 8. What remains open

* **BLK-02A-16 (new)** — `nullspace.null_vector` and `_rel_residual` remain in
  the tree and are now unused by `assembly`/`modal`. They have not been deleted
  because other code still imports them. Their metric is unsound; marking them
  deprecated is safer than a silent removal.
* **BLK-02A-17 (new)** — the `φ`/`D_z` removal rests on `e31 = e33 = e15 = 0`
  for the material in question, not on normal incidence. A configuration with a
  piezoelectric *and* a non-piezoelectric half-space now behaves correctly
  (phi row kept), but the reduction reason string still says "sagittal
  decoupling", which is the wrong cause. Wording, not mathematics.
* **G4 / G5 remain OPEN.** Nothing here implements or claims `R + T + A = 1`.
  Incident–reflected cross terms are preserved. Absorption is not computed as
  `1 − R − T`.
* **BLK-02A-01/-02/-03/-05/-13/-15** unchanged.

## 9. Validity of the Task 2 numbers

Task 2's normal-incidence scattering amplitudes are **INVALID / SUPERSEDED**.
They were computed from a system with the shear rows deleted. The oblique
results (78 of 80 configurations) never exercised the defect — the trivial-row
reduction only fired at exactly `k = 0` — and are unaffected; they are
**RETAINED**. Regeneration is carried out under experiment ID
`PHASE_02A_TASK2R` in the companion CSV.
