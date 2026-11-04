# Phase 2A Task 1R — Modal Inventory Data-Integrity Repair

**Branch:** `arena/a3ede4fe-negative-reflection`
**Base commit:** `3cee0741d32f99d83f38d4076aa147aeccee54c7`
**Supersedes experiment:** `PHASE-02A-MODAL-CORRECTED-1` (INVALID)
**Corrected experiment:** `PHASE-02A-MODAL-CORRECTED-2`
**Scope:** Task 1R only. Task 2 NOT started.

---

## 1. What was wrong

The supervisor found that the Task 1 exporter published incorrect per-branch
values. This was confirmed before any repair was attempted.

### P0-1 — cross-branch contamination

In `code/scattering/export_inventory.py`, `build_rows()` computed each branch
into loop-local variables `xi`, `jz`:

```python
for br, xi in enumerate(roots):          # xi, jz bound here
    ...
    jz = fc.flux_z_mp(m, k, xi, w, a)
    modes.append(md)                     # only the CLASSIFICATION was kept
...
for br, md in enumerate(modes):          # second loop, after the first ended
    "Re_xi": _f(mp.re(xi)),              # <-- xi is the LAST branch's root
    "J_z": _f(jz),                       # <-- jz is the LAST branch's flux
```

`xi` and `jz` survived the first loop holding the final iteration's values, so
**every row in a material/angle block carried branch 7's root and flux.**

Confirmed in the published file — PZT-5H at 30°, `validation/phase_02a_modal_inventory.csv`:

| branch | exported Re ξ | true Re ξ | exported J_z |
|---|---|---|---|
| 0 | +2.5694e6 | **−2.5694e6** | +1.957e27 |
| 1 | +2.5694e6 | −3497.7 | +1.957e27 |
| … | +2.5694e6 | … | +1.957e27 |
| 7 | +2.5694e6 | +2.5694e6 | +1.957e27 |

Branch 0's sign was inverted. This is not a rounding error and not a
presentational issue: the published artifact asserted that a mode decaying
toward −z has a positive real wavenumber.

**Why no test caught it.** Task 1 had no test that compared an exported value
against the calculation that produced it. The artifact was internally consistent
— every row agreed with every other row — which is exactly what a contaminated
block looks like.

### P0-3 — incomplete radiation condition

Separately, `screen_halfspaces` assigned admissibility from `sign(Im ξ)` alone
and filed every branch with `Im ξ = 0` under `non_decaying`. Al6061 branches 1
and 6 have `Im ξ` exactly zero and carry `⟨J_z⟩ = ∓1.815e20 W/m²` — genuine
lossless propagating waves, one incoming and one outgoing. Filing both under
`non_decaying` discarded the only information that distinguishes them, and the
classifier called them `WEAKLY_DAMPED` despite a damping of exactly zero.

---

## 2. The repair

### 2.1 Structured branch records (P0-1)

`_Branch` is a frozen dataclass holding, for ONE branch:

`branch` index · complex `xi` · `eigenvector` · `jz` · `classification` ·
`nullspace` diagnostics · `admissibility` · `outgoing_status`

`build_branch_records()` returns a list of these, index-aligned with the root
bank. The half-space screen needs the whole set, so it is applied in a second
pass that rebuilds each record with `dataclasses.replace` — the record already
owns its own root and flux, so there is no window in which a value can be read
from a neighbour.

`build_rows()` reads fields from the record by name. Nothing is recovered by
parsing the `reasons` prose; that text is generated *from* the record, and
`INV01` asserts the two cannot contradict.

### 2.2 Radiation condition (P0-3)

New module `code/scattering/radiation.py`. The test, in order:

1. **GAUGE** → `GAUGE_EXCLUDED`. Bounded, but no mechanical motion, so not a
   scattering channel. Boundedness is **not** admissibility.
2. **`Im ξ ≠ 0`** → the mode is localised; it is admissible in the half-space it
   decays into, and `OUTGOING` there.
3. **`Im ξ = 0`** → the mode does not decay anywhere and is bounded on **both**
   sides. The direction is carried by the axial energy flux. For an isotropic
   elastic medium `⟨J_z⟩ = ±cE` with `c = √(μ/ρ) > 0`, so `sign(⟨J_z⟩)` is the
   axial phase/energy-velocity direction. Admissible on both sides; `OUTGOING`
   in the side the flux flows into, `INCOMING` in the other.
4. **`Im ξ = 0` with unresolvable flux** → `UNRESOLVED`, reported as such. A
   bounded mode is not thereby an outgoing mode.

The limiting-absorption continuation justifies step 3: continuing `Im ω → 0⁺`
gives a flux direction that agrees with `sign(⟨J_z⟩)` for a lossless mode, and
removes the tie for a nearly-lossless one. Where the continuation gives no
direction (flux zero at every order) the code refuses to invent one.

`Im ξ = 0` is therefore **not** treated as inadmissible, and a gauge is **not**
treated as outgoing. Both are now covered by tests.

### 2.3 Normal incidence and degeneracy

At `k = 0` the physical scale `|Im ξ|/|k|` does not exist; the ratio falls back
to `|Im ξ|/|ξ_th|` and the basis is written into every row. Two further issues
found here:

- The thermal tolerance was calibrated at oblique incidence, where thermal
  branches land at 1.000–1.016. At `k = 0` the ratio genuinely degrades to
  `1/√2 = 0.707`, i.e. the tolerance is being asked to work outside its measured
  range. **All** `k = 0` classifications are now capped at `MEDIUM` confidence.
  The cap is applied once, globally — an earlier attempt capped it inside three
  of the five classification arms and six rows escaped reporting `HIGH`.
- The two electric wire gauges sit at `ξ = 0` exactly at normal incidence, so
  the dispersion polynomial has a **genuine double root**. This is real
  degeneracy, not contamination, and the integrity test distinguishes the two by
  consulting the independent root bank rather than by counting duplicates.

---

## 3. Independent verification (P0-2)

`validation/tests/test_13_inventory_integrity.py` — 7 cases.

The oracle does **not** re-read the CSV. For every row it recomputes:

| check | independent oracle |
|---|---|
| 1. branch index ↔ root | index against the independently rebuilt root bank |
| 2. `Re_xi`, `Im_xi` | (a) evaluated directly in the dispersion polynomial — proves it is a root at all; (b) compared with bank entry `[branch]` — proves it is the *right* root |
| 3. `abs_xi` | recomputed from the exported components |
| 4. `J_z` | flux recomputed for that row's **own** eigenvector; contamination moves it by orders of magnitude |
| 5. spatial class | re-derived from the root |
| 6. radiation dir | `sign(Im ξ)` when non-zero, else from the flux |
| 7. admissibility | full two-condition screen |
| 8. outgoing channels | gauge and zero-flux modes rejected; **and** `outgoing_status` recomputed per side |
| 9. confidence | checked against the diagnostics that are supposed to support it |
| 10. prose | `reasons` must not contradict the numeric fields |

All materials, all five angles, all eight branches: **120 rows, 0 defects**.

`INV03` adds a physical identity that is independent of the flux code: for a
lossless mode `sign(⟨J_z⟩)` must equal `sign(Re ξ)` with `c = √(μ/ρ) > 0`.

### Proof the verifier catches the real defect

Run against the **published Task 1 artifact**, the same validator reports
**684 defects**, including:

```
PZT-5H@0.0 br0: Re_xi 2.5693678583e+06 != root[0] -2569367.9
PZT-5H@0.0 br0: J_z 1.9568958497e+27 != recomputed -1.9568958e+27
PZT-5H@0.0 br0: J_z_sign + != -
```

### Row-by-row comparison report

`reports/PHASE_02A_TASK1R_ROW_COMPARISON.md` — all 120 rows, generated by
`python3 -m scattering.compare_report`, showing MATCH/FAIL for every field
alongside what the invalidated file contained for the same row.

---

## 4. Mutation protection (P0-4)

Every mutation is addressed **positionally**. An earlier version looked rows up
again by the `branch` value it had just modified, so the index-shift mutation
matched nothing and silently became a no-op — a mutation that "passes" while
testing nothing.

| # | injected defect | detected | problems raised |
|---|---|---|---|
| 1 | last-branch root repeated across rows | **YES** | 300 |
| 2 | last-branch flux repeated across rows | **YES** | 168 |
| 3 | branch-index shift | **YES** | 689 |
| 4 | opposite radiation direction | **YES** | 114 |
| 5 | missing outgoing propagating mode | **YES** | 2 |
| 6 | gauge accepted as outgoing | **YES** | 42 |
| 7 | inconsistent `reasons` vs numeric fields | **YES** | 14 |
| 8 | broken manifest checksum (`INV05`) | **YES** | detected |

**8/8 detected.**

Mutation 5 initially escaped. That was a **real validator hole**, not a bad
mutation: the checks rejected a gauge being called `OUTGOING` but never noticed
the opposite failure — an outgoing channel being declared non-transporting,
which silently removes a physical mode from the scattering problem. Check 8b now
recomputes `outgoing_status` per side from the root and flux and compares it
field by field.

Task 1's own protections were re-verified after the semantics changed, so they
were not silently weakened:

| mutation | caught by |
|---|---|
| tautological energy test | MAC06 **DETECTED** |
| no gauge exclusion in screening | MAC06 **DETECTED** |
| radiation by `sign(Re ξ)` | MAC03 **DETECTED** |
| `Im ξ = 0` ⇒ radiation `NONE` (P0-3 regression) | MAC03/MAC07 **DETECTED** |

---

## 5. Evidence preservation

`validation/phase_02a_modal_inventory.csv` is **not overwritten and not
deleted**. It stays on disk as the record of what was published and what was
wrong with it.

New `evidence/artifact_registry.json`:

| artifact | status | sha256 (16) |
|---|---|---|
| `phase_01d_independent_comparison.csv` | HISTORICAL | `76c7174c08923ded` |
| `phase_02a_modal_inventory.csv` | **INVALID** | `1f09b475d7133dc6` |
| `phase_02a_modal_inventory_v2.csv` | ACTIVE | `a3b2965319a49fb2` |

`INV06` enforces this: the invalid artifact must be marked `INVALID`, must still
exist, must be byte-identical, must name its replacement, and the Phase 1D
evidence must be unmodified.

---

## 6. Corrected results

`PHASE-02A-MODAL-CORRECTED-2` — 120 rows, sha256 `a3b2965319a49fb21bd0d045840cd8810ee2562b15100767a1fef7d4cefd5f46`.

| `spatial_class` | count |
|---|---|
| `WEAKLY_DAMPED` | 54 |
| `THERMAL_DIFFUSIVE` | 30 |
| `PROPAGATING_NON_DECAYING` | 14 |
| `GAUGE` | 14 |
| `EVANESCENT` | 8 |
| `PROPAGATING` | 0 (retired) |

Confidence: 68 `HIGH`, 46 `MEDIUM`, 6 `LOW`. Labels changed vs Phase 1D: 78.

The 14 `PROPAGATING_NON_DECAYING` branches are the lossless modes that were
previously buried in `non_decaying` with no direction. Each now carries a
role — 7 incoming, 7 outgoing — which is information the Task 1 artifact did not
contain at all.

---

## 7. Validation executed

| level | cases | result |
|---|---|---|
| FAST | 72 | 71 PASS / 0 FAIL / 1 UNRESOLVED |
| STANDARD | 72 | 71 PASS / 0 FAIL / 1 UNRESOLVED |
| FULL | 72 | 71 PASS / 0 FAIL / 1 UNRESOLVED |

Registry state after this task: **43 equations** (was 41) and **49 claims**
(was 45). The four new claims are `CLM-02A-01..04`; three are `VERIFIED_STRONG`
and one (`CLM-02A-03`, lossless-mode roles) is deliberately `VERIFIED_WEAK`
because it has not been confirmed against an independent elastic solver.

`SCT02` is the single UNRESOLVED case. It is pre-existing, deliberate, and
listed in `validation/known_open.json` under blocker `BLK-G5`: there is no
coupled scattering solver, so R/T/A cannot be evaluated. It reports UNRESOLVED
rather than passing, and it is **not** a regression from this task.

---

## 8. Unresolved issues

1. **`kernel_scaled` is still the production null-vector route.** Phase 2A Task 1
   showed it can return the wrong branch depending on working precision (Al6061
   br3 at dps 90/110/120/135). Task 1R did not change it, because the Phase 1D
   evidence depends on it. The corrected inventory uses
   `scattering.nullspace.null_vector` instead. This remains open.
2. **The Phase 1D BLK-U3 closure did not test precision invariance.** It
   established agreement at one precision. Still open as stated in the Task 1
   audit.
3. **`jz_scale` relative test is effectively unreachable.** The absolute
   `1e-40` test fires first for every real case, so the relative arm never
   executes. It is dead code carrying a live-looking name; it should either be
   exercised or removed.
4. **The `PROPAGATING_NON_DECAYING` classification is new and untested against
   an independent solver.** It rests on `Im ξ = 0` exactly plus a resolvable
   flux. The elastic reference solver (`code/elastic_psv.py`) has not yet been
   used to confirm the lossless-branch roles.

**Task 2 is not started and awaits supervisor approval of this inventory.**
