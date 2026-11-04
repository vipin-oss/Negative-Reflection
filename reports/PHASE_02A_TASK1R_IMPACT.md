# Phase 2A Task 1R — Impact Analysis

**Branch:** `arena/a3ede4fe-negative-reflection`
**Base commit:** `3cee0741d32f99d83f38d4076aa147aeccee54c7`
**Corrected experiment:** `PHASE-02A-MODAL-CORRECTED-2`

Per the Phase 2A rule that an impact analysis must precede changes to
foundational equations: **no foundational equation was changed by task 1R.**
This analysis covers what the defect and its repair actually touch.

---

## 1. What was changed

| file | change | nature |
|---|---|---|
| `code/scattering/export_inventory.py` | rewritten around frozen `_Branch` records | defect repair |
| `code/scattering/radiation.py` | **new** — radiation conditions | defect repair (P0-3) |
| `code/scattering/modal.py` | radiation direction, lossless class, `k=0` confidence cap | defect repair (P0-3) |
| `code/scattering/compare_report.py` | **new** — generates the row comparison | verification tooling |
| `validation/tests/test_13_inventory_integrity.py` | **new** — 7 cases | verification + mutation |
| `validation/tests/test_12_modal_classification.py` | MAC03, MAC07 updated | see §4 |
| `evidence/artifact_registry.json` | **new** | evidence governance |

Not changed: `code/foundation_checks.py`, `theory/*.tex`, the equation registry,
the material table, `code/elastic_psv.py`.

---

## 2. Claim exposure

The invalidated artifact `ART-02A-INVENTORY-V1` was published at `3cee074` and
is the only artifact carrying contaminated numbers. No claim registry entry was
bound to it: task 1 was committed without adding claims, so **no claim in
`evidence/claim_registry.csv` currently cites any Phase 2A inventory number.**

This is a narrow escape, and it is worth naming plainly: the containment came
from the claims not having been written yet, **not** from any check that would
have caught a contaminated claim. Had task 1's numbers been promoted into
claims before the supervisor's inspection, those claims would now require
withdrawal.

### Claims that remain unaffected

| claim area | why unaffected |
|---|---|
| Phase 1C/1D equation and claim registers | **extended, not altered**: 41 → 43 equations and 45 → 49 claims. The 4 new claims are `CLM-02A-01..04`; the 2 new equations are `EQ-02A-RAD-01` and `EQ-02A-ADM-01`. No existing row was modified, so every Phase 1C/1D binding is unchanged. |
| `BASELINE-01D-1` (CANDIDATE) | derived from Phase 1D artifacts, not from task 1 |
| energy-flux claims (`sec:ls-closure`, `prop:relax`, `rem:twofluxes`) | depend on `flux_z_mp`, which task 1R does not modify |

### Claims that would have been affected

The contaminated columns were `Re_xi`, `Im_xi`, `abs_xi`, `J_z`. Any claim
about **modal energy flux**, **decay direction**, or **mode admissibility**
derived from the inventory would have been wrong for 105 of 120 rows (all rows
except the last of each block, which were correct by accident). Specifically the
inverted sign on branch 0 in every block would have asserted that the strongly
decaying thermal branch propagates toward `+z`.

---

## 3. Dependency impact

| dependent | impact |
|---|---|
| `screen_halfspaces` | signature extended with optional `jzs`; old call sites still work |
| `classify_mode` | new keys `radiation_source`, new enum values; existing keys unchanged |
| `SPATIAL_CLASSES` | gained `PROPAGATING_NON_DECAYING`; `PROPAGATING` remains retired |
| `RADIATION_DIRS` | gained `TRAVELS_TOWARD_PLUS_Z` / `TRAVELS_TOWARD_MINUS_Z` |
| any consumer switching on `radiation_dir` | **must** handle the `TRAVELS_*` values; a consumer matching only `DECAYS_*` will silently drop the 14 lossless branches |
| Task 2 (not started) | unaffected — no Task 2 code exists |

The `RADIATION_DIRS` change is the one real forward hazard. A scattering
assembly that selects outgoing modes by testing
`radiation_dir == "DECAYS_TOWARD_PLUS_Z"` would now exclude the lossless
outgoing branches, losing Al6061 br6 and its counterparts. This is flagged in
the blockers report as a Task 2 requirement: **select on `outgoing_status`, not
on the decay enum.**

---

## 4. Test changes, and why they are not weakening

Two Task 1 tests were updated. Both encoded a model that task 1R supersedes,
and in both cases the corrected test is **stricter**, not looser.

### MAC03 — radiation direction

Was: `radiation_dir` must equal the `sign(Im ξ)` rule, else `NONE`.
Now: non-zero `Im ξ` → decay direction; zero `Im ξ` → direction must follow
`sign(⟨J_z⟩)`, and the label must be `TRAVELS_TOWARD_*`.

The old rule asserted that a lossless mode has no direction — the P0-3 defect.
The new rule additionally requires the flux-derived direction to be *consistent
with the flux*, which the old test never checked. Reintroducing `Im ξ = 0 ⇒
NONE` is caught (**verified by mutation**).

### MAC07 — admissibility

Was: each branch must land in exactly one of `z>0` / `z<0` / `non_decaying` /
`gauge_excluded`.
Now: gauges appear only in `gauge_excluded`; decaying modes in exactly one
half-space; lossless modes admissible in **both** halves **and** listed as
non-decaying **and** carrying a flux-derived direction.

The "exactly one category" rule was satisfied by the old code precisely because
lossless modes were dumped into `non_decaying` and nothing else. The new test
requires three things of such a mode where the old required one.

### Verification that the mutation protection still fires

| mutation | caught |
|---|---|
| tautological energy test | MAC06 **DETECTED** |
| no gauge exclusion in screening | MAC06 **DETECTED** |
| radiation by `sign(Re ξ)` | MAC03 **DETECTED** |
| `Im ξ = 0` ⇒ `NONE` | MAC03/MAC07 **DETECTED** |

No test was removed. No tolerance was widened. No case was withdrawn.

---

## 5. Numerical changes vs the invalidated artifact

| quantity | invalidated v1 | corrected v2 | why |
|---|---|---|---|
| `WEAKLY_DAMPED` | 68 | 54 | 14 lossless modes reclassified |
| `PROPAGATING_NON_DECAYING` | — | 14 | new class for `Im ξ = 0` |
| `EVANESCENT` | 8 | 8 | unchanged |
| `THERMAL_DIFFUSIVE` | 30 | 30 | unchanged |
| `GAUGE` | 14 | 14 | unchanged |
| `HIGH` confidence | 68 | 68 | unchanged |
| `MEDIUM` confidence | 34 | 46 | 12 more from the global `k=0` cap |
| `LOW` confidence | 0 | 6 | gauges at `k=0` with no resolvable flux |
| labels changed vs Phase 1D | 78 | 78 | unchanged |

The class counts by category are stable; what changed is the *identity* of the
rows — v1 assigned 120 correct-looking labels to 120 rows of which 105 carried
the wrong root. **The counts agreeing between v1 and v2 is not evidence that v1
was right**; only the row-level comparison establishes that.

---

## 6. Blockers affected

None closed. None opened for the defect itself, because it is repaired and
verified. Four items remain open and are carried to
`reports/PHASE_02A_OPEN_BLOCKERS.md`:

1. `kernel_scaled` precision instability (pre-existing, unmodified by 1R).
2. Phase 1D BLK-U3 closure did not test precision invariance.
3. `jz_scale` relative test is unreachable dead code.
4. `PROPAGATING_NON_DECAYING` roles not yet confirmed against an independent
   elastic solver.

---

## 7. Assessment

The defect was real, it was in published evidence, and no test would have caught
it. The repair is structural rather than a patch to the two offending lines: the
data model now makes cross-branch reads impossible, and the new verifier
recomputes every number from an independent route instead of checking that the
file agrees with itself.

Running the verifier against the published artifact produces **684 defects**,
which is the strongest available statement that the verifier is not vacuous.

What this does **not** establish: that the physics of the corrected inventory is
right. It establishes that every exported number corresponds to the correct
underlying computation. That is the stated objective of task 1R, and it is a
necessary but not a sufficient condition for the science.
