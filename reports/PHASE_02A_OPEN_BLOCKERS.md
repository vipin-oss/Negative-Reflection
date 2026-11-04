# Phase 2A — Open Blockers

**As of commit:** `3cee0741d32f99d83f38d4076aa147aeccee54c7` + task 1R working tree
**Scope:** Phase 2A. Task 1R complete; Task 2 not started.

Nothing here is closed by a small residual, and nothing is closed because a test
turned green. Each entry states what evidence would close it.

---

## Pre-existing, inherited from Phase 1D

### BLK-G4 — energy accounting (P0, OPEN)

Gate G4 is not passed. The least-squares energy identity and the
cycle-averaged stored-energy derivative do not yet agree to the required
standard across all interface cases.

*Would close with:* a derivation that fixes the functional **before** the
numbers are examined, then a test that the same functional satisfies the
identity for independent mode constructions. Choosing the functional that makes
the residual small does not close this.

### BLK-G5 — coupled scattering energy balance (P0, OPEN)

No coupled piezothermoelastic R/T/A exists, so R + T + A = 1 cannot be
evaluated. `SCT02` reports `UNRESOLVED` for exactly this reason.

*Would close with:* Task 2–5 complete and an independent energy audit.

### BLK-G0 / BLK-G1 — foundation gates (P0, OPEN)

Unchanged from Phase 1D. Not touched by task 1R.

### BLK-U5 — (P1, OPEN)

Unchanged from Phase 1D. Not touched by task 1R.

### BLK-LS-STORAGE — (P1, OPEN)

Unchanged from Phase 1D. Not touched by task 1R.

### BLK-U3 — **CLOSED in Phase 1D, but the closure is weaker than it appeared**

Phase 1D closed BLK-U3 by showing the production null vector agreed with an
independent construction to ~1e-58. That agreement was real. But both routes
were evaluated at a **single working precision**, so neither exposed the
precision-dependent branch selection found in Phase 2A:

```
Al6061, 30 deg, branch 3
  dps  60..85  -> correct thermal branch
  dps  90      -> WRONG: wire gauge
  dps 110,120  -> WRONG
  dps 135      -> WRONG
```

The closure stands as written. It is not treated as covering precision
invariance, and it should not be cited as if it does.

---

## Opened or surfaced in Phase 2A

### BLK-02A-01 — `kernel_scaled` is precision-dependent (P1, OPEN)

`code/foundation_checks.py::kernel_scaled` guards column equilibration with
`if mcol > 0`. For a centrosymmetric medium the electric-potential column is
**exactly zero**, so whether the guard fires depends on the root's last digits
and therefore on working precision. When it does not fire, a spurious 3×3 minor
wins the free-variable search and the routine returns a **different null
vector** — the pure wire gauge instead of a thermal branch — with no error.

*Repaired downstream, not in place.* `code/scattering/nullspace.py` provides a
precision-stable replacement and the corrected inventory uses it exclusively.
`kernel_scaled` itself is unchanged because the Phase 1D evidence depends on it.

*Would close with:* either retrofitting `kernel_scaled` and regenerating the
Phase 1D evidence under a new experiment ID, or formally deprecating it and
migrating every call site. Neither is done.

### BLK-02A-02 — gauge classification rests on a measured gap (P2, OPEN)

`GAUGE_MECHANICAL = 1e-40` is set inside a measured gap: true gauges are exactly
`0.0` relative to scale, and the smallest physical mode observed is `4.01e-13`.
This is a **measured** separation over 4 materials × 5 angles, not a proof.

*A risk:* a material with weaker piezoelectric coupling could place a physical
mode below `4.01e-13`, where the test would still call it physical. The
threshold is far enough below that gap (`1e-40`) that this is comfortable but
unproven.

*Would close with:* an analytic bound on the minimum mechanical content of a
physical piezoelectric branch for the material class.

### BLK-02A-03 — `jz_scale` relative test is unreachable (P3, OPEN)

In `classify_mode`, the absolute test `|jz| <= 1e-40` fires for every real case,
so the relative arm
`|jz| <= 1e-20 * jz_scale` never executes. It is dead code carrying a
live-looking parameter. Either exercise it or remove it; leaving an untested
branch in a classifier invites a future false PASS.

### BLK-02A-04 — lossless branch roles unconfirmed independently (P1, OPEN)

The new `PROPAGATING_NON_DECAYING` classification identifies 14 branches
(`Im ξ = 0` exactly with resolvable flux) and assigns each a role from
`sign(⟨J_z⟩)`, justified by `⟨J_z⟩ = ±cE` with `c = √(μ/ρ) > 0`.

The identity has been checked (`INV03`: `sign(⟨J_z⟩) == sign(Re ξ)`, `c > 0`),
but the roles themselves are **not** yet confirmed against the independent
elastic solver `code/elastic_psv.py`. Since `code/elastic_psv.py` was written
independently of the coupled path, it can settle this cheaply.

*Would close with:* a test that `elastic_psv` reproduces the same
incoming/outgoing assignment for Al6061 br1/br6 at each angle. **This is a Task
3 deliverable and is required before any coupled result uses these branches.**

### BLK-02A-05 — `radiation_dir` enum grew; Task 2 must not select on it (P1, OPEN)

`RADIATION_DIRS` gained `TRAVELS_TOWARD_PLUS_Z` / `TRAVELS_TOWARD_MINUS_Z`. A
scattering assembly that selects outgoing modes by matching
`DECAYS_TOWARD_*` will silently drop the 14 lossless branches, losing the
outgoing Al6061 br6 and its counterparts.

*Requirement for Task 2:* select on `outgoing_status` from
`scattering.radiation`, never on the decay enum.

---

## Closed in Phase 2A Task 1R

### BLK-02A-10 — CLOSED: cross-branch value contamination in the exported inventory

Exported `Re_xi`, `Im_xi`, `abs_xi`, `J_z` were read from a loop variable after
iteration, so every row in a block carried the last branch's values (sign
inverted on branch 0). Repaired by generating rows from frozen per-branch
records; `ART-02A-INVENTORY-V1` marked INVALID and preserved.

*Closed by:* `INV01` (120/120 rows verified against independent recomputation),
`INV02`, and 7/7 mutation detection. Verifier produces 684 defects on the
invalid artifact, proving it is not vacuous.

### BLK-02A-11 — CLOSED: lossless modes discarded by the radiation screen

`Im ξ = 0` was treated as having no direction and filed under `non_decaying`.
Now admissibility and outgoing/incoming are separate questions; lossless modes
are admissible on both sides and carry a flux-derived role.

*Closed by:* `radiation.py` + `MAC07` + the P0-3 mutation test.

### BLK-02A-12 — CLOSED: `k = 0` classifications could report HIGH confidence

The `k = 0` fallback rests on the thermal scale, and the thermal tolerance was
calibrated at oblique incidence only. All `k = 0` classifications are now capped
at `MEDIUM`, applied globally.

*Closed by:* `MAC08` + the global cap verified across all three materials.

---

## Summary

| id | severity | status |
|---|---|---|
| BLK-G4 | P0 | OPEN |
| BLK-G5 | P0 | OPEN |
| BLK-G0, BLK-G1 | P0 | OPEN |
| BLK-U5, BLK-LS-STORAGE | P1 | OPEN |
| BLK-02A-01 kernel_scaled precision | P1 | OPEN (worked around) |
| BLK-02A-02 gauge threshold is empirical | P2 | OPEN |
| BLK-02A-03 unreachable `jz_scale` arm | P3 | OPEN |
| BLK-02A-04 lossless roles unconfirmed | P1 | OPEN (Task 3) |
| BLK-02A-05 `radiation_dir` enum hazard | P1 | OPEN (Task 2 requirement) |
| BLK-02A-10 cross-branch contamination | — | **CLOSED** |
| BLK-02A-11 lossless modes discarded | — | **CLOSED** |
| BLK-02A-12 `k=0` HIGH confidence | — | **CLOSED** |

**Public verdict remains CONDITIONAL PASS.** Task 1R repaired evidence; it did
not pass G4 or G5, and it does not authorise any reflection-coefficient result.
