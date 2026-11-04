# Phase 2A Task 2 — Open Blockers

This file records what is **not** resolved. Per the brief, a converged boundary
solve is not evidence of a validated physical result, and no gate is closed
without sufficient evidence.

---

## 1. Gate status summary

| gate | status | evidence gap |
|---|---|---|
| **G4** — coupled piezothermoelastic energy conservation | **OPEN** | no independent coupled energy gate passed; the classical energy checks in `classical_benchmarks.py` do not carry over to the coupled problem |
| **G5** — coupled thermodynamic consistency | **OPEN** | no independent coupled thermodynamic gate passed; the elastic-limit and electrostatic/thermal-limit checks in `state.py` validate the *state definition*, not the *scattering result* |
| G0, G1 | **OPEN** (inherited) | unchanged from Task 1R |
| `BLK-02A-04` | **CLOSED** | 14/14 lossless branches, independent oracle |
| `BLK-02A-10`, `-11`, `-12` | **CLOSED** | inherited |

**Consequence.** Every reflection/transmission/absorption amplitude produced by
`scattering.solver` is labelled **`DIAGNOSTIC`**. That label is carried in the
returned record (`"label": "DIAGNOSTIC"`) and is not cleared by convergence,
by a small interface residual, or by rescaling invariance.

---

## 2. `BLK-02A-13` — NEW: coupled energy accounting not yet gated

**Severity** P0. **Status** OPEN.

The coupled solver produces complex amplitudes and closes the interface to
`2.7e-51`, but **no independent energy balance has been established** for the
coupled system. The classical `R + T = 1` identity does not transfer directly:

1. The piezothermoelastic flux `<J_z>` includes electric (`iωD_z φ*`) and
   thermal (`q_z ϑ*/T₀`) contributions whose signs and normalisations are
   convention-dependent. `state.py` implements one convention and validates it
   in three limits, but the *sum over modes* has not been checked against an
   independent total-field computation.
2. Four of the eight channels are diffusive (thermal) rather than propagating,
   so "reflected" and "transmitted" do not partition the flux the way they do
   for two acoustic branches.
3. Energy absorbed by the thermal and electrical channels has no independent
   check yet — absorption `A` cannot be claimed.

**What would close it:** a total-field flux computation
`⟨J_z⟩ = −½Re[σ_xz(iωu_x)* + σ_zz(iωu_z)* + iωD_z φ*] + ½Re[q_z ϑ*]/T₀`
evaluated on both sides of the interface from the reconstructed physical fields
and compared against the incident flux, plus a passivity bound
`0 ≤ A ≤ 1` verified independently. This is designed but **not implemented**.

---

## 3. `BLK-02A-14` — NEW: Al6061 | Al6061 at normal incidence is under-determined

**Severity** P1. **Status** OPEN.

At normal incidence for elastically identical half-spaces the assembled system
is 4 equations × 6 unknowns, rank 4, equilibrated condition number `1.1e21`. The
`u_x`, `φ` and `D_z` continuity rows are identically zero because the sagittal
geometry decouples them and the medium has no piezoelectric coupling.

Five of six coefficients are therefore not individually determined. This is
**physical degeneracy** (an elastically and electrically identical interface is
invisible), not a numerical failure — but the code reports the configuration as
unresolved rather than pretending otherwise.

**What would close it:** an explicit reduced formulation for the identical-medium
case that fixes the gauge freedom by an additional stated condition (for example
fixing the reflected branch to zero by symmetry), or exclusion of the degenerate
case from reported results with the exclusion documented. Not done: choosing a
regularisation would be exactly the "arbitrary regularisation" the brief
prohibits.

---

## 4. `BLK-02A-15` — NEW: two classification rules changed after Task 1R

**Severity** P1. **Status** OPEN (measurement complete, decision pending).

Two fixes in this task change modal labels relative to the published Task 1R
inventory:

| change | rows | why |
|---|---|---|
| `GAUGE` → `ELECTROSTATIC` | 4 | at `k ≠ 0` the `ξ = ∓ik` mode has `D_z = −ε₃₃k ≠ 0` and satisfies `D_z` continuity; calling it a gauge deleted a physical channel |
| `EVANESCENT_NO_TRANSPORT` → `INCOMING` | 42 | attenuation `6.1e-14 … 4.6e-7` per wavelength, i.e. decay lengths `8.4e5 … 1.6e7` m — propagation, not localisation |
| `GAUGE_EXCLUDED` → `OUTGOING` / `EVANESCENT_NO_TRANSPORT` | 8 | consequence of the same `GAUGE` → `ELECTROSTATIC` reclassification |

**The published artifact has not been modified.** `phase_02a_modal_inventory_v2.csv`
is byte-identical to its pinned sha256 (`a3b2965319a49fb2…`) and `INV00`
continues to enforce that. The legacy rules remain reachable through
`classify_mode(..., legacy_gauge=True)` and `classify_role(..., strict_legacy=True)`,
so the artifact stays reproducible, and the divergence is measured in
`validation/phase_02a_task2_convention_divergence.json`.

**Decision required from the supervisor:** whether to (a) issue a v3 inventory
under the new convention with the divergence recorded as a superseding change,
or (b) keep v2 as the published record and treat the new rules as
assembly-only. Option (b) is what is implemented now; it means the inventory and
the assembly use different conventions, which is defensible but must be stated
in any later citation of the inventory's role labels.

---

## 5. `BLK-02A-01` — `kernel_scaled` precision instability

**Severity** P0. **Status** OPEN (worked around).

The active numerical path still contains the legacy `kernel_scaled`. The new
`code/scattering/*` modules do **not** use it: `state.py`,
`elastic_reference.py`, `assembly.py` and `solver.py` all call
`foundation_checks.M_mp` / `flux_z_mp` directly. The legacy file is retained and
`B10` pins its defective behaviour so the historical evidence is not lost.

**Not done:** replacing `kernel_scaled` inside the legacy path. That would alter
Phase 1C/1D evidence and requires a separate superseding decision.

---

## 6. `BLK-02A-02`, `-03`, `-05` — inherited, unchanged

| blocker | status | note |
|---|---|---|
| `BLK-02A-02` gauge threshold empirical | OPEN | `GAUGE_MECHANICAL = 1e-40` is a measured machine-zero, not derived. The new `D_z` test removes the dependence for the *electrostatic* distinction, but the mechanical-zero test itself is still empirical |
| `BLK-02A-03` unreachable `jz_scale` arm | OPEN | still unreachable; not touched |
| `BLK-02A-05` `radiation_dir` enum hazard | OPEN | selection now goes through `outgoing_status` everywhere in `scattering/*`; the enum remains present in the inventory and in `radiation.py`'s reason strings. The hazard is mitigated, not removed |

### 6.1 A new threshold that must be justified in review

`ATTENUATION_PROPAGATING = 1e-4` per wavelength is a **stated physical
criterion**, not a fitted tolerance:

- the largest case that must propagate: `3e-13` per wavelength (PZT-5H SV, 30°)
- the smallest case that must stay localised: `6.7e-11` per wavelength
- the threshold sits ten orders above the first and four below the second

It is nonetheless a **new empirical constant**. It should be reviewed as such.
The `strict_legacy=True` flag makes its effect exactly reversible.

---

## 7. Numerical limits observed

| observation | value | implication |
|---|---|---|
| over-determined coupled path | `cond ≈ 2e4`, normal equations square it to `4e8` | residual floor `~1e-37`, not `1e-51`. Stated in the tolerance with the reason |
| square coupled path | `cond ≈ 1e4–3e4` | residual `2.7e-51`; trustworthy |
| identical-media interfaces | `σ_min ≈ 1e-5`, `cond 1.8e5` | near-singular by physics; reported, not regularised |
| Al\|Al normal incidence | `cond 1.1e21` | under-determined; `BLK-02A-14` |

`COND_TRUST = 1e8` is the threshold above which a coupled result is reported
but **not trusted**. It is not widened to obtain agreement.

---

## 8. What is explicitly *not* claimed

- No coupled `R`, `T` or `A` value is claimed as a validated piezothermoelastic
  result.
- No energy conservation statement is made for the coupled problem.
- No finite-interphase calculation, parameter optimisation, anti-reflection
  design or manuscript material has been produced — all out of scope.
- The classical solver is validated for **isotropic welded** contact only; the
  spring-contact path in the legacy file is untested and not part of this task.

---

## 9. Recommended next actions

1. Implement the coupled total-field flux gate (`BLK-02A-13`) — this is the
   shortest route to closing G4.
2. Decide v2-vs-v3 for the inventory convention (`BLK-02A-15`).
3. Either reduce the identical-medium normal-incidence case analytically or
   exclude it explicitly (`BLK-02A-14`).
4. Review `ATTENUATION_PROPAGATING` as a new empirical constant.


---

## TASK 2R VALIDITY MARKING

**Added by Task 2R.** The normal-incidence scattering amplitudes reported by
Task 2 are marked **INVALID / SUPERSEDED**.

**Reason.** The interface reduction removed a row when it was zero across the
admitted basis, without requiring it to be zero on the incident right-hand side.
At exactly `k = 0` that deleted the `u_x` and `txz` rows, which are the only
equations distinguishing a shear mode from a longitudinal one. Measured at
Al6061 `k = 0`: the shear state is `(u_x, u_z) = (1, 0)` with
`sigma_xz = 5.3147726e13j`.

**Additionally superseded:** the `phi` continuity row was a gauge, not a
boundary condition, at `k = 0`. The electric governing equation reduces to
`e33 u_z - eps33 phi = 0`, which forces `D_z = 0` identically on every mode while
`phi` remains fixed only up to a constant.

**Scope of the invalidation.** Only `k = 0` results are affected. The 78 oblique
configurations never exercised the defect -- the trivial-row reduction fires only
at exactly `k = 0` -- and their numbers are **RETAINED**.

**Regenerated under** experiment ID `PHASE_02A_TASK2R`:
`validation/phase_02a_task2r_independent_comparison.csv` and
`validation/phase_02a_task2r_test_results.json`.

**The original Task 2 reports and data are retained unchanged as historical
evidence.** See `reports/PHASE_02A_TASK2R_NORMAL_INCIDENCE_AUDIT.md` section 9.

**BLK-02A-14 is CLOSED as an artifact.** It recorded `Al|Al` normal incidence as
under-determined on the strength of `cond = 1.10e21` and a rescaling verdict of
`DEFERRED_TO_CONDITIONING`. With the row truncation repaired the system is 6x6
rank 6 for both polarizations, the reflected field is `1.1e-37`, and continuity
is exact.
