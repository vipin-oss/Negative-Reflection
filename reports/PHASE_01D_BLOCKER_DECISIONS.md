# Phase 1D — Blocker Decisions

Every decision below records the criterion, the measurement, and the evidence path.
**No blocker was closed because a suite ran green.** Two were left open deliberately
when the brief offered an easy path to closing them.

---

## BLK-U3 — CLOSED

**Severity:** P0
**Title (original):** Coupled scattering amplitude normalisation / null-vector
direction accuracy

### The recorded diagnosis, and its status

> "The null-vector DIRECTION for the 6 mm characteristic matrix (cond ~1e22) is
> accurate only to ~1e-5 … It is the root cause of every coupled R/T/A result being
> unavailable."

**Status: WRONG. Retracted with reasons.** The direction is accurate to better than
1e-40. The residual came from two defects in `flux_z_mp`, both now fixed.

### Closure criteria (all five required by the brief)

#### 1. Production reproduces the independent formulation within justified tolerance

Two genuinely different constructions:

* **(a) production** — column equilibration + 3×3 Cramer solve for the free component
  (`kernel_scaled`, `code/foundation_checks.py`);
* **(b) reference** — two-sided row+column equilibration, eigenvector of the
  equilibrated Gram matrix `CᴴC` at **dps 120**, root Newton-refined at **dps 160**
  from the exact coefficient polynomial.

They share no code path beyond `M_mp`.

| quantity | value |
|---|---|
| worst direction difference, 92 coupled configurations | **6.39e-58** |
| worst energy-balance residual, same configurations | **3.25e-56** |
| declared tolerance (unchanged) | **1e-8** |

Evidence: `validation/phase_01d_independent_comparison.csv`

#### 2. Stable under increased precision

The residual must shrink with precision, not plateau. It does:

| dps | 60 | 80 | 100 | 140 |
|---|---|---|---|---|
| residual | 3.44e-62 | 3.06e-81 | 5.51e-102 | **8.38e-142** |

**This is the diagnostic that distinguishes a correct formulation from a structural
error.** The *original* defect was precision-independent at 1.8602e-05 from dps 30
through 120 — a plateau, which is the signature of a wrong formula rather than a
roundoff limit. Test `U3D` now asserts the shrinking behaviour.

#### 3. Physical boundary-state quantities agree

The balance is stated, and tested, in **traction** `σ_xz`, **electric displacement**
`D_z` and **heat flux** `q_z` — not in eigenvalues. The independent two-route flux
comparison agrees to **1.95e-16** (check `D1`).

> A small eigenvalue residual is explicitly *not* accepted as sufficient here. The
> original defect had a characteristic residual of 2.2e-27 while being 1.86e-5 wrong
> in the energy balance, which is precisely why the brief forbids relying on it.

#### 4. Degenerate and near-degenerate configurations handled explicitly

At the test wavenumber, branches 3/4 are an exact conjugate pair, and 1↔6 and 2↔5 are
conjugate branch pairs. Each branch is solved for its **own** null vector and measured
independently in `U3C`; none are mixed. 20 further configurations are classified
`DEGENERATE_EXCLUDED` in the comparison CSV (the `u_x` column collapsing to 1e-104 at
normal incidence) and 8 as `ELASTIC_DECOUPLED` (`θ ≡ 0` for Al6061) — excluded by
classification, not silently inverted.

#### 5. Executed regression tests reproduce it

`U3A`, `U3B`, `U3C`, `U3D` all PASS. Mutations `M12` (direction), `M13` (scaling),
`M14` (gauge removal), `M15` (relaxation operator) are all DETECTED.

### What the closure does *not* claim

Closing BLK-U3 removes an obstacle. It does **not** produce a coupled scattering
solution and does **not** pass G4 or G5. `SCT02` remains UNRESOLVED with zero cases
executed.

---

## BLK-LS-STORAGE — RETAINED OPEN

**Severity:** P1
**Status:** `MATHEMATICAL_PART_CLOSED_PHYSICAL_PART_OPEN`

The brief says: *"Do not infer a stored-energy law solely from a vanishing harmonic
average"* and *"Preserve the LS vanishing-time-derivative vs vanishing-stored-energy
question as open until a complete thermodynamic derivation settles it."*

### What the derivation settles

`theory/energy_flux_derivation.tex`, Proposition `prop:relax`:

```
⟨q·q̇⟩ = 0            identically, for every monochromatic field
⟨q·q⟩ = ½|Q̂|² ≠ 0   always
```

Hence the EIT stored-energy density `τ_q q·q/(2κT₀)` is genuinely non-zero, and the
earlier inference ("the relaxation term averages to zero, therefore there is no
stored energy") is invalid. This is now a permanent test (`LS01`), which fails if
either half is withdrawn.

### Why the blocker is not closed

The mathematical question is not the physical question. Two flux functionals exist and
must not be interchanged:

* the **physical first-order flux** `J = σu̇ + φḊ − (θ/T₀)q`, linear in `q`;
* the **cycle-averaged effective wave-energy flux**
  `⟨J_z⟩ = −½Re[σ_xz v_x* + σ_zz v_z* + iωD_zφ*] + ½Re[q_zθ*]/T₀`, bilinear.

The relaxation enters neither, but it *does* enter the energy density, and hence any
definition of absorbed energy built on a stored-energy change. **Which functional an
experiment measures requires a complete scattering solution.** That does not exist.

`promotion_blocked` remains `true`. Retaining it is a deliberate choice made against
the easier option of closing it on the strength of §3.1 alone.

---

## BLK-G4 — NOT PASSED (unchanged)

The brief: *"Do not assert G4/G5 PASS because a special elastic limit closes."*

The energy account now has an internally consistent derivation, positivity of bulk and
interface dissipation, and a production balance that closes to 3.25e-56 on propagating
modes. **None of that passes G4.** G4 requires the coupled energy partition to close
as a scattering statement, and no scattering solution exists. The elastic reference
(`SCT03`, Fresnel to 1e-15) is a *limit* and is explicitly not admissible as evidence
for G4.

---

## BLK-G5 — NOT PASSED (unchanged)

G5 requires a complete scattering solution with interface terms evaluated on the
actual fields. No such solution exists: `SCT02` is UNRESOLVED with `n_exec = 0`.
Readiness (`SCT03`–`SCT06`) is not a solution. G5 may remain open past Phase 1D; the
brief permits this explicitly.

---

## BLK-G0, BLK-G1 — NOT PASSED (unchanged)

Not in scope for Phase 1D, not touched, not promoted. Recorded as unchanged so that
the transition from "open" to "closed" cannot happen by omission.

---

## BLK-DOC-01 — CLOSED

The withdrawn `τ_q` causality floor was still presented as a boxed
`\textsc{Derived}` result in `theory/governing_equations.tex` (line ~328) two phases
after `CLM-H1-02` was marked INVALID. Found by test `THY02`, written during this
phase.

Now:

* the active inequality is removed and replaced with an explicit retraction;
* `CLM-H1-02` remains in `evidence/claim_registry.csv` as `INVALID` with its original
  statement intact (`THY02` fails if it is deleted);
* `H1`'s label no longer says "LS causality bounds";
* `M19` re-injects the floor and is DETECTED.

---

## BLK-U5 — OPEN (unchanged)

P1, out of scope for the Phase 1D brief. Not touched.

---

## Deliberately *not* closed

Recorded because each was a plausible route to a greener report:

| Temptation | Why refused |
|---|---|
| Close BLK-LS-STORAGE on Proposition `prop:relax` | It settles mathematics, not physics |
| Pass G4 because the elastic limit closes | Forbidden explicitly; a limit is not the coupled account |
| Pass G5 because readiness tests pass | Readiness is not a solution; `SCT02` executed 0 cases |
| Relax the `D2` tolerance to make the old residual pass | The identity was made to close by correcting physics; tolerance untouched at 1e-8 |
| Weaken `M06` when it survived | The survivor was a real validation gap (`DAT01` became a tautology); closed with a new independent check `PRV06` |
| Accept `M03`'s first result | The mutation targeted code the detector never ran; the survivor was a genuine finding |

---

## Gate state after Phase 1D

| Gate | State | Basis |
|---|---|---|
| G0 | **NOT PASSED** | unchanged |
| G1 | **NOT PASSED** | unchanged |
| G2 | PASSED | unchanged |
| G3 | PARTIAL | unchanged |
| G4 | **NOT PASSED** | no coupled scattering energy partition |
| G5 | **NOT PASSED** | `SCT02` UNRESOLVED, `n_exec = 0` |

`ci_gate.py` exits **0** — "no contradictions, every unknown is a registered, owned,
open item. This is NOT a statement that the research claims are accepted."
