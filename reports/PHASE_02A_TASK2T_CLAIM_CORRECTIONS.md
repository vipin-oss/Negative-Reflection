# PHASE 02A — TASK 2T
# Claim Corrections: Four Scientific Claims Audited

**Branch** `arena/a3ede4fe-negative-reflection`
**Starting commit** `006d059164d7f8910d432932fa70b7411e821ad9`

This report audits the four claims named in the Task 2T brief. Each is assigned
**SUPPORTED**, **CONDITIONAL**, **INCORRECT** or **UNRESOLVED**, with the measurement behind the
verdict. Where a claim is INCORRECT the original wording is preserved here as historical evidence and
the superseding statement is given.

---

## Claim A

> "`D_z` continuity at `k = 0` imposes insulation."

### Verdict: **INCORRECT**

### Why

`D_z` continuity is `n · (D2 − D1) = sigma_free` with `sigma_free = 0`. That says **no free charge is
stored on the interface**. It does not say the interface is insulated, and insulation is not the
distinguishing feature.

An ideal dielectric contact and an open circuit satisfy the *same* local condition. They differ in a
**global** property — whether the total charge on the isolated body is constrained — which a single
local row cannot express. Measured: at `k/k_ac = 1e-3` the `CONDUCTING_ELECTRODE`, `SHORT_CIRCUIT`
and `OPEN_CIRCUIT` models return numerically **identical** systems
(`validation/phase_02a_task2t_boundary_comparison.csv`). If `D_z` continuity identified insulation,
these three would not be degenerate.

**Additionally**, at `k = 0` the condition is not merely uninformative about insulation — it is
**vacuous**. Gauss's law reduces to `i·xi·D_z = 0`, giving `D_z = 0` on both sides, so the jump is
`0 − 0 = 0` whatever the field. `EQ-02A-MAXW-02` and the derivation in the boundary audit, §3.2.

### Superseded by

`CLM-02A-20`: "For a dielectric contact with no free surface charge, continuity of potential AND
continuity of `D_z` hold simultaneously; `D_z` continuity alone does not identify an insulated
interface."

### Consequence for the codebase

`BLK-02A-23` was re-derived in Task 2T. The reference configuration is a **dielectric contact**,
stated as such, with five attached assumptions that would each have to be re-examined for a different
experiment.

---

## Claim B

> "The bonded model imposes potential continuity but not displacement continuity."

### Verdict: **INCORRECT**

### Why

The two conditions are not alternatives, and the system does not choose between them. They descend
from two different Maxwell equations:

- potential continuity from `curl E = 0` (tangential `E`), and
- normal-`D` continuity from `div D = 0` (normal flux).

Neither implies the other, and there is no regime in which a dielectric contact has one without the
other. The claim is not a simplification — it describes a boundary condition that no interface in
this problem satisfies.

### Measurement

Solving the full eight-row dielectric-contact system in an adaptive basis at `k/k_ac = 1e-3`:

| quantity | residual |
|---|---|
| potential jump | `9.237e-38` |
| normal-`D` jump | `2.032e-42` |
| linear-solve residual | `2.976e-29` |
| condition number | `3.969e9` |

**Both hold simultaneously, to 1e-38 or better, in a well-conditioned system.** There is no conflict
to resolve.

### Superseded by

`CLM-02A-20`, as for Claim A.

### Note on how this claim arose

It is the natural conclusion to draw if one imposes both rows and finds the system unsolvable — which
is exactly what happened in Task 2S. But the unsolvability was a **basis defect**, not a physical
incompatibility: the unit-flux basis drops the near-gauge mode that carries the potential jump. See
`PHASE_02A_TASK2T_NORMAL_INCIDENCE_LIMIT.md`, §5.2 and §4.4.

---

## Claim C

> "A modal potential diverging as `1/k` proves `phi` continuity cannot be imposed at `k = 0`."

### Verdict: **INCORRECT** — with a CORRECT sub-observation inside it

### The part that is wrong

A modal eigenvector has arbitrary normalisation, so "the potential diverges" is a statement about the
scaling convention and not about the field. The production convention (`roots._normalise`, pinning
`u_x = 1`) pins a component that **itself vanishes like `k`**, which manufactures a `1/k` in every
other component.

Measured three independent ways on PZT-5H branch 5 over `k/k_ac` from `1e-1` to `1e-6`:

| normalisation | `phi` over six decades |
|---|---|
| `\|u_z\| = 1` | spread `2.78e-5` (bounded) |
| displacement magnitude | spread `2.33e-3` (bounded) |
| `\|u_x\| = 1` (production) | grows by `1.002e5` |

The mechanism is confirmed independently: `u_x/u_z = 1.046e-6 · (k/k_ac)`, exactly linear.

**And the conclusion is refuted directly.** Potential continuity is imposed at `k = 0` — in the
reduction it is not imposed only because it is vacuous there, and in the full system it holds to
`9.2e-38`. It was never impossible.

### The part that is correct

There *is* a genuine singularity at `k = 0`, and it is worth separating carefully because it is the
grain of truth that made the wrong claim persuasive:

> A physical mode degenerates into the gauge mode as `k -> 0`: its `xi ∝ k`, the condition number
> grows like `k^-2` (`3.949e16` at `k/k_ac = 1e-7`), and its modal coefficient must grow like `1/k`
> while the total field stays finite.

That is real, is measured (`EQ-02A-COND-T1`), and is recorded as `BLK-02A-25`. It means the `k = 0`
answer must be **extrapolated**, not evaluated. It does **not** mean any boundary condition is
unimposable.

### Superseded by

`CLM-02A-21` (the `1/k` is a normalisation artefact) and `CLM-02A-24` (the singularity is a
degeneration, not an obstruction).

### Why this claim was especially dangerous

It converts a **numerical conditioning problem** into a **physical impossibility**. That inversion
would have justified permanently discarding a genuine boundary condition on the grounds that physics
forbade it, and the recorded reason would have been unfalsifiable by construction — any failure to
impose the row would appear to confirm it.

---

## Claim D

> "Full-versus-reduced agreement at nonzero `k` proves normal-incidence row removal."

### Verdict: **INCORRECT**

### Why

Two systems implementing the same boundary equations at `k > 0` will agree at `k > 0`. That is a
consistency check on the shared equations. It is silent about the limit, because the limit is a
different operator: rows that carry information at `k > 0` become vacuous at `k = 0`, and whether the
reduced system still lands on the correct limit is precisely what is in question.

The brief states this directly: "Agreement between two systems implementing the same boundary
equations at `k>0` is not sufficient proof of the exact `k=0` reduction."

### What was done instead

The reduction was validated against the **`k -> 0` limit of the full system**, obtained by:

1. solving the full eight-row system at a sequence of decreasing `k` in a `k`-independent basis,
2. classifying each component's order in `k` by measured log-log slope,
3. extrapolating only `CONSTANT` components, and
4. comparing against the production `k = 0` reduction.

All seven material pairs, both polarisations: worst `CONSTANT`-component disagreement `1.6e-7` to
`1.5e-4` of the field scale, each within 100× the limit's own measured convergence.

### Two things this corrected

1. **The reduction is not exact.** Potential continuity is approached at `O(k)`, worst `1.1e-5` of
   field scale, because that row is vacuous at `k = 0` and is therefore not imposed. The reduction is
   the limit, not the operator. Recorded as `BLK-02A-24`.
2. **The naive extrapolation was wrong.** Applying Richardson extrapolation to an unclassified
   component produced a spurious `5e-4` residual on Al6061|Al6061 P `t_xz`, whose true form is
   `2.8668193·k` — a cancellation between sides, limit zero, and production's `0` was right. The
   order must be classified **before** extrapolating.

### Superseded by

`CLM-02A-23` (the reduction is accurate to `O(k)`, not exact) and `EQ-02A-LIM-02` (agreement is
established against the limit, not against a finite-`k` sibling).

---

## Summary table

| claim | verdict | key measurement |
|---|---|---|
| **A** `D_z` continuity at `k = 0` imposes insulation | **INCORRECT** | three "free" models numerically identical; `D_z` jump `≤ 3.4e-71` is a cancellation identity, not physical content |
| **B** bonded model imposes `phi` but not `D_z` continuity | **INCORRECT** | both hold to `9.2e-38` / `2.0e-42` simultaneously, cond `3.969e9` |
| **C** `1/k` potential proves `phi` continuity impossible | **INCORRECT** (contains a correct sub-observation about `cond ∝ k^-2`) | `\|u_z\|=1` bounded (spread `2.78e-5`); `\|u_x\|=1` grows `1.002e5`; `u_x/u_z = 1.046e-6·k` |
| **D** finite-`k` agreement proves the removal | **INCORRECT** | agreement established against the extrapolated limit, `1.6e-7`–`1.5e-4` |

**None of the four is SUPPORTED.** All four are recorded here as historical evidence and superseded
in the registries.

---

## Registry changes

### `evidence/claim_registry.csv` — 64 claims

- `CLM-02A-13` **rewritten**. Old statement: *"A potential row is a gauge, not a boundary condition,
  at normal incidence; posing it as continuity forces the transmitted potential to zero and makes the
  mixed system unsolvable."* New statement records the degeneration of `k[phi] = 0` instead, and the
  label is annotated *"CORRECTED in Task 2T; the earlier gauge wording is superseded"*.
- `CLM-02A-20` **added** — both conditions hold simultaneously; `D_z` continuity is not insulation.
- `CLM-02A-21` **added** — the `1/k` is a normalisation artefact; the total field is finite.
- `CLM-02A-22` **added** — both rows are identically satisfied at `k = 0` and neither is a model
  choice.
- `CLM-02A-23` **added** — the reduction is accurate to `O(k)`, not exact.
- `CLM-02A-24` **added** — the `k -> 0` problem is singular by mode degeneration, `cond ∝ k^-2`.

`CLM-02A-14` remains **unbound**, as required.

### `evidence/equation_registry.csv` — 64 equations

- `EQ-02A-GAUGE-01` **rewritten and renamed**. The old name described the potential row as a gauge.
  The new name records the degeneration, the `latex_summary` gives the derivation from `E_x = -i k phi`,
  and the `notes` field records the conditional promotion explicitly:

  > Task 2S gave two different wrong reasons in succession: first a vacuous rank-gain argument, then a
  > `1/k` divergence argument drawn from a single eigenvector normalised to `u_x = 1`. Both are
  > withdrawn. The correct statement is the one recorded here.

- `EQ-02A-MAXW-01` **added** — the five electrical models, from Maxwell.
- `EQ-02A-LIM-01` **added** — the `1/k` is an artefact; this entry also records the `k_order` bug.
- `EQ-02A-LIM-02` **added** — the reduction agrees with the limit on classified components.
- `EQ-02A-COND-T1` **added** — the degeneration of the near-gauge mode, `cond ∝ k^-2`.

### `reports/PHASE_01C_OPEN_BLOCKERS.json` — 18 blockers

- `BLK-02A-19` **re-closed** with a correction notice. The *dimensional-scaling* content stands
  unchanged and is retained; the *electrical interpretation* attached to it in Task 2S is withdrawn.
  The closure record now separates the two.

  This is deliberate: the blocker's original defect was real and its fix is retained, but a closure
  record must not be allowed to launder an invalid interpretation into permanence.

- `BLK-02A-23` **re-derived** — the reference configuration is a dielectric contact, now derived
  rather than asserted.
- `BLK-02A-24` **added** — the reduction enforces potential continuity only to `O(k)`.
- `BLK-02A-25` **added** — the `k -> 0` system is numerically singular; the limit is extrapolated.
- `BLK-02A-22` **added** — oblique incidence not analysed for the limit.

---

## What was preserved and why

The brief requires invalidated statements to be **preserved as historical evidence** while being
superseded. Three mechanisms are used:

1. **In-place rewrite with an explicit correction note.** `CLM-02A-13` and `EQ-02A-GAUGE-01` keep
   their IDs and their old wording is quoted in this report and in the registry `notes` field.
2. **Superseding claims rather than edits.** `CLM-02A-20` through `CLM-02A-24` are new rows, so the
   trail from the old claim to the new one is visible in the registry itself.
3. **Blockers re-closed, not deleted.** `BLK-02A-19` records both what was right and what was wrong
   about its Task 2S closure.

No claim, equation or blocker was deleted. The Task 2S report files
(`PHASE_02A_TASK2S_ELECTRICAL_GAUGE_AUDIT.md`, `PHASE_02A_TASK2S_FULL_VS_REDUCED.md`) are **retained
unmodified** and are marked as superseded by reference from this file.

---

## Appendix: `GAUGE06` was re-derived, not deleted

`validation/tests/test_16_scattering_task2s.py::test_GAUGE06` asserted, in its Task 2S form:

> "`D_z` MUST be removed: at `k = 0` it is Gauss's law `D_z = 0`, which is the **INSULATED**
> condition. Imposing it on a bonded contact is a different physical model."

That is Claim A, and it is incorrect. The row is still correctly removed, but for the vacuity reason
in §3.2 of the boundary audit, not because continuity is the wrong model.

**The test was re-derived rather than deleted**, because the part of it that does real work is
structural and remains valid: no mechanical or thermal row may be dropped at `k = 0`; every removed
row must carry a recorded reason; the kept plus removed counts must account for all eight. That is
the protection against quiet, convenience-driven row removal, and deleting the test would have
removed it.

Mutation `MR09` (which makes production keep the `D_z` row at `k = 0`) is still DETECTED by the
re-derived test. All 7 Task 2S tests pass.

The withdrawn physics now lives only in `test_17_scattering_task2t.py`, where both conditions are
shown to hold simultaneously.

---

## Task 2T test and mutation inventory

| test id | what it establishes | status |
|---|---|---|
| `MAXW01` | the five models are distinct and correctly specified | PASS |
| `MAXW02` | `D_z` continuity is not insulation; both conditions hold | PASS |
| `MAXW03` | the tangential condition vanishes at `k = 0` | PASS |
| `LIM01` | the `1/k` is a normalisation artefact | PASS |
| `LIM02` | the total field has a finite limit (incl. positive controls) | PASS |
| `LIM03` | the reduction matches the limit on classified components | PASS |
| `LIM04` | the reduction does not enforce potential continuity exactly | PASS |
| `CONDT1` | the `k -> 0` problem is singular, `cond ∝ k^-2` | PASS |

| mutation | what it does | status |
|---|---|---|
| `MR12` | dielectric contact replaced by a conducting electrode | DETECTED |
| `MR13` | `D_z` continuity discarded with no derivation | DETECTED |
| `MR14` | zero free charge imposed when a charge is prescribed | DETECTED |
| `MR15` | global gauge confused with a relative potential jump | DETECTED |
| `MR16` | eigenvector rescaling treated as physical divergence | DETECTED |
| `MR17` | `k = 0` limit validated from a single finite-`k` sample | DETECTED |

Suite total: **40/40 mutations detected.**
