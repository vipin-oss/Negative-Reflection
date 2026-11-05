# PHASE 02A — TASK 2T FINAL ACCEPTANCE CHECK

**Branch** `arena/a3ede4fe-negative-reflection`
**Starting commit** `929ad6b7d49a9202fda079fb054bf97d8b7ce340`

---

## 0. Verdict

**CONDITIONAL ACCEPT.**

The Maxwell derivation, the dimensional-scaling fixes and the recovery of the missing near-gauge
channel are confirmed correct and are retained unchanged. **But this checkpoint found a defect that
supersedes my own Task 2T conclusion:** the production `k = 0` reduction discards the potential row
by an *incident-amplitude* test, and for every metal|piezoelectric configuration that returns a
potential that is wrong by O(1) in a normalisation-invariant measure. My Task 2T report presented
this as an O(k) accuracy limitation of `1.1e-5`. That was wrong, and the reason I missed it is
recorded in §2.4.

Three questions resolved:

| Q | subject | outcome |
|---|---|---|
| 1 | exact normal-incidence electrical contact | **RESOLVED** — the finite jump is a normal dipole layer; the paper forbids it; it is **not** admissible and must be imposed, not inferred |
| 2 | numerical limit and error analysis | **RESOLVED** — order 2 demonstrated for PZT\|PZT; the metal\|piezo failures are **not** a tolerance issue, they are O(1) and structural |
| 3 | near-gauge modal normalisation | **RESOLVED** — `hybrid_mode` is defective in two independent ways; a threshold-free basis fixes the mode loss and improves Al\|Al by 6 orders |

---

## 1. Q1 — Exact normal-incidence electrical contact

### 1.1 The tangential condition is the *only* potential condition Maxwell provides

Maxwell gives exactly two interface conditions here:

```
n . (D2 - D1) = sigma_free          from div D = rho_free
n x (E2 - E1) = 0                   from curl E = 0
```

In the sagittal plane, with `E_x = -i k phi`, the second reads

```
k [phi] = 0                                                    (T)
```

**"Continuity of potential" is not a separate assumption.** It is what (T) says wherever `k != 0`.
There is no additional independent condition available, and no choice to make between them.

### 1.2 What a finite potential jump physically is

A surface dipole layer with moment per unit area **p** normal to the interface shifts the potential
across the interface by `p / eps0` while leaving the tangential electric field continuous — because
`E_x = -d(phi)/dx` differentiates *along* the surface, and a constant jump cancels in the derivative.

Therefore, in this geometry:

```
[phi] != 0     <=>    a normal surface dipole layer is present
```

This is why (T) carries a factor `k` and the dipole condition does not: at `k = 0` the tangential
derivative that would detect the layer vanishes identically, so **Maxwell alone cannot constrain a
normal dipole layer at exactly normal incidence.** `k [phi] = 0` at `k = 0` is `0 = 0`.

This is the precise sense in which the two are distinct, and it is the answer to the supervisor's
first question: the tangential-E condition and potential continuity are the *same* condition for
`k != 0`, and they *diverge* at `k = 0` — not because potential is redundant, but because the
detector degenerates while the physical quantity remains.

### 1.3 Does the paper permit dipole layers? No.

The blueprint specifies the electrical contact explicitly (lines 279–282):

> **Electrical contact.** Choose one: continuous potential benchmark, or charge balance
> `D_z^(2) - D_z^(1) = sigma_free` with a passive capacitive gap `D_z = C_s [[phi]]`
> (`C_s >= 0` in F/m²). State whether free charge, storage or interfacial current is present.

Two options, and the paper takes the first: **continuous potential**. The capacitive-gap option is
the only one that admits `[[phi]] != 0`, and it is not adopted. The blueprint's `C_s = 0` limit of the
second option is exactly the first option.

**Therefore the paper forbids interfacial dipole layers. `[[phi]] = 0` is a required physical
condition at every `k`, including `k = 0`.**

### 1.4 Is the finite jump in the `k = 0` reduced solution admissible? No.

Measured, production `k = 0` against the `k -> 0` limit of the full eight-row system, reported as
`|phi| / field_scale` so the comparison is normalisation-invariant:

| configuration | production `k=0` | true limit | status |
|---|---|---|---|
| Al6061\|Al6061 P | `0.0` | `0.0` | agree |
| Al6061\|Al6061 SV | `0.0` | `0.0` | agree |
| PZT-5H\|PZT-5H P | `8.254e-6` | `8.254e-6` | agree to `7e-15` |
| PZT-5H\|PZT-5H SV | `0.0` | `2.896e-13` | limit is itself O(k) |
| **Al6061\|PZT-5H P** | **`8.254e-6`** | **`8.249e-6`** | agree to `6.3e-4` rel |
| **Al6061\|PZT-5H SV** | **`0.0`** | **`1.462e-5`** | **wrong by O(1)** |
| **PZT-5H\|Al6061 P** | **`0.0`** | **`1.644e-5`** | **wrong by O(1)** |
| **PZT-5H\|Al6061 SV** | **`0.0`** | **`2.285e-5`** | **wrong by O(1)** |

The jump in the transmitted potential is a *model output*, not an imposed quantity, and in four of
eight configurations it is **zero where the limit is not, or nonzero where the limit is `8.249e-6`**.
That is a dipole layer appearing accidentally. It is **not admissible** under a continuous-potential
benchmark.

**The `1.1e-5` figure I reported in Task 2T is not an accuracy limit. It is the size of an
unimposed physical condition, and the worst case is not 1.1e-5 but a complete failure.**

### 1.5 The potential condition is not redundant, and is not "automatically satisfied"

To be explicit, because Task 2T's wording invited this reading: `k [phi] = 0` degenerating at `k = 0`
does not make `[[]phi[]] = 0` vacuous. It makes *that equation* unable to express it. The condition
must be imposed by another route at `k = 0` — shown in §2.3 to be a genuine, independent constraint.

---

## 2. Q2 — Numerical limit and error analysis

### 2.1 Demonstrable asymptotic order

Al6061|PZT-5H, P incidence, `|phi/ref - production|` against `k/k_ac`, 70-digit arithmetic, with
`ref` a normalisation-invariant component ratio:

| `k/k_ac` | `|phi/ref - prod|` | local slope |
|---|---|---|
| `1e-2` | `5.01011e-9` | — |
| `3e-3` | `5.20191e-9` | `-0.03` |
| `1e-3` | `5.21877e-9` | `-0.003` |
| `1e-4` | `5.22086e-9` | `-3.0e-5` |
| `1e-6` | `5.22089e-9` | `-4.7e-9` |

The difference saturates at `5.22089e-9`. **The slope is zero: there is no asymptotic order for this
quantity, because production is not converging to the limit at all.** It lands on a different answer
and stays there.

**PZT-5H|PZT-5H, by contrast, shows clean order 2:**

| `k/k_ac` | `|phi/ref - prod|` | local slope |
|---|---|---|
| `1e-2` | `6.00751e-12` | — |
| `1e-3` | `6.00382e-14` | `2.0` |
| `1e-5` | `6.00378e-18` | `2.0` |
| `1e-6` | `6.00378e-20` | `2.0` |

Six exact decades at slope `2.0`. **This is the demonstrable order the supervisor asked for, and it
exists exactly where the reduction is correct.** Where it fails, the failure is not a tolerance.

### 2.2 Precision convergence

`phi_t` for Al6061|PZT-5H P at `k/k_ac = 1e-5`, identical to all printed digits at 60, 80, 100 and
140 digits, with solve residuals `1.3e-47`, `2.1e-67`, `4.8e-87`, `3.5e-127`. **The defect is not
numerical.** It is in the formulation.

### 2.3 Why: the reduction drops the potential row by an *incident-amplitude* test

`assembly.py` removes a row when the incident wave carries no content in it. The comment at that
site is explicit and, for the `u_x` case it was written for, correct. But the same rule is applied to
`phi`, and **Al6061 has `e = 0`, so an incident P or SV wave in aluminium has `phi = 0` exactly**.

Consequently production removes the `phi` row (and almost always the `Dz` row with it) whenever the
incident medium is non-piezoelectric:

| configuration | production keeps | drops | `phi` result |
|---|---|---|---|
| Al6061\|Al6061 P | 6 mechanical/thermal | `phi`, `Dz` | `0` — correct, both sides non-piezo |
| PZT-5H\|PZT-5H P | **plus `phi`** | `Dz` | `8.254e-6` — correct |
| Al6061\|PZT-5H P | 6 mechanical/thermal | `phi`, `Dz` | `8.254e-6` — `6.3e-4` off |
| PZT-5H\|Al6061 P | 6 mechanical/thermal | `phi`, `Dz` | `0` — wrong by O(1) |
| all SV in metal\|piezo | 6 mechanical/thermal | `phi`, `Dz` | `0` — wrong by O(1) |

**The rule conflates "the incident datum is empty" with "the transmitted field is unconstrained."**
For Al6061|Al6061 those coincide. For a metal|piezoelectric pair they do not: the incident side
carries no potential, but the transmitted piezoelectric side does, and its value is physical.

The empty-`rhs` and empty-`basis` tests are defensible — a row whose right-hand side and columns are
both zero genuinely adds nothing. The `inc_empty` test is not: it drops a row whose `rhs` is zero but
whose *columns are not*, which deletes a real constraint.

### 2.4 Why Task 2T missed this — a metric error, recorded

Task 2T measured the discrepancy with `asymptotic.dielectric_contact_system`, whose `flux`
normalisation **silently drops the non-piezoelectric modes** (`sqrt(0) = 0` -> `None`). For Al6061|Al6061
that is 2 of 8 modes. The resulting system cannot satisfy the mechanical rows, so it returned a
spurious `0.03234` residual and a wrong reference. **I compared production against a defective
reference, and the "6.342e-4" I reported was an artefact of that reference, not a property of
production.** This checkpoint's numbers use a basis that retains every admitted mode (§3).

### 2.5 Exact equivalence, asymptotic convergence, finite approximation — distinguished

The supervisor asked for these to be kept apart. They are three different statements and only the
third is ever available here:

| statement | status |
|---|---|
| **Exact mathematical equivalence** of the reduced and full `k=0` systems | **FALSE.** They are different operators. The reduced system imposes 6 or 7 conditions where the interface has 8. |
| **Asymptotic convergence**, `reduction -> limit` at a stated order in `k` | **TRUE for PZT-5H\|PZT-5H (order 2, six decades). FALSE for metal\|piezo**, where the slope is zero. |
| **Finite numerical approximation** at a stated precision | **TRUE everywhere**; residuals `1e-47` to `1e-127` at 60–140 digits. |

The `100x`-drift acceptance tolerance I used in Task 2T is **withdrawn**: it was a permissive
device that could not distinguish order-2 convergence from saturation. A tolerance is not evidence
of an order. The tests now assert the measured order instead.

---

## 3. Q3 — Near-gauge modal normalisation

### 3.1 `hybrid_mode` is defective in two independent ways

**(a) Its fallback is dead code for the modes it was written for.** The guard is `abs(jz) > 0`. For a
piezoelectric mode `|jz|` is never exactly zero — measured, it *grows* as `k` falls (`1.71e+38` at
`1e-2`, `1.71e+54` at `1e-10`) — so the branch is always taken and the potential fallback never runs.

**(b) For a piezoelectric mode the flux scale is exactly proportional to the potential scale.** This
is the decisive measurement, and it was not made in Task 2T:

| `k/k_ac` | `|jz|` | `|phi|` | `|jz| / |phi|^2` |
|---|---|---|---|
| `1e-2` | `1.71e+38` | `4.428e+14` | `8.722e+8` |
| `1e-5` | `1.71e+44` | `4.428e+17` | `8.722e+8` |
| `1e-10` | `1.71e+54` | `4.428e+22` | `8.722e+8` |

Constant to ten decades. So `sqrt(|jz|) ∝ |phi|`, and the "adaptive" flux branch is a potential
normalisation in disguise. `hybrid_mode` is not adaptive; it was a threshold-shaped patch.

### 3.2 The scale was never the real problem — mode *retention* was

`unit_flux_mode` returns `None` when `sqrt(|jz|) = 0`, which happens for every purely mechanical mode
in a non-piezoelectric material. Measured across all four pairs and both polarisations:

| pair | pol | admitted | `flux` | `hybrid` | **`max`** | `displacement` |
|---|---|---|---|---|---|---|
| Al6061\|Al6061 | P, SV | 8 | **6** | 8 | **8** | 6 |
| PZT-5H\|PZT-5H | P, SV | 8 | 8 | 8 | **8** | 8 |
| Al6061\|PZT-5H | P, SV | 8 | **7** | 8 | **8** | 7 |
| PZT-5H\|Al6061 | P, SV | 8 | **7** | 8 | **8** | 7 |

`hybrid_mode` does retain all 8 — through the *first* loop, because the discarded modes come from
`unit_flux_mode` inside the `flux` path, not from `hybrid`. So `hybrid` is not the mode-loser. It is
the wrong *reason*: it keeps them by accident of scale choice, not by a rule that guarantees it.

### 3.3 The replacement: `max_component_mode`, with no threshold at all

```
scale = max(|component|) over the eight state components
```

No physical judgement, no cutoff, no fitted parameter. It is nonzero for any nonzero modal vector by
construction, so it cannot discard a mode, and the per-column scale is a free constant the coefficient
absorbs — which is why the reconstructed field is scale-independent (§3.5).

### 3.4 Effect on the acceptance criteria

| criterion | `flux` | `hybrid` | **`max`** |
|---|---|---|---|
| modes retained, all pairs | 6–8 of 8 | 8 of 8 | **8 of 8** |
| Al6061\|PZT-5H P solve residual | `0.0323` | `1.3e-47` | **`1.6e-61`** |
| Al6061\|PZT-5H P boundary residuals | `4e-6` | `4e-54` | **`1e-66`** |
| k-order classifier usable | no | yes | **yes** |
| Al6061\|Al6061 P, reduction vs limit | n/a | `4.9e-13` | **`4.9e-13`** |
| PZT-5H\|PZT-5H P, same | n/a | — | **`2.1e-13`** |
| threshold fitted to a result | — | implicit `>0` | **none** |

Al6061|Al6061 improves from the `6.3e-4` order of magnitude reported in Task 2T to `4.9e-13` when
measured with a basis that keeps its modes — a factor of `1e9`, and it shows that the apparent
disagreement was the reference's fault.

### 3.5 Cross-basis comparison, done correctly

Absolute totals **cannot** be compared across bases: the basis sets the incident amplitude, and the
total scales with it. Measured — the `hybrid`/`max` ratio of every component of the transmitted field
equals the incident scaling factor `5884.71837416` with spread `3.2e-67`.

The correct comparison is a component **ratio**, which is invariant. Measured, `hybrid` vs `max` on
Al6061|PZT-5H P, every component normalised by `tzz`:

| component | hybrid | max | relative difference |
|---|---|---|---|
| `phi` | `4.363809215e-11 - 8.248585726e-6j` | *identical* | `7.9e-46` |
| `tzz` | `1.0` | `1.0` | `0.0` |
| `qz` | `-0.01830478142 + 0.01835276979j` | *identical* | `2.0e-51` |

All to `1e-46` or better. **The two bases describe the same field; only the amplitude label differs.**
My earlier claim that "the three bases give wildly different totals" was measuring the incident
normalisation, not the physics — an error worth recording, because it briefly looked like a
contradiction in the solver.

---

## 4. What must change, and what does not

**Does not change.** The Maxwell derivation (P0), the model table, the five-model comparison, the
degeneration analysis of both rows at `k = 0`, the `1/k` artefact finding, the `cond ∝ k^-2`
singularity, and every Task 2S scaling fix. All stand.

**Must change.**

1. `hybrid_mode` becomes `max_component_mode` throughout; `hybrid` is retained only as a documented
   counter-example.
2. The `inc_empty` row-removal rule must not be applied to `phi`. Whether the incident datum is empty
   is irrelevant to whether the transmitted potential is constrained.
3. The `k = 0` reduced system for a metal|piezoelectric pair must impose `[[phi]] = 0` explicitly, and
   if that over-determines the mode set the resolution is to re-derive the `k = 0` mode count — not to
   drop the condition. Measured: with an explicit `phi` row the system is `7x6` with residual `0.47`,
   so the mode set genuinely cannot supply 7 conditions. That discrepancy is a real, open derivation
   question and is recorded as a blocker rather than patched.
4. The `100x`-drift tolerance is withdrawn in favour of asserted asymptotic order.

---

## 5. Conclusion

**CONDITIONAL ACCEPT.**

Accepted: the corrected Maxwell equations, the dimensional-scaling fixes, the recovery of the missing
near-gauge channel, and the `1/k` artefact finding.

Conditional because: the production `k = 0` reduction is **correct for piezoelectric|piezoelectric and
non-piezoelectric|non-piezoelectric pairs** (order-2 convergence demonstrated) and **wrong by O(1) for
metal|piezoelectric pairs**, for a reason now identified exactly — an incident-amplitude test applied
to a transmitted-side condition. No claim about normal-incidence amplitudes across a metal|piezo
interface may be made until that is repaired.

Coupled energy conservation is **not** begun, and no `R + T + A` claim is made or implied.

---

## Appendix A — Evidence inventory

| artefact | content |
|---|---|
| `reports/PHASE_02A_TASK2T_FINAL_ACCEPTANCE.md` | this report: the three questions and the verdict |
| `validation/tests/test_18_scattering_task2t_accept.py` | `ACPT01`–`ACPT06`, six independent tests executing the production formulation |
| `validation/phase_02a_task2t_acceptance_results.json` | machine-readable form of every number quoted here |
| `code/scattering/asymptotic.py` | `max_component_mode` added; `hybrid_mode` retained as a documented counter-example |
| `validation/mutation_tests/run_mutations.py` | `MR18`–`MR21`; suite 43 mutations |
| `reports/PHASE_01C_OPEN_BLOCKERS.json` | `BLK-02A-26`..`29` added; `BLK-02A-24` superseded with closure evidence |
| `evidence/claim_registry.csv` | `CLM-02A-25`..`27` added, 67 claims |
| `evidence/equation_registry.csv` | `EQ-02A-ACPT-01`..`02` added, 66 equations |

**Suite:** 104 PASS, 1 FAIL (`ACPT02`, the open defect), 1 UNRESOLVED (`SCT02`, no coupled solver).
**Test results:** `ACPT01`, `ACPT03`, `ACPT04`, `ACPT05`, `ACPT06` **PASS**.
`ACPT02` **FAILS**, deliberately: it is the measurement of the open defect, and it is left failing
rather than weakened to green. The CI gate correspondingly reports **GATE NOT SATISFIED** on
`ACPT02` — which is the correct behaviour for an open P0 blocker, and is recorded as such rather
than being silenced.

**Mutations:** suite total **43/43 DETECTED**. `MR18` (finite jump declared admissible), `MR19`
(basis chosen by a fitted threshold) and `MR20` (reported incident inconsistent with the columns
used) are all **DETECTED**. `MR21` is a
*repair probe* rather than a defect: it removes the incident-amplitude test for `phi` and is
recorded to show whether the acceptance test is sensitive to the actual defect. Measured — the
single-line repair does **not** clear `ACPT02`, because the `phi` row is also removed by the
basis-empty test. That is why the defect is logged as needing a derivation rather than a patch.

**A detector error found and fixed.** `MR18`'s first detector included `ACPT02`, which already fails
on a pristine tree; the mutation harness correctly reported `INVALID_DETECTOR` because the detector
fired before the mutation was applied. Pointing it at `ACPT01` alone then made it *survive*, because
`ACPT01` did not read the model table. `ACPT01` now asserts the active declaration
(`model_spec(PAPER_MODEL)["phi"] == "CONTINUOUS"`) as well as measuring the identity, and `MR18` is
detected. Both intermediate states are recorded because a detector that passes for the wrong reason
is the failure mode this suite exists to prevent.

---

## Appendix B — Two errors made during this checkpoint, recorded

Consistent with the brief's requirement to report discrepancies honestly, including where this
work was itself wrong.

**B.1 Cross-basis totals compared directly.** I initially concluded that `hybrid` and `max` "give
wildly different totals" (`0.0323` / `0.0646` / `1.1e-5` for the same configuration) and briefly
treated this as a solver contradiction. It is not: a basis sets the incident amplitude, so absolute
totals *must* differ. The ratio of the two bases' totals equals the incident scaling factor
`5884.71837416` to `3.2e-67`, and every normalisation-invariant component ratio agrees to `1e-46`.
`ACPT05` now asserts this so the mistake cannot recur silently.

**B.2 The 5e-3 "residual" on Al|Al `t_xz`.** An unclassified Richardson extrapolation applied to a
quantity whose true form is `2.8668193·k` produced a spurious non-zero limit. The order must be
classified before extrapolating — the same lesson as Task 2T's `t_xz` error, arrived at again from a
different direction.

---

## Appendix C — What this checkpoint did NOT do

- No coupled energy conservation, no `R + T + A` claim, explicitly including the metal|piezo case.
- No anti-reflection optimisation, no finite-interphase modelling.
- No repair of the `phi` row removal. The defect is **diagnosed and recorded**, not patched: the
  correct treatment requires re-deriving the `k = 0` mode count (`BLK-02A-27`), and a patch applied
  before that derivation would be exactly the "remove a condition for numerical convenience" failure
  the brief forbids.
