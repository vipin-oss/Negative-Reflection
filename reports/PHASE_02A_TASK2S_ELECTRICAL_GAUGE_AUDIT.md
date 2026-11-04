# Phase 2A Task 2S — Electrical gauge versus interface-continuity audit

**Commit at generation:** see `validation/phase_02a_task2s_results.json` (`commit_at_generation`).
**Supersedes:** the Task 2R justification for the normal-incidence row removal, and
the Task 2R normal-incidence amplitudes.
**Gates:** G0, G1, G4, G5 remain **OPEN**. No `R + T + A = 1` is claimed anywhere in
this document, and none may be inferred from it.

---

## 0. The question, and the answer

The supervisor asked one question: **is the Task 2R reduction — which removes the
`phi` and `Dz` interface rows at normal incidence — equivalent to a properly
gauge-fixed full formulation, or does it merely happen to solve?**

**Answer: it is equivalent. But the reason Task 2R gave for it was wrong, and the
production assembly had two defects that a successful solve concealed.**

Three separate things were found, in increasing order of severity:

| # | Finding | Status |
|---|---|---|
| 1 | The rank-gain argument that licensed the row removal is **vacuous** — it cannot distinguish a redundant row from an inconsistent one | Argument withdrawn, replaced by a limit argument (§4) |
| 2 | Continuity rows were built from **per-side dimensionless** states, so for dissimilar media the "continuity" condition was between two different unit systems | Fixed; the physical violation was 0.2616 (§5) |
| 3 | Column normalisation by the **Euclidean norm over mixed physical units** annihilated the displacement content of every column | Fixed; distinct columns were scoring as proportional (§6) |

Item 2 is the one that matters most for the record. The production system solved
its own matrix to **2.3e-41** while violating physical traction continuity by
**26%**. A small residual is not evidence of a correct formulation. The supervisor
warned about exactly this; the warning was correct.

---

## 1. Normal-incidence electrical equations from the 6 mm constitutive relations

For a 6 mm (poled-ceramic) medium the constitutive relations in the reduced
form the code uses are

```
sigma_xx = c11 exx + c13 ezz - e31 Ez + (thermal)
sigma_zz = c13 exx + c33 ezz - e33 Ez + (thermal)
tau_xz   = c44 gamma_xz          - e15 Ex
D_x      = e15 gamma_xz + eps11 Ex
D_z      = e31 exx + e33 ezz + eps33 Ez + p theta
```

with the poling axis along `z`, `E = -grad(phi)`, and `p` the pyroelectric
coefficient. The state vector is

```
STATE_FIELDS = (ux, uz, phi, theta, txz, tzz, Dz, qz)
```

and for a plane wave `exp(i k x + i xi z - i w t)` the strain and field
components are

```
exx = i k ux        ezz = i xi uz        gamma_xz = i xi ux + i k uz
Ex  = -i k phi      Ez  = -i xi phi
```

(verified against `code/scattering/state.py`: `Ez = -I*xi*phi`).

Gauss's law with no free volume charge is `div D = 0`:

```
i k D_x + i xi D_z = 0
```

**At normal incidence `k = 0` this collapses to**

```
i xi D_z = 0     =>     D_z = 0     for every xi != 0.
```

This is not a modelling choice; it is the only content of Gauss's law at `k = 0`.
It was checked directly against the assembly's second row and against the state
vectors (`EQ-02A-GAUGE-A1`): `|i xi D_z|` is at the roundoff floor for every
admitted mode of both materials.

### 1.1 `phi` is slaved, not free

Substituting the constitutive law into `D_z = 0`:

```
e33 (i xi uz) + eps33 (-i xi phi) + p theta = 0
=>  -xi^2 (e33 uz - eps33 phi) + i p xi theta = 0
=>  phi = (e33/eps33) uz  -  i (p/eps33) theta / xi
```

`phi` is a **function of `(uz, theta)`** — it carries no independent information.
Substituting this expression back into `D_z` returns

```
D_z = -p theta + p theta = 0
```

identically. So `D_z = 0` is a **cancellation identity, not an independent
constraint**. Its numerical accuracy is therefore a property of the null vector,
not of the formulation: measured **2.1e-17** relative to the largest cancelled
term, against terms of order `1e7`. This floor is reported in
`EQ-02A-GAUGE-A1` rather than hidden behind an absolute bound.

### 1.2 The `xi = 0` sector

The `xi = 0` branch exists in every medium and its mode vector is

```
(ux, uz, phi, theta, txz, tzz, Dz, qz) = (0, 0, 1, 0, 0, 0, 0, 0)
```

— **a pure additive potential**, with no strain, no stress, no heat flux and no
charge. It is the electrostatic/gauge sector, and it is not a physical wave. The
distinction the brief asks for is exactly this: `xi != 0` modes are physical;
the `xi = 0` mode is a constant offset of the potential.

---

## 2. Additive potential constants, and the gauge count

Each half-space admits its own `xi = 0` branch, so each contributes an additive
constant:

```
phi_1(z) = phi_1^physical(z) + C_1
phi_2(z) = phi_2^physical(z) + C_2
```

The physical fields are unchanged by

```
(C_1, C_2) -> (C_1 + c, C_2 + c)      for any constant c,
```

because only `E = -grad(phi)` enters the constitutive laws. That is **one global
gauge freedom**, and it requires exactly **one fixing condition**.

The **difference `C_1 - C_2` is not gauge**: it is a physical jump in potential
across the interface, and for a bonded contact it is fixed to zero by potential
continuity. Confusing "one global gauge" with "a relative offset" is the error the
brief warns against, and it is the error `MR08` injects.

**Gauge count, as required by the acceptance gate:**

| quantity | value |
|---|---|
| additive potential constants per interface | **2** |
| global gauge freedoms | **1** |
| fixing conditions required | **1** |
| free relative offsets permitted for a bonded contact | **0** |

Verified in `EQ-02A-GAUGE-A2` for all four material pairs and both
polarisations: `n_global_gauge_conditions == 1` in every configuration, and the
`xi = 0` vector is a pure potential to `1e-40`.

---

## 3. Interface electrical boundary models — kept apart

The brief forbids silently interchanging them. The three models, stated
explicitly:

| model | mechanical | potential | normal `D` | physical picture |
|---|---|---|---|---|
| **Bonded (shorted)** — *the model this system imposes* | `u_i`, `t_iz` continuous | `phi` continuous | **free** — a surface charge may sit on the interface | welded contact, electrodes not insulated from each other |
| **Insulated (open)** | `u_i`, `t_iz` continuous | free to jump | `D_z` continuous | a thin insulating film between the media |
| **Coated / other** | as bonded | `phi` jumps by a prescribed function | `D_z` jumps by the coating charge | not modelled here |

**The system imposes the bonded model.** That is why `phi` continuity is the
electrical condition and `D_z` continuity is *not*.

At `k = 0` the picture sharpens considerably:

* the `D_z` row is Gauss's law, i.e. `D_z = 0`. Imposing it *across the
  interface* is the **insulated** condition. It is not a duplicate of the bonded
  condition — it is a **different physical model**;
* the `phi` row is a condition on a quantity that has no finite `k = 0` value
  (§4), so it cannot be imposed at `k = 0` at all.

Both rows are therefore removed at `k = 0`, and for **two different, separately
justified reasons**. `EQ-02A-GAUGE-A6` asserts this classification: at `k = 0`
`Dz` *must* appear among the removed rows (keeping it would impose insulation on a
bonded contact), and no mechanical or thermal row may ever be dropped. `MR09`
keeps the `Dz` row and is caught.

---

## 4. Why the removal is legitimate — and why Task 2R's reason was not

### 4.1 Task 2R's argument, and why it fails

Task 2R justified the removal by rank: adding the `phi` row yields **zero rank
gain**, so the row is implied.

**This is vacuous.** Once a system has at least as many rows as columns and the
columns are independent, the rank equals the column count no matter what the
extra row says. `EQ-02A-GAUGE-A3` demonstrates it directly: appending an
*arbitrary* row of constants — a row that certainly cannot be implied by the
interface equations — also produces **zero rank gain**. A test that cannot fail
is not a proof.

The measured zero is real (`phi`-row gain **0** on the mode-only subsystem for all
8 configurations, rank 6 -> 6), but it carries no evidential weight. It is
reported for completeness and explicitly **withdrawn as a justification**.

### 4.2 The valid argument: `k = 0` is a limit, not an evaluation

Polynomial mode vectors are amplitudes multiplying `exp(i k x + i xi z)`. As
`k -> 0` the **potential entry** of the piezoelectric branches grows like `1/k`:

| `k/k_ac` | `phi` of PZT branch 5 |
|---|---|
| `1e-2` | `2.63e11` |
| `1e-4` | `2.63e13` |
| `1e-6` | `2.63e15` |
| `1e-8` | `2.63e17` |

A clean factor of 100 per decade of `k`. There is therefore **no `k = 0` value of
the potential to impose continuity on**. The `k = 0` potential condition is the
limit of a condition, not the condition itself.

### 4.3 The consequence, measured

Assembling the full system at exactly `k = 0` anyway — the additive constants
absorb the divergent potential, so the system still "solves" — gives, for
Al6061|PZT-5H P-incidence, a transmitted/incident `t_zz` of

```
k = 0 (full)      1.3316249147 + 0.00115753873 j     <-- wrong in the 4th digit, invents an imaginary part
k = 1e-4          1.3317978949 + 0.00000483... j
k = 1e-5          1.33179789532 + 0.000004839... j
k = 0 (reduced)   1.33179789533 + 0.00000484... j     <-- the limit
```

**The reduced `k = 0` answer is the limit; the full `k = 0` answer is not.** This
was a defect in *my* full formulation, not in the reduction, and it is now fixed
by refusing `k = 0` in `full_system` with an explanatory error. `EQ-02A-GAUGE-A5`
enforces the limit property by extrapolation in `k`, an oracle independent of the
row-removal rule.

---

## 5. Defect 1 — continuity built in the wrong units

`Side.states` divides each field by **that medium's own** reference scales. For two
different media the scales differ, and subtracting the two sides is not a
continuity condition. Measured scale ratios, Al6061 vs PZT-5H at `k = 0`:

| field | ratio |
|---|---|
| `ux`, `uz` | 1.0 |
| `theta` | **2.46881** |
| `txz` | **0.738449** |
| `tzz` | **0.597113** |
| `qz` | **178.461** |
| `Dz` | **2.46731e-13** |
| `phi` | **5.58798e-10** |

`states1 - states2 = -incident` therefore compared `u/scale1` with `u/scale2`.

**Measured consequence, Al6061|PZT-5H, SV, `k = 0`:**

```
production coefficients  [0, -0.160877, 0, 0, 2.35032, 0]
  matrix residual        2.296e-41        <-- converged
  PHYSICAL continuity    0.2616           <-- 26% traction jump

raw-state coefficients   [0, -0.219604, 0, 0, 0.780396, 0]
  matrix residual        8.078e-28
  PHYSICAL continuity    2.492e-41
```

The coefficients differ in the **third digit**. The defect was invisible because
the residual was measured on the matrix the defective system had built for itself.

**Fix:** `assembly.py` builds columns and right-hand side from `raw_states`;
`solver.py` reconstructs the interface fields from `raw_states` too.

---

## 6. Defect 2 — Euclidean column norms over mixed units

With physical columns restored, the columns now span `1e-37` (`D_z`) to `1e22`
(`q_z`) and the system is unsolvable as written (`cond = 1.6e62` at 30°;
`mp.lu_solve` raises `matrix is numerically singular`). The obvious fix — divide
each column by its Euclidean norm — is **wrong**, because the norm is dominated by
the loudest component:

```
Al6061|PZT-5H, SV: column q_z ~ 1e19, column u_x = 1
after Euclidean normalisation:
  Al  shear column  (1.88e-14, -1j)      u_x information gone
  PZT shear column  (1.20e-14, -1j)
  Cauchy-Schwarz overlap = 1.0           "duplicate column" -- false
```

The two columns are *not* duplicates; they differ by the shear-modulus ratio
`sqrt(C44_Al / C44_PZT) = 1.066`. The metric destroyed the evidence.

**Fix — diagonal equilibration, both factors invertible, so the solution set is
unchanged:**

```
row scale    s_i = max_j |state_i| over admitted modes, per field
column scale n_j = max_i |state_i / s_i| over the fields of mode j
```

Effect: `cond` falls from **1.6e62 to 7.7e9** at 30°, and every oblique residual
lands at the arithmetic floor. Amplitudes are recovered as `a_j = x_j / n_j`.

---

## 7. Equivalence result

Full formulation: all eight interface rows, the `xi = 0` sector kept as explicit
constants `C_1`, `C_2`, and **one** gauge-fixing row (`C_1 = 0`).
Reduced formulation: the production `assemble`, with the rows it removes at
`k = 0` and no constants.

Because the full formulation has no valid `k = 0` evaluation, equivalence is
tested at `k > 0`, where both are well posed, and the `k = 0` reduction is then
connected to it by the limit test of §4.

**Result — 4 material pairs x {P, SV} x `k/k_ac` in {1e-2, 1e-4}, all
gauge-invariant components (`ux, uz, theta, txz, tzz, Dz, qz`):**

```
worst relative disagreement between full and reduced:  1.707e-48
  located at PZT-5H|Al6061, SV, k/k_ac = 1e-4, component q_z
full system rank:     full for every configuration (8 of 8 unknowns, 9 rows)
reduced system rank:  full for every configuration
```

`phi` is **excluded** from the comparison by construction, and this is the point
the brief permits: the two formulations differ from `phi` by a constant on each
side, and a common gauge constant is allowed. The measured offsets were exactly
constant on each side (all seven other components agreeing to `1e-48`), which is
what "differs by a gauge constant" means.

Detail in `reports/PHASE_02A_TASK2S_FULL_VS_REDUCED.md` and
`validation/phase_02a_task2s_gauge_comparison.csv`.

---

## 8. Precision sensitivity

The brief forbids inferring rank from Gram eigenvalues without reporting precision
dependence. All five working precisions were run (`EQ-02A-COND-A1`):

| dps | rank (Al\|PZT, PZT\|Al, `k = 0`) | solve residual | cond | null-vector worst residual |
|---|---|---|---|---|
| 30 | 6 / 6 | 0.0 | 1.05e3 | 1.26e-35 |
| 40 | 6 / 6 | 1.1e-41 | 1.05e3 | 6.80e-46 |
| 60 | 6 / 6 | 7.8e-62 | 1.05e3 | 1.84e-66 |
| 80 | 6 / 6 | 0.0 | 1.05e3 | 8.50e-86 |
| 100 | 6 / 6 | 0.0 | 1.05e3 | 1.98e-106 |

**The rank does not move at any precision, and the residual falls as digits are
added.** A genuine rank deficiency would persist; a numerical artefact would
improve with precision. Neither happens, so full rank is a property of the system,
not of the arithmetic. Condition numbers were measured with `mp.eighe` on the
Hermitian part of `A^H A`, because `mp.svd` and `mp.cond` both raise on complex
input.

---

## 9. Mutation testing — including the electrical interface conditions

`validation/mutation_tests/run_mutations.py` now runs **34** mutations, all
detected (**34/34**). The detectors execute the **production** `assemble`; the
four added here inject deliberately incorrect *gauge treatments*, as the brief
requires:

| id | mutation | detector | result |
|---|---|---|---|
| MR08 | relative potential offset treated as free at every angle | `EQ-02A-GAUGE-A6`, `04` | DETECTED |
| MR09 | `D_z` continuity kept at `k = 0`, imposing insulation on a bonded joint | `EQ-02A-GAUGE-A6` | DETECTED |
| MR10 | global gauge condition removed from the full system | `EQ-02A-GAUGE-A2` | DETECTED |
| MR11 | per-field unit normalisation removed (raw SI columns span 1e59) | `EQ-02A-GAUGE-A5`, `COND-01` | DETECTED |

### 9.1 Two mutations were initially survivors — both findings are kept

**MR09 initially survived**, and it was a real gap: `EQ-02A-GAUGE-A6` only
checked that no *mechanical* row was dropped, so production could keep the `Dz`
row — silently imposing the insulated model — and still pass. The test now
asserts that `Dz` **must** be removed at `k = 0` and states why. This is exactly
the failure mode the brief names: silently interchanging electrical boundary
models.

**MR08 was initially a no-op** — the mutation inserted a `pass` and changed
nothing, so its "survived" verdict was an artefact of the mutation, not a test
gap. It now genuinely removes the `phi` row at every angle.

**MR11 was initially mis-specified.** The first version restored the Euclidean
column norm on the theory that it was the production defect. It is not: with the
per-field row scale already applied, a Euclidean norm is one scalar per column — a
legitimate diagonal rescaling that cannot change rank or solution set — and the
probe correctly reported it undetectable. The real defect was the *absence* of
per-field normalisation, so that is what is now injected. Detecting it required
moving the detector from a rank test to a conditioning test.

---

## 10. Acceptance gates

| requirement | status | evidence |
|---|---|---|
| interface electrical boundary model explicitly specified | **MET** | §3, `EQ-02A-GAUGE-A6` |
| gauge count mathematically justified | **MET** | §2, `EQ-02A-GAUGE-A2` |
| full and reduced agree in physical observables | **MET** | §7, 1.707e-48, `EQ-02A-GAUGE-A4` |
| no required interface equation simply discarded | **MET** | §4 — each removal has a stated reason; `Dz` removal is *required*, not tolerated |
| independent tests detect an incorrect gauge treatment | **MET** | §9, MR08–MR11, 34/34 |
| rank **not** inferred from Gram eigenvalues alone | **MET** | §8, five precisions |
| no `R + T + A = 1` claimed | **MET** | not claimed; G0/G1/G4/G5 open |

---

## 11. What changed in the code

| file | change |
|---|---|
| `code/scattering/assembly.py` | columns and rhs from `raw_states`; per-row triviality tolerance (a matrix-wide scale misfired once the rows are SI and span 1e-37..1e19); diagonal row+column equilibration; `row_scale` exported |
| `code/scattering/solver.py` | `boundary_residuals` and `_mode_fields` reconstruct from `raw_states` |
| `code/scattering/gauge.py` | `solve_path` gains column equilibration; `full_system` refuses `k = 0` with a derivation-backed message |
| `validation/tests/test_16_scattering_task2s.py` | **new**, 7 cases |
| `validation/tests/test_15_scattering_task2r.py` | SCTR04 rebuilt on rank (leave-one-out), Cauchy–Schwarz on mixed-unit columns withdrawn |
| `validation/tests/test_14_scattering_task2.py` | SCTA05 residual made **relative**, justified in the record |
| `validation/scattering/assembly_validation.py` | enforced-row residual made **relative** |
| `validation/mutation_tests/run_mutations.py` | MR08–MR11 added |
| `evidence/equation_registry.csv` | 7 entries added (`EQ-02A-GAUGE-A1..A6`, `EQ-02A-COND-A1`); `EQ-02A-GAUGE-01` **corrected** |
| `reports/PHASE_01C_OPEN_BLOCKERS.json` | BLK-02A-19/20/21 closed with closure evidence; BLK-02A-22 (oblique gauge) and BLK-02A-23 (bonded vs insulated) opened |

### 11.1 An ID collision, and why it matters

`evidence/equation_registry.csv` already contained an entry `EQ-02A-GAUGE-01` for a
*different* claim, so the Task 2S checks were renamed into the `EQ-02A-GAUGE-A`
(audit) series rather than shadowing it.

That existing entry then had to be corrected anyway, and the correction is
instructive: it stated the removal reason as *"at k=0 the electric governing
equation forces `D_z = 0` identically, so `phi` is fixed only through its
gradient"* — which repeats the withdrawn reasoning. Its name, summary, oracle,
evidence, failure mode and notes now carry the Task 2S derivation and the measured
`1.3316249147 + 0.00115753873j` versus `1.3317978953` comparison. The
`verification_status` stays `VERIFIED_STRONG`, because the *conclusion* it records
was never wrong — only the argument offered for it.

Two tolerance changes deserve their own note, because changing a tolerance is
where a task can quietly turn a failure into a pass:

* **SCTA05** compared an *absolute* interface residual against `1e-40` while the
  fields it reconstructs reach `1e20` (`q_z`). That bound was measuring the
  magnitude of the fields, not the quality of the solve. It is now relative, and
  the measured value is `1.76e-34` against fields of order `1e20` — a relative
  `1e-54`. **The accuracy is unchanged; only the units are.**
* **assembly validation** C4 had the same shape, `1e-24` absolute against
  `1e20` fields, relative `1e-47`. Same treatment.

Neither change loosens a physical claim. Both are recorded in the test
justifications with the measured numbers.

---

## 12. Superseded

* **Task 2R's rank-gain justification** — withdrawn; vacuous (§4.1).
* **Task 2R normal-incidence amplitudes** — superseded; computed with the per-side
  dimensionless continuity defect of §5.
* **Task 2R's `SCTR04` duplicate-column evidence** — withdrawn; it used a
  Cauchy–Schwarz overlap on mixed-unit columns, which reported `1.0` for genuinely
  distinct columns.

These are preserved as superseded rather than deleted.

---

## 13. What this task does **not** establish

* No energy conservation result. **`R + T + A = 1` is not claimed and must not be
  inferred.** G4 and G5 remain open.
* No claim that the full formulation is redundant. It is the reference the
  reduction was checked against, and its `k > 0` behaviour is what makes the
  reduction meaningful.
* No oblique-incidence gauge analysis. Equivalence is demonstrated at `k > 0`
  approaching normal incidence and at the `k = 0` limit; the oblique case inherits
  the corrected scaling but was not the subject of the derivation.

**Task 2S is complete. Work stops here pending supervisor review of this audit.**
