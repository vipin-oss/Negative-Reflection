# PHASE 02A — TASK 2T
# Maxwell Interface Conditions and the Electrical Boundary-Value Problem

**Branch** `arena/a3ede4fe-negative-reflection`
**Starting commit** `006d059164d7f8910d432932fa70b7411e821ad9` (Task 2S)
**Scope** Derive the electrical interface conditions this model actually imposes, from Maxwell's
equations, independently of any numerical convenience. Correct the electrical-boundary
interpretation of Task 2S. **Not** a broad rewrite: the dimensional-scaling and
numerical-conditioning fixes of Task 2S are retained unchanged.

---

## 0. Summary of what changed and what did not

**Retained from Task 2S, unmodified.**

- Per-side dimensionless `raw_states` alongside the dimensionless `states` (BLK-02A-19, dimensional
  half). The defect it fixed is real and was measured: building cross-interface rows from `.states`
  produced a matrix residual of `2.296e-41` while the *physical* traction continuity was violated by
  `0.2616`, i.e. 26%. That is a scaling bug, not an interpretation.
- Per-field row scaling, per-field column equilibration, `TRIVIAL_ROW_TOL`, and the
  `"row_scale"` diagnostics.

**Withdrawn from Task 2S.**

- The claim that the potential row is "a gauge, not a boundary condition."
- The claim that potential continuity "forces the transmitted potential to zero and makes the mixed
  system unsolvable."
- The claim that a modal potential diverging as `1/k` proves potential continuity cannot be imposed.

All three are **empirically false** and are corrected below (§3, §5) and in
`reports/PHASE_02A_TASK2T_CLAIM_CORRECTIONS.md`.

**New.** `code/scattering/maxwell.py` is the executable statement of the model table. It is the
authority for what the system claims to impose; no numerical routine determines the physics.

---

## 1. Derivation

### 1.1 The electroquasistatic reduction actually used by this model

The model works at `f = 1 MHz`, where the acoustic wavelength is of order millimetres and the
electromagnetic wavelength is of order `300 m`. Two reductions follow, and both are *approximations
that must be stated* rather than facts:

1. **Neglect magnetic induction.** `curl E = -dB/dt` becomes `curl E = 0`, so `E = -grad(phi)`. The
   neglected term is of order `(|v|/c)^2 ~ 1e-10` for the mechanical velocities involved.
2. **Neglect free volume charge in the bulk.** `div D = rho_free` becomes `div D = 0` in each
   material. Any free charge in this problem sits on the interface, which is where the interface
   condition enters.

Under these two assumptions the field equations in each bulk are

```
div D = 0 ,   curl E = 0 ,   E = -grad(phi) ,   D = eps E + e : grad_s u      (piezoelectric)
```

with the mechanical equations coupled in the usual way. Nothing below depends on the piezoelectric
constitutive law; it enters only through the relation between `D`, `E` and the displacement field.

### 1.2 The interface conditions

For a stationary or moving interface with unit normal `n` pointing from side 1 into side 2, the exact
Maxwell conditions are

```
n . (D2 - D1) = sigma_free          (normal D -- from div D = rho, integrated across the surface)
n x (E2 - E1) = 0                   (tangential E -- from curl E = 0)
n . (B2 - B1) = 0                   (identically satisfied: B is neglected)
n x (H2 - H1) = K_free              (decoupled from the problem: H is neglected)
```

Two rows survive, and they are **independent**. They come from two different Maxwell equations. One
is a flux condition; the other is a circulation condition. Neither implies the other, and no
combination of them is redundant.

### 1.3 What the two rows look like in this geometry

The problem is set in a sagittal plane with the interface at `z = 0`, normal `n = z-hat`, and a
horizontal wavevector `k` along `x`. Write the horizontal dependence as `exp(i k x)`. Then

```
E_x = -d phi/dx = -i k phi
```

and the tangential condition `n x (E2 - E1) = 0` becomes, in the sagittal plane,

```
E_x2 - E_x1 = 0     =>     -i k (phi_2 - phi_1) = 0     =>     k [phi] = 0
```

while the normal condition is

```
D_z2 - D_z1 = sigma_free
```

**These two equations are the whole electrical content of the problem.** Everything that follows is
consequences of them.

---

## 2. The five electrical models

The brief requires five cases to be stated separately. `code/scattering/maxwell.py` records each as
a row of `MODEL_SPEC`, with the equations that are actually imposed, so the table can be executed and
audited rather than read.

| model | potential | normal `D` | `sigma_free` | equations imposed |
|---|---|---|---|---|
| `DIELECTRIC_CONTACT` | CONTINUOUS | CONTINUOUS | ZERO | `phi`, `Dz` |
| `CONDUCTING_ELECTRODE` | EQUIPOTENTIAL | TERMINATES | FREE | `phi` only |
| `SHORT_CIRCUIT` | FIXED_ZERO | TERMINATES | FREE | `phi = 0` |
| `OPEN_CIRCUIT` | FLOATS | TOTAL_CHARGE_ZERO | FREE_ZERO_MEAN | `phi`, with a global charge constraint |
| `PRESCRIBED_CHARGE` | CONTINUOUS | JUMPS_BY_GIVEN | GIVEN | `phi`, `Dz - sigma` |

### 2.1 (1) Ideal dielectric–dielectric interface, no free surface charge

`sigma_free = 0`, so `D_z2 = D_z1` and `phi_2 = phi_1`. **Both** hold. This is the model the paper
studies, and it is the model `PAPER_MODEL` names.

### 2.2 (2) Electrically perfect contact between two dielectric solids

Identical to 2.1 **as an electrical statement**. "Perfect contact" is a mechanical statement
(displacements and tractions continuous) and adds no new electrical equation. The two must not be
conflated: mechanical continuity does not imply electrical continuity, and vice versa.

### 2.3 (3) Ideal conducting electrode / short-circuited boundary

A conductor cannot support a tangential electric field, so the surface of a conductor is an
equipotential. It also carries whatever free charge is needed to terminate the normal displacement
field. Therefore:

- `phi` is **prescribed** on the conductor (`phi = 0` if it is grounded), and
- `D_z` **terminates**: it is not continuous across it, and the jump *is* the induced surface charge.
- Tangential `E` vanishes on the conductor surface, which is the physical meaning of equipotential.

Imposing `D_z` continuity on a conductor is **wrong**, because the whole point of a conductor is that
it can supply free charge.

### 2.4 (4) Insulated / open-circuit boundary

`sigma_free = 0` locally, so the normal `D` is continuous. But an insulated interface has an
additional global property: the **total** charge on an isolated conductor is its specified value
(usually zero), which is one scalar constraint on the field over the whole surface. At a single
planar interface under a single plane wave this constraint is not independent of the local condition,
which is why the numerical comparison cannot distinguish the two at fixed `k` (see §4.3 of
`PHASE_02A_TASK2T_NORMAL_INCIDENCE_LIMIT.md`).

**The brief's warning is exactly right here: continuity of `D_z` does *not* uniquely identify an
insulated interface.** `D_z` continuity says only that no free charge is stored. An idealized
dielectric contact (2.1) and an open circuit (2.4) satisfy the same local condition and differ only
in the global charge constraint and in whether the potential is free to float.

### 2.5 (5) Interface carrying a specified surface charge

`D_z2 - D_z1 = sigma_free(x)`. The potential remains continuous, **because the tangential condition
is unaffected by `sigma_free`**. A surface charge changes the normal flux but not the tangential
circulation. Imposing both continuity of potential and a charge jump is the correct and
*non*-redundant combination.

---

## 3. Why the two rows both vanish at `k = 0`, and why removing them is legitimate

This is the crux of the task, and the Task 2S reasoning was wrong here.

### 3.1 Tangential condition at `k = 0`

The condition is `k [phi] = 0`. At `k = 0` — normal incidence — this is

```
0 . [phi] = 0        for ANY finite [phi]
```

It is **identically satisfied**, regardless of the value of the potential jump. This is a
**degeneration**, not a gauge condition. The equation has lost the ability to carry information about
`[phi]`; it has not declared `[phi]` meaningless.

The distinction matters. A gauge freedom means *different potentials describe the same physics*. A
degenerate equation means *this equation says nothing here, but the quantity it would have
constrained is still physical*. `[phi]` at `k = 0` is still physical: a potential jump across the
interface is a surface dipole layer, and a dipole layer is a real charge distribution. It is
measurable, it stores energy, and it is not removable by a gauge transformation.

### 3.2 Normal condition at `k = 0`

Gauss's law in each bulk, `div D = 0`, with `exp(i k x)` dependence, reads

```
i k D_x + dD_z/dz = 0
```

At `k = 0` this reduces to `dD_z/dz = 0`, i.e. `D_z` is constant through each layer, and with
radiation/decay conditions it gives `D_z = 0` on *both* sides. Hence

```
D_z2 - D_z1 = 0 - 0 = 0
```

automatically, for `sigma_free = 0`. The normal condition is also **identically satisfied** at
`k = 0`.

### 3.3 The consequence

Both electrical rows are satisfied *whatever* the field configuration, at `k = 0`. Attempting to
impose them does not over-determine the problem in a physical sense — it over-determines the
*linear system*, because it adds two rows that every candidate solution satisfies equally.

**Removing them at `k = 0` is therefore legitimate, and for a reason the brief explicitly allows:
the rows are identically satisfied there.** The brief forbids removing a physical condition "solely
to get a consistent/square matrix." That is not what is happening. The rows are empty equations at
`k = 0`, and empty equations belong in the system only if one wishes to state them.

This is a **different and stronger justification** than the one Task 2S gave. Task 2S argued that
the rows were gauge/unimposable, which is false. The correct argument is that they are vacuous at
this one point in `k` and carry information at every other point.

### 3.4 What this does *not* license

It does not license removing the potential row at `k > 0`. At any `k != 0` the condition `k[phi] = 0`
is equivalent to `[phi] = 0` and does real work. `EQ-02A-MAXW-03` pins the `k = 0` behaviour; the
`k > 0` system keeps both rows.

---

## 4. Which model the paper studies, and the assumptions attached to it

**The paper studies the dielectric contact (model 1), with `sigma_free = 0`.**

The justification is not that this makes the matrix square. It is that the configuration is two
elastic solids bonded at a planar interface with no electrode present: there is no conductor to
supply free charge, and no external circuit to impose a charge. The interface therefore carries no
free charge, which by §1.2 forces both `phi` and `D_z` continuity.

**Assumptions this record depends on, each of which could fail in a different experiment:**

1. **No electrode at the interface.** If the two solids are bonded with a metallic interlayer, the
   correct model is (3) or a thin-conductor limit of it. The paper's configuration has no such layer.
2. **No trapped interfacial charge.** A real bonded interface can trap charge from processing or from
   prior electrical history. Modelled by (5) with `sigma_free != 0`. The paper assumes zero.
3. **Electroquasistatic reduction valid.** `(|v|/c)^2 ~ 1e-10` and `div D = 0` in the bulk. Both hold
   at 1 MHz.
4. **Interface is sharp.** No interphase layer. Finite-interphase modelling is explicitly out of
   scope for Task 2T.
5. **No free carrier transport across the interface.** In a real dielectric contact at high field
   there can be hopping or trap-assisted conduction. Assumed absent.

All five are recorded in `BLK-02A-23`, which stays open precisely so that adopting a different
electrical model requires an explicit, re-derived decision.

---

## 5. The numerical check that the derivation survives contact with the code

The derivation predicts that for a dielectric contact, **both** conditions hold. That is testable
against the production states without any reference to the reduction.

At `k/k_ac = 1e-3`, solving the full eight-row system with an adaptive basis:

| quantity | residual |
|---|---|
| potential jump | `9.237e-38` |
| normal `D` jump | `2.032e-42` |
| linear-solve residual | `2.976e-29` |
| condition number | `3.969e9` |

Both conditions hold simultaneously. **This alone refutes the Task 2S claim that imposing potential
continuity forces the transmitted potential to zero**: the system is perfectly solvable with both
rows imposed, once the basis carries all admitted modes.

### 5.1 A real bug found and fixed while doing this

The first version of the comparison fed the **scattered** field to the electrical residual function
instead of the **total** field. Every jump was therefore shifted by the incident amplitude, and the
`phi` row appeared to carry a residual of `0.03234` that was entirely an artefact.

The fix has three parts: `asymptotic.dielectric_contact_system` now returns `f1_scattered` and the
total `f1` separately; `maxwell.electrical_residuals` documents the total-field contract at the top
of its docstring; and `EQ-02A-MAXW-01` fails loudly if given a scattered field.

A second, related defect: the basis built from unit-flux modes **dropped** the near-gauge mode,
because that mode carries no flux and `sqrt(|jz|) -> 0`. A seven-column basis for eight admitted
modes produced a `2.121e41` condition number and a residual with the same `0.03234` magnitude as the
first bug — which is why the two were initially confused. `hybrid_mode` scales by flux when there is
flux and by the potential otherwise, and the count of columns against admitted modes is now asserted
in the test.

---

## 6. Verdicts for the P0/P1 questions

| question from the brief | verdict |
|---|---|
| Does `D_z` continuity uniquely identify an insulated interface? | **No.** It identifies the no-free-charge condition. Insulation adds a global constraint. |
| Does potential continuity replace `D_z` continuity? | **No.** They come from `curl E = 0` and `div D = 0` respectively and are independent. |
| Which model does the paper intend? | The **dielectric contact**, no free charge. Stated, assumed, and recorded with five attached assumptions. |
| Was the model chosen to make the matrix square? | **No.** The model was chosen from the physics; the matrix is square at `k > 0` with both rows and the `k = 0` rows are vacuous, not merely inconvenient. |
| Are both electrical rows automatically satisfied at `k = 0`? | **Yes.** Tangential: `k[phi] = 0` degenerates. Normal: `D_z = 0` on both sides by Gauss. |
| Is removing them legitimate? | **Yes**, because they are identically satisfied — the reason the brief permits. |
| Is a relative potential offset a gauge freedom? | **No.** It is a physical surface dipole layer. |

---

## 7. Evidence

| artefact | content |
|---|---|
| `code/scattering/maxwell.py` | the five models as executable `MODEL_SPEC`; `PAPER_MODEL = DIELECTRIC_CONTACT`; `electrical_residuals`; `distinguishes` |
| `validation/phase_02a_task2t_boundary_comparison.csv` | 70 rows: 7 material pairs x 2 polarisations x 5 models, with `phi_jump`, `D_jump`, `tangential_E`, solve residual, `sigma_free` |
| `validation/tests/test_17_scattering_task2t.py` | `EQ-02A-MAXW-01/02/03` — 8 tests, all PASS |
| `validation/mutation_tests/run_mutations.py` | `MR12`–`MR17`, six negative controls, all DETECTED; suite now 40/40 |
| `reports/PHASE_02A_TASK2T_CLAIM_CORRECTIONS.md` | the four claims audited in P4 |
| `reports/PHASE_02A_TASK2T_NORMAL_INCIDENCE_LIMIT.md` | the `k -> 0` limit analysis, P2/P3 |
