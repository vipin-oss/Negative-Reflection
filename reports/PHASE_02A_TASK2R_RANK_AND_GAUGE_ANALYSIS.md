# PHASE 2A — TASK 2R: RANK AND GAUGE ANALYSIS

**Companion to `PHASE_02A_TASK2R_NORMAL_INCIDENCE_AUDIT.md`.**

This report answers the question Task 2 got wrong: when the assembled system has
fewer independent rows than unknowns, what has actually happened?

The short answer: **four distinguishable things can produce that symptom, and
none of them is "the physics is non-unique".** This report separates them and
shows which one occurred.

---

## 1. The four possible causes

| # | cause | what it looks like | how to tell |
|---|---|---|---|
| 1 | **true gauge freedom** | a direction with zero physical content (no traction, no flux) | the state vector is zero except for a potential; it couples to no equation |
| 2 | **algebraic overparameterisation** | two columns proportional | normalised column overlap `\|aᴴb\|/(\|a\|\|b\|) = 1` |
| 3 | **degenerate modal basis** | two roots of the dispersion relation coincide, giving one mode twice | the modal `(ξ, state)` pairs are identical |
| 4 | **actual physical non-uniqueness** | two different *physical fields* solve the same boundary-value problem | only after 1–3 are excluded, and only if distinct physical fields survive |

Task 2 reported the symptom — `rows < unknowns`, `cond = 1.10e21`, rescaling
`DEFERRED_TO_CONDITIONING` — and filed it as cause 4. It had not tested 1, 2
or 3. In fact the cause was **none of the four**: an interface row carrying real
information had been deleted, so the row count was wrong rather than the physics.

---

## 2. What was actually removed at normal incidence

At Al6061, `k = 0`, the reduction removed exactly two rows, `φ` and `D_z`:

```
kept     ux, uz, theta, txz, tzz, qz
removed  phi, Dz
```

Justification for each removed row, established by measurement rather than by
argument:

| row | zero for every admitted mode? | zero on the incident RHS? | verdict |
|---|---|---|---|
| `phi` | yes — 0.0 on all six admitted modes | yes — Al6061 has `e31 = e33 = e15 = 0`, and the incident P/SV wave carries no potential | legitimately removable |
| `Dz` | yes — 0.0 on all six | yes | legitimately removable |

Both conditions now have to hold. Before the repair only the first was tested,
and `rhs_is_zero` was computed into the record without influencing the decision.

The critical consequence: the `u_x` and `t_xz` rows **are not in that table**,
because they are not zero. They carry `1.0` and `5.31e13j` on the shear branches.
Under the old logic they were removed anyway, because the incident *P* wave has
`u_x = 0` — a property of one wave generalised to the whole basis.

### Cross-check on a piezoelectric medium

The same reduction on PZT-5H, where `φ` is **not** zero:

| configuration | n_eq | n_unk | removed |
|---|---|---|---|
| Al \| Al | 6 | 6 | `phi`, `Dz` |
| PZT \| PZT | 7 | 6 | `Dz` |
| Al \| PZT | 7 | 6 | `Dz` |
| PZT \| Al | 7 | 6 | `Dz` |

The `φ` row is **kept** wherever a piezoelectric half-space is present and
removed only for the non-piezoelectric pair. The reduction tracks the material,
not the angle. That is the signature of a correct rule.

---

## 3. The three candidate causes, tested and excluded

### 3.1 True gauge freedom (cause 1) — present, and correctly identified

Branches 3 and 4 have `ξ = 0` and state `(0, 0, 1, 0, 0, 0, 0, 0)`: a pure
electric potential, zero displacement, zero traction, zero flux. This is genuine
gauge — it couples to nothing and carries no energy. It is excluded by two
independent gates (`SPATIAL_CLASS == GAUGE` in `Side.admissible`, and role
`GAUGE_EXCLUDED`).

This is real, and it is **not** what made Al|Al look under-determined.

### 3.2 Algebraic overparameterisation (cause 2) — absent

`SCTR04` scans every column pair of the reduced matrix on the Cauchy–Schwarz
equality case:

```
worst normalised overlap = 1.2e-1   (over all pairs, Al|Al 0°, Al|Al 30°, PZT|PZT 0°)
```

Nothing is near 1. No duplicate columns.

### 3.3 Degenerate basis (cause 3) — absent

The eight roots at `k = 0` are distinct: `0, ±979.77, ±2005.57, ±(216714 ±
216712j)`. The two `ξ = 0` roots share a value but are the *same* gauge mode,
counted twice by the polynomial, and both are excluded. No physical mode is
duplicated.

### 3.4 Physical non-uniqueness (cause 4) — absent

After the repair, Al|Al at normal incidence is **square and full rank**:

| incidence | n_eq | n_unk | rank |
|---|---|---|---|
| P | 6 | 6 | 6 |
| SV | 6 | 6 | 6 |

There is nothing left to be non-unique about. The solution is unique, the
reflected field is zero to `1e-41`, the transmitted field equals the incident
field to `0.0` in every component.

---

## 4. Why `rows < unknowns` was never evidence

The reasoning that produced BLK-02A-14 was:

> the system has fewer independent equations than unknowns, therefore the
> physical solution is not determined.

That inference is invalid in general and was invalid here. Two independent
reasons:

**A row can be dropped and the system stay determined.** Removing a genuinely
zero equation removes no constraint. If 7 rows constrain 7 unknowns and the 8th
row reads `0 = 0`, the solution is unique *provided the 7 remaining rows are
independent*. Counting rows before checking which were removed, and why, is not
a determination test.

**A gauge direction is not resolvable because it has no physical content.**
`φ = const` with everything else zero is invisible to every boundary condition
and to the energy flux. Its coefficient is arbitrary by construction. That makes
*the coefficient* non-unique while leaving *the physical field* unique. This
distinction is the one Task 2 collapsed, and it is the distinction that matters
for any subsequent energy audit: an energy computation depends on the physical
field, not on a gauge coefficient.

The correct test, applied here, is: exhibit the physical field, show it is
unique, and classify each deficient direction by what it physically is.

---

## 5. P3 — numerical solver validity

### 5.1 Equation counts and rank

Reported per configuration in `validation/phase_02a_task2r_independent_comparison.csv`.
The exact count, the numerical rank, and the nullity are all carried, along with
the tolerance used, so a rank that depends on its tolerance is visible.

### 5.2 Tolerance sensitivity

`matrix_rank_nullity` returns the tolerance with the rank. The tolerance is
relative to the largest singular value of the **row-equilibrated** matrix. It is
set at `10^-(dps-10)`, i.e. ten orders above the arithmetic floor at 40 digits.

The rank is not sensitive in the range tested: at Al|Al `k = 0` the smallest
surviving singular value is separated from the largest by more than 20 orders,
so any tolerance between `1e-30` and `1e-12` returns rank 6.

### 5.3 Consistency residual and whether the RHS is in range

For every configuration the RHS is checked to lie in the range of the operator:
the reduced matrix is square with full rank, so the system is solvable for any
RHS, and the residual `‖B a − rhs‖` is reported. Measured `0.0` for Al|Al at
`k = 0` in both polarizations.

**This is deliberately not treated as proof of a correct solution.** See §6.

### 5.4 Coefficients versus physical fields

Two distinct uniqueness questions, answered separately:

* **Coefficient uniqueness.** With no gauge columns admitted and no proportional
  columns, the coefficient vector is unique. Verified by the full-rank check.
* **Physical-state uniqueness.** Verified independently by reconstructing the
  eight-component fields on both sides and checking continuity. This is the
  answer that matters, and it is the one that is reported.

No configuration was found where the physical field is unique but the modal
coefficients are not, once the gauge directions are excluded.

### 5.5 Precision convergence

Roots were verified at 30, 40, 50 and 60 decimal digits. `xi_roots` matches
`w/cT` to 0.0 or ≤2.7e-48 at every level; the null-vector residual stays at the
arithmetic floor and the row-0 roundoff tracks `10^-dps` as it should. No
quantity was found that converges to a nonzero limit — which is what a genuine
numerical artifact would look like.

---

## 6. A converged solve is not a correct solve

Task 2's own history already contains this failure once: the assembly was built
`[B1 +B2]` while the docstring claimed `[B1 −B2]`, and the linear-algebra
residual was `1e-37` while the interface field was wrong by exactly twice the
incident state. It was caught only by recomputing the interface field
independently.

Task 2R adds a second instance of the same class of error: the matrix residual
for the wrongly-reduced Al|Al system was also essentially zero, because the
system it solved was square and consistent — it was just the wrong system. A
residual measured against a matrix cannot detect a missing equation. That is why
`SCTR07` re-derives the boundary conditions from the modal states and refuses a
reduction whose removed rows are not zero in the incident data as well as the
basis.

---

## 7. Declarations

* **Energy conservation is not implemented and not claimed.** G4 and G5 remain
  OPEN. No `R + T + A = 1` statement appears anywhere in this repair.
* **Cross terms are preserved.** The interface and total-field outputs retain
  incident–reflected interference. Absorption is **not** computed as
  `1 − R − T`.
* **Task 2 normal-incidence numbers are INVALID / SUPERSEDED.** Task 2 oblique
  numbers are RETAINED and unaffected; the defect only fires at exactly `k = 0`.
* **BLK-02A-14 is CLOSED** as an artifact of the truncation, not a property of
  the physics.
