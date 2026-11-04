# PHASE 02A — TASK 2T
# The Normal-Incidence Limit: Is the Field Finite, and Does the Reduction Reach It?

**Branch** `arena/a3ede4fe-negative-reflection`
**Starting commit** `006d059164d7f8910d432932fa70b7411e821ad9`

---

## 0. The two questions

1. **P2.** Is the reported `1/k` growth of the potential a physical divergence, or an artefact of
   how a modal eigenvector is normalised?
2. **P3.** Does the production `k = 0` reduction actually reach the `k -> 0` limit of the full
   dielectric-contact system — measured independently, not by asserting agreement at `k > 0`?

Both are answered by measurement. Neither can be answered by inspecting one eigenvector.

---

## 1. Why a single eigenvector cannot answer P2

An eigenvector is defined up to a scalar. If `v` solves the modal problem then so does `alpha v` for
any nonzero complex `alpha`, and every statement about "the potential of the mode" changes with
`alpha` while the physics does not. A claim that the potential "diverges like `1/k`" is therefore a
claim about a *particular* `alpha`, and is meaningless until that `alpha` is named.

The production code names it in `roots._normalise`, which divides every eigenvector by its **largest
mechanical component** — in practice pinning `u_x = 1`. That is a perfectly reasonable choice for
root-finding. It is not a physical statement, and it must not be read as one.

---

## 2. P2: the `1/k` is an artefact — measured three ways

Take PZT-5H, branch 5 (the mode that carries the large potential), and track the same mode direction
as `k` decreases from `1e-1` to `1e-6` of `k_ac`. Three independent normalisations:

| normalisation | behaviour of `phi` over six decades | verdict |
|---|---|---|
| `\|u_z\| = 1` | spread `2.78e-5` — **bounded** | finite |
| displacement magnitude `sqrt(\|u_x\|^2 + \|u_z\|^2) = 1` | spread `2.33e-3` — **bounded** | finite |
| `\|u_x\| = 1` (production) | grows by a factor of **`1.002e5`** | divergent |

The divergence exists in exactly one normalisation, and that normalisation pins a component which
itself vanishes like `k`. For a quasi-longitudinal mode `u_x/u_z` is proportional to `k`; dividing by
`u_x` therefore multiplies the whole vector by `1/k`. The `1/k` in the potential is the `1/k` in
`u_z/u_x`, seen from the other side of the fraction.

**The mechanism is confirmed by the ratio itself:** `u_x/u_z = 1.046e-6 * (k/k_ac)`, exactly linear,
measured over the same range. The normalisation is not merely suspicious; the quantity it divides by
is `k`-proportional by construction.

### 2.1 An honest limitation of the artefact claim

The hybrid basis (flux where there is flux, potential otherwise) does **not** remove the ratio
`phi/u_x` — it grows by `1e5` there too. This is unsurprising once stated: rescaling multiplies every
component by the same number, so a *ratio* of two components is invariant. Two intermediate versions
of the check compared hybrid magnitudes and then hybrid ratios and both reported `1.0x` growth, i.e.
they measured nothing.

The correct statement is narrower and is the one recorded in `CLM-02A-21`:

> The modal ratio `phi/u_x` is an intrinsic property of the mode and does diverge. What is finite is
> the **total field**, because the coefficient multiplying this mode vanishes at the same rate. The
> artefact is in reading a component ratio off a single basis vector and calling it the field.

The three normalisations establish the artefact in the sense that matters: a physical quantity
(the potential of the mode up to scale) does not have a normalisation-independent `1/k`; it is the
*choice of pinned component* that produces it. The field-level question is settled in §3.

---

## 3. The total field has a finite limit

This is the normalisation-invariant statement the brief demands. In a **unit-z-flux basis**
(`|j_z| = 1` per mode, which is `k`-independent for a propagating mode), every configuration was
examined at `k/k_ac` from `1e-2` to `1e-5` and each component classified by its measured log-log
order before any extrapolation.

**No component in any configuration diverges.** Across 7 material pairs × {P, SV}, every component
is either `CONSTANT` (finite nonzero limit) or `VANISHING` (limit zero, approached like a power of
`k`). `k_order` reports a slope, so this is a measurement and not a classification asserted by
default — and it had to be verified that the classifier *can* detect divergence, see §3.2.

Selected limit values (unit-z-flux basis, transmitted side):

| pair | pol | `phi` limit | `t_zz` limit |
|---|---|---|---|
| Al6061\|Al6061 | P | `0.0` | `5884.7184` |
| Al6061\|Al6061 | SV | `0.0` | `0.44906` |
| PZT-5H\|PZT-5H | P | `0.068572202` | `8307.9396` |
| PZT-5H\|PZT-5H | SV | `1.6379065e-6` | `0.30511` |
| Al6061\|PZT-5H | P | `0.064646274` | `7837.2555` |
| Al6061\|PZT-5H | SV | `0.073347938` | `0.64801` |
| PZT-5H\|Al6061 | P | `0.091267105` | `5551.3213` |
| PZT-5H\|Al6061 | SV | `0.091693763` | `0.30163` |

**Verdict on P2: the physical divergence claim is dead.** The total field is finite in every
configuration tested.

### 3.1 A taxonomy that the earlier "k-orderings" wording obscured

Not all vanishing components vanish for the same reason, and conflating them produced a real error
in this task. Two categories must be reported separately:

**(a) Indeterminate 0/0 with a nonzero limit.** `t_xz` at P incidence. Side 1 and side 2 each carry
`t_xz = C·k` with the same `C`, so their difference is `O(k^2)` — a cancellation, not a vanishing.
Measured on Al6061|Al6061 P: `t_xz = 2.8668193·k` exactly, constant to seven digits over
`k/k_ac` from `1e-3` to `1e-5`. Production returns `0`, and `0` is the correct limit.

**(b) Pure `O(k)`, true zero limit.** `t_xz` at SV incidence; `u_z` and `theta` at P incidence.
Measured: `t_zz` sequence `5.5858264, 0.55858265, 0.055858265` at `k/k_ac = 1e-3, 1e-4, 1e-5` —
exactly a factor of ten per decade.

The error this caused: an unclassified quadratic (Richardson) extrapolation applied to a category-(b)
component reports a spurious residual of `5e-4` on a quantity whose limit is genuinely zero. The rule
adopted is **classify the order first, extrapolate only `CONSTANT` components**.

### 3.2 The order classifier had a real bug, caught by its own positive control

`k_order` initially tested `p < 0.05` first and returned `CONSTANT`. A clean `1/k` sequence has
`p = -1.0`, which is also less than `0.05`, so **a divergence would have been classified as
constant** and "no divergence found" would have been vacuous.

Caught by adding a positive control that feeds the classifier a synthetic `1/k` and a synthetic
constant sequence and requires the right verdict for each. The fix tests the divergence branch first.
Verified:

```
1/k    -> order -1.0 -> DIVERGENT
const  -> order  0.0 -> CONSTANT
k      -> order  1.0 -> VANISHING
k^2    -> order  2.0 -> VANISHING
```

Both positive controls now live inside `EQ-02A-LIM-02` and are what mutation `MR17` protects.

---

## 4. P3: independent validation of the reduction

### 4.1 What "independent" means here

Agreement between two systems implementing the same boundary equations at `k > 0` proves nothing
about `k = 0`. The check performed is therefore:

1. Solve the **full eight-row dielectric-contact system** at a sequence of small `k`, in a
   `k`-independent basis, with 60-digit arithmetic.
2. Classify each component's order in `k` and extrapolate only the `CONSTANT` ones.
3. Compare against the **production `k = 0` reduced system**, which removes the two vacuous rows.
4. Require agreement on every `CONSTANT` component, to within the limit's *own* measured convergence.

### 4.2 Results

| pair | pol | worst `CONSTANT` disagreement (relative to field scale) | limit drift |
|---|---|---|---|
| Al6061\|Al6061 | P | `1.6e-7` | ~`1e-9` |
| Al6061\|Al6061 | SV | `4.8e-9` | ~`1e-10` |
| PZT-5H\|PZT-5H | P (`phi`) | `1.2e-9` | ~`1e-11` |
| PZT-5H\|PZT-5H | SV | `1.9e-9` | ~`1e-10` |
| Al6061\|PZT-5H | P (`phi`) | `6.342e-4` | `6.3e-6` |
| PZT-5H\|Al6061 | P | `1.64e-5` | `1.6e-7` |
| Al6061\|PZT-5H | SV | `9.46e-5` | `9.4e-7` |
| PZT-5H\|Al6061 | SV | `1.48e-4` | `1.5e-6` |

All seven pairs and both polarisations agree within 100× the limit's own convergence, which is the
tolerance the test asserts and justifies. **The reduction reaches the limit.**

### 4.3 A finding about the five models

Solving all five models at `k/k_ac = 1e-3` in the same basis
(`validation/phase_02a_task2t_boundary_comparison.csv`, 70 rows):

- `DIELECTRIC_CONTACT` and `PRESCRIBED_CHARGE` (with `sigma = 1e-9`) satisfy their conditions to
  `1e-50`, with `phi_jump` and `tangential_E` at `1e-56`—`1e-63`.
- `CONDUCTING_ELECTRODE`, `SHORT_CIRCUIT` and `OPEN_CIRCUIT` satisfy their own conditions to `1e-50`
  but show `phi_jump` from `0.119` up to `4.958e4` of the field scale.
- **The three "free" models are numerically identical** at fixed `k`.

The third point is itself a finding: at a single planar interface under a single plane wave, open and
short circuit differ only in a global charge constraint that the local system cannot see. This is
why `D_z` continuity cannot be used to identify an insulated interface (§2.4 of the boundary audit),
and why the open-circuit model is listed as requiring an explicit global constraint rather than a
local row.

### 4.4 The reduction does *not* enforce potential continuity exactly

The production `k = 0` solve was checked for continuity of the quantities it does not impose:

| configuration | `phi_jump / field_scale` | `D_z_jump / field_scale` |
|---|---|---|
| Al6061\|PZT-5H P | `8.25e-6` | `≤ 3.4e-71` |
| PZT-5H\|Al6061 P | `1.1e-5` | `≤ 3.4e-71` |
| all others | `0.0` | `≤ 3.4e-71` |

The normal-`D` condition is satisfied to machine noise, because at `k = 0` it is a cancellation
identity (`D_z = 0` on both sides) and any decent solve reproduces a zero as a zero. That is **not**
evidence of physical content — it is exactly the vacuity described in §3.2 of the boundary audit.

The potential jump is **not** exactly zero where it matters, because the row is not imposed. It is
approached at first order in `k`: `8e-6` to `1e-5` of the field scale. This is a genuine, measured
limitation of the reduction, recorded as `BLK-02A-24` and pinned by `EQ-02A-LIM-04`.

**It does not refute the reduction**, which is accurate to `O(k)` by construction, but it does mean
the reduction must never be described as *exactly* equivalent to the full dielectric contact. The
`k = 0` answer is the limit; it is not the same operator.

---

## 5. The structure of the `k -> 0` singularity

The system is mathematically non-singular at `k = 0` but numerically singular in its `k > 0`
representation, and the two must not be confused.

| quantity | `k > 0` | `k = 0` |
|---|---|---|
| admitted branches per side | 4 | 3 |
| matrix rank | 8 | 6 |
| columns | 8 | 6 |

The branch that disappears **is** the `xi = 0` gauge vector: it degenerates into it, with
`xi ∝ k` and state `(0, 0, 1, 0, ...) + O(k)`. Consequently its `D_z` flux vanishes like `k`, and
its coefficient must grow like `1/k` to produce a finite contribution to the total field.

Measured: `cond ∝ k^-2.0`, reaching `3.949e16` at `k/k_ac = 1e-7`; `xi ∝ k^1.0`. Both slopes are
recovered to within `0.1` of the integer.

Two practical consequences, both recorded as `BLK-02A-25`:

1. **The `k = 0` answer cannot be obtained by evaluating the `k > 0` system at small `k`.** The
   individual coefficients are unbounded even though their sum is finite. It must be extrapolated.
2. **Any basis used for the limit must carry this mode.** A unit-flux basis silently drops it
   (`sqrt(|jz|) -> 0`), and the resulting seven-column system has condition number `2.121e41` and a
   `0.03234` residual concentrated on the `phi` row. That residual is *not* the `1/k` problem; it is
   the near-gauge mode being missing, and it is what made potential continuity look impossible.

---

## 6. P3 scope: what was not done

- **Finite-interphase models** are explicitly out of scope for Task 2T and were not implemented.
- **Oblique incidence** was not re-analysed for the limit; recorded as `BLK-02A-22`. The boundary
  derivation applies at any `k`, but the numerical limit analysis is normal-incidence only.
- **Only two piezoelectric materials** (PZT-5H and BaTiO3, the latter with assumed `beta`/`p`
  constants) and CdSe (assumed `kappa` constants) were exercised. Admissible source-backed constants
  for a second independent piezoelectric are not in the workspace.

---

## 7. Verdicts

| question | verdict |
|---|---|
| Is the `1/k` potential divergence physical? | **No.** Artefact of pinning `u_x = 1`, a `k`-proportional component. |
| Is the total field finite at `k = 0`? | **Yes**, in every component of every configuration tested, under a `k`-independent basis. |
| Does the production reduction reach the limit? | **Yes**, for `CONSTANT` components, within 100× the limit's own convergence. |
| Does it do so exactly? | **No.** Potential continuity is approached at `O(k)`; worst `1.1e-5` of field scale. |
| Is the `k -> 0` problem numerically singular? | **Yes**, `cond ∝ k^-2`, so the limit is extrapolated, not evaluated. |

---

## 8. Evidence

| artefact | content |
|---|---|
| `code/scattering/asymptotic.py` | `hybrid_mode`, `unit_flux_mode`, `unit_displacement_mode`, `dielectric_contact_system(rows=, normalisation=)`, `k_order` |
| `validation/tests/test_17_scattering_task2t.py` | `EQ-02A-LIM-01/02/03/04`, `EQ-02A-COND-T1` — all PASS |
| `validation/mutation_tests/run_mutations.py` | `MR16` (rescaling treated as divergence), `MR17` (single-sample limit) — both DETECTED |
| `validation/phase_02a_task2t_boundary_comparison.csv` | 70 rows, five models, seven pairs, two polarisations |
| `validation/phase_02a_task2t_results.json` | machine-readable form of every number in this report |
| `reports/PHASE_02A_TASK2T_MAXWELL_BOUNDARY_AUDIT.md` | the derivation this rests on |
