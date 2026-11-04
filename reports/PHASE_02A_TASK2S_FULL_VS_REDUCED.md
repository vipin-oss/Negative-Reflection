# Phase 2A Task 2S — Full versus reduced formulation, independent benchmark

Companion to `reports/PHASE_02A_TASK2S_ELECTRICAL_GAUGE_AUDIT.md`.
Machine-readable data: `validation/phase_02a_task2s_gauge_comparison.csv`,
`validation/phase_02a_task2s_results.json`.

**Gates G0, G1, G4, G5 remain OPEN. No `R + T + A = 1` is claimed.**

---

## 1. What is being compared

Three formulations, all at normal incidence, all built from **raw physical**
states:

**(a) Full system with explicit gauge variables.**
All eight interface rows. The `xi = 0` branch of each half-space is carried as an
explicit unknown `C_1`, `C_2`. One global gauge condition (`C_1 = 0`) is supplied.
Two constants, one invariance, one fixing row.

**(b) Properly derived gauge-invariant reduced system.**
The eight rows with the `xi = 0` sector projected out, so no constants appear and
no fixing row is needed. This is the reference the reduction is measured against.

**(c) The existing Task 2R production reduction.**
`code/scattering/assembly.py` `assemble()` at `k = 0`, which removes the rows it
classifies as gauge or empty and solves the remainder.

The benchmark is (a) versus (c) on physical observables, with (b) establishing
what "the physical answer" is.

---

## 2. Why the comparison is run at `k > 0` and then at the limit

The full formulation **cannot be evaluated at exactly `k = 0`**. The potential
entry of a piezoelectric mode grows like `1/k`:

| `k/k_ac` | `phi`, PZT branch 5 |
|---|---|
| `1e-2` | `2.63e11` |
| `1e-4` | `2.63e13` |
| `1e-6` | `2.63e15` |
| `1e-8` | `2.63e17` |

The additive constants absorb this divergence, so the system still *solves* at
`k = 0` — and produces an answer that is **not** the limit. Measured,
Al6061|PZT-5H P-incidence, transmitted/incident `t_zz`:

```
k = 0   full     1.3316249147 + 0.00115753873 j     <- 4th digit wrong, spurious imaginary part
k = 1e-4         1.3317978949 + 0.00000483... j
k = 1e-5         1.33179789532 + 0.000004839... j
k = 0   reduced  1.33179789533 + 0.00000484... j     <- the limit
```

The reduced answer is the limit; the full answer at `k = 0` is not. `full_system`
therefore now rejects `k = 0` with an error explaining why, and the benchmark is:

1. compare (a) and (c) at `k/k_ac = 1e-2` and `1e-4`, where both are well posed;
2. confirm separately that (c) at `k = 0` is the limit of (c) at small `k`.

Step 2 is what makes the `k = 0` reduction meaningful; step 1 is what shows it
agrees with a properly gauge-fixed formulation.

---

## 3. Method

For each configuration:

* build both systems independently;
* solve each by row- and column-equilibrated least squares;
* reconstruct the eight interface fields on the transmitted side;
* form the **transfer ratio** `field_component / incident_component`, choosing the
  incident component of largest magnitude as the reference.

The transfer ratio is the only physically meaningful comparison: the polynomial
mode normalisation is `k`-dependent, so the raw field magnitudes are not
comparable between formulations or between wavenumbers. Reported differences are
relative to `max(|full|, |reduced|, 1)`.

**`phi` is excluded from the agreement metric, by design.** The two formulations
differ from `phi` by a constant on each side, and the brief explicitly allows
potential values differing by a common gauge constant. The test that this really
*is* a constant is that all seven other components agree while the `phi` offset is
the same on every mode — verified.

---

## 4. Result

**Al6061|Al6061, PZT-5H|PZT-5H, Al6061|PZT-5H, PZT-5H|Al6061 × {P, SV} ×
`k/k_ac` ∈ {1e-2, 1e-4} — 16 configurations, 112 gauge-invariant comparisons.**

```
worst relative disagreement between full and reduced:   1.707e-48
at:  PZT-5H|Al6061, SV, k/k_ac = 1e-4, component q_z
```

Representative transfers (real part; the imaginary parts agree to the same
precision):

| pair | pol | `k/k_ac` | full | reduced | rel. diff |
|---|---|---|---|---|---|
| Al6061\|Al6061 | P | 1e-2 | 1.0000000000000 | 1.0000000000000 | 3.1e-127 |
| PZT-5H\|PZT-5H | SV | 1e-2 | 1.0000000000000 | 1.0000000000000 | 0.0 |
| Al6061\|PZT-5H | P | 1e-2 | 1.3317935982779 + 4.275e-6j | 1.3317935982779 + 4.275e-6j | 1.2e-60 |
| Al6061\|PZT-5H | SV | 1e-4 | 1.4430177914333 − 2.865e-6j | 1.4430177914333 − 2.865e-6j | 1.9e-52 |
| PZT-5H\|Al6061 | P | 1e-4 | 0.66819471015765 + 2.524e-6j | 0.66819471015765 + 2.524e-6j | 1.6e-61 |
| PZT-5H\|Al6061 | SV | 1e-4 | 0.98859117888161 − 1.154e-5j | 0.98859117888161 − 1.154e-5j | 4.6e-53 |

Identical media return a transfer of exactly `1.0` at every wavenumber, as they
must — no reflection, no asymmetry. That is a check on the whole pipeline
independent of the gauge question.

Ranks:

| | full | reduced |
|---|---|---|
| unknowns | 8 (6 modes + 2 constants) | 6 (6 modes) |
| rows | 9 (8 conditions + 1 gauge fixing) | 6 at `k = 0`, 8 at `k > 0` |
| rank | 8 in every configuration | full in every configuration |

---

## 5. The `k = 0` limit test

The reduction removes rows *because* `k = 0` is a limit. So the reduction must
actually reproduce the limit. Tested by extrapolation in `k`, an oracle that does
not use the `k = 0` row-removal rule at all:

**Al6061|PZT-5H, P-incidence** — relative change in the transfer ratio against the
`k = 0` reduced answer:

| `k/k_ac` | deviation |
|---|---|
| 1e-2 | 3.254e-6 |
| 1e-3 | 5.412e-8 |
| 1e-4 | 4.359e-9 |

Decade ratios 60 and 12 — **first order, converging on the `k = 0` answer.** The
same first-order behaviour (ratio ≈ 10 per decade) holds for every pair and both
polarisations where the deviation exceeds the roundoff floor (`EQ-02A-GAUGE-A5`).
Identical media are exact at every `k`, as they must be.

This is the substantive justification for the row removal. At `k = 0`:

* the **`Dz` row** is Gauss's law `D_z = 0`, which is the **insulated** interface
  condition — a different physical model from the bonded contact the rest of the
  system imposes. Keeping it would silently switch models;
* the **`phi` row** is a condition on a quantity with no finite `k = 0` value.

Both are removed, for two separately justified reasons. Neither is "discarded
because the system was over-determined".

---

## 6. Consistency with the earlier `k = 0` measurements

An intermediate version of this benchmark compared the two formulations at exactly
`k = 0` and reported disagreements of order `1.0` in `ux`, `theta`, `tzz` and `qz`
for Al|PZT and PZT|Al. Those numbers are **withdrawn**: they were produced by the
defective full system described in §2, which had no valid `k = 0` evaluation.

The corrected picture is:

* at `k > 0`, full and reduced agree to `1.7e-48`;
* at `k = 0`, the reduced answer is correct and the full answer does not exist as
  an evaluation — only as a limit.

Reporting the earlier disagreement as "the reduction is broken" would have been
wrong, and reporting the agreement without the limit test would have been
unfounded. Both are kept on the record.

---

## 7. Reproduce

```
cd validation
python3 -m pytest tests/test_16_scattering_task2s.py -q      # 7 cases, all PASS
python3 mutation_tests/run_mutations.py                      # 34/34 DETECTED
python3 scattering/assembly_validation.py                    # 120/120
```

The comparison table is regenerated by the procedure in §3; the committed copy is
`validation/phase_02a_task2s_gauge_comparison.csv` (16 rows, one per
configuration).

**G4 and G5 remain OPEN. No energy-conservation result is claimed. No
`R + T + A = 1`.**
