# Phase 1D — Foundation Repair

**Author:** Vipin Gupta, ORCID 0000-0002-5238-1963
**Branch:** `arena/a3ede4fe-negative-reflection`
**Start commit:** `5bf01f10ba66c638a57f89cc73656c6538dd50eb`
**Delivered commit:** `30ea2cc60d25e2e488930287fba88d8c65afee3b`
**CI:** run 41825495880 — completed / success
**Date:** 4 November 2026

---

## 0. Verdict

**CONDITIONAL PASS.** One P0 blocker (**BLK-U3**) is genuinely closed, with per-criterion
evidence. Two further defects were found and fixed that no prior phase had detected.
The Lord–Shulman stored-energy question is resolved *mathematically* and deliberately
left **open physically**. **BLK-G4 and BLK-G5 remain NOT PASSED**, and the energy gate
was **not** closed on a special elastic limit.

Nothing in this phase is a promotion. The promotion state is unchanged: `G0`, `G1`,
`G4`, `G5` not passed; `G2` passed; `G3` partial; coupled solver not ready.

---

## 1. What was actually wrong — and what was *not*

The single most important result of Phase 1D is that **the recorded diagnosis of the
P0 blocker was wrong**, and the way it was wrong is instructive.

### 1.1 The recorded claim (Phases 1B and 1C)

> "The null-vector DIRECTION for the 6 mm characteristic matrix (cond ~1e22) is
> accurate only to ~1e-5. Because the energy balance identity is scale-invariant in
> the null vector, this residual cannot be a magnitude error."

### 1.2 What was actually true

Two independent defects existed in **one function**, `flux_z_mp`, and neither had
anything to do with the eigenvector:

| # | Defect | Size of effect | How it was found |
|---|---|---|---|
| **U3-A** | `D_z` was assembled with `eps11 * E_z` instead of `eps33 * E_z` | **1.8555e-5** — exactly the reported residual | `D_z` is `e31·e_xx + e33·e_zz + eps_zz·E_z`; for a z-poled medium `eps_zz = eps33`. Every other `D_z` in the file (lines 139, 317, 373, 1240, 1554) already used `eps33`; line 525 alone did not |
| **U3-B** | the transverse **dissipation** channel `<D_x>` was omitted | **4.67e-8** = exactly `(k/ξ)²` at 30° | measured 4.663e-8 against the predicted 4.672e-8 |

The eigenvector direction was measured directly and is accurate to **better than
1e-40** (test `U3C`), i.e. some 35 orders of magnitude more accurate than the
residual it was being blamed for.

### 1.3 Why the wrong diagnosis survived two phases

Three structural reasons, each now closed by a test:

1. **The residual was used as a proxy for the direction.** A test that measures a
   residual cannot measure a direction. `U3C` now measures the direction against an
   independent construction; `U3A`/`U3B`/`U3D` measure the residual. They are
   different measurements and can no longer be substituted for one another.
2. **The residual normalisation hid the defect.** `eigen_residual_mp` divides by the
   *global max entry*. With a 23-decade spread in column norms, that denominator
   makes any directional error look like 1e-26. Componentwise per-row residuals —
   `|Σ_j A_ij a_j| / Σ_j |A_ij||a_j|` — are the only honest measure, and are what
   `U3D` uses.
3. **Every energy test called the same function it was testing.** `check_D2_balance`
   compared `flux_z_mp` against `_dissipation_terms_mp`; both inherited the `eps11`
   error, so the comparison was self-consistent and passed for the *wrong* reason.
   A residual test is only meaningful when pinned to an independently *derived*
   expression, not merely an independently *computed* number.

---

## 2. BLK-U3 — closure against the stated criteria

The brief requires **all** of the following before BLK-U3 may close. Each is listed
with the measurement that satisfies it.

| Criterion | Measurement | Verdict |
|---|---|---|
| Production reproduces an independent formulation within justified tolerance | production balance closes at **3.45e-42** (PZT-5H), **3.12e-63** (CdSe) against a tolerance of **1e-8** | ✅ |
| Stable under increased precision | residual **shrinks** with precision: 3.44e-62 (dps 60) → 3.06e-81 → 5.51e-102 → **8.38e-142** (dps 140) | ✅ |
| Physical boundary-state quantities agree | the balance is stated in **traction** `σ_xz`, **electric displacement** `D_z` and **heat flux** `q_z`; the two-route flux comparison agrees to 1.95e-16 | ✅ |
| Degenerate / near-degenerate configurations handled explicitly | all 8 branches carry their own null vector and are each measured in `U3C`; 3/4 are an exact conjugate pair, 1↔6 and 2↔5 are conjugate pairs, none are mixed | ✅ |
| Executed regression tests reproduce it | `U3A`, `U3B`, `U3C`, `U3D` all PASS; mutations `M12`–`M15` all DETECTED | ✅ |

### 2.1 Two genuinely different formulations

Independence was established by construction, not by assertion:

* **(a) Production** — column-equilibrated Gaussian elimination with a 3×3 Cramer
  solve for the free component (`kernel_scaled`).
* **(b) Reference** — two-sided (row *and* column) equilibration, then the eigenvector
  of the equilibrated Gram matrix `CᴴC` at **dps 120**, with the root Newton-refined
  at **dps 160** from the exact coefficient polynomial.

These share no code path beyond `M_mp`. **Result: they agree to 1e-55 … 1e-58** after
a common physical normalisation, across 92 coupled configurations
(`validation/phase_01d_independent_comparison.csv`).

> **Note on what was *not* selected.** My first attempt at (b) used the eigenvector of
> `CᴴC` at the *same* dps=50 root that production uses. That agreed to 1e-25 and I
> nearly reported it as confirmation. It was not independent in the relevant sense:
> the root itself carries only ~50 digits, so both sides inherited the same error.
> The reference must refine the root at higher precision, or it is measuring the same
> thing twice.

### 2.2 Numerical conditioning (measured, not estimated)

| Quantity | Value |
|---|---|
| Column-norm spread of `M` | **23 decades** (3.08e23, 1.55e24, 3.08e14, 3.63e12) |
| Raw double-precision condition number | `inf` (singular in float) |
| Column-equilibrated singular values | `[1.73, 1.00, 2.87e-3, 4.06e-19]` |
| Smallest elimination pivot at a *true* root | **1.92e-125** |
| Perturbation between production and reference direction | **6.4e-58** |

The 1.92e-125 pivot confirms the matrix is numerically rank-deficient at working
precision even though its exact rank is 4 — which is precisely why cofactor/Cramer
methods must carry high precision and why a bare `np.linalg.svd` on the float matrix
(what `Medium6mm.mode_vector` does) is meaningless. That routine is **not** on the
production path; `kernel_scaled` is.

> **A correction I had to make to my own work.** I initially reported a
> "76 %–141 %" directional disagreement between production and an exact Cramer null
> vector. That number was an artefact of my own probe: I compared vectors normalised
> by two *different* conventions. Under a common physical normalisation the
> disagreement is 1e-58. The superseded figure is retained in
> `reports/PHASE_01D_PRE_REPAIR_INPUT_CHECKSUMS.json` and in this note so the error
> is on record rather than quietly dropped.

---

## 3. The Lord–Shulman derivation — and the one thing it does *not* settle

`theory/energy_flux_derivation.tex` gained a full section (`sec:ls-closure`)
covering all ten items the brief requires: potential and conjugate variables,
reversible energy, the relaxation contribution, the local balance, entropy
production and positivity, the exact harmonic average, transverse conduction,
interface dissipation, factor-of-two conventions, and boundary terms for
discontinuous `φ` or `θ`.

### 3.1 The resolution, and why it is precise

The question carried open since Phase 1 was whether the Lord–Shulman relaxation
stores energy. **Proposition (prop:relax) settles the mathematics:**

```
⟨q·q̇⟩  = 0                     identically, for every monochromatic field
⟨q·q⟩  = ½|Q̂|²  ≠ 0           always
```

*Proof:* write `q = ½(Q̂e^{iωt} + Q̂*e^{−iωt})`. Multiplying by `q̇` and averaging
annihilates every `e^{±2iωt}` term, leaving `⟨q·q̇⟩ = (iω/4)(Q̂·Q̂* − Q̂*·Q̂) = 0`
because the two terms are complex conjugates. The second average is
`¼(Q̂·Q̂* + Q̂*·Q̂) = ½|Q̂|²`.

So the *extended-irreversible-thermodynamics* stored-energy density
`τ_q q·q / (2κT₀)` is **genuinely non-zero**, and an earlier reading of this project —
which concluded there was no stored energy because `τ_q q·q̇/κ` averages to zero —
was an inference from a vanishing harmonic average. The brief forbids exactly that
inference, and it was being made.

### 3.2 Why this does **not** close BLK-LS-STORAGE

Two distinct objects must be kept apart, and the note now states this as
Remark `rem:twofluxes`:

* the **physical first-order flux** `J = σu̇ + φḊ − (θ/T₀)q`, in which `q` enters
  *linearly*;
* the **cycle-averaged effective wave-energy flux**
  `⟨J_z⟩ = −½Re[σ_xz v_x* + σ_zz v_z* + iωD_zφ*] + ½Re[q_zθ*]/T₀`, in which the
  thermal part is *bilinear*.

The relaxation affects the second **not at all**, and the first only at instant level.
It *does* affect the energy **density**, and therefore any definition of absorbed
energy based on a stored-energy change. Which functional an **experiment** measures
requires a complete scattering solution. **BLK-LS-STORAGE is therefore retained
OPEN** with status `MATHEMATICAL_PART_CLOSED_PHYSICAL_PART_OPEN`.

### 3.3 What the derivation establishes, and what it refuses to

Established: internal consistency of the account; positivity of bulk dissipation
(`⟨𝒟⟩ = ⟨q·∇θ⟩/T₀ ≥ 0` by positive-definite `κ`) and of the interface form
(`eq:Dint`, which vanishes identically for a welded, perfectly conducting bond);
the correct flux including both conduction directions; dimension checks.

**Refused, explicitly:** the derivation does **not** establish G4 or G5. The elastic
limit closing is a special case and cannot certify the coupled account. The note says
so in `§Status of gates G4 and G5`, and test `THY04` fails if either gate is ever
marked passing.

### 3.4 A transverse-conduction trap I fell into and had to undo

While deriving the flux I added a term `+½Re[q_xθ*]/T₀` to `⟨J_z⟩`, reasoning that a
mode at angle θ carries `∂θ/∂x = ikθ` too, and that omitting it "leaves exactly
`(k/ξ)²`". The arithmetic looked right and `check_D2_balance` did improve.

**It was wrong.** `⟨J_z⟩` is the *z*-component of the flux and therefore carries
`q_z` only. The transverse gradient contributes to the **dissipation**
(`q·∇θ/T₀` includes `q_x ∂θ/∂x`), not to `J_z`. A 2×2 table made this unambiguous:

| | `D = D_z` only | `D = D_z + D_x` |
|---|---|---|
| `J = mech + elec + q_z` | 4.67e-8 | **6.85e-72** ← correct pairing |
| `J = mech + elec + q_z + q_x` | 4.60e-8 | 7.11e-10 |

The correct pairing cancels to machine precision. Adding `q_x` to the flux *destroys*
an identity that otherwise closes exactly. The code now carries a comment recording
this as a **deliberate non-change**, and test `U3B` fails if anyone "fixes" it back.

---

## 4. Stale documentation (BLK-DOC-01)

Test `THY02`, written for this phase, immediately found that **BLK-DOC-01 was still
live**: `theory/governing_equations.tex` line 328 still presented

```
τ_q ≥ κ/(ρ c_E c_L²)        [boxed, labelled \textsc{Derived}]
```

as an active result, with per-material minimum values — two phases after
`CLM-H1-02` had been marked `INVALID`. This is precisely the failure mode the brief
names, and it had survived because nothing compared the theory text against the claim
registry.

It is now replaced with an explicit retraction that states the withdrawn inequality,
explains why it does not follow, and carries no `\textsc{Derived}` label. `H1`'s check
label was likewise changed from `"thermodynamic stability + LS causality bounds"` to
`"thermodynamic stability (universal conditions only)"`; its *code* was already
correct and reported the ratio as `OPTIONAL (not enforced)`.

The active theory now separates four categories explicitly, as required:

* **(A) active scientific assumptions** — quasi-static electric, small strain,
  linearisation about `T₀`, full coupling retained;
* **(B) withdrawn claims** — the `τ_q` causality floor, the linear thermal flux;
* **(C) optional modelling restrictions** — `k = 0`, `p = 0`, `β = 0`, plane strain
  (flagged as diagnostics, never as piezothermoelastic results);
* **(D) unresolved physical questions** — the stored-energy *physical* question, the
  `ASSUMED` thermoelastic/pyroelectric constants, `τ_q` provenance.

`CLM-H1-02` remains in the registry as `INVALID` with its original statement intact —
the history is preserved, which `THY02` also checks.

---

## 5. Material-data provenance

### 5.1 The defect

`build_materials()` hard-coded every constant in Python while
`literature/material_constants.csv` — a versioned table with units, sources and
status codes for 80 values — was *documentation that nothing read*. Two independent
sources for the same numbers, free to diverge, with no guarantee a reported result
used the versioned one.

### 5.2 The fix

`load_material_table()` is now the authoritative path and **validates** rather than
merely parses:

* every required symbol must be present, or `ProvenanceError`;
* every value must carry a unit, a source, and a status from a controlled vocabulary;
* a missing critical constant **raises** — it is never defaulted, because a default is
  an invented constant;
* `ASSUMED` values are returned **flagged**, so a quantitative claim can refuse to
  rest on them;
* `build_materials()` loads the table, with the Phase-1C literals retained only as a
  fallback (`REQUIRE_TABLE` forces strictness).

### 5.3 What the validator found on its first run

It rejected the real table three times before accepting it, each rejection a genuine
labelling issue or a vocabulary gap:

1. `ZERO(centrosym)` (Al6061 — symmetry-forced zeros, *stronger* than VERIFIED) was
   not in the vocabulary;
2. `VERIFIED_TO_ORDER_OF_MAGNITUDE` was not in the vocabulary;
3. `p3`, `beta1`, `beta3` were absent from `REQUIRED_SYMBOLS` — so a table with the
   pyroelectric coefficient **deleted** validated cleanly. That hole was found by
   test `PRV03` and is now closed.

### 5.4 The provenance guard is now a tautology, and that is reported

`DAT01` (added in Phase 1C) compared the code against the table. That was a real check
*while the code hard-coded literals*. **After this change production loads the table,
so `DAT01` compares an object against the rows it was built from — it is now a
tautology and can no longer detect a corrupted table.**

Mutation `M06` (a 1000× unit slip in `C11`) demonstrated this by *surviving*. The
honest response was not to weaken `M06` but to add an independent check: `PRV06`
tests each tabulated value against **externally stated decade-wide physical bands**
(handbook ranges, not derived from the code). It catches the unit-slip class and
nothing finer — a limitation stated in the test rather than papered over.

### 5.5 No constant was invented

Every `ASSUMED` value in the table remains `ASSUMED` and is now flagged on the medium
object. `beta1`, `beta3`, `p` for PZT-5H and BaTiO₃, and `kappa11`, `kappa33` for CdSe,
are assumptions and are labelled as such. Their presence is a blocker on quantitative
claims about the materials that need them, exactly as the brief requires.

---

## 6. Reduced-scattering readiness

Task 6 asks for a *readiness* test, not measured coefficients. Four cases were added,
and **none of them reports a piezothermoelastic R/T/A.**

| Test | What it establishes |
|---|---|
| `SCT03` | the classical reduced elastic P–SV substrate reproduces Fresnel to 1e-15 at normal incidence and closes energy to <1e-12 at 0–60° |
| `SCT04` | all 8 branches classify into acoustic/gauge families under **one common physical normalisation** (`u_z = 1`), with no gauge mode counted as acoustic |
| `SCT05` | the incoming/outgoing classification is invariant under arbitrary complex rescaling of each eigenvector (40 rescalings, 0 changes) |
| `SCT06` | the electric gauge degeneracy of a centrosymmetric medium is **detected** rather than mistaken for an accuracy failure |

### 6.1 A correction to my own test, and the finding underneath it

`SCT04` first reported that **every** PZT-5H branch was "electrically dominated" with
mechanical content at the 1e-9…1e-13 level, which would have been a catastrophic
defect. It was **my test that was wrong**: I classified each branch by its components
normalised by *that vector's own largest entry*. That is a gauge, not a physical
scale — it differs per branch, and it makes volts (order 1e9) swamp metres (order 1).

Under one common physical normalisation (`u_z = 1`, applied identically to every
branch and every material) the components carry their real units and the picture is
sane:

| branch | `\|u_x\|` | `\|u_z\|` | `\|φ\|` (V) | `\|θ\|` (K) |
|---|---|---|---|---|
| P-fast | 0.60 | 1 | 1.77e9 | 1.58e5 |
| S-fast | 1.08 | 1 | 4.68e10 | 9.47e4 |
| thermal | 2.7e-4 | 1 | 1.50e9 | 5.52e11 |

`|φ| ~ 1e9` is **correct**: an acoustic wave in a piezoelectric carries kilovolts to
megavolts of potential. The production `kernel_scaled` reproduces this reference to
**1e-58** under that normalisation.

### 6.2 The genuine readiness finding: centrosymmetric gauge degeneracy

For Al6061 the piezoelectric coupling is exactly zero by symmetry (`e = 0`, `p = 0`),
so the electric-displacement equation decouples. At 30° the `φ` row of `M` has
magnitude **8.0e-87** against mechanical rows of order **1.5e17** — 104 orders down.
The potential is then a true gauge: *any* `φ` satisfies that row to working precision.

This matters for readiness, not for physics: a componentwise residual on those
branches measures **4.6e-3** — not because the mode is wrong, but because the
denominator is a numerically zero row. A coupled solver that treats that residual as
an accuracy metric would reject correct modes; one that inverts an 8×8 operator
without rank detection would produce garbage. `SCT06` therefore asserts the
**detection** of the degeneracy (a ~100-order separation between Al6061 and PZT-5H)
rather than a small residual.

The same mechanism accounts for the 20 `DEGENERATE_EXCLUDED` rows in the comparison
CSV: at normal incidence the `u_x` column collapses to 1e-104, because for a *pure*
`u_z` wave `∂u_x/∂x ≡ 0`, so `u_x` drops out of the constitutive relations entirely.
That is a correct exclusion, not an error.

---

## 7. Adversarial mutation results

**19 mutations, 19 detected, 0 survivors, 0 invalid detectors.**
(`validation/mutation_tests/last_mutation_run.json`)

Eight are new in Phase 1D, one per defect the brief names:

| ID | Injected defect | Detector | Result |
|---|---|---|---|
| `M12` | inaccurate eigenvector **direction** despite a tiny characteristic residual | `U3C` | DETECTED |
| `M13` | wrong displacement-vs-potential scaling | energy balance (residual 120.2) | DETECTED |
| `M14` | incorrectly removed gauge mode | root count `[4,4]` ≠ `[8,8]` | DETECTED |
| `M15` | missing Lord–Shulman relaxation operator | energy balance (residual 2.3e-6) | DETECTED |
| `M16` | invalid material-data provenance (literal drifts from table) | `PRV02` | DETECTED |
| `M17` | stale registry after a theory correction | frozen-input check names the edited file | DETECTED |
| `M18` | promotion with an open mandatory blocker / non-computing test | claim binding | DETECTED |
| `M19` | withdrawn `τ_q` causality floor restored to active theory | `THY02` | DETECTED |

### 7.1 Three survivors I had to fix — reported, not hidden

All three were **defects in the validation, not in the theory**, and each is now
closed with the reason recorded in the mutation file:

1. **`M03` (re-injected linear thermal term) SURVIVED.** The mutation patched
   `flux_z`, the *float* routine, while `check_D2_balance` calls `flux_z_mp`. The
   mutation was applied to code the detector never executes. Fixing the target
   exposed a second, deeper problem: `check_D2_balance` reported `UNRESOLVED` for
   every non-closure, on the (then-justified) theory that a residual might reflect
   an unknown numerical limit. With BLK-U3 closed that semantics is wrong — a
   residual of 1.00 is a **defect**, not an unknown. `FAIL` is now returned when a
   case executes and exceeds tolerance; `UNRESOLVED` only when nothing executed.
2. **`M06` (1000× unit slip) SURVIVED** — via the `DAT01` tautology described in
   §5.4. Closed by adding `PRV06`.
3. **`M18` (promotion with a non-computing test) SURVIVED.** `BIND01` trusted the
   evidence map entirely and never compared the map's `check` field against the
   registry's own `test_id`, so rewriting a claim's `test_id` was invisible.
   Cross-checking them immediately found a **pre-existing** ambiguity (`MAT` vs
   `MAT01`/`MAT02`), now disambiguated in the registry.

A separate collision I caused and fixed: my first provenance test reused the id
`DAT01`, shadowing the pre-existing cross-source check. Renumbered to `PRV01`–`PRV06`.

---

## 8. Evidence integrity during the repair

The frozen-input and manifest checks **fired correctly** on every legitimate change I
made, which is the system working:

```
EXP-FOUND-01B: 1/2 output hashes match after re-run
   HASH MISMATCH: the recorded output is not what the code produces now
2/3 manifests verified
```

Nothing was overwritten silently. The procedure followed the documented policy in
`evidence/approved_baselines.json`:

1. `reports/PHASE_01D_PRE_REPAIR_INPUT_CHECKSUMS.json` — every pre-repair checksum;
2. `reports/phase_01d_pre_repair_evidence/` — pre-repair manifests, raw outputs and
   stdout for all three experiments, verbatim;
3. `reports/PHASE_01D_MANIFEST_VERIFY_BEFORE.txt` — the verify output *before*
   re-issue, including the failure;
4. `BASELINE-01D-1` recorded as `CANDIDATE`, **not** `APPROVED`, with promotion still
   blocked;
5. `regen_manifests.py --rewrite`, which archives the superseded manifests as
   `manifest.superseded.<sha>.json` rather than deleting them.

**No tolerance was relaxed anywhere.** `check_D2_balance` still declares `1e-8`; the
identity now closes at 3.45e-42 because the physics was corrected. The registry entries
carry explicit retractions of the earlier misdiagnosis rather than silent edits.

---

## 9. Test execution

| Suite | Result |
|---|---|
| **L1 FAST** | **55 PASS, 0 FAIL, 1 UNRESOLVED (SCT02), 0 WITHDRAWN, 0 NOT_APPLICABLE — 56 total, 9449 executed cases**, exit 2 |
| L1 baseline before Phase 1D | 35 PASS, 0 FAIL, 1 UNRESOLVED — 36 total, 9292 cases |
| L1B foundation suite | **16 of 16 PASS** (D2 improved 1.86e-5 → 3.45e-42) |
| L2 STANDARD | exit 2, same single UNRESOLVED |
| Mutations | **19/19 detected** |
| Manifest verify | 3/3 after deliberate re-issue |
| CI gate | exit 0 — "no contradictions, every unknown registered and owned" |

`SCT02` remains `UNRESOLVED` with `n_exec = 0` because **no coupled solver exists**.
It is not counted as passing, and it cannot be: the harness converts a zero-case PASS
to `UNRESOLVED` by construction.

Registries after Phase 1D:

* claims — 45 total: **29 VERIFIED_STRONG**, 9 VERIFIED_WEAK, 5 UNVERIFIED,
  2 INVALID, **0 UNRESOLVED** (was 1);
* equations — 41 total: **20 VERIFIED_STRONG** (was 17), 6 WEAK, 7 UNRESOLVED,
  3 DEFINED_NOT_TESTED, 2 PARTIAL, 3 INVALID (unchanged).

Every changed claim has an executed test and a code binding (`BIND01`).

---

## 10. Blocker decisions

See `reports/PHASE_01D_BLOCKER_DECISIONS.md` for the per-criterion record. Summary:

| Blocker | Sev | Decision |
|---|---|---|
| **BLK-U3** | P0 | **CLOSED** — all five criteria met, per-criterion evidence recorded |
| BLK-LS-STORAGE | P1 | **OPEN** — mathematical half resolved; physical half needs a scattering solution. The question is deliberately *not* closed |
| BLK-G4 | P0 | **NOT PASSED** — unchanged. An elastic limit cannot certify the coupled account |
| BLK-G5 | P0 | **NOT PASSED** — unchanged. No complete scattering solution exists |
| BLK-G0 | P0 | **NOT PASSED** — unchanged |
| BLK-G1 | P0 | **NOT PASSED** — unchanged |
| BLK-U5 | P1 | OPEN — unchanged |
| BLK-DOC-01 | — | **CLOSED** — the withdrawn floor is out of the active theory (`THY02`, `M19`) |

---

## 11. What remains open, stated plainly

1. **No coupled scattering solution exists.** `SCT02` is UNRESOLVED with zero cases
   executed. No piezothermoelastic R/T/A may be quoted.
2. **G4 and G5 are not passed.** The energy gate has an internally consistent
   derivation and a production balance that closes to 3e-56 on propagating modes, but
   a *scattering-level* balance has never been evaluated. G5 stays open.
3. **The stored-energy question is open physically.** §3.1 settles the mathematics;
   which functional an experiment measures does not follow from it.
4. **Thermoelastic and pyroelectric constants for PZT-5H and BaTiO₃ are ASSUMED.**
   They block quantitative claims about those materials.
5. **The centrosymmetric gauge degeneracy is detected but not solved.** A coupled
   solver must handle it explicitly; `SCT06` currently only proves the detection works.
6. **Four branches at normal incidence are genuinely excluded** (`u_x` column at
   1e-104). Any solver must exclude rather than invert them.
7. **CI was inspected, after the push, on the actual commit.** Run **41825495880**
   on `30ea2cc` — `completed` / `success` (job "Level 1 (fast) validation: success").
   Recorded in `validation/phase_01d_test_results.json` under `ci`. The prior green
   run 41819120262 belongs to the Phase 1C commit and is superseded, not inherited.
   The CI *artifact* was not downloaded (the sandbox cannot reach the blob host), and
   no artifact content is claimed.

---

## 12. Reproduce

```bash
# baseline integrity
python3 validation/regen_manifests.py            # VERIFY mode, never writes
python3 validation/freeze_inputs.py --check 2>/dev/null || true

# foundation suite (16 checks)
python3 code/foundation_checks.py

# L1 / L2
cd validation && python3 run_fast.py             # 55/0/1, exit 2
cd validation && python3 run_standard.py

# mutation suite (19 mutations, each with an unmutated control)
python3 validation/mutation_tests/run_mutations.py

# evidence gate
cd validation && python3 ci_gate.py              # exit 0

# impact analysis
python3 evidence/impact_analysis.py --changed EQ-ENG-02
python3 evidence/impact_analysis.py --check
```

Component evidence:

* `validation/phase_01d_test_results.json` — executed counts, tallies, conditioning
* `validation/phase_01d_independent_comparison.csv` — 120 configurations, two
  formulations, per-class classification
* `reports/PHASE_01D_CHANGE_IMPACT.md` — dependency traversal
* `reports/PHASE_01D_BLOCKER_DECISIONS.md` — closure criteria

---

## 13. Recommendation

**CONDITIONAL PASS.**

The P0 blocker that has gated this programme since Phase 1 is closed, and it was
closed by finding the actual defect rather than by relaxing a criterion. Three real
defects were found in Phase 1D that every earlier phase had missed — the `eps11`
index, the transverse dissipation pairing, and the stale causality floor — and each
now has a regression test and a mutation that proves the test can fail.

The condition is unchanged and non-negotiable: **no energy-ratio claim, no
piezothermoelastic reflection coefficient, and no promotion of G4 or G5** until a
complete scattering solution exists and the material-provenance blockers on the
thermoelastic and pyroelectric constants are cleared.

The most valuable thing this phase produced is not the closed blocker. It is that the
recorded diagnosis of that blocker was confidently wrong for two phases, and the
reason it survived was structural: a residual was used as a proxy for a direction, a
global normalisation hid the error, and every test compared a function with itself.
All three are now closed by construction rather than by vigilance.
