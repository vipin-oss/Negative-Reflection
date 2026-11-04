# Phase 01B — Response to Supervisor

**To:** Supervisor
**From:** Arena (Phase 1B agent)
**Date:** 4 November 2026
**Branch:** `arena/a3ede4fe-negative-reflection`
**Subject:** Response to the Phase 1B mandatory mathematical repair instructions

---

## 1. Executive summary

**All three P0 items were substantiated.** They were not cosmetic, and two of them
invalidated Phase-1 deliverables outright. All three are now corrected, and each
correction is backed by an executed test whose reference is external to the code
being tested.

| Item | Your finding | Verdict | Current state |
|---|---|---|---|
| P0-1 parity | `S = diag(1,1,1,-1)` inconsistent | **CONFIRMED** | corrected to `diag(1,-1,-1,1)`; verified 16/16 entrywise |
| P0-2 energy flux | `+½Re[q_z]` not a valid cycle-averaged term | **CONFIRMED** | withdrawn; bilinear form derived and verified |
| P0-3 test credibility | zero-case PASS, unverified claims | **CONFIRMED — 7 of 7 requirements violated** | accounting rebuilt; fake claim replaced by a real test |
| P1 causality | speed-ordering rule unsupported | **CONFIRMED** | rule withdrawn; a parameter restored |

Two things I want to state plainly rather than bury:

**First, the parity defect was compounded by a false claim.** The Phase-1 document
stated that check C1 "evaluates all 16 entries ... symbolically and obtains zero."
It did no such thing. That is a worse failure than the wrong matrix, because the
wrong matrix was discoverable by exactly the test that was claimed to have been run.

**Second, the energy-flux defect was 34 orders of magnitude, not a factor of two.**
On the thermal branch the Phase-1 thermal term was 1.99 × 10³¹ W m⁻² against a
correct value of −1.30 × 10³ W m⁻². Had that formula reached a scattering
calculation it would not have produced a slightly wrong reflectance; it would have
produced a meaningless one.

**Gates G4 and G5 are NOT PASSED** and are now marked as such.

---

## 2. Item-by-item response

### 2.1 P0-1 — Parities

Your proposed `S = diag(1,-1,-1,1)` is **correct**. It satisfies
`s_i s_j = +1` on the eight ξ-even entries and `−1` on the eight ξ-odd entries.

The Phase-1 `S = diag(1,1,1,-1)` fails on **six** entries — `(1,2) (1,3) (1,4)
(2,1) (3,1) (4,1)` — precisely row 1 and column 1. The structural reason is that
row 1 carries ξ-odd entries in columns 2, 3 *and* 4, so a pattern of the form
`s₁s_j` can make all three negative only if `s₂ = s₃ = s₄ = −s₁`, which
`diag(1,1,1,-1)` does not do.

An exhaustive search over all 2⁴ diagonal sign patterns returns exactly
`(1,-1,-1,1)` and `(-1,1,1,-1)` — one symmetry up to overall sign.

The determinant parity `det M(-ξ) = det M(ξ)` follows from `det(S)² = 1`. **This is
why the error survived:** `det(S)² = 1` holds for *any* sign pattern, so a
determinant-only test cannot distinguish the two candidates. Both pass. Phase 1's
check was structurally incapable of detecting the defect it was nominally guarding
against.

Retested: roots, branches and residuals unchanged (C2 residual 2.2e−27; I1 external
closed forms 1.65e−15 over 24 cases). The correction changes a proof, not a result.

### 2.2 P0-2 — Energy identity

Re-derived from the same Ψ as `theory/governing_equations.tex`, against your
eight-point checklist. Result:

```
<J_z> = -½ Re[ σ_xz v_x* + σ_zz v_z* + iω D_z φ* ]  +  ½ Re[ q_z Θ* ] / T0
```

Checklist disposition:

| # | Item | Finding |
|---|---|---|
| 1 | Kinetic factors | unchanged — Phase 1 was correct here |
| 2 | Electric energy and potential-jump terms | unchanged in the bulk; the interface form needs the jump counted, already done in `interface_conditions.tex` |
| 3 | Heat-flux and entropy-production **sign** | **corrected** — sign now fixed by `⟨D⟩ ≥ 0`, not by residual size |
| 4 | τ_q-dependent stored energy | **corrected** — the Phase-1 "relaxation term" is identically **zero** for a monochromatic field and cannot be storage or dissipation |
| 5 | Correct quadratic thermal contribution | **corrected** — `q_z θ/T0`, not `q_z` |
| 6 | Bulk and interface positivity | bulk proved and tested; interface still carries the ½-bookkeeping item U-5 |
| 7 | All factors of ½ | **corrected** — every term now carries exactly one ½; the Phase-1 claim of a difference of *origin* was a misdiagnosis |
| 8 | Dimensions of every term | verified to W m⁻² for all three terms |

**One thing I got wrong while doing this**, recorded because it bears on method: my
first corrected attempt appeared to fail by a factor matching a missing term. Rather
than adjust a sign until it vanished, I varied τ_q over six decades — the residual
did not move, which excluded relaxation effects — and then tested the transverse
conduction term. It closed to 7.1e−17. The missing piece was `q_x ∂_xθ/T0`, a
geometric term of relative size (k/ξ)² = 1.2e−7, **not** a sign error. Had I tuned
a sign to make the number go away I would have destroyed the physics.

### 2.3 P0-3 — Test credibility

Audited against your seven requirements: **all seven were violated** in Phase 1.
The table is in §3.1 of `PHASE_01B_CRITICAL_REPAIR.md`. Fixes:

- `n_exec` is now **mandatory**; a missing count is a **FAIL**, not a default pass.
- `n_exec == 0` with `PASS` is automatically demoted to **UNRESOLVED**.
- Four categories — PASS / FAIL / UNRESOLVED / WITHDRAWN — all counted separately,
  every non-PASS check named in the summary. Nothing can drop out of the total.
- Exceptions are FAILs. There is no skip path.
- The script exits **1** while anything is unresolved.
- The false C1 claim is replaced by a real symbolic test on an independent
  construction (`symbolic_M()`), not the numeric tables and not the cached
  determinant.
- D1's reference route no longer shares an assumption with the implementation: the
  cycle average is built from separately evaluated **real** fields.
- The decisive flux test (D2) uses a dissipation functional sharing **no amplitude**
  with the mechanical or electrical flux terms.

**The guard proved itself during this work.** When `n_exec` was added to `record()`
but not yet to the individual checks, twelve previously-green checks immediately went
red with `FAIL / n_exec=0`. That is the mechanism working.

### 2.4 P1 — Causality

Your reassessment is correct and I withdraw the rule.

`v_T = √(κ/(ρc_Eτ_q))` is the speed of a **constitutive relaxation front** in a
truncated Taylor expansion. It is not a signal speed. Lord–Shulman theory contains
no mechanism bounding it by `c_L`, and `c_L` is not even the relevant bound — the
relativistic one is `c`, which `v_T ~ 10³ m/s` sits nine orders of magnitude below.
Rejecting sourced parameters on this basis was unsupported.

Retained as **universal** (enforced): ρ>0, c_E>0, T_0>0, C>0, ε>0, κ>0, τ_q≥0.
The last is genuinely required: Cattaneo must relax *toward* Fourier, not away.
Reported only: v_T ≤ c_L, and the frequency ceiling derived from it.

**Consequence — a parameter is restored.** Phase 1 had changed CdSe τ_q from
1.0e−13 s to 5.0e−13 s *purely* to clear the withdrawn rule. It is now back at
1.0e−13 s, with κ corrected to 6 W/(m·K) (midpoint of the sourced 5–10 bulk range).
CdSe now correctly reports v_T/c_L = 1.76 and is correctly **accepted**. The
material table no longer contains a value tuned to satisfy a rule that does not exist.

### 2.5 P1 — Coupled scattering preparation

Your caution was right and Phase 1's proposed remedy was wrong. Established two ways:

**Numerically.** The 8×8 coupled interface operator for Al6061|Al6061 at 30° has
rank **4**, not 8, *even at full coupling*, with smallest singular values
3.4e−20 and 9.4e−38. As λ→0 the outgoing mode set cannot be completed at all.

**Exactly** (check J1, residuals 0):

```
e = p = 0      :  det M = (ε33ξ² + ε11k²) · det M₃      rank 4 → 3
                  the φ factor is a GAUGE quadratic, not a wave branch
+b = 0         :  det M₃ = det(u_x,u_z) · M_θ            rank 3 → 2
isotropic      :  det(u_x,u_z) = C11 C44 (ξ²−ξ_P²)(ξ²−ξ_S²)
```

So the decoupled limit is singular because the formulation **retains two
zero-energy variables and their 0 = 0 conditions**. The remedy is **reduction**,
not λ→0 continuation. Your instruction is followed exactly.

**Built and validated:** `code/elastic_psv.py`, a reduced elastic P–SV solver
independent of the piezothermoelastic path:

| | |
|---|---|
| Identical media, 6 angles | R = 0 exactly |
| Normal incidence vs Fresnel | **1.67e−16** |
| Oblique, 7 angles, energy closure | ≤ 1e−15 |
| Spring → welded as K→∞ | monotone |

**Is the coupled solver ready? No.** The blocking item is now precisely localised:
the null-vector **direction** is good only to ~1e−5 at condition number 10²², and
because the balance identity is *scale-invariant* in the null vector, that residual
cannot be a magnitude error. This is U-3, unchanged in cause, now much better
characterised. No interphase parameter optimisation was started.

---

## 3. Current test status

```
foundation_checks: 15 PASS, 0 FAIL, 1 UNRESOLVED, 0 WITHDRAWN, 16 total
executed test cases across all checks: 151
  UNRESOLVED: D2  local balance d<J_z>/dz = -<D> (thermal branch)
```

Exit status **1**. Machine-readable: `validation/phase_01b_test_results.json`.

| Check | Status | n_exec | What it establishes |
|---|---|---|---|
| A1 | PASS | 4 | material tensor symmetry, stability |
| A2 | PASS | 26 | SI dimensions of every `M` entry |
| B1 | PASS | 6 | constitutive set from one potential |
| B2 | PASS | 4 | LS equation from the entropy balance |
| **C1** | **PASS** | **16** | **parity, entry-by-entry, both candidate S** |
| C2 | PASS | 8 | root branches, eigen residual 2.2e−27 |
| **D1** | **PASS** | **12** | **corrected flux: 3 external references** |
| **D2** | **UNRESOLVED** | **2** | **local balance — blocked by U-3** |
| E1 | PASS | 10 | flux z-conservation, corrected functional |
| F1 | PASS | 3 | interface passivity |
| G1 | PASS | 6 | Fresnel, Prop. 1(a)(b)(c), quarter-wave |
| **H1** | **PASS** | **4** | **universal only; speed rule withdrawn** |
| I1 | PASS | 24 | roots vs external closed forms, 1.65e−15 |
| I2 | PASS | 3 | interface rank 8 |
| **J1** | **PASS** | **5** | **reduction chain, all exact** |
| **J2** | **PASS** | **18** | **reduced elastic solver vs Fresnel, 1e−16** |

New checks this phase: **C1 (rewritten), D2, J1, J2**. Withdrawn or replaced:
the Phase-1 D1 and E1 flux paths, and the Phase-1 H1 pass condition.

---

## 4. Gate status — explicit

| Gate | Status | Reason |
|---|---|---|
| G0 prior art | **NOT PASSED** | full text still not accessed (U-1) |
| G1 constants | **NOT PASSED** | PZT-5H / BaTiO₃ β, p still ASSUMED (U-4); τ_q(CdSe) order-of-magnitude |
| G2 characteristic matrix | **PASSED** | parity now correctly proved (C1); roots validated externally (I1) |
| G3 validation | **PARTIAL** | reduced elastic reference done (J2); coupled routes unavailable |
| **G4 LS energy identity** | **NOT PASSED** | flux now correctly derived and verified, but the coupled balance does not close through the production code path (D2 UNRESOLVED) |
| **G5 interface dissipation + numerical balance** | **NOT PASSED** | requires the coupled solver; U-5 bookkeeping also open |

I have not upgraded G4 on the strength of the independent dps = 60 evaluation. The
identity closes there, but the gate asks for it to hold in the solver that will be
used, and it does not yet.

---

## 5. Open items

| ID | Item | Blocks | Status |
|---|---|---|---|
| U-1 | Gate G0: full-text prior art | any novelty claim | unchanged |
| U-3 | Amplitude normalisation of the coupled solve | G4, G5, all R/T/A | **root cause now localised** (direction, not magnitude) |
| U-4 | Assumed β, p for PZT-5H / BaTiO₃ | quantitative thermal-channel results | unchanged |
| U-5 | ½-bookkeeping between mechanical and thermal interface forms | E2 | **re-characterised** as convention, not physics |
| U-6 | `h → 0` reduction only to first order | interphase analysis | unchanged |
| U-7 | Thermal flux at realistic planes argued, not demonstrated | evaluation-plane choice | unchanged |
| U-8 | Interphase configuration not analysed | H2 | unchanged |
| U-9 | Bulk absorption not separated from reflection in practice | A vs R | unchanged |
| **U-11** | **Process: audit every Phase-1 "verified in check X" claim** | confidence in Phase-1 prose | **NEW** |
| **U-12** | **Case counts are not a uniform unit across checks** | precision of the 151 figure | **NEW** |

**On U-11.** P0-1 was not primarily an algebra error; it was a document asserting a
verification that had not been performed. I have re-checked the claims in the
material I rewrote, and I have checked the ones the corrections touch. I have **not**
systematically re-audited every "checked", "verified" and "tested" statement in all
Phase-1 prose. That audit is the natural first task of any further review, and I
recommend it before Phase 2 begins.

---

## 6. What I did not do

- **No Phase 2 work.** No coupled scattering coefficient, no `R`, `T` or `A`.
- **No manuscript results.**
- **No interphase parameter optimisation** — explicitly instructed, and not started.
- **No tolerance loosening to convert an UNRESOLVED into a PASS.** D2 is reported
  unresolved with its number attached.
- **No repair by sign-tuning.** Where a correction did not close, the cause was
  investigated (τ_q sweep, transverse-term test) rather than absorbed into a sign.
- **No `.tex` compilation claim.** No TeX engine is available; all five documents
  are lint-checked as text only (environment nesting, macro definitions, math-mode
  discipline, dangling cross-references — all clean).

---

## 7. Recommendation

The mathematical foundation is now sound where it has been repaired: parity is
proved correctly, the energy functional is correct and verified against the Second
Law, and the decoupled limit is understood as a reduction with an exact
elastic reference built to match.

The blocker for Phase 2 is **U-3 and nothing else**. It is not a modelling question
— it is the accuracy of a null-vector direction at condition number 10²². The
concrete next step, if you approve Phase 2, is to attack that numerically
(structured scaling of the state unknowns by their physical magnitudes, then a
rank-revealing decomposition rather than a free-column Cramer solve), with the
newly validated reduced elastic solver as the acceptance reference.

I do not recommend proceeding to any reflection-suppression claim until G4 and G5
are passed. That is the position this repair supports.

**Stopped after Phase 1B.**
