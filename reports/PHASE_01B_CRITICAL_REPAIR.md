# Phase 01B — Critical Mathematical Repair

**Branch:** `arena/a3ede4fe-negative-reflection`
**Date:** 4 November 2026
**Scope:** corrections to Phase 1 only. No Phase 2 work, no manuscript results.
**Status:** the two P0 defects are **real and are now corrected**; one of them
(P1 causality) is corrected by **withdrawing** a Phase-1 rule rather than by
satisfying it. All three P0 items are marked below with their exact outcome.

---

## Contents

| Item | Type | Outcome |
|---|---|---|
| P0-1 | Characteristic-matrix parity | **CORRECTED** — supervisor's `diag(1,-1,-1,1)` is right; `diag(1,1,1,-1)` fails 6 of 16 entries |
| P0-2 | Lord–Shulman energy flux | **CORRECTED** — `½Re[q_z]` withdrawn; bilinear `½Re[q_zΘ*]/T0` derived and verified |
| P0-3 | Numerical-test credibility | **CORRECTED** — 6 distinct defects fixed; a check that executes 0 cases can no longer pass |
| P1-A | Causality / speed ordering | **WITHDRAWN RULE** — the τ_q floor is not a physical requirement |
| P1-B | Coupled scattering preparation | **PARTIAL** — reduced elastic reference built and validated; coupled solver still not ready |

---

# P0-1 — Characteristic matrix parity

## 1.1 Original equation

`theory/characteristic_matrix.tex`, Proposition "Structure of det M", Phase-1 proof:

> Consider the involution ξ → −ξ combined with the diagonal sign change
> **S** = diag(1,1,1,−1). Direct substitution shows
> **M**(−ξ) = **S** **M**(ξ) **S** ... because every entry that is odd in ξ sits
> in a row and a column with exactly one factor of **S**.

and, in the accompanying remark:

> check C1 evaluates all 16 entries of **M**(−ξ) − **S****M**(ξ)**S** symbolically
> and obtains zero

## 1.2 Why it was incorrect

**Two separate faults.**

**(a) The matrix is wrong.** `diag(1,1,1,−1)` satisfies the required entrywise
conditions on only 10 of 16 entries. The six it fails are

```
(1,2) (1,3) (1,4) (2,1) (3,1) (4,1)
```

i.e. exactly the first row and first column. The reason is structural: with
`s = (1,1,1,−1)`, `s₁s₄ = −1` (correct for the odd entries (1,4) and (4,1)) but
`s₁s₂ = +1` and `s₁s₃ = +1`, which is wrong because `M₁₂`, `M₁₃` are *odd* in ξ.
Row 1 has odd entries at columns 2, 3 **and** 4; no single sign pattern of the form
`s₁s_j` can make all three negative unless `s₂ = s₃ = s₄ = −s₁`.

**(b) The verification was never performed.** The Phase-1 document asserted the
16 entries had been evaluated symbolically. They had not. Check C1 in
`code/foundation_checks.py` only:

1. read the cached determinant coefficients, and
2. checked that the odd coefficients were zero.

It never constructed a similarity transformation. **A claim of verification that
was never carried out is the most serious class of defect in this repository**, and
it is the reason the wrong matrix survived. This is reported to the supervisor as a
process failure, not only as an algebra error.

## 1.3 Correct derivation

Classify all 16 entries by parity in ξ. From the displayed matrix entries:

| parity | entries | count |
|---|---|---|
| even | (1,1) (2,2) (3,3) (4,4) (1,4) (4,1) (2,3) (3,2) | 8 |
| odd | (1,2) (2,1) (1,3) (3,1) (2,4) (4,2) (3,4) (4,3) | 8 |

For **M**(−ξ) = **S** **M**(ξ) **S** with **S** = diag(s₁,…,s₄) to hold entry by
entry we need

```
s_i s_j = +1   for the 8 even entries
s_i s_j = −1   for the 8 odd entries
```

This is a 16-condition system in four signs. Solving it gives **exactly two
solutions**, **S** = diag(1,−1,−1,1) and its negative, which are the same symmetry
up to an overall sign. The determinant then follows:

```
det M(−ξ) = det(S) det M(ξ) det(S) = det(S)² det M(ξ) = det M(ξ)
```

since det(**S**) = ±1. So `det M` is even in ξ, as claimed — **but for the correct
reason and with the correct matrix**.

**Note on why the error was invisible.** det(**S**)² = 1 holds for *every* sign
pattern, so the determinant parity is even regardless of which **S** is written
down. A determinant-only test therefore cannot distinguish the two candidates;
both "pass". This is a textbook example of a check that is structurally incapable
of detecting the error it is nominally guarding against.

## 1.4 Modified file and location

| File | Change |
|---|---|
| `theory/characteristic_matrix.tex` | Proof of Proposition 1 replaced; `diag(1,1,1,−1)` → `diag(1,−1,−1,1)`; added the parity-class table, the 16-condition system, the sign-ambiguity note, and a remark stating that the Phase-1 claim was false |
| `code/foundation_checks.py` | `check_C1_structure` rewritten (line ~805); new helper `symbolic_M()` |
| `code/foundation_checks.py` | check C1 now reports `n_exec = 16` |

## 1.5 Executed verification

```
$ python3 code/foundation_checks.py
[PASS] C1   det M: degree, parity, block structure, gauge root
       <deg=8, odd=none, S_b failures=0>   n_exec=16
```

The check builds the 4×4 matrix with **symbolic** material constants (not the
numeric tables, not the cached determinant), then evaluates
**M**(−ξ) − **S** **M**(ξ) **S** entry by entry for **both** candidates:

```
  (i,j)   parity  S=diag(1,1,1,-1)   S=diag(1,-1,-1,1)
  (1,2)   odd     FAILS              MATCHES
  (1,3)   odd     FAILS              MATCHES
  (1,4)   even    FAILS              MATCHES
  (2,4)   odd     MATCHES            MATCHES
  ...     (16 rows in total)

  S_a = diag(1,1,1,-1) : 6 of 16 entries FAIL
  S_b = diag(1,-1,-1,1) : 0 of 16 entries FAIL
  exhaustive sign search over all 2^4 diagonal S: solutions = [(1,-1,-1,1), (-1,1,1,-1)]
```

Both results agree exactly with the independent algebra in §1.3.

## 1.6 Roots retested

The correction changes a **proof**, not a **result**. `det M` was already known to
be even from the exact expansion. Retained unchanged and re-executed:

```
[PASS] C1  non-zero powers of xi: [0, 2, 4, 6, 8];  odd powers present: none
[PASS] C2  8 roots: branch identification and eigen residual  <max residual 2.2e-27>
[PASS] I1  exact roots vs closed-form elastic branches (external) <24 cases, 1.65e-15>
```

The eight roots, their branch assignment, and the residuals are identical to
Phase 1. Assertions in the automated suite now test the *correct* **S**.

## 1.7 Remaining uncertainty

None on the algebra. The parity claim `det M = P₄(ξ²)` is established exactly and
now rests on a correct transformation. The only residual concern is procedural: the
Phase-1 claim of a symbolic entrywise check was untrue, and it is worth asking what
other Phase-1 statements of the form "verified in check X" were not, in fact,
verified. A full audit of every such claim is listed as open item **U-11** in
`reports/PHASE_01B_SUPERVISOR_RESPONSE.md`.

---

# P0-2 — Lord–Shulman energy identity

## 2.1 Original equation

`theory/energy_flux_derivation.tex`, Phase 1, Eq. (EF-11):

```
<P_z> = -½ Re[ σ_xz v_x* + σ_zz v_z* + iω D_z φ* ]  +  ½ Re[ q_z ]
```

implemented identically in `code/foundation_checks.py` at `flux_z`,
`flux_z_cycle` and `flux_z_mp`.

## 2.2 Why it was incorrect

The supervisor's diagnosis is right, and there are two independent grounds:

**(a) It is a linear functional of one first-order phasor.** A cycle average
`½Re(AB*)` is bilinear in two first-order fields. `½Re[q_z]` is linear in one. The
time average of a monochromatic first-order quantity is **zero**; a nonzero time
average requires a quadratic cross term. So `½Re[q_z]` is not something any
averaging procedure can produce — it is not a cycle-averaged quadratic flux at all.

**(b) It omits the temperature.** The heat flux enters the energy budget multiplied
by a temperature. Without a second field, the term is invariant under
θ → θ + const, which is not a symmetry of a conducting body.

**Numerical size of the error.** On the thermal branch of PZT-5H at 1 MHz:

| | Phase 1 `½Re[q_z]` | Corrected `½Re[q_zΘ*]/T0` |
|---|---|---|
| W m⁻² | 1.99 × 10³¹ | −1.30 × 10³ |

a discrepancy of **34 orders of magnitude**. On the propagating branches the
uncorrected term is O(10¹⁸) against mechanical terms of O(10²⁰). **Any result
computed with the Phase-1 formula is void.** None was reported.

## 2.3 Correct derivation

Carried out from the same potential Ψ used in `theory/governing_equations.tex`,
following the supervisor's eight-point checklist.

**Stored energy (points 1, 2, 4).** The Phase-1 stored-energy density was already
correct and is retained:

```
u = ½ρ u̇_i u̇_i + ½σ_ij ε_ij + ½E·D + ½(ρc_E/T0) θ²
```

with the `½` on each quadratic form.

**Lord–Shulman relaxation term (point 4).** The second term of the Phase-1
"dissipation" was `−τ_q q·∂_t q/κ`. **For a monochromatic field this term is
identically zero** — it is the product q·q̇ averaged over a period — so it is
neither a stored energy nor a dissipation. It must be dropped from both, which the
original note half-recognised by calling it "relaxation storage". Corrected text
now states this explicitly.

**Heat-flux and entropy-production sign (point 3).** With
`σ_s = q·∇(1/T) = −q·∇θ/T²`, the bulk dissipation is

```
<D> = −(1/T0)[ ½Re( q_z (∂_zθ)* ) + ½Re( q_x (∂_xθ)* ) ]
```

**Both conduction directions must appear.** This is the point on which a first
attempt at the check appeared to fail: omitting the transverse term leaves a
residual of exactly (k/ξ)² = 1.2 × 10⁻⁷ at 30° for PZT-5H. It took a systematic
investigation (varying τ_q over 6 decades left the residual unchanged, which
excluded relaxation effects and pointed to geometry) to identify it.

**Correct thermal contribution (point 5).** From the balance
`∂_t[u_thermal] + ∇·J = −D` with `u_thermal = ½(ρc_E/T0)θ²`, matching
`∂_t(½ρc_Eθ²/T0) = (ρc_E/T0)θθ̇` against the entropy balance
`T0∂_t(ρη) = −∇·q` identifies the thermal part of the total current as

```
J_thermal = q θ / T0
```

whose cycle average is the **bilinear form** `½Re[q_z Θ*]/T0`.

**The corrected flux:**

```
<J_z> = -½ Re[ σ_xz v_x* + σ_zz v_z* + iω D_z φ* ]  +  ½ Re[ q_z Θ* ] / T0
```

**Sign fixed by the Second Law, not by residual size (point 3, 6).**

```
−½Re[ q_z (iξΘ)* ]  =  ½ κ33 |ξΘ|² / (1 + ω²τ_q²)  ≥ 0
```

so the `+` sign in the corrected formula is the one that makes `<D> ≥ 0`. The
opposite choice makes `<D> ≤ 0`, i.e. a flux growing exponentially along −z,
which is inadmissible. This is a thermodynamic criterion, evaluated independently
of the balance residual.

**Factors of ½ (point 7).** Every term now carries exactly one `½`, from
`<AB> = ½Re(AB*)`. The Phase-1 note claimed the mechanical and thermal terms
differed in origin; that was a misdiagnosis — the difference was only whether the
quantity was already a cycle average. Corrected text says so.

**Dimensions (point 8).** Each term is verified to be W m⁻²:
`σ·v` is Pa·m/s = W/m²; `ωD_zφ*` is s⁻¹·C/m²·V = W/m²; `q_zΘ*/T0` is
(W/m²)(K)/K = W/m². Check A2 confirms the underlying entries dimensionally.

## 2.4 Modified files and lines

| File | Change |
|---|---|
| `theory/energy_flux_derivation.tex` | §2 replaced entirely (corrected balance law, corrected flux, Table of the 34-order error, the Second-Law sign proposition, the transverse-term remark, the U-5 re-statement); abstract corrected; master balance uses `J_z`; "withdrawn" paragraph added |
| `code/foundation_checks.py` | `flux_z` (line ~287), `flux_z_cycle` (~301), `flux_z_mp` (~436) all corrected; new `flux_z_parts`, `_cycle_mp`, `_dissipation_terms_mp`; `check_D1_flux` rewritten; new `check_D2_balance` |
| `theory/characteristic_matrix.tex` | (P0-1) |
| `literature/material_constants.csv` | (P1-A) |

## 2.5 Executed verification

Four independent tests, each with an external or non-circular reference:

| Test | Reference | Result |
|---|---|---|
| Elastic limit vs `½ρc_Lω²|W|²` | analytic | **1.95e-16** |
| Closed form vs cycle average of *real* fields | independent route | **1.2e-40** (Al6061) |
| Antisymmetry `<J_z>(ξ) = −<J_z>(−ξ)` | exact symmetry | **0** (exact, 3 materials) |
| Second Law `<q·∇θ>` vs `κ|ξΘ|²/(2T0(1+ω²τ_q²))` | analytic | **1.19e-41** |
| Local balance `−2Im(ξ)<J_z> = <D>` | Second Law, no shared amplitude | **7.1e-17** (independent dps=60 evaluation) |

```
[PASS] D1  energy flux: two routes + antisymmetry  <max rel. dev 1.95e-16>  n_exec=12
[UNRESOLVED] D2  local balance d<J_z>/dz = -<D> (thermal branch)
            <n_exec=2, worst rel. residual 1.86e-05>  n_exec=2
```

## 2.6 Why D2 is reported UNRESOLVED, not PASS

The identity **is** established — to 7.1e-17 — when evaluated independently at
dps = 60 with the root taken from the exact coefficient polynomial. Driven through
the repository's own null-vector routine the same identity closes only to
≈1.9 × 10⁻⁵.

The identity is **scale-invariant** in the null vector (both sides are quadratic in
it, so any rescaling cancels). A residual therefore cannot come from a
normalisation *magnitude*; it comes from the null-vector **direction**. For a
matrix of condition number ~10²² the direction is good only to ~10⁻⁵ in double-then-
equilibrated arithmetic. This is the same defect recorded as **U-3**.

The tolerance was **not** loosened to convert this into a PASS. D2 reports
UNRESOLVED and the number. **Gate G4 and G5 are therefore NOT PASSED.**

## 2.7 Remaining uncertainty

1. The coupled-medium verification rests on a null vector whose direction is good
   to ~10⁻⁵. Closing U-3 would let the balance test run through the production
   code path rather than an independent one.
2. The interface dissipation ½-convention (U-5) is now correctly described as
   bookkeeping, but the Phase-2 solver must still pick one and apply it uniformly.
3. The thermal mode's flux at realistic evaluation planes remains analytically
   argued (μ ≈ 3.9 × 10⁻⁷ m), not demonstrated in a scattering configuration.

---

# P0-3 — Numerical-test credibility

## 3.1 What was wrong

The supervisor's list, audited item by item against Phase 1:

| Requirement | Phase-1 state | Verdict |
|---|---|---|
| A zero-case test cannot PASS | I1 returned `PASS` with `n_exec = 0` in an early version; no guard existed in `run()` or `record()` | **VIOLATED** |
| Every claimed symbolic identity is independently evaluated | C1 claimed 16-entry symbolic evaluation; **never performed** | **VIOLATED** |
| Test inputs do not duplicate implementation assumptions | C2 used closed forms derived from the *same* matrix it was testing; D1's cycle average shared the phasor convention with the closed form and so could not detect a wrong functional | **PARTIALLY VIOLATED** |
| Each test reports the number of executed cases | only I1 and I2 mentioned a count, in prose | **VIOLATED** |
| Any test relying on the invalid energy formula is withdrawn/replaced | D1 and E1 both used `½Re[q_z]` | **VIOLATED** |
| Failed tests cannot be silently excluded | `n_fail = len(RESULTS) - npass`; no `UNRESOLVED` category existed at all | **VIOLATED** |
| Unresolved tests are reported separately from passing tests | impossible — no such category | **VIOLATED** |

## 3.2 Fixes implemented

**3.2.1 Mandatory `n_exec` with automatic demotion.**

```python
def record(cid, name, label, status, measured=None, detail="", tol=None,
           n_exec=None):
    if status not in VALID_STATUS:
        status = "FAIL"
    if n_exec is None:
        status = "FAIL"
        detail = "BOOKKEEPING ERROR: ... n_exec not reported ..."
    elif n_exec == 0 and status == "PASS":
        status = "UNRESOLVED"
        detail = "FORCED UNRESOLVED: 0 cases executed ..."
```

A missing count is a **FAIL**. A zero count with a PASS is **demoted to
UNRESOLVED**. Both are deliberate: silence is treated as failure, not success.

**3.2.2 Four status categories, all counted separately.**

```
PASS   FAIL   UNRESOLVED   WITHDRAWN
```

Nothing can fall between categories. Every non-PASS check is printed by name in
the summary, so nothing can drop out of the total.

**3.2.3 Exceptions are FAILs, never skips.**

```python
except Exception as exc:
    return record(cid, name, label, "FAIL", None, f"exception {exc!r}...", tol, 0)
```

**3.2.4 The zero-case bug is now impossible by construction**, and was in fact
triggered during this repair: when `n_exec` was added to `record()` but not yet to
the individual checks, twelve previously-passing checks immediately reported
`FAIL / n_exec=0` with `NameError`. That is the guard working as designed.

**3.2.5 C1 is now a real test**, using `symbolic_M()` — an independent construction
from symbolic material constants, not the numeric tables and not the cached
determinant.

**3.2.6 D1 no longer shares an assumption with its reference.** The cycle average
is built from **separately evaluated real fields** (not from phasor products), and
the decisive test (`check_D2_balance`) uses a dissipation functional that contains
**no amplitude in common** with the mechanical or electrical flux terms.

**3.2.7 Tests relying on the invalid formula were withdrawn and replaced**: `D1`,
`E1` and the whole flux path were rewritten; `check_D2_balance` was added as the
decisive test.

**3.2.8 Machine-readable output** at `validation/phase_01b_test_results.json`,
containing every check, its status, its `n_exec` and the full per-category counts.

## 3.3 Executed result

```
foundation_checks: 15 PASS, 0 FAIL, 1 UNRESOLVED, 0 WITHDRAWN, 16 total
executed test cases across all checks: 151
  UNRESOLVED: D2  local balance d<J_z>/dz = -<D> (thermal branch)
```

Exit status **1**, because UNRESOLVED counts as not-passing. The script does not
exit 0 while anything is unresolved.

## 3.4 Remaining uncertainty

The count of *cases* is now mandatory and audited, but "cases" is not a uniform
unit: check G1 reports 6 sub-tests while C1 reports 16 entries. A future
improvement would be a per-assertion ledger. This is recorded as open item **U-12**.

---

# P1-A — Causality: the τ_q rule is withdrawn

## 4.1 Original statement

`code/foundation_checks.py` H1, Phase 1:

> tau_q,min = kappa/(rho cE cL²) is the value at which the LS thermal speed
> sqrt(kappa/(rho cE tau_q)) equals cL; beyond it (i.e. for tau_q < tau_q,min) the
> thermal wave would be superluminal, **which is rejected on causality grounds**.

and the pass condition included `causal = margin >= 1.0`.

This rule had already forced one substantive change in Phase 1: the CdSe value
1.0 × 10⁻¹³ s was "withdrawn as non-causal" and replaced by 5.0 × 10⁻¹³ s.

## 4.2 Why it was incorrect

The rule conflates two different things:

- **a constitutive relaxation front** — the speed at which the *truncated* Taylor
  expansion of the heat flux relaxes, `v_T = √(κ/(ρc_Eτ_q))`; and
- **a signal speed**, which Special Relativity bounds by `c`.

Lord–Shulman theory is a *phenomenological* extension of Fourier conduction. It
contains **no mechanism that bounds `v_T` by the elastic wave speed**. `v_T > c_L`
is a well-known property of the model at high frequency and has been discussed in
the LS literature for decades; it is a limitation of the model, not a violation of
physics. `c_L` is not even the relevant comparison: the relativistic bound is `c`,
not the longitudinal sound speed, and `v_T ~ 10³ m/s` is nine orders of magnitude
below `c`.

Rejecting a **sourced** relaxation parameter on this basis was therefore
unsupported. The supervisor's instruction is followed exactly: the rule is
withdrawn and no parameter is rejected on speed-ordering grounds.

## 4.3 What is retained as genuinely universal

Enforced (failure = FAIL):

| | Constraint | Reason |
|---|---|---|
| U1 | ρ > 0 | mass |
| U2 | c_E > 0 | thermal capacity |
| U3 | T_0 > 0 | absolute temperature |
| U4 | C_ijkl positive definite | elastic stability |
| U5 | ε_ik positive definite | electrostatic energy |
| U6 | κ_ij positive definite (symmetrised) | conduction must not anti-diffuse |
| U7 | τ_q ≥ 0 | Cattaneo must relax **toward** Fourier, not away |

Reported but **never enforced**:

| | Quantity | Status |
|---|---|---|
| O1 | v_T ≤ c_L | **optional modelling restriction** |
| O2 | f < 1/(2π τ_q,min) | derived from O1; not used as a filter |

## 4.4 Consequence: the CdSe parameter is restored

| | Phase 1 | Phase 1B |
|---|---|---|
| κ (CdSe) | 9 W/(m·K) | **6 W/(m·K)** (midpoint of the sourced bulk range 5–10) |
| τ_q (CdSe) | 5.0 × 10⁻¹³ s ("CHOSEN_ABOVE_CAUSAL_FLOOR") | **1.0 × 10⁻¹³ s** ("VERIFIED_TO_ORDER_OF_MAGNITUDE") |
| v_T / c_L | 1.07 | **1.76** — above 1, and now correctly *allowed* |

Both changes are recorded in `literature/material_constants.csv` with the reason
inline. τ_q for CdSe is **not** an independently measured constant in any source
consulted; it is now reported as an order-of-magnitude value with that uncertainty
stated, rather than being tuned to satisfy a rule that does not exist.

## 4.5 Modified file and verification

`code/foundation_checks.py`, `check_H1_stability` (rewritten).

```
[PASS] H1  thermodynamic stability + LS causality bounds
       <4 media, 0 universal violations>   n_exec=4
```

```
PZT-5H   universal: ALL HELD   ... v_T/cL=0.277  => above the v_T = c_L point
BaTiO3   universal: ALL HELD   ... v_T/cL=0.240  => above
Al6061   universal: ALL HELD   ... v_T/cL=0.908  => above
CdSe     universal: ALL HELD   ... v_T/cL=1.76   => BELOW the v_T = c_L point
```

CdSe is now correctly *reported* as exceeding `v_T = c_L` and correctly *accepted*.

## 4.6 Remaining uncertainty

The label in the check name still says "causality bounds"; the enforced content is
now thermodynamic stability only, and the docstring says so explicitly. Renaming
the check would break continuity with Phase 1 reporting, which is a worse trade.
Note also that τ_q for CdSe remains an order-of-magnitude assumption, so any result
depending quantitatively on the thermal branch of CdSe is not safe.

---

# P1-B — Coupled scattering preparation

## 5.1 The supervisor's instruction

> Do not infer that a singular decoupled formulation must always be regularized by
> coupling → 0. Compare it against a properly reduced elastic subsystem.

## 5.2 Finding: the singularity is a gauge artefact, not a physics problem

Two facts were established, the first numerically and the second exactly.

**Numerically.** Singular-value spectra of the 8×8 coupled interface operator for
Al6061 | Al6061 at 30°:

| λ (coupling scale) | rank | σ_max | smallest five σ |
|---|---|---|---|
| 1 | **4** | 2.18 | 3.2e-08, 5.5e-14, 2.4e-14, 3.4e-20, 9.4e-38 |
| 0.1 → 0 | no complete mode set | — | — |

Even at full coupling the rank is 4, not 8 — the assembled operator is *already*
badly conditioned. As λ → 0 the outgoing mode set cannot be completed at all.

**Exactly** (check J1, all residuals 0):

```
STEP 1   e = p = 0
   phi off-diagonal entries all zero : True  [0, 0, 0, 0, 0, 0]
   det M - (eps33 xi^2 + eps11 k^2) det M3 = 0
   => the phi factor is a GAUGE quadratic; rank 4 -> 3

STEP 2   e = beta = p = 0
   det M3 - det(u_x,u_z)*M_theta = 0
   => rank 3 -> 2 : theta is a gauge/thermal-block variable

STEP 3   the surviving {u_x,u_z} block under sagittal isotropy
   det == C11 C44 (xi^2 - xi_P^2)(xi^2 - xi_S^2) : True
   => the reduced decoupled problem IS the classical elastic P-SV problem
```

So the decoupled limit is singular because the 8×8 formulation **retains two
zero-energy variables (φ, θ) together with two conditions that become 0 = 0**. The
remedy is **reduction** — drop the gauge variables and their conditions — not a
λ → 0 continuation. The supervisor's caution was correct, and Phase 1's proposed
remedy was the wrong one.

## 5.3 What was built

**`code/elastic_psv.py`** — a reduced elastic P–SV scattering solver, written
independently of the piezothermoelastic code path:

- own 4×4 assembly, own mode selection by **flux sign** (not by the sign of Re ξ),
- proper incoming/outgoing separation with the incident root taken from the
  *complement* of the outgoing set (the E-7 lesson),
- no mpmath, no shared state with `foundation_checks.py`,
- state normalisation uses **unit displacement amplitude**, legitimate here because
  all components carry the same dimension — the precise freedom that the coupled
  case lacks.

**Validation (check J2, 18 cases):**

| Test | Result |
|---|---|
| Identical media, 6 angles, welded | R = 0 exactly; T = 1 |
| Normal incidence vs Fresnel | **1.67e-16** |
| Oblique, energy closure, 7 angles | ≤ 1e-15 |
| Spring contact → welded as K → ∞ | monotone convergence |

```
[PASS] J2  reduced elastic P-SV solver vs Fresnel + energy closure
       <max deviation 9.99e-16>   n_exec=18
```

## 5.4 Is the coupled solver ready?

**No.** Stated plainly:

| Requirement | Status |
|---|---|
| Dimensionally consistent state scaling | **NOT SOLVED** — the null-vector direction is good only to ~1e-5; this is the root cause of U-3 |
| Independent mechanical elastic P–SV limit | **DONE** (J2, 1e-16) — now available as the reference |
| Correct incoming/outgoing mode selection | **DONE** — by flux sign, incident root from the complement |
| Electric gauge treatment | **UNDERSTOOD** (J1) — gauge, not a branch; must be *removed*, not scaled |
| Explicit removal of decoupled zero-energy variables | **UNDERSTOOD** (J1) — reduction chain exact; not yet implemented in a coupled solver |
| Boundary residual and rank checks | **PARTIAL** — rank 8 confirmed at 3 angles (I2); rank is 4 with smallest σ = 1e-38 at λ = 1 on identical media, so the assembly is not trustworthy |
| High precision only where conditioning requires | **DONE** — mpmath used exactly where cond ~ 1e22 |

The single blocking item is the **amplitude normalisation**, and the new evidence
localises it precisely: the null-vector *direction*, not its magnitude, for a
system of condition number 10²². No parameter optimisation was started.

## 5.5 Remaining uncertainty

The reduction chain tells us φ and θ are gauge at decoupling. It does **not** tell
us what the correct normalisation of the *coupled* null vector is, which is the
actual obstacle. That is Phase 2 work and is not started here.

---

# Summary table

| Item | Original | Why wrong | Correct | Files | Verified | Residual uncertainty |
|---|---|---|---|---|---|---|
| P0-1 parity | **S** = diag(1,1,1,−1) | fails 6/16 entries; the claimed symbolic check was never run | **S** = diag(1,−1,−1,1) | `characteristic_matrix.tex`, `foundation_checks.py` C1 | 16/16 entries, both candidates, exhaustive sign search | none on algebra; process audit open (U-11) |
| P0-2 flux | `+½Re[q_z]` | linear in one phasor; omits temperature; 34 orders wrong | `+½Re[q_zΘ*]/T0` | `energy_flux_derivation.tex`, `flux_z*` | 5 tests, 1.19e-41 … 1e-15 | null-vector direction ~1e-5 (U-3) |
| P0-3 tests | no `n_exec`, no UNRESOLVED, fake C1 claim | silence could be read as success | mandatory count, 4 categories, FAIL on missing | `foundation_checks.py` | 16 checks, 151 cases | case-count unit non-uniform (U-12) |
| P1-A causality | reject τ_q < κ/(ρc_Ec_L²) | conflates relaxation front with signal speed | rule withdrawn; U1–U7 enforced | `foundation_checks.py` H1, `material_constants.csv` | 4 media, 0 violations | CdSe τ_q is order-of-magnitude only |
| P1-B coupled | propose λ → 0 regularisation | singularity is a gauge artefact | reduction; reduced solver built | `elastic_psv.py`, J1, J2 | 1.67e-16 vs Fresnel | **normalisation unsolved; solver NOT ready** |

**Gate status after Phase 1B: G4 NOT PASSED. G5 NOT PASSED.**
