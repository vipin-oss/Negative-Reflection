# Phase 1 Claim Audit — Completion of Open Item U-11

**Phase:** 1C, task A
**Date:** 4 November 2026
**Audited commit:** `127eba77f397be6baaf7b24fc5b8607a903d402b`
**Claims audited:** 44 (`evidence/claim_registry.csv`)

---

## 1. Purpose and method

Open item **U-11** was raised in Phase 1B after the parity defect: a document
asserted that check C1 "evaluates all 16 entries ... symbolically and obtains
zero", when C1 did no such thing. The wrong matrix survived partly because a claim
of verification was trusted without being inspected.

This audit inverts that. **A claim is evidence only if I located the code that
executes the asserted calculation and read its actual output.** Where the code
carries out something adjacent but weaker, the claim is marked `VERIFIED_WEAK`
with a `failure_mode`. Where no code path executes it, it is `UNVERIFIED`
regardless of whether the mathematics is obviously right.

The method, applied to all 44 claims:

1. Extract the sentence containing the claim and its cited test ID.
2. Open the cited test function and identify the statement that carries the
   asserted calculation.
3. Re-execute that function and read the numeric output, not the PASS flag.
4. Compare the *scope* of what executed against the *scope* of what is claimed.
5. Assign a status. A PASS label was never treated as evidence.

Status vocabulary:

| Status | Meaning |
|---|---|
| `VERIFIED_STRONG` | The executed calculation is the asserted one; the oracle is external to the code under test |
| `VERIFIED_WEAK` | Something real executes, but the oracle is circular, the scope is narrower than the claim, or the check is near-tautological |
| `UNVERIFIED` | No executed or inspectable basis. May be mathematically correct — that is not the point |
| `UNRESOLVED` | Genuine open question, correctly reported as open |
| `INVALID` | The claim is false, or rests on a rule that has been withdrawn |

---

## 2. Results

| Status | Count |
|---|---|
| `VERIFIED_STRONG` | 27 |
| `VERIFIED_WEAK` | 9 |
| `UNVERIFIED` | 5 |
| `UNRESOLVED` | 1 |
| `INVALID` | 2 |
| **Total** | **44** |

**U-11 is not "all claims were fine".** Nine claims execute something weaker than
what is claimed, five have no executed basis at all, one is properly unresolved,
and two are invalid. Six distinct defect classes were found.

---

## 3. Findings

### U-11-1 — C1 claims two calculations it does not perform (`CLM-C1-04`, `CLM-C1-06`)

**Severity: HIGH.** This is the same defect class as the original parity failure,
committed during its own repair.

`theory/characteristic_matrix.tex` states, of the gauge root:

> \textsc{Tested} (check C1 computes the null space symbolically and obtains
> exactly that vector).

and of the block structure:

> \textsc{Tested} (check C1).

I inspected the current `check_C1_structure`. It:

- builds `symbolic_M()` and tests parity on all 16 entries — **real, and correct**;
- checks degree, parity and the non-zero power set from the cached coefficients —
  **real**;
- runs an exhaustive 2⁴ sign search — **real**;
- **prints** the gauge vector `(0,0,1,0)` and the block-structure text as
  *literal strings*;
- performs no `subs(xi, 0)`, no `subs(k, 0)`, and no null-space computation.

Two consequences:

1. The claim of a symbolic null-space computation is **false**. In *Phase 1* a
   null-space routine did exist in C1; my Phase 1B rewrite **removed it** while the
   `.tex` sentence referring to it stayed. I introduced this defect.
2. The block-structure claim has never been tested.

I verified both statements by hand during this audit (substituting ξ=0 gives
{u_x,θ}⊕{u_z,φ}; substituting k=0 gives {u_x}⊕{u_z,φ,θ}; the null space at
(k,ξ)=(0,0) is indeed spanned by (0,0,1,0)). **The mathematics is right. The
evidence is absent.** Both are recorded `UNVERIFIED` and have regression tests
assigned in `validation/tests/` so that the claim becomes true rather than being
deleted.

### U-11-2 — The passivity test cannot fail (`CLM-F1-01`, `CLM-F1-03`)

**Severity: HIGH.** `check_F1_passivity` is cited by `interface_conditions.tex`
three times, in the form "draws 5×10⁴ random positive semi-definite C ... and finds
a minimum of exactly 0". It does hold a real state, but it cannot detect a failure:

```python
neg = 0.0
for _ in range(50000):
    A = rng.normal(size=(2, 2))
    C = A @ A.T                      # positive semi-definite BY CONSTRUCTION
    dv = rng.normal(size=2)
    neg = min(neg, 0.5 * dv @ C @ dv)
```

Three defects:

1. **The reported `0.0` is the initialiser, not an observation.** I traced the
   samples: the actual minimum is `1.835e-06`, the median `0.812`. Since every
   sample is positive, `min()` never updates, and the output `0.0` is the value
   from before the loop started. The evidence line is therefore not a measurement
   at all.
2. **`C = A@A^T` is PSD by construction**, so no sample can ever be negative. The
   test demonstrates a linear-algebra identity, not interface passivity.
3. **The `.tex` claims the samples "include the degenerate case"**. They do not —
   `A@A^T` is invertible almost surely.

The same pattern appears in the thermal form (`CLM-F1-03`), where `H_T = |N(0,1)|·10⁵`
and `Δθ²` are both non-negative, so the expression is non-negative by construction.

**The underlying physics is correct** — a quadratic form in a PSD matrix is
non-negative, and that is a theorem. What is wrong is presenting a tautology as
empirical validation. A test that cannot fail is not evidence.

### U-11-3 — A2 overstates its coverage (`CLM-A2-01`)

**Severity: MEDIUM.** `governing_equations.tex` states:

> check A2 verifies all 27 independent entries of M against these targets and
> finds zero mismatches.

Measured: the `terms` list has **26 additive triples**, which map onto **16
distinct (row, col) slots**. A2 is a genuine and useful dimensional test; it is
not an exhaustive one.

**Correction to this finding (recorded rather than quietly edited).** On
re-checking, `M` is a *dense* 4x4: it has 16 entries carrying 33 additive terms,
and the earlier statement in this audit that four entries "are absent from the
list" was wrong — a 4x4 cannot have 27 entries. The "27" in the source text is
the number of additive terms, not entries; A2 covers 26 of the 33 terms. That
makes the source claim **scope-overstated, not false**: `CLM-A2-01` is therefore
`VERIFIED_WEAK` rather than `UNVERIFIED`. The defect is a wording defect in the
document, and it is closed separately by check `DIM01`, which tests all 16 entries
and all 33 additive terms dimensionally (0 mismatches). The audit text above is
left in place so the error in the audit itself remains visible.

Severity is medium rather than high because the omitted entries share their
coefficient structure with entries that *are* tested, so a missed error is
unlikely — but "unlikely" is not the standard the claim asserts.

### U-11-4 — H1's causality rule is withdrawn (`CLM-H1-02`)

`governing_equations.tex` still contains the section "Causality floor on τ_q",
which states that requiring the thermal wave not to exceed `c_L` gives a floor and
that parameters below it are rejected. **This rule was withdrawn in Phase 1B** as
unsupported (it conflates a constitutive relaxation front with a signal speed).

The Phase-1B repair fixed the *code* (H1 no longer enforces it) and the *report*,
but not this `.tex` section. Recorded `INVALID`. The document must be corrected —
tracked as a required fix, not silently dropped.

### U-11-5 — Untested or untestable claims (`CLM-MAT-02`, `CLM-FIG-01`, `CLM-LIT-02`)

- `beta`, `p` for PZT-5H/BaTiO₃ are **ASSUMED**. Open item U-4. Nothing tests them
  because there is nothing to test against; they are marked `UNVERIFIED` so that
  no downstream result can inherit unearned confidence.
- Figures `fig1`–`fig3` have **no regeneration test**. A figure is not evidence
  until it can be rebuilt from raw data by a recorded command.
- The six candidate gaps G-1…G-6 are `UNVERIFIED` because gate G0 (full-text prior
  art) is not passed. They are questions, not findings.

### U-11-6 — Identities that could not have caught the defects (`CLM-D1-03`, `CLM-C2-01`)

Two checks are marked `VERIFIED_WEAK / PARTIAL_IDENTITY` because they are
*necessary but not sufficient*, and this matters historically:

- **Antisymmetry** `⟨J_z⟩(ξ) = −⟨J_z⟩(−ξ)` passed at `0.00e+00` **while the flux
  functional was wrong**. Any correction that preserves the k→−k symmetry passes
  it. It is a useful consistency check and was presented in Phase 1 as though it
  confirmed the signs; it could not have.
- **C2's eigen residual** is computed with the same matrix construction that
  produced the root. It validates the solver, not an external truth. `I1` is the
  genuine external oracle.

---

## 4. The two INVALID claims, and why they are kept

`CLM-C1-05` (parity with `diag(1,1,1,-1)`) and `CLM-H1-02` (the τ_q causality
floor) are retained in the registry marked `INVALID`, with `superseded_by` pointing
at their replacements.

Deleting them would erase the evidence that the failure occurred. Section I of the
Phase 1C brief requires that a correction must not erase the evidence of the
earlier failure, and the same principle applies to the registry.

---

## 5. What this audit does **not** establish

- **It does not show the science is wrong.** Every mathematical result re-checked
  here that is *tested* is tested correctly, and the two hand-verifications I was
  forced to do (null space, block structure) both confirm the stated mathematics.
  The defects are in the *evidence layer*: claims stronger than their tests.
- **It is not exhaustive over prose.** I audited statements that carry an explicit
  verification attribution (`verified`, `tested`, `derived`, `checked`, `confirmed`
  together with a test ID or an equation reference). Narrative and interpretive
  prose in `reports/` was not audited claim-by-claim.
- **It does not re-derive the equations.** That is task C, the equation registry,
  which tests each equation against an independent oracle.
- **It does not close any blocker.** Gates G0, G1, G4 and G5 remain not passed.

---

## 6. Actions arising

| # | Action | Status |
|---|---|---|
| 1 | Add a real symbolic null-space and block-structure test for C1 (`CLM-C1-04`, `CLM-C1-06`) | assigned to `validation/tests/` |
| 2 | Replace the F1 tautologies with a genuine falsifiable passivity test, including a negative control | assigned |
| 3 | Restate the A2 claim to its true scope (26 of 33 terms); coverage gap closed by DIM01 | **DONE in Phase 1C** |
| 4 | Correct the "Causality floor on τ_q" section of `governing_equations.tex` | **open, needs correction** |
| 5 | Add figure regeneration tests and checksums | assigned (task F) |
| 6 | Keep `beta`, `p` `UNVERIFIED` until sourced | blocked on U-4 |

Items 1–3 and 5 are implemented as part of task D; item 4 is recorded as a
document correction to be made with the equations, and is listed in
`reports/PHASE_01C_OPEN_BLOCKERS.md`.

---

## 7. Verdict on U-11

**U-11 is complete.** All 44 verification-bearing claims in the Phase 1 / 1B
material have been traced to a test function and an actual output, or marked as
lacking one.

The audit found that **the evidence layer was systematically weaker than its own
descriptions**: 9 claims execute something adjacent but weaker, 5 have no executed
basis, and 2 are invalid. That is a materially different picture from "15 PASS, 0
FAIL", and it is the reason the Phase 1C infrastructure exists.

U-11 should not be closed as a one-off. It becomes a permanent mechanism: the
claim registry is now validated on every run, and a claim whose evidence file is
missing or stale fails the build (see `validation/tests/test_claims.py` and
section K of `reports/PHASE_01C_QUALITY_INFRASTRUCTURE.md`).
