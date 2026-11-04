# Negative Reflection of Elastic and Coupled Waves — Research Discovery Dossier

> **STATUS: AWAITING INDEPENDENT RESEARCH REVIEW — NOT APPROVED**
>
> This is a research-discovery exercise, not a manuscript. Nothing here has been approved, funded or
> started, and no experimental work is claimed.

Prepared for **Dr. Vipin Gupta**, Department of Mathematics, Gurugram University, Gurugram, India
(ORCID [0000-0002-5238-1963](https://orcid.org/0000-0002-5238-1963)) · compiled 2026-11-04.

---

## The ten deliverables

| # | Document | Phase | What it establishes |
|---|----------|-------|---------------------|
| 01 | [`01_Researcher_Profile.md`](reports/01_Researcher_Profile.md) | 1 | Verified 42-work record, homonym disambiguation, and an explicit map of what is *established* versus *absent* (no guided-wave edge scattering, no ZGV, no experiments). |
| 02 | [`02_Verified_Literature.xlsx`](reports/02_Verified_Literature.xlsx) | 2 | 60 DOI-verified records in 6 classes, each with a relevance note; plus an excluded-records sheet and the search-method statement. |
| 03 | [`03_Negative_Reflection_Taxonomy.md`](reports/03_Negative_Reflection_Taxonomy.md) | 2 | Ten phenomena separated by operational tests; four distinct NR mechanisms (M1–M4). Core thesis: *a negative reflection coefficient is neither necessary nor sufficient for negative reflection.* |
| 04 | [`04_Research_Gap_Evidence.md`](reports/04_Research_Gap_Evidence.md) | 2 | Four graded gaps (GAP-1…4) plus one lemma gap, each with what is established, what is not, how it was searched, and its strength. |
| 05 | [`05_Physical_Feasibility_Foundation.md`](reports/05_Physical_Feasibility_Foundation.md) | 3 | The mechanism, the flux identity (verified to 5.3 × 10⁻⁷), and three computed premise results in a coupled plate. |
| 06 | [`06_Candidate_Gap_Matrix.md`](reports/06_Candidate_Gap_Matrix.md) | 4 | Candidates A–J plus four new ones (K–N) scored on nine questions and labelled VERIFIED / PLAUSIBLE / UNVERIFIED / ALREADY INVESTIGATED. |
| 07 | [`07_Five_Research_Problems.md`](reports/07_Five_Research_Problems.md) | 5 | Five mechanism-distinct problems, P1–P5, each against all 18 required elements. |
| 08 | [`08_Mathematical_Feasibility_Audit.md`](reports/08_Mathematical_Feasibility_Audit.md) | 6 | 15 checks × 5 problems. No outright rejections; two problems carry a standing auto-reject gate. |
| 09 | [`09_Comparative_Research_Ranking.md`](reports/09_Comparative_Research_Ranking.md) | 7 | Weighted 0–10 ranking with published arithmetic, a weight-sensitivity test, and the three named awards. |
| 10 | [`10_RESEARCH_PROPOSAL_FOR_APPROVAL.md`](reports/10_RESEARCH_PROPOSAL_FOR_APPROVAL.md) | 8 | Full design of P1 + P4, P3 as escalation path, and the decision requested from the reviewer. |

## Headline outcome

| Problem | Score | Rank | Recommendation |
|---|---|---|---|
| **P1** Bias-tuned NR in a polarized piezoelectric plate | **7.80** | 1 | **Primary — execute** |
| **P4** NR at a spring-type joint | 7.68 | 2 | **Companion — execute** |
| P2 Grading as mode-shape control | 7.38 | 3 | Fallback if bias data fail |
| P3 Higher-order BCs in a strain-gradient plate | 7.30 | 4 | Escalation path, not now |
| P5 Anisotropic NR criterion in LiNbO₃ | 6.88 | 5 | Not recommended now |

**Awards:** BEST OVERALL **P1** · BEST HIGH-RISK/HIGH-REWARD **P3** · BEST MATHEMATICALLY TRACTABLE **P2** (P4 tied).

## Three computed results that shape everything else

These were obtained in this project, not assumed, and each one closed or opened a research direction:

1. **A forward/backward fold exists in a coupled plate** at f = 1.12975 MHz, k = 809.0 m⁻¹
   (PZT-5H, 1.5 mm) — the premise of the whole programme, quantified.
2. **Piezoelectric coupling shifts it by −35.4 %, identically at 0.05, 0.5 and 1.5 mm.** Since
   Π₁ = e₃₁²/(ε₃₃c₁₁) contains no length, "size-dependent piezoelectric negative reflection" is
   dimensionally inadmissible. That framing was removed.
3. **The electrical open/short switch shifts this fold by −0.0006 %** — so "switch the electrode to
   tune the reflection" is a dead end, killed before a paper was written about it. This is why P1 is
   built on **DC bias** (a stiffness control) instead.

A fourth result — that the mode's potential is *antisymmetric* through the thickness
(`reports/figures/fig3b`) — replaced an earlier, wrong explanation ("the potential is nearly uniform")
that had been written into the working notes. The corrected statement is the one carried forward.

## Reproducing the numbers

```bash
cd Negative-Reflection
python3 code/gen_references.py            # -> data/references.csv, data/excluded_records.csv
python3 code/check_01_elastic_validation.py   # -> data/check_01_elastic_validation.json
python3 code/check_02_flux_and_zgv.py     # -> data/check_02_flux_and_zgv.json  (~50 s)
python3 code/make_figures.py              # -> reports/figures/fig1..fig3 (.pdf/.png)
```

Validation tolerances achieved: 6.3 × 10⁻⁷ (cutoffs vs exact), 5.3 × 10⁻⁷ (dispersion vs exact
Rayleigh–Lamb), 9.1 × 10⁻⁷ (ZGV vs independent tracker), 8.0 × 10⁻⁹ (mesh convergence),
5.3 × 10⁻⁷ median (flux identity v_e = P/(2⟨T⟩) = dω/dk). Three real bugs found by these checks are
recorded in `05` §7.

## Figures

| Figure | Content |
|---|---|
| `reports/figures/fig1_dispersion_ve.*` | The S₁/S₂ᵦ fold with the sign of the energy flux; the flux changes sign exactly at the ZGV. No branch tracking used. |
| `reports/figures/fig2_validation_premise.*` | (a) solver validation against exact theory; (b) f_ZGV ∝ 1/h across four configurations. |
| `reports/figures/fig3_premise_limits.*` | (a) per-branch electrical sensitivity — most branches are strongly active, this one is not; (b) the mode's antisymmetric potential. |

## Integrity notes

- No paper, title, author, DOI or journal in this dossier was invented. Every DOI was resolved, and
  one secondary-source attribution was found wrong and **excluded rather than silently corrected**
  (`data/excluded_records.csv`: `10.1121/1.2932071`).
- Nothing is claimed to be "unexplored" in an absolute sense. Scopus and Web of Science were not
  reachable from this environment, so every gap is graded as a *degree of evidence*. The investigator
  must repeat the search on the subscription indexes before committing.
- No Q-rank, impact factor or acceptance probability is claimed anywhere.
- Negative reflection is **not** claimed to exist in any coupled material. Its existence in isotropic
  plates is established in the literature; the existence of the *fold* in a coupled plate is computed
  here; whether that fold yields usable negative reflection is the open question the proposal is
  designed to answer, starting with a cheap go/no-go gate.

---

## Audit of Blueprint V2.0 (4 Nov 2026)

Referee audit of the author's own `FINAL RESEARCH BLUEPRINT — V2.0`
(piezothermoelastic anti-reflection). Verdict: **REDESIGN**.

- `ARENA_V2_CRITICAL_AUDIT.md` — full referee report, 4 critical / 6 major / 6 minor
  findings with exact corrections, plus a staged redesign.
- `ARENA_V2_CLAIM_EVIDENCE.csv` — 23 evidence rows, one per checkable claim.
- `code/check_04_blueprint_audit.py` — all supporting numbers (Parts A–H).
- `data/check_04_blueprint_audit.json`, `data/check_04_audit_log.txt` — results, log.

Run: `python3 code/check_04_blueprint_audit.py` (~14 s).
Three self-checks must pass on every run: uniform-medium limit (|R| = 3e-16),
lossless energy balance (residual <= 7e-16), matched-layer quarter-wave limit
(R_E = 6e-10).

---

## Controlled audit V3 (4 Nov 2026)

`CONTROLLED_BLUEPRINT_AUDIT_V3.md` — re-audit of Blueprint V2 after supervisor
dispute. Verdict: **CONDITIONAL REVISION** (prior audit's REDESIGN downgraded).

Five prior-audit claims are corrected here, two of them rejected outright:
- **D2** "with loss only gain closes the reflection floor" — rejected as a
  general claim; `Re Z_in < Z₁` in 198/366 configurations once the interphase
  impedance is free.
- **D3** "the thermal field is inert" — rejected; it is a statement about τ_q.
  Replaced by the constraint `v_T/c_L = √(2π f_on τ_min)` with causality
  `τ_q ≥ κ/c_L²`.
- **D9** an "exact reflection zero at k = 2890.43" found during this audit was
  **retracted** — two compounding sign/term errors. Corrected: no zero exists.
- **D10** *(new)* piezoelectric coupling **destroys** the elastic SH reflection
  minimum (elastic: |R| = 8.9e-7; piezo: strictly monotone, min = bulk mismatch).
  Consequence: reflection minima must be a P–SV phenomenon.

Recommended: **Option A** (two half-spaces + imperfect interface), 8.52/10.

- `code/check_05_dispute_resolution.py` — Parts Q1/Q2/Q4/Q5/Q6, 8 self-checks
- `data/check_05_dispute_resolution.json`, `data/check_05_log.txt`

Run: `python3 code/check_05_dispute_resolution.py` (~7 s).

---

# Phase 1 (Blueprint V5.0) — Foundation deliverables and reproduction

> **STATUS: PHASE 1 COMPLETE AND PHASE 1B-REPAIRED, NOT APPROVED.** No reflection,
> transmission or absorption coefficient is reported, and **no novelty claim is
> made**. Gates **G4 and G5 are NOT PASSED**.
>
> Three Phase-1 defects were confirmed by supervisor review and corrected; see
> [`reports/PHASE_01B_CRITICAL_REPAIR.md`](reports/PHASE_01B_CRITICAL_REPAIR.md)
> and [`reports/PHASE_01B_SUPERVISOR_RESPONSE.md`](reports/PHASE_01B_SUPERVISOR_RESPONSE.md).
> A fourth change was a **withdrawal** rather than a correction: the τ_q causality
> rule was wrong and is gone.

Phase 1 implements the theoretical foundation demanded by Blueprint V5.0 §3–§7.
Phase 1B repairs it. Both stop before any numerical scattering experiment.

## What was produced

| Path | Content |
|---|---|
| [`blueprint/Piezothermoelastic_AntiReflection_Master_Blueprint_V5.tex`](blueprint/Piezothermoelastic_AntiReflection_Master_Blueprint_V5.tex) | Faithful LaTeX transcription of the V5.0 blueprint (source text only — **not compiled**, no TeX engine available) |
| [`literature/verified_references.csv`](literature/verified_references.csv) | 10 Crossref-verified records, R02–R11, 15 columns. Every row `read_status = metadata only` |
| [`literature/material_constants.csv`](literature/material_constants.csv) | 80 constitutive constants, 4 materials, each tagged VERIFIED / ASSUMED / DERIVED / ZERO / VERIFIED_TO_ORDER_OF_MAGNITUDE |
| [`literature/research_gap_matrix.md`](literature/research_gap_matrix.md) | 3 direct competitors, 6 candidate gaps (all UNRESOLVED), 7 closed derived results |
| [`theory/governing_equations.tex`](theory/governing_equations.tex) | Single-potential 6 mm piezothermoelastic formulation, Lord–Shulman conduction, stability conditions |
| [`theory/interface_conditions.tex`](theory/interface_conditions.tex) | Spring–dashpot, capacitive-gap and thermal-conductance contact laws; dissipations; the `h → 0` interphase map |
| [`theory/characteristic_matrix.tex`](theory/characteristic_matrix.tex) | The 4×4 `M(ξ,k,ω)`, its exact `det M`, parity, block structure, gauge roots, branch inventory |
| [`theory/energy_flux_derivation.tex`](theory/energy_flux_derivation.tex) | Normal total-field flux, bulk positivity, interface dissipation, the composite power balance from which `R + T + A = 1` *follows* |
| [`validation/benchmark_plan.md`](validation/benchmark_plan.md) | All 19 benchmarks + 3 energy gates, each with expected value, source, acceptance rule and honest status |
| [`reports/UNRESOLVED_ISSUES.md`](reports/UNRESOLVED_ISSUES.md) | 9 self-corrected errors and 10 open items (partly superseded by Phase 1B) |
| [`reports/PHASE_01B_CRITICAL_REPAIR.md`](reports/PHASE_01B_CRITICAL_REPAIR.md) | Every Phase-1B issue: original statement, why it was wrong, correct derivation, file and line, executed verification, remaining uncertainty |
| [`reports/PHASE_01B_SUPERVISOR_RESPONSE.md`](reports/PHASE_01B_SUPERVISOR_RESPONSE.md) | Point-by-point response, gate status, open items |
| [`code/elastic_psv.py`](code/elastic_psv.py) | Reduced elastic P–SV scattering solver — the validated reference for the coupled problem |
| [`validation/phase_01b_test_results.json`](validation/phase_01b_test_results.json) | Machine-readable results for all 16 checks |
| [`code/foundation_checks.py`](code/foundation_checks.py) | The 13-check automated suite (output in [`data/foundation_checks.json`](data/foundation_checks.json)) |

## Reproducing the foundation checks

```bash
pip install --break-system-packages numpy sympy mpmath
python3 code/foundation_checks.py          # ~6 s
```

Expected final lines:

```
foundation_checks: 15 PASS, 0 FAIL, 1 UNRESOLVED, 0 WITHDRAWN, 16 total
executed test cases across all checks: 151
  UNRESOLVED: D2  local balance d<J_z>/dz = -<D> (thermal branch)
```

**Exit status 1**, because one check is unresolved. That is intentional: the
script does not report success while anything is unresolved, and a check that
executes zero cases can never be reported as PASS.

Machine-readable results: [`validation/phase_01b_test_results.json`](validation/phase_01b_test_results.json).

### What the 16 checks test

| Check | Status | n_exec | Tests | Independent reference? |
|---|---|---|---|---|
| A1 | PASS | 4 | Material tensor symmetry and stability | no — internal |
| A2 | PASS | 26 | SI dimensional consistency of every entry of `M` | yes — dimensional analysis |
| B1 | PASS | 6 | Constitutive set derives from **one** potential | yes — potential identity |
| B2 | PASS | 4 | Lord–Shulman heat equation from the entropy balance | yes — entropy balance |
| C1 | PASS | 16 | `det M` degree, **parity entry-by-entry for both candidate S**, gauge root | yes — symbolic structure |
| C2 | PASS | 8 | 8 roots: branch identification, eigen residual 2.2e-27 | no — internal |
| D1 | PASS | 12 | Corrected flux: analytic limit, independent cycle average, antisymmetry | yes — ½ρc_Lω²\|W\|² |
| **D2** | **UNRESOLVED** | 2 | Local balance `d<J_z>/dz = −<D>` | yes — Second Law |
| E1 | PASS | 10 | Flux conservation along z | yes — analytic decay |
| F1 | PASS | 3 | Interface passivity | yes — quadratic forms |
| G1 | PASS | 6 | Fresnel, Proposition 1(a)(b)(c), quarter-wave | **yes — closed forms** |
| H1 | PASS | 4 | Universal stability only; speed-ordering rule withdrawn | yes — thermodynamic |
| I1 | PASS | 24 | Exact roots vs closed-form elastic branches | **yes — external, 1.65e-15** |
| I2 | PASS | 3 | Interface rank | no — internal |
| J1 | PASS | 5 | Decoupled limit is a **reduction**, not coupling→0 | yes — symbolic factorisation |
| J2 | PASS | 18 | Reduced elastic P–SV solver vs Fresnel + energy closure | **yes — 1e-16** |

Seven checks compare against references external to the code.

## Two results worth naming

**The exact determinant.** For 6 mm symmetry the characteristic determinant is
`det M = P₄(ξ²)`, with powers `{0,2,4,6,8}` and leading coefficient

```
P₈ = −C₄₄ · κ₃₃ · (C₃₃ · ε₃₃ + e₃₃²)
```

so the piezoelectrically stiffened 1-D modulus `C₃₃ + e₃₃²/ε₃₃` appears in the
highest-order coefficient. Verified symbolically, not asserted.

**A causality rule that was wrong, and has been withdrawn.** Phase 1 imposed
`τ_q ≥ κ/(ρ c_E c_L²)` and, on that basis, changed the CdSe working value from
1.0e-13 s to 5.0e-13 s. The rule conflated a **constitutive relaxation front** with
a **signal speed**: Lord–Shulman theory contains no mechanism bounding the thermal
wave speed by the elastic speed, and `v_T ~ 10³ m/s` sits nine orders of magnitude
below `c`. The rule is **withdrawn**, the sourced value is **restored**, and the
material table no longer contains a value tuned to satisfy a rule that does not
exist. What is retained as universal is only ρ>0, c_E>0, T_0>0, C>0, ε>0, κ>0,
τ_q≥0 — the last because Cattaneo must relax *toward* Fourier, not away.

## What Phase 1 did not do

Stated plainly, because the gates demand it:

- **No coupled scattering solve.** No `R`, `T` or `A` was computed for any
  configuration. Gate G5/E3 is `UNAVAILABLE`. The amplitude normalisation of the
  mode-matching system could not be validated, so the attempt was withdrawn rather
  than tuned (U-3).
- **No full-text prior-art check.** All ten references are metadata-only, so gate
  G0 is not passed and the "candidate gaps" are questions, not findings.
- **No novelty claim.**
- **The `.tex` files were not compiled** — no TeX engine is installed. They are
  text-verified only.
- **No anti-reflection curve, parameter search, or optimised contact.**

## A note on the nine corrected errors

`reports/UNRESOLVED_ISSUES.md` Part 1 lists nine real defects in this phase's own
work, each found by running the checks rather than by review. Three were the
dangerous kind that produce a *silently wrong* answer instead of a visible failure:
a validation check that reported `PASS` after executing zero cases, an eigenvector
normalisation that mixed metres with volts, and a causality constraint that was
computed but never tested. All three are now guarded in code.
