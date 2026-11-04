# Phase 1C — Research Evidence, Validation and Regression Infrastructure

**Repository:** `https://github.com/vipin-oss/Negative-Reflection`
**Branch:** `arena/a3ede4fe-negative-reflection`
**Base commit at phase start (verified against the remote):** `127eba77f397be6baaf7b24fc5b8607a903d402b`
**Author:** Vipin Gupta, ORCID 0000-0002-5238-1963
**Scope:** infrastructure only. No new simulations, no reflection optimisation, no
manuscript text. Phase 2 is not started.

---

## 1. What this phase was for

Phase 1B repaired two serious defects — a wrong parity transformation and a
thermal flux term that was 34 orders of magnitude wrong on the thermal branch.
Both had passed review. The reason they passed is worth stating precisely,
because it determines what this phase had to build:

> The documents cited tests. The tests ran. The tests passed. **But the tests did
> not execute the calculations the documents cited them for.**

`check_C1_structure` printed the phrase "compute null space" as a string literal;
it had not computed one since a Phase-1B rewrite deleted the routine. The claim
above it continued to say the null space had been obtained symbolically. Nobody
was lying; there was simply no mechanism connecting a claim to the code that was
supposed to support it.

Phase 1C builds that mechanism, plus the surrounding apparatus — dependency
graph, equation registry, mutation suite, experiment manifests, baseline history —
and audits every existing claim against it.

**The headline result is uncomfortable and should be read first:** of 44 audited
claims, 27 are `VERIFIED_STRONG`, 9 `VERIFIED_WEAK`, 5 `UNVERIFIED`, 1
`UNRESOLVED` and **2 `INVALID`**. Two statements in the committed theory documents
are false. The infrastructure built here now makes that impossible to overlook.

---

## 2. Deliverables

All paths are relative to the repository root.

| # | deliverable | status |
|---|---|---|
| A | `evidence/claim_registry.csv`, `evidence/phase_01_claim_audit.md` | complete |
| B | `evidence/dependency_graph.json`, `evidence/impact_analysis.py` | complete |
| C | `evidence/equation_registry.csv` | complete |
| D | `validation/tests/` (5 modules) + `validation/harness.py` | complete |
| E | `validation/mutation_tests/` (11 mutations) | complete |
| F | `validation/evidence/*/manifest.json`, `validation/make_manifests.py` | 3 experiments |
| G | `validation/run_fast.py`, `run_standard.py`, `run_full.py` | complete |
| H | `.github/workflows/research-validation.yml` | **green on the remote**: run 41818854579, all 15 steps |
| I | `evidence/approved_baselines.json` | 3 baselines, 2 retained as rejected |
| J | `evidence/manuscript_evidence_map.json`, `test_05_manuscript.py` | complete |
| K | `reports/PHASE_01C_OPEN_BLOCKERS.{md,json}` | 8 blockers, 7 blocking |
| L | all items above plus `requirements-research.txt`, this report | complete |
| M | `validation/demonstrate.py` (8 demonstrations) | executed, see §7 |
| N | commit, push, remote SHA verification | see §10 |

---

## 3. Task A — the claim audit (U-11)

**44 claims, each with:** a stable id (`CLM-<AREA>-<nn>`), the exact source file
and location, the named test and test function, whether that test *actually
executes the asserted calculation*, the input fixture, the actual output, a
status, and a failure mode.

### 3.1 Status tally

| status | count | meaning |
|---|---|---|
| `VERIFIED_STRONG` | 27 | executed, independent oracle, justified tolerance |
| `VERIFIED_WEAK` | 9 | executed, but the oracle is circular or the scope is narrower than the claim |
| `UNVERIFIED` | 5 | no executed or inspectable basis — plausibility is not evidence |
| `UNRESOLVED` | 1 | an execution exists and does not close |
| `INVALID` | 2 | **the statement is false** |

### 3.2 The six findings

**U-11-1 (HIGH) — a claim to a calculation that no longer exists.**
`characteristic_matrix.tex` states that check C1 "computes the null space
symbolically and obtains exactly that vector" and tests the block structure. The
Phase-1B rewrite removed the null-space routine; C1 now emits literal strings.
The mathematics is right — I confirmed by hand that the gauge vector is `(0,0,1,0)`
at `ξ=0` and that `{u_x,θ} ⊕ {u_z,φ}` is the correct block structure — but
**nothing executes it**. `CLM-C1-04` and `CLM-C1-06` are `UNVERIFIED`.
*Closed in this phase* by `test_EIG03` and `test_EIG04`, which do the computation
for real. The original claims stay `UNVERIFIED`, because the *original function*
still does not do it.

**U-11-2 (HIGH) — a test that cannot fail.**
`check_F1_passivity` forms `C = A @ Aᵀ`, which is positive semi-definite by
construction, then reports `neg = 0.0` as "minimum observed sample". Both are
wrong: the 0.0 is the **initialiser**, never overwritten because the loop's `min`
is only reached when a negative value appears, and no negative value can appear.
Measured truth: the actual minimum sample is `1.835e-06`, median `0.812`. The
document's claim that the samples "include the degenerate case" is false.
`CLM-F1-01` and `CLM-F1-03` are `VERIFIED_WEAK` with failure mode
`CIRCULAR_ORACLE`.
*Closed* by `test_IFC01`, which includes a **negative control**: an indefinite
matrix for which the form must go negative. The test fails if the control is not
detected. That is what makes it falsifiable.

**U-11-3 (MEDIUM) — scope overstated.**
`governing_equations.tex` says check A2 "verifies all 27 independent entries of M
… and finds zero mismatches". A2 executes 26 additive terms. **Correction recorded
during this phase:** M is a dense 4×4, so it has 16 entries carrying 33 additive
terms — the "27" is a term count, not an entry count, and an earlier draft of this
audit wrongly said four entries were untested. A 4×4 cannot have 27 entries. The
claim is therefore `VERIFIED_WEAK` (`SCOPE_OVERSTATED`), not `UNVERIFIED`, and the
audit text retains the original error so the error in the audit stays visible.
*Closed* by `test_DIM01`, which tests all 16 entries and all 33 terms
dimensionally with 0 mismatches.

**U-11-4 (MEDIUM) — a withdrawn rule still printed.**
`governing_equations.tex` still contains the "Causality floor on τ_q" section
asserting a lower bound on τ_q from requiring the thermal wave not to outrun `c_L`.
That rule was withdrawn in Phase 1B as non-physical. The code was fixed, the
report was fixed, the document was not. `CLM-H1-02` is `INVALID`.
**Still open** as `BLK-DOC-01`: the document correction is owed and is not done in
this phase, because silently patching a `.tex` while auditing it would destroy the
audit's own evidence.

**U-11-5 — assumed constants and unregenerable figures.**
`beta`, `p` are `ASSUMED`; figures `fig1`–`fig3` have no regeneration test and no
raw data. `G-1…G-6` are `UNVERIFIED` because gate G0 has not been passed.

**U-11-6 (INFO) — an identity that passes no matter what.**
The antisymmetry check `⟨J_z⟩(ξ) = −⟨J_z⟩(−ξ)` passed at `0.00e+00` while the flux
was **wrong**, because any correction symmetric under `k → −k` satisfies it. C2's
eigenvalue residual reuses the same construction that produced the root, so it is
not independent; I1 (passivity of every root) is the genuine external oracle.

### 3.3 Why `INVALID` claims are kept

`CLM-C1-05` (the false parity `S = diag(1,1,1,−1)`) and `CLM-H1-02` (the withdrawn
causality floor) remain in the registry permanently. A correction that erases the
record of the failure is worse than the failure: it removes the only evidence that
the review process missed something. The same principle applies to
`EQ-GOV-12`, `EQ-CM-06` and `EQ-ENG-09` in the equation registry.

---

## 4. Tasks B and C — dependency graph and equation registry

### 4.1 The graph

48 nodes, 56 edges, 0 fatal integrity problems, 0 warnings:

```
source → constant → equation → module → test → claim → experiment → figure → manuscript
```

`evidence/impact_analysis.py` answers the operational question. Changing the
Lord–Shulman heat equation:

```
$ python3 evidence/impact_analysis.py --changed EQ-GOV-08
affected downstream nodes: 13
  d1  EQ-CM-01    equation  VERIFIED_WEAK    via EQ-GOV-08
  d1  EQ-ENG-03   equation  VERIFIED_STRONG  via EQ-GOV-08
  d1  EQ-GOV-10   equation  VERIFIED_STRONG  via EQ-GOV-08
  d2  EQ-CM-02    equation  VERIFIED_STRONG  via EQ-CM-01
  d2  EQ-ENG-04   equation  UNRESOLVED       via EQ-ENG-03
  d2  EQ-SCT-01   equation  NOT_IMPLEMENTED  via EQ-CM-01
  d3  BLK-G4      blocker   OPEN             via EQ-ENG-04
  d3  BLK-G5      blocker   OPEN             via EQ-SCT-01
  d3  EQ-ENG-08   equation  UNRESOLVED       via EQ-ENG-04
  d3  EQ-SCT-02   equation  UNRESOLVED       via EQ-SCT-01
  d4  EQ-ENG-10   equation  DEFINED_NOT_TESTED via EQ-ENG-08
```

Two design decisions matter:

1. **Non-propagating edges.** `EQ-CM-06 (INVALID) --correction--> EQ-CM-03` exists
   so the replacement is traceable, but impact analysis must not follow it,
   otherwise a correction would invalidate the thing that corrects it. Verified:
   querying `EQ-CM-06` returns an empty downstream set.
2. **Documentation is not modelled.** `README.md` is not a node, and a query on it
   exits non-zero. That is the mechanism by which "an unrelated documentation
   change does not invalidate numerical results" is *true by construction* rather
   than by assertion.

### 4.2 The equation registry

41 equations (40 at audit, plus `EQ-EIG-01` added when the root-set was registered). Each carries assumptions, variables with dimensions, derivation
status, the symbolic test, an **independent** oracle description, downstream
dependencies and a verification status.

| status | count |
|---|---|
| `VERIFIED_STRONG` | 17 |
| `VERIFIED_WEAK` | 7 |
| `UNRESOLVED` | 9 |
| `DEFINED_NOT_TESTED` | 3 |
| `PARTIAL` | 2 |
| `INVALID` | 3 |

The registry enforces the independence rule explicitly. `EQ-GOV-11`'s oracle is
marked `PARTIAL` because I1 verifies energy positivity (external) while H1 reuses
`eigvalsh` on the same matrix. Two functions calling the same routine are **not**
independent verification, and the registry records that rather than hiding it.

---

## 5. Task D — the automated suite

### 5.1 The harness makes un-evidenced tests unwritable

`validation/harness.py` refuses to construct a case that does not declare: a
stable test id, the scientific purpose, the oracle **and whether it is
independent**, the input fixture, the expected result, the tolerance **and its
justification**, and the level. A case with `tolerance=None` must justify why.

Three rules are enforced in code, not by convention:

- `n_exec == 0` can never be `PASS` or a successful `NOT_APPLICABLE` — it is
  forcibly converted to `UNRESOLVED`;
- `FAIL`/`UNRESOLVED`/`WITHDRAWN` are never collapsed into `PASS`;
- an exception inside a case body is a `FAIL`, never a skip.

`finish()` then turns any non-`PASS` into a pytest failure, so the runner and the
harness can never disagree.

### 5.2 Coverage of the 12 required categories

| # | category | tests |
|---|---|---|
| 1 | source/material-data integrity | `MAT01`, `MAT02`, `DAT01` |
| 2 | dimensional consistency | `DIM01` |
| 3 | constitutive symmetry / thermodynamics | `CONS01` |
| 4 | characteristic-matrix identities | `MATX01`, `MATX02` |
| 5 | eigenvalues / eigenvectors | `EIG03`, `EIG04` |
| 6 | interface BCs | `IFC01`, `IFC02` |
| 7 | energy identity / passivity | `ENG01`, `ENG02` |
| 8 | analytical limiting cases | `LIM01`, `LIM02` |
| 9 | scattering validation | `SCT01`, `SCT02` |
| 10 | numerical convergence | `CONV01` |
| 11 | reproducibility / artifact integrity | `ART01`, `ART02`, `REP01`, `REP02`, `BIND08` |
| 12 | research-claim traceability | `CLM01`–`CLM03`, `EQ01`, `BIND01`–`BIND07`, `DEP01`, `DEP02` |

### 5.3 `code/foundation_checks.py` is integrated, not discarded

`REP01` re-runs the entire Phase-1B suite in a clean subprocess and compares its
check count, status tally and executed-case total against the committed
`validation/phase_01b_test_results.json`. The 16 checks remain the primary evidence
for the constitutive, determinant, flux and reduction results; the new suite adds
the oracles they lacked. Nothing was deleted to make the new numbers look better.

---

## 6. Task E — the mutation suite

**11 mutations, each in an isolated temporary copy, each with a control run.**

The suite refuses to count "any non-zero exit" as detection. Each mutation names
the detector that must trip, and — critically — **every detector is also run
against an unmutated copy.** A detector that fires there would "detect"
everything, so an `INVALID_DETECTOR` counts against the suite exactly like a
survivor.

That control run found **five real defects in my own first version**, all of which
would have produced false confidence:

| detector | what was wrong | lesson |
|---|---|---|
| M01 | mutated the *harness guard*, i.e. the detector itself | a degenerate mutation proves nothing; inject into the producer |
| M02 | swapped a sign in a symmetric off-diagonal block | that is invisible to a ξ-parity identity; the parity test catches ξ-parity errors only |
| M06 | asserted `C44=(C11−C12)/2` | **PZT-5H violates it by 5.0e7 Pa** — it is transversely, not sagittally, isotropic |
| M07 | tested the *premise* (PZT isn't isotropic) | true of the pristine tree too; detect the misuse, not the material |
| M09/M10 | read a report file, not the test's exit code | stopped and read the code |

The final state:

| id | mutation | detected by |
|---|---|---|
| M01 | zero executed cases reported as PASS | harness forces `UNRESOLVED`, pytest fails |
| M02 | wrong characteristic-matrix parity | `M(−ξ) = S M S` breaks on 1 entry |
| M03 | invalid heat-flux expression | Second-Law balance oracle (`D2`) |
| M04 | wrong polynomial coefficient ordering | relative polynomial residual = 1.000 |
| M05 | wrong incident/outgoing mode selection | 4 growing modes selected |
| M06 | wrong material units | cross-source `DAT01` comparison |
| M07 | isotropic benchmark on an anisotropic medium | probe present + closed form evaluated |
| M08 | incorrect interface sign | rank is insensitive; sign defect detected directly |
| M09 | missing mandatory evidence artifact | `ART01` |
| M10 | wrong source/data checksum | `ART02` |
| M11 | claim promoted to VERIFIED without execution | `BIND01` token binding |

**Detection rate 11/11 = 100%, no survivors, no invalid detectors, origin tree
unmodified** (`origin_untouched: true`, hashes of all 7 frozen inputs re-verified
before and after).

Note on M08: the mutation is genuinely detected, but the detection rests on the
sign defect being *present in the assembled operator*, not on a physical
discriminator. Rank is provably insensitive to that defect (it stays 8), and a
true welded-limit test (identical media must give R = 0) would be the right
detector — but it cannot be built until the coupled solver exists. This is
recorded rather than papered over: **M08 is the weakest of the eleven.**

---

## 7. Task M — the eight mandatory demonstrations, executed

`python3 validation/demonstrate.py` runs each one for real and records the command,
output, case count and exit code. The JSON is at
`validation/evidence/demonstration.json`.

| # | demonstration | command | result |
|---|---|---|---|
| 1 | a previously false PASS is rejected | `test_IFC01` + the old construction re-run | PASS: the old form's 50 000 samples produce a minimum of exactly 0.0 (the initialiser) and can never be negative; the new test's negative control must be detected |
| 2 | incorrect parity detected | `run_mutations.py --only M02` | `DETECTED`, 1 entry broken |
| 3 | invalid thermal flux rejected by an independent oracle | `run_mutations.py --only M03` | `DETECTED`, `D2` → `UNRESOLVED` |
| 4 | an altered equation invalidates dependents | `impact_analysis.py --changed EQ-GOV-08` | 13 nodes reachable |
| 5 | an unrelated doc change does not | `impact_analysis.py --changed README.md` | exit 2, "unknown node" |
| 6 | missing evidence prevents promotion | `_demo_missing_evidence.py` | `PROMOTION BLOCKED`, exit 1 |
| 7 | a correct run is reproducible | `run_fast.py` twice | identical tallies |
| 8 | every VERIFIED claim has an inspectable basis | `test_BIND01` | PASS: all 36 VERIFIED claims resolve to a function containing the asserted tokens |

---

## 8. Tasks F, I, J — manifests, baselines, manuscript traceability

**Manifests.** Three experiments have a full manifest: `EXP-ELASTIC-REF`
(`PASS`), `EXP-PARITY-FLUX` (`PASS`), `EXP-FOUND-01B` (`UNRESOLVED`, by design —
it records `process_exit_code: 1` honestly). `validation/regen_manifests.py`
re-runs each recorded command and compares output hashes; it **verifies by
default and refuses to overwrite**. Result: 3/3 manifests reproducible. Re-issuing
a manifest requires `--rewrite`, which archives the old one and appends to
`MANIFEST_CHANGE_LOG.json`.

**Baselines.** Three, append-only. `BASELINE-01B-1` records the Phase-1 rejected
state, with the exact SHA marked **unrecoverable** rather than guessed — the parity
and flux corrections were committed together with the Phase 1B report, so the
pre-repair tree does not exist in this history. `BASELINE-01B-2` records the
repaired state at `127eba77…`, `CANDIDATE`, never authorised. `BASELINE-01C-1` is
this phase, also `CANDIDATE`.

Promotion requires **both** zero non-passing cases **and** no open
promotion-blocking blocker. Today it is refused on both counts.

**Manuscript.** No manuscript text is written. Six stable evidence ids exist
(`MS-EQ-001…003`, `MS-FIG-001`, `MS-RES-001`, `MS-CON-001`), each bound to
equation and claim ids. `MS-CON-001` ("the framework is validated") is registered
`FORBIDDEN` with seven blockers attached, so that if anyone drafts that sentence
in Phase 2 the map already says it is unsupported.

---

## 9. Task K — what has *not* been closed

This is the part to read sceptically. The suite is green in the places where the
science is sound, and **not one scientific gate moved**.

| blocker | severity | status |
|---|---|---|
| BLK-U3 amplitude normalisation | P0 | open |
| BLK-G4 LS energy identity | P0 | open |
| BLK-G5 interface / numerical balance | P0 | open |
| BLK-G0 full-text prior art | P0 | open |
| BLK-G1 sourced constants | P0 | open |
| BLK-U5 half-factor bookkeeping | P1 | open |
| **BLK-LS-STORAGE** | P1 | **open, preserved as a scientific question** |
| BLK-DOC-01 stale `.tex` section | P1 | open, non-blocking |

### The Lord–Shulman stored-energy question

Phase 1 called `τ_q q·∂_t q/κ` a "relaxation storage" term, then a dissipation.
For a monochromatic field the cycle average of `q·∂_t q` is **identically zero** —
measured `3.27e-6` against a `1e34` dissipation scale. A term that vanishes
identically is neither storage nor dissipation *in the harmonic regime*.

What that does **not** establish: it says nothing about the internal energy,
nothing about the transient regime, and nothing about whether LS has a genuine
out-of-phase stored-energy contribution in some other form. **A vanishing
cycle-averaged time derivative and a vanishing stored-energy term are different
statements**, and the second does not follow from the first. Settling it needs a
complete thermodynamic derivation, not a numerical test. It is preserved as a
question and no test in this phase pretends to answer it.

### Two new defects found while building the infrastructure

1. **The code did not read its own data table.** `build_materials()` hard-codes
   every constant in Python; `literature/material_constants.csv` was decorative.
   Consequence: an `ASSUMED` label could not influence or be caught by any
   computation, and a 1000× unit slip in the table was invisible to the whole
   suite. Fixed by `test_DAT01` (72 constants compared cross-source), which
   immediately found `Al6061.beta` differing between the two by 2.3e-4.
2. **`tauq` was still justified by the withdrawn rule.** Three of four `DERIVED`
   rows sourced the relaxation time to "the causal floor from check H1" — the rule
   withdrawn in Phase 1B. Values unchanged; provenance corrected, because leaving
   them would have preserved the withdrawn reasoning inside the data file.

---

## 10. Tasks G, H, N — regression levels, CI, delivery

**Three levels.** Fast (every commit, ~11 s, 35 cases), Standard (pre-submission,
adds the mutation suite and the full foundation suite), Full (milestone, adds
manifest regeneration, sensitivity, precision study, figure rebuild). Plus targeted
dependency-chain mode (`impact_analysis.py`), which does **not** substitute for the
Level 1 gate.

**Exit codes distinguish three scientific states,** which was a deliberate choice:
`0` all pass, `1` something is *wrong*, `2` something is *unknown*, `3` nothing
executed. Collapsing 1 and 2 would let "we have no evidence" masquerade as a
failure to investigate, and `--promotion` adds a fourth: promotion refused while a
blocker is open, even with zero failures.

**CI — observed, and it took three attempts.** `.github/workflows/research-validation.yml`
pins every dependency. The remote results, read from actual workflow runs rather
than assumed:

| run | result | cause | fix |
|---|---|---|---|
| 41818228896 | **failure** | `actions/setup-python` with `cache: pip` looks for `requirements.txt`/`pyproject.toml`; this project has only `requirements-research.txt`, so the step errored before anything ran | explicit `cache-dependency-path`; `3.11` instead of the exact `3.11.9`, which is not guaranteed to exist as a prebuilt toolcache image |
| 41818366238 | **failure** | setup-python fixed; the Level 1 step then failed because `run_fast.py` correctly exits 2 for "unknown" and a CI step reads any non-zero exit as failure. With `SCT02` permanently UNRESOLVED, the build would have been red forever | **not** suppressed. Added `validation/known_open.json` (enumerated tolerated cases, each with a reason and an owning blocker) plus `validation/ci_gate.py` |
| 41818854579 | **success** | — | all 15 steps green, 53 s |

The gate is the interesting part, because the obvious fix — making unknowns exit
0 — would have destroyed the property the whole phase is about. Instead:

- any `FAIL` fails the build, always;
- any `UNRESOLVED`/`WITHDRAWN` **not** pre-registered fails the build;
- a registered case that did not run fails the build (an invisible unknown is not
  tolerated, it is hidden);
- a register entry whose blocker has since closed fails the build, so the register
  cannot outlive its justification.

Verified in both directions with injected controls: real FAIL → exit 1,
unregistered unknown → exit 2, registered test absent → exit 2, clean report →
exit 0. A tolerated unknown is still reported as `UNRESOLVED` everywhere in every
report; it is tolerated by decision, never converted into a pass.

**Caveat, stated rather than glossed:** the uploaded CI artifact could not be
downloaded inside this sandbox (egress restriction), so the remote *tallies* are
inferred from the gate's contract — the gate step passed, which is only possible
with zero failures and every unknown registered — rather than read out of the
artifact. The job and all fifteen steps are confirmed directly from the run view.

---

## 11. What I would tell a reviewer

The infrastructure works and I trust it in one specific sense: it now fails when
it should, which is the only property that distinguishes validation from theatre.
The mutation suite proves that by attacking itself, and the control runs caught
five defects in my own work.

It does **not** mean the physics is accepted. Seven P0/P1 blockers are open, no
coupled scattering coefficient has ever been computed, and two statements in the
committed theory documents are false. The registry says so in the same file that
says 27 claims are strongly verified, which is the point: the evidence system is
only useful if it is allowed to be unflattering.

The single most valuable thing this phase produced is not a passing test. It is
`evidence/claim_evidence_map.json` — the mechanism that downgrades a claim
automatically when the code supporting it stops executing. That is what was
missing in Phase 1, and its absence is why a 34-orders-of-magnitude error cleared
review.

---

## 12. Reproducing this phase

```bash
pip install -r requirements-research.txt

python3 evidence/impact_analysis.py --check        # graph integrity
python3 validation/freeze_inputs.py                # only if inputs changed deliberately
python3 validation/run_fast.py                     # Level 1
python3 validation/run_standard.py                 # Level 2, includes mutations
python3 validation/run_full.py                     # Level 3, includes manifests
python3 validation/mutation_tests/run_mutations.py # mutations alone
python3 validation/regen_manifests.py              # verify manifests, never overwrite
python3 validation/demonstrate.py                  # the 8 demonstrations
```

---

## Appendix — machine-readable results at the time of writing

```json
{
  "tests":            {"PASS": 35, "FAIL": 0, "UNRESOLVED": 1, "WITHDRAWN": 0,
                       "NOT_APPLICABLE": 0, "total": 36, "executed_cases": 9292},
  "foundation_1b":    {"PASS": 15, "UNRESOLVED": 1, "total": 16, "executed_cases": 151},
  "mutations":        {"detected": 11, "total": 11, "survivors": [], "invalid_detectors": []},
  "demonstrations":   {"satisfied": 9, "total": 9, "unresolved": []},
  "claims":           {"VERIFIED_STRONG": 27, "VERIFIED_WEAK": 9, "UNVERIFIED": 5,
                       "UNRESOLVED": 1, "INVALID": 2},
  "equations":        {"VERIFIED_STRONG": 17, "VERIFIED_WEAK": 7, "UNRESOLVED": 9,
                       "DEFINED_NOT_TESTED": 3, "PARTIAL": 2, "INVALID": 3},
  "manifests":        {"reproducible": 3, "total": 3},
  "promotion_allowed": false,
  "ci":               "VERIFIED_GREEN"
}
```

The single `UNRESOLVED` test is `SCT02`, which exists to state that no coupled
scattering solver is implemented. It is not a failure and it is not a pass; the
harness refuses to let a zero-execution case be recorded as either.
