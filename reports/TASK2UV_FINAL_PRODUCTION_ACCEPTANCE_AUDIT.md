# TASK 2U-V — Final production acceptance audit (evidence-only)

**Audited commit:** `8b3bb9d87da6cb18aa3a2b479a47aec81ad080fb` on
`arena/1aa803b8-negative-reflection`. Verified: local HEAD ==
`origin/arena/1aa803b8-negative-reflection` == that SHA, working tree clean.
`origin/main` = `07b5ae5c...` is an ancestor of HEAD — **main untouched**.
Handoff `handoff/task2uv-20261105` @ `1183d34` incorporated by merge; its
evidence was re-executed (oracle 10/10 checks + 8/8 detectors at 40 and 70
digits; finite-k mutations 4/4 at both; root certificate 27/27).
**CI run** [41994618431](https://github.com/vipin-oss/Negative-Reflection/actions/runs/41994618431):
conclusion **success**; verified steps include Level-1 suite, the
contradictions/unregistered-unknowns gate, manifest verification, mutation
spot-checks, graph integrity. Machine-readable evidence:
`validation/task2uv_final_acceptance_evidence.json`.

## 1. Files changed by the milestone

Production: `code/scattering/assembly.py`, `code/scattering/gauge.py`.
Tests/validators: `validation/tests/test_15_scattering_task2r.py`,
`validation/tests/test_18_scattering_task2t_accept.py`,
`validation/tests/test_19_scattering_task2u.py` (new),
`validation/mutation_tests/run_mutations.py`.
Registries: `reports/PHASE_01C_OPEN_BLOCKERS.json`, `validation/known_open.json`.
Evidence/report: `validation/evidence/last_run_1.json`,
`validation/mutation_tests/last_mutation_run.json`,
`reports/PHASE_02A_TASK2U_PRODUCTION_REPAIR.md`.

## 2. Exact-axis (k = 0) implementation — source-verified

| Requirement | Location | Measured |
|---|---|---|
| 6 mech/thermal rows kept | `assembly.py:472-508` route only `phi`; kept = ux,uz,theta,txz,tzz,qz | 6×6, rank 6, square, all 5 pairs × P/SV, dps 40 & 60 |
| Incident sign | `assembly.py:329` `rhs_full[i] = -st[i]/row_scale[i]`; side-2 negated `:321` | MR25 (sign flip) DETECTED |
| Outgoing basis | `assembly.py:80` `BASIS_ROLES=("OUTGOING",)` via `radiation.classify_role` | MR26 (swap) DETECTED |
| One relative constant | `gauge.py:317-334` `reconstruct_potential`: `c2 = f1[phi]-f2[phi]` on TOTAL fields (incident included, `gauge.py:305-311`) | C1=0, C2−C1 fixed, recorded in `gauge_constants` |
| [φ] from fields, not assumed | ACPT02/ACPT08 recompute `abs(f1tot[phi]-f2[phi])/FS` | 0.0 .. 2.1e-70 |
| One global gauge remains | C1=0 fixing; only C1+C2 free | stated + recorded |
| D_z excluded only on-shell | `assembly.py:433-470`: per-column check vs constitutive scale | worst rel 1e-62 ≤ 1e-25 |
| No electrostatic channel lost | at k>0 all 8 rows kept (verified k=1e-2,1e-3,1e-4: eq=8, `rows_gauge_only=[]`) | — |
| Rescale invariance | `solver.py:150` `rescale_invariance` | ok, worst 2.7e-51; manual mech check 0.0 |

Solve residuals: 0 .. 7.8e-62 (dps 60), 0 .. 1.1e-41 (dps 40), 10 cases.
Gauge-aligned φ vs k→0 limit: rel 1.09e-61 (dps 60) / 1.6e-41 (dps 40).

## 3. k ≠ 0 — verified independently

All eight ideal dielectric-contact rows remain at every nonzero k tested; no
incident/basis-emptiness heuristic fires (`potential_rows=[]`,
`rows_gauge_only=[]`). Production vs independent finite-k solver
(`asymptotic.dielectric_contact_system`) at k/k_ac=1e-4, dps=50, per field,
both sides: ux 1.8e-42, uz 2.5e-43, phi 3.4e-42, theta 3.8e-39, txz 5.8e-43,
tzz 2.3e-43, Dz 3.4e-42/5.8e-39, qz 4.0e-38. Convergence to the k=0 6×6 is
measured per component in ACPT07 (PASS), not by a single matrix residual.

## 4. ACPT02 / ACPT07 / ACPT08 review

All three call the production path (`gauge.reduced_system` →
`solve_path` → `physical_fields_reduced`) and compare against the
independent finite-k assembly. Tolerances (1e-6) are justified by the measured
order-~2 convergence; the old O(1) defect is far outside it. Old-defect
negative controls: **MR21** (incident-dependent φ removal re-injected) and
**MR22** (no reconstruction) both DETECTED — ACPT02 fails on the mutated
tree ("1 failed") and passes on the pristine control ("1 passed"). The two
comparison paths do not share the defect: MR21/MR22 mutate only production
code.

## 5. Mutations MR21–MR27 (from `last_mutation_run.json`, origin_untouched=true)

| ID | Mutated expression | Detector | Pristine | Mutated |
|---|---|---|---|---|
| MR21 | `assembly.py` phi branch gated on `inc` empty | ACPT02 | pass | **fail** |
| MR22 | `potential_rows: []` (no [[φ]]=0 reconstruction) | ACPT02 | pass | **fail** |
| MR23 | D_z on-shell row kept at k=0 (7×6) | ACPT02/ACPT07 | pass | **fail** |
| MR24 | txz row never kept | ACPT07/ACPT08 | pass | **fail** |
| MR25 | `rhs_full[i] = +st[i]/...` | ACPT07/ACPT08 | pass | **fail** |
| MR26 | `radiation.py` OUTGOING/INCOMING swapped | ACPT07/ACPT08 | pass | **fail** |
| MR27 | `foundation_checks.py` qz ∝ g instead of κ/g, synthetic ωτ≈2 | LS identity probe | 5.6e-52 | **7.13** |

Suite totals: 50 listed, **44/44 executed DETECTED**, 0 survivors, 0
inconclusive. MR01/02/05/07/08/09 are SUPERSEDED with reasons: the code they
mutated was removed by this repair, and each behavior is covered by MR22–MR24
(phi reconstruction, D_z row, shear row).

## 6. Registry integrity

BLK-02A-26: status CLOSED with non-empty `closure_evidence` (repair, tests,
mutations, residuals) and `closed_by`; THY04 PASS confirms every CLOSED
blocker is evidenced and G0/G1/G4/G5 remain unpassed. `known_open.json` now
lists only SCT02. BLK-02A-27 remains **OPEN** pending supervisor acceptance
of the conditional lemma. SCT02, G0, G1, G4, G5 unchanged. No previously
failed observation was removed: the ACPT02 UNRESOLVED measurement history is
preserved in `reports/PHASE_02A_TASK2T_FINAL_ACCEPTANCE.md` and git history.

## 7. Verdicts

- **A. Production implementation correctness: PASS**
- **B. BLK-02A-26 closure: ACCEPT**
- **C. Readiness to address BLK-02A-27: READY** (evidence delivered; formal
  closure requires supervisor acceptance of the H1–H3 conditional lemma)

Green CI was not treated as automatic acceptance; every verdict above rests
on the executed measurements cited. No energy, R/T/A, conductor-contact, or
new-phase work was performed.
