# Phase 2A Task 2U — Production scattering repair (targeted), 2026-11-05

**Branch:** `arena/1aa803b8-negative-reflection` (this work). Handoff branch
`handoff/task2uv-20261105` was merged in read-only fashion; `main` untouched.
**Verdict:** the known production defect BLK-02A-26 is repaired and ACPT02 now
PASSES on measured, gauge-aligned physical observables. BLK-02A-27, SCT02, G4,
G5 remain OPEN. **No energy-conservation, R/T/A, conductor-contact or
global-in-frequency claim is made.**

## 1. What was reproduced before any edit

Handoff evidence re-executed on this machine (CPython 3.11, mpmath 1.3.0,
pinned `requirements-research.txt`):

| Script | Result | Shipped evidence |
|---|---|---|
| `task2uv_analytic_oracle.py --digits 40 70` | 10/10 analytic checks, 8/8 mutation detectors at each precision | `analytic_oracle_results.json` — matches |
| `task2uv_finitek_mutations.py --digits 40 70` | 4/4 interface-row mutations detected at each precision | `finitek_mutation_results{,_70}.json` — matches |
| `task2uv_root_gap_certificate.py` | 27/27 sampled spectra separated (min P/thermal gap 1.54704e-4, SV gap 8.00250e-5) | `root_separation_certificate.json` — matches |
| `validation/run_fast.py` (baseline, pre-repair) | 104 PASS, 0 FAIL, 2 UNRESOLVED (SCT02, ACPT02) | matches Task 2T acceptance record |

## 2. The repair (production path only)

`code/scattering/assembly.py`:

* **Removed** the rule that deleted the `phi`/`Dz` rows at k = 0 when a
  half-space basis or the incident wave had no content in the row. That rule
  was BLK-02A-26: Al6061 has e = 0, so an incident wave in it has phi = 0
  exactly, and the rule fired for every metal|piezo pair.
* **At k = 0:** the `Dz` row is removed only after verifying the on-shell
  identity D_z = 0 (Gauss, i xi D_z = 0) per column against each column's own
  constitutive scale |e33 xi u_z| + |eps33 xi phi| + |p theta| (measured worst
  relative residual ~1e-62). The `phi` row is routed — whenever it carries
  content, independent of basis/incident emptiness — to the potential-gauge
  handler. The kept rows form the square 6x6 mechanical/thermal system
  (ux, uz, theta, txz, tzz, qz; 6 unknowns = 3 outgoing modes per side).
* **At k != 0:** nothing is removed; all eight ideal dielectric contact rows
  ([phi] = [Dz] = 0) are kept.

`code/scattering/gauge.py`: `physical_fields_reduced` now calls
`reconstruct_potential` at k = 0: phi_i = phi_i^modes + C_i with the relative
offset fixed by [[phi]] = 0 on the TOTAL fields (incident included on side 1),
C_2 - C_1 = phi_1^tot - phi_2^modes, and one global gauge fixed by C_1 = 0 and
recorded (`gauge_constants`).

## 3. Measured results after the repair (DPS = 60, W = 2 pi 1e6)

* All five pairs (Al6061|Al6061, PZT-5H|PZT-5H, Al6061|PZT-5H, PZT-5H|Al6061,
  PZT-5H|BaTiO3) x {P, SV}: k = 0 system is square 6x6, solve residual
  0 .. 7.8e-62, [[phi]]/FS = 0 (to 2.1e-70 worst) after reconstruction.
* ACPT02 PASS: gauge-aligned transmitted potential agrees with the k -> 0
  limit of the full eight-row finite-k system to < 1e-6; measured phi k-slope
  1.76 .. 2.00 over k/k_ac 1e-2..1e-7 (order ~2); all six mechanical/thermal
  components classified by measured k-order agree with the limit (CONSTANT
  components to < 1e-6 relative, VANISHING components to zero).
* ACPT07 PASS (new, `test_19_scattering_task2u.py`): finite-k 8x8 vs k=0 6x6
  per-component convergence with measured orders, all pairs, P and SV.
* ACPT08 PASS (new): every interface contact residual reported per field on
  its own physical scale; no unscaled mixed-SI norm anywhere.

## 4. Tests and mutations

* `validation/tests/test_18...` ACPT02 redefined to gauge-invariant /
  gauge-aligned observables (contact residual, per-component mech/thermal
  limit agreement, gauge-aligned phi and its k-slope) — PASS.
  ACPT03 updated to the gauge-invariant order measurement — PASS.
  test_15 SCTR07 updated to accept the Task 2U potential-gauge route — PASS.
* New mutations MR21–MR27 on the production path (phi reconstruction disabled,
  D_z row kept at k=0, shear row dropped, RHS sign flip, outgoing/incoming
  swap, LS g misplacement probed at synthetic omega*tau ~ 2): **7/7
  DETECTED**, each with a pristine-tree control that does not fire.
* Full Level-1 suite after repair: see `validation/evidence/last_run*.json` and
  the CI run on the pushed branch.

## 5. Blocker states

* **BLK-02A-26: CLOSED** (production defect repaired; ACPT02 PASS; registry
  entry removed from `validation/known_open.json`).
* **BLK-02A-27: OPEN.** The formulation-level answer is delivered — 3+3
  outgoing modes, 6x6, [[phi]]=0 as relative-gauge reconstruction — supported
  by the conditional normal-incidence lemma (H1–H3, not a global theorem) and
  by ACPT07's measured convergence. Formal closure awaits supervisor
  acceptance of the lemma's hypotheses.
* **SCT02, G4, G5: OPEN, unchanged.** No coupled energy balance, no R/T/A, no
  finite-interphase work was started.

## 6. Honest limits

* Al6061 remains a dielectric surrogate (eps ~ eps0), not a conductor; no
  metal-electrode contact claim.
* p3, beta1, beta3, cE ambiguities for PZT-5H/BaTiO3 remain ASSUMED in the
  material table; 40/70-digit agreement is arithmetic, not provenance.
* The k -> 0 route is the oblique-incidence (normal-incidence) limit under
  H1–H3; the uniform-D0 affine branch is excluded only as the bounded
  k -> 0+ selected family, exactly as the lemma states.
