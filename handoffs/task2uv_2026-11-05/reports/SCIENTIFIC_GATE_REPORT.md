# Phase 2A Task 2U-V — Adversarial Constitutive / Mutation Scientific Gate

**Date:** 2026-11-05  
**Repository baseline (untouched):** `8713aea631d58bbeff608012ae4f20bfd2e9bdfd`  
**GitHub push:** NONE  
**Verdict:** **CONDITIONAL scientific milestone; NOT production-ready**.

## Motivation

Claude's second referee report showed previous `negative_controls.py` tests could pass tautologically and showed that comparing a mutated finite-k operator to its own mutated k=0 limit cannot validate the original constitutive laws. This milestone includes an independent *analytic* 1D Lord–Shulman thermoelastic benchmark plus **actual altered calculations** with measurements against a fixed pristine reference.

## Actual executed results

From `evidence/analytic_oracle_results.json`:

- Two arbitrary-precision runs, 40 and 70 digits; **10/10 analytic checks** at each precision: one classical no-coupling elastic-interface reflection benchmark and nine LS dispersion comparisons (three material models times synthetic ωτ=10⁻⁵,0.1,1).
- **8/8 perturbation detectors** per precision in a reduced 1D physical interface / spectral test harness. The three material perturbations were applied to the *input material coefficients*, yielding new interface fields, not just comparison with a self-mutated limit.
- Recorded 40-digit relative transmission-temperature changes for PZT-5H | BaTiO3: p3 sign 0.0397925; beta3 sign 0.890042; e33 sign 0.0397925. The equality of the p3/e33 sensitivities is consistent with beta* depending on their product; this is not two independent constitutive sign anchors.
- Mutating the LS `g=1-iωτ` placement shifts the root set by 0.1 and 0.999829 at synthetic ωτ=0.1 and 1; these are **algorithmic stress tests, not physical relaxation-time measurements**.
- An actual 1D assembled 4-row contact system is re-solved with the incident RHS replaced by zero; recomputed original physical interface jumps detect the omission.
- Reversing the left outgoing selection in the 1D solve admits the incident channel as outgoing and yields different reflected fields; the detector checks both the wrong channel and the changed observable.
- At exact normal incidence, dropping SV traction from the *physical two-row shear system* leaves the equation underdetermined; a concrete remaining-row-satisfying candidate violates omitted traction by 0.15562 at the recorded normalized scale. **This last item is a structural counterexample, not a unique result of a square mutant solve.**

From `evidence/finitek_mutation_results.json` and `_70.json`:

- At 40 and 70 digits, **4/4 injected 8×8 interface-row replacements** detected for PZT-5H|BaTiO3 P incidence, `k/k_ref=10^-4`.
- Explicit altered solves replace respectively physical phi, Dz, ux and shear traction continuity by a coefficient constraint; after solving the *mutant* 8×8 system, original removed-row **uncancelled-term-normalized residuals** are approximately 1.0 while the remaining physical interface rows close. These demonstrate that the imposed rows are nonredundant for this full finite-k configuration.
- The 8×8 source comes from the previous **independent experimental** bulk state-space implementation; it is not the repository production `assembly.py`. The row-constraint replacement is an adversarial model mutation, not the exact behavior of removing a row from an underdetermined 7×8 system.

From `evidence/root_separation_certificate.json`:

- **27/27 sampled spectra** have separate longitudinal/thermal roots and distinct SV wave at k=0 (three models × three frequencies × three tau conditions).
- Minimum sampled normalized P/thermal root separation ~1.54704e-4; minimum SV-to-other separation ~8.00250e-5. These are positive but do not prove simplicity outside sampled regimes or k>0 analytic continuation.

## Independent analytic physics anchors

1. For zero thermo/piezo coupling, displacement reflection is `(Z_L-Z_R)/(Z_L+Z_R)` and transmission is `2 Z_L/(Z_L+Z_R)` with `Z=√(rho C33)`. An independently built 2×2 displacement/stress system matches to ~1e-41 at 40-digit precision.
2. Isotropic/1D LS thermoelastic dispersion, from the two governing PDEs, is
   `(C λ²+ρ ω²)(κ λ²/g+iωρc)+iω T0 β² λ²=0` with `g=1-iωτq`. For `D_z=0` piezoelectric contact, the identical algebraic structure uses `C*, β*, R*` as derived in the accompanying lemma. The quartic is independently checked against a direct 1D first-order differential operator.

## Threats to validity and what we do NOT claim

- 1D analytic LS oracle **does not verify** the entire 8-state finite-k operator, nor the original research production scattering path.
- Material `p3`, `beta1`, `beta3` are marked ASSUMED; their sensitivity is not proof of their sourced values.
- Al6061 is modeled electrically as an epsilon0 *dielectric surrogate*, not a conductor or electrode.
- The earlier negative-control suite is **not rehabilitated** by these new tests; it should be replaced.
- Finite-k PZT-5H|BaTiO3 uses sample CSV parameters whose provenance must be audited before material-specific publication claims.
- The analytic lemma proves uniform invertibility **conditionally**; the sampled frequency certificate does not prove the hypotheses for all frequencies or material combinations.
- Independent coupled energy dissipation, G4/G5 and SCT02 remain unresolved.
- BLK-02A-26 (production defect) and BLK-02A-27 (outgoing and gauge admissibility acceptance) stay OPEN.
- The PZT/BaTiO3 sign comparisons do not independently verify Maxwell reciprocity or boundary conventions for a real conducting electrode.

## Prioritized next action after Claude review

1. Referee checks the analytic LS quartic and `g` placement against the PDEs; inspect the full mutation harness for hidden tautologies or algebraic re-use.
2. Verify radiation selection and analytic continuation across the relevant k–frequency domain; consider independently deriving nonzero-k closed-form special limits.
3. If accepted, implement the **reduced 6×6 exact-axis production repair**, with phi relative-gauge postprocessing and correct no-charge dielectric assumptions. Test `ACPT02` using invariant physical observables, rather than an arbitrary raw potential magnitude.
4. Add production-path mutations and CI, then adjudicate BLK-02A-26/27. Do not start G4/G5 until the scattering baseline is accepted.

**STOP:** this archive is a review-only scientific checkpoint. No production files changed; do not promote to baseline based on its check counts alone.
