# CLAUDE — Task 2U Finite-k / Limiting-Absorption Milestone: Adversarial Referee Report
5 Nov 2026 · Package `6cd91c57-…FINITE_K_FOR_CLAUDE.zip` + two copies of `SCIENTIFIC_MILESTONE_REPORT.md`. Repository not modified, nothing pushed, no new phase.

## 0. Separate verdicts

| Component | Verdict |
|---|---|
| Bulk 8-state operator | **PASS** (executed independent re-derivation). Physical constants/conventions remain model assumptions (§4). |
| k>0 radiation selection | **CONDITIONAL PASS.** Correct for the tested materials at one frequency; rule is brittle and not flux-based; LS factor g is effectively untested. |
| Limiting-absorption (LAP) argument | **CONDITIONAL PASS, upgraded from my earlier state.** The finite-dimensional block-triangular argument is sound and numerically supported; it is not yet a written proof (3 explicit hypotheses, §3). |
| Exact-axis 6×6 (+ gauge constant) | **PASS as the k→0⁺ (grazing-angle) limit of the mechanical/thermal fields**, conditional on LAP hypotheses H1–H3 and on the dielectric-surrogate reading of Al6061. |
| Production readiness | **FAIL.** Unchanged: production `assembly.py` untouched, BLK-02A-26/27 open, SCT02/G4/G5 open. |

Prior P0-1 (wrong D0-exclusion argument): **substantively answered** by the milestone, pending the written lemma. Prior P0-2 (Al6061 not a conductor): **acknowledged and correctly disclosed**; remains a physical-interpretation limit.

## 1. Evidence ledger — what I actually did

**Executed by me**
1. SHA-256 check of all package files: all OK.
2. Read all six scripts, `material_constants.csv`, both prior-review files, `LAP_CERTIFICATE.json`, `production_crosscheck.json`, the milestone report.
3. **Re-ran** `independent_finite_k_state_space.py` (reproduces eigen-residuals ~1e-42/1e-72; small roots/k = ±1.0 (Al) and ±1.0744 (PZT)).
4. **Re-ran** `independent_scattering_limit.py` at 40 digits (k/k_ref = 1e-2,1e-4,1e-6; 24 cases, 0 exceptions). All 24 `max_mech_thermal_diff_fixed_scale` values equal the shipped JSON (relative difference 0). I did not re-run the 70-digit set or `validate_lap_evidence.py`/`crosscheck_production_characteristic.py`.
5. **Re-ran** `negative_controls.py` (PASS at 40 and 70 digits) — and found it largely tautological (§5).
6. **My own independent bulk derivation** (`scratchpad/mychecks.py`): coded the 4×4 plane-wave matrix M(k,λ) directly from the second-order PDEs (stress, Gauss, LS balance), without the first-order reduction. For both materials and k/k_ref = 1e-2, 1e-3, 0.3, **every one of the 8 eigenvalues of the package operator makes my M singular** (min/max singular value ≤ 7e-30 at 40 digits). Conclusion: the first-order reduction is algebraically correct.
7. **Thermodynamic sign check** (derivation, not numerics): with σ = Cε − eE − βθ and D = eε + εE + pθ, the Maxwell relations require ρη = βε + pE + ρcθ/T0; so the LS source T0(β ε̇ + p Ė), with E_z = −φ′, i.e. `-p3*php`, is the thermodynamically consistent sign. The package's heat row `q_z' = −κ11k²θ/g + iω[ρcθ + T0(ikβ1ux + β3uz' − pφ')]` with q = −κθ′/g is exactly LS: (1−iωτ)q = −κ∇θ, ∇·q = iωS.
8. **Electrostatic pair, independent closed form.** For λ→0 inertia dominates, so the electrostatic mode is mechanically clamped (u_z ≈ O(λ), strain ≈ 0) and thermally adiabatic (S = 0). That gives ε33,eff = ε33 − T0p²/(ρc_E), ε11,eff = ε11, s_e = √(ε11/ε33,eff), ε_eff = √(ε11 ε33,eff). I computed: Al6061 s_e = 1.0, ε_eff = 8.8541878128e-12; **PZT-5H s_e = 1.074422244, ε_eff = 1.39889136589e-8 — identical to the package values.** (My first attempt, stress-free quasi-statics, gave 1.246 and 2.20e-8: wrong, because it ignored the inertial term; recorded as my error, corrected.) The scaled 2×2 determinant (ε_L+ε_R)/ε_ref reproduces −2 and −1.00063 (= 1 + ε_Al/ε_PZT).
9. **Radiation direction by energy flux and group velocity** (`scratchpad/flux.py`): for Al6061 and PZT-5H at k/k_ref = 0, 1e-3, 0.3, for every propagating root (P, SV) the sign of Im λ, of the group velocity dω/dκ (finite difference in ω), and of the mechanical+electrical time-averaged flux −½Re(σ·v*) + ½Re(φ(−iωD_z)*) all agree. The diffusive LS thermal roots (Re λ ≈ −Im λ) carry thermal flux ½Re(θq*)/T0 with the sign of Im λ and decay in the direction of travel. Electrostatic roots ±s_e k are real and reactive (flux ~1e-24): decay sign is the only valid selector, and it is unambiguous.
10. **Component-wise convergence with natural scales** (not the package's arbitrary fixed scales): e.g. PZT|Al P: |Δu_z|/|u_z| = 1.4e-8 → 1.4e-10, |Δθ|/|θ| = 1.6e-8 → 1.3e-10 (k/k_ref 1e-4 → 1e-6): slope 1.0, small constants. Quantities that vanish at k=0 (u_x, τ_xz for P; u_z, θ, τ_zz, q_z for SV) are genuinely O(k) amplitudes.
11. **Mutation test on the independent path**: flipping sign of p3 changes the exact-axis transmitted fields by 1.0 %; flipping β1,β3 by 201 % (i.e. factor 3) — yet the package's finite-k-vs-exact-axis convergence is **unchanged to all printed digits (8.68e-9, 8.66e-11)** in both mutants.

**Derived by hand, not computed** [D]: the Jordan-block/Schur-complement structure of the (φ,D_z) subsystem; the O(k²) form of the Jordan perturbation; the LAP block structure; remarks on energy flux.

**Not checked**: 70-digit and 100-digit re-runs; frequency other than 1 MHz; frequencies near 1/τ_q; full field profiles in z; prior small-k datasets from the first package; baseline repository beyond what was reviewed earlier; the production cross-check script.

## 2. Mathematically supported results

**(M1) Operator.** The 8-state operator is the correct first-order form of the stated model (executed, item 6).

**(M2) Electrostatic pair.** λ_e = ±s_e k + o(k), s_e = √(ε11/ε33,eff); D_z/(kφ) → ∓ε_eff (derived + executed, item 8). The near-zero pair arises from perturbing the 2×2 Jordan block [[0, −1/ε33,eff],[−k²ε11,0]] (D_z′ = −k²ε11φ, φ′ = −D_z/ε33,eff), whose perturbation is exactly O(k²); hence the eigenvalues ±s_e k are analytic in k and the potential-normalized eigenvectors are analytic.

**(M3) Block structure ⇒ LAP.** Order the scattering system with columns [6 mechanical/thermal | 2 electrostatic] and rows [6 mech/thermal | φ | D_z/k]. (i) The 6 mechanical/thermal roots are simple at k=0 (P4 has four distinct roots, plus ±shear), so their eigenvalues/eigenvectors are analytic in k and converge to the exact-axis ones. (ii) Electrostatic columns have mechanical/thermal/traction entries O(k) and φ = 1 (potential-normalized). (iii) After dividing the D_z row by k, mechanical columns contribute D_z^{(1)}/k = O(1) and electrostatic columns contribute ∓ε_eff(1+o(1)). Hence the limit matrix is block lower-triangular with diagonal blocks (a) the exact-axis 6×6 and (b) E0 = [[1,−1],[−ε_L,−ε_R]], det = −(ε_L+ε_R) ≠ 0 **for any positive permittivities**. Both blocks are invertible ⇒ the inverse is uniformly bounded in a punctured neighbourhood of k=0 ⇒ electrostatic amplitudes are O(1), D_z = O(k), and the mechanical/thermal amplitudes converge to the exact-axis solution. A nonzero D0 would require unbounded electrostatic amplitude ~1/k, which is excluded by this bound. This is exactly the argument I requested in P0-1 and it is mathematically correct given H1–H3.

**(M4) Interpretation.** The selector is a k→0⁺ (grazing-incidence, k = ω sinθ/c) regularization of a genuinely non-unique exact-axis problem; it is not a literal ω→ω+i0 limiting-absorption principle. I recommend calling it "oblique-incidence/k→0⁺ limit" in the manuscript.

## 3. Unresolved scientific assumptions (hypotheses H1–H3)

- **H1 (simple mechanical roots).** The six mechanical/thermal roots are simple at k=0 and stay on their outgoing branch; verified only at 1 MHz for two materials. Needs a frequency sweep including ωτ_q ~ 1 where the thermal mode may become wave-like.
- **H2 (Jordan perturbation order).** The O(k²) claim for the (φ,D_z) block is derived by me for the model with ε11 and ε33,eff; the written lemma should give it explicitly including shear and thermal terms.
- **H3 (uniform boundedness).** Numerics sample 3 k-values per precision; boundedness is proved only via M3 for the limiting matrix, not uniformly for all k in (0,k0). A one-line perturbation bound (Neumann series) is enough, but must be written.
- **H4 (physical):** the regularization selected (k→0⁺) is the physically relevant one for the intended experiment (finite lateral aperture). This is a modelling statement, not a theorem.
- **H5 (interface electrical law):** ideal dielectric contact, [φ]=0, [D_z]=0, no free charge, no dipole layer.

## 4. Defects and findings (with locations)

### P0
- **None new that invalidates the exact-axis formulation.** The previously open P0-1 is answered in substance by M3; closure requires only the written lemma (W1 below). **Al6061 as ε0 dielectric (`material_constants.csv`, Al6061 eps rows, note "the metal's own permittivity is infinite; vacuum value used…") is not a conductor. This does not invalidate the 6×6 mechanical/thermal results, but it does invalidate any manuscript statement of "metal|piezo electrical contact".** Treat the Al6061|PZT-5H results as "elastic non-piezoelectric dielectric | piezoelectric". A conductor needs a different electrical sector (φ = const, E = 0 in the metal, free surface charge σ_f = [D_z]); my expectation is that the mechanical/thermal 6×6 is unchanged but C, charge and the validity of the D0=0 argument on the piezo side must be re-derived.

### P1
- **P1-1 `negative_controls.py` is largely tautological (executed + read).**
  - L137: `shear_detected=all((7!=B.cols) for _ in [4])` evaluates `7 != 8` and is always `True`; no row is ever removed.
  - L139-140: `reversed_R=[m for m in m2 if m['left_out']]; all(not m['right_out'] …)`: `left_out` and `right_out` are mutually exclusive by construction (`independent_scattering_limit.py` L47-48), so this is always `True` (and vacuously so if the list is empty).
  - L129-134: "omit incident", "wrong φ RHS", "wrong D_z RHS" solve a square nonsingular system with a perturbed RHS and then compare against the original RHS: the residual equals the injected perturbation by linear algebra. They test nothing about physics or radiation.
  Replace by constitutive and radiation mutations whose effect is measured on independent observables (items 11, §6).
- **P1-2 Finite-k-vs-exact-axis convergence is blind to constitutive errors (executed, item 11).** Both paths use the same `operator()`. Flipping p3 (1 % effect on transmitted fields via β* = β3 − e33p/ε33) or β (200 %) leaves the convergence figures unchanged. Convergence therefore certifies the limiting procedure, not the physics. The only physics anchors are (a) the cross-check against production M, which shares the same authors' conventions, (b) my thermodynamic consistency argument (item 7). Needed: a published/analytic benchmark (e.g. isotropic LS thermoelastic half-space dispersion, closed form; PZT limiting case with κ→∞/τ→0), and a symmetry/reciprocity test.
- **P1-3 The Lord–Shulman factor g is practically unexercised.** `tauq` in the CSV is 3.7e-13 (PZT), 2.1e-12 (Al) s; at ω = 2π·10⁶ rad/s, ωτ ≈ 2e-6…1.3e-5, so g−1 ~ 1e-6. A wrong g placement (κ/g vs κg, g on the source) would change results at ~1e-6, below every reported tolerance except the slope tests (which are blind, P1-2). Test with an artificially large ωτ (e.g. τ×10⁶) against the analytic LS dispersion relation of an isotropic thermoelastic solid.
- **P1-4 Pyroelectric p3 and β are `ASSUMED` in the CSV (PZT-5H β1, β3, p3; BaTiO3 likewise).** Their effect on exact-axis results is 1 % (p) and O(1) (β). Any quantitative R/T claim carries this uncertainty. Also `cE` must be specific heat at constant strain and field (c_ε) for the LS energy equation; the CSV does not state this (Al6061 handbook c_p vs c_v differs by a few %).
- **P1-5 Half of the scattering test matrix is trivial.** `PZT|PZT` and `Al|Al` are homogeneous (no scattering; reflected amplitudes ≡ 0, electrostatic amplitudes 0). The slope-2 D_z result for PZT|PZT is a bulk property. Only the two mixed pairs test the LAP block with nonzero electrostatic amplitudes. Add genuinely different piezo|piezo pairs (PZT-5H|BaTiO3, rotated poling) to exercise E0 with two non-trivial permittivities and e15/e31 couplings.
- **P1-6 Fixed convergence scales are arbitrary** (`independent_scattering_limit.py` L145: θ 1e6, τ 1e14, q 1e9, u 1). They divide O(k)-vanishing quantities (e.g. q_z in SV, whose exact-axis value is ~1e-29) by an unrelated scale; the headline "9.61e-6 at 1e-8" is dominated by q_z of SV, an O(k) amplitude whose prefactor is ~680 relative to that scale. Report per-field errors against the natural magnitude (item 10), and report O(k) coefficients separately.
- **P1-7 Radiation selector is a hard-coded phase/decay rule** (`group`, L43-50; threshold `|Re λ/Im λ| < 1e-4` and `|Im λ| > 1e-3 k_ref`; `electric` flag `|a| ≤ 5k/k_ref` L23; exact-axis zero filter `|a| > 1e-15` L86-87). Correct here (items 9), but: (a) a lossy mode with a ratio above 1e-4 is switched to the decay rule (still correct for passive media, but the discontinuity is undocumented); (b) incident wave search `incident()` L52-57 raises if an attenuated P exceeds the threshold; (c) the absolute `1e-15` zero-eigenvalue filter relies on defective-eigenvalue perturbation ~√eps (1e-20 at 40 digits) and breaks if working precision is lowered or the operator is rescaled; (d) only one frequency tested. Recommended: classify by Re λ alone when Re λ ≠ 0 at working precision, use flux/group velocity only for the lossless Re λ = 0 case, and assert agreement of both.

### P2
- `crosscheck_production_characteristic.py` L192-196: the revised error denominator uses amplitude × physical-scale products (`amp*phys_scales[j]`) rather than the actual |M_ij v_j| terms; it can inflate denominators when components are far below their nominal scale, and the pass threshold (1e-16) is much looser than what is achieved (1e-33). My inspection: the Al failure (raw metric 1.0) is a legitimate cancellation artefact in the decoupled ε0(k²+ξ²)φ row; PZT's raw metric already passes at 1e-21…1e-25, so the coupled constitutive structure is not masked. Preferred metric: Σ|M_ij v_j| (term-wise backward error), keep the raw metric visible (it is, which is good practice). Note also that this test shows production M and the independent operator agree; they share the same stated model.
- `independent_scattering_limit.py` L103-107: exact-axis φ constant `C = fL[2]-fR[2]; fR[2]+=C` is the gauge row by construction; `phi_jump_over_fixed_scale` in finite-k rows (L152) is a residual of an imposed equation, not a validation of [φ]=0. State this in the report (it already warns "do not conflate common offsets with jumps").
- `validate_lap_evidence.py` L34-38: acceptance slope window 0.85–1.15, D_z window 0.8–2.2 and amplitude cap 100 are loose; tighten to the observed ±1e-6 and justify the cap physically (amplitudes are O(1) in units of 1e9 V).
- Equal-material scattered amplitudes of exactly 0 are reported as bounded electrostatic amplitudes (certificate shows 0.0); exclude from LAP evidence.
- 40- vs 70-digit comparison of convergence diagnostics (`cross_precision`) is a determinism check, not a correctness check.
- Missing script referenced in prior reconciliation (`independent_jordan_oracle.py`, `independent_jordan_results.json`) is not in this bundle; the Jordan-chain claim (M2/§2) cannot be rerun from the package. I reproduced the structure myself (toy in the previous round; here via the eigenvalue pair ±s_e k).

## 5. Answers to the eight mandatory tests
1. **Operator derivation**: done independently (items 6, 7); correct including pyroelectric sign and LS source.
2. **Electrostatic pair and determinant**: verified independently in closed form (item 8). Note ε_eff is the clamped-adiabatic value, not the stress-free value.
3. **Block-triangular inverse / mode switching**: argument valid (M3); mode switching is excluded for simple roots by analyticity; the one delicate branch is the coalescing pair, which is analytic in k (M2). Numerical samples alone would not suffice; the structure does.
4. **Radiation direction**: verified by flux and group velocity for P, SV, thermal, electrostatic (item 9). Phase-sign rule adequate here; see P1-7.
5. **Recompute 40/70/100 and field profiles**: 40-digit set reproduced exactly; 70/100-digit and z-profiles **not** done. Common potential offsets are not conflated in my comparisons (φ excluded; I used u, θ, σ, q).
6. **Al6061 surrogate**: see P0 note. It does not invalidate the mechanical/thermal 6×6 results nor the LAP block argument (it only needs ε>0), but it does invalidate a conductor-contact interpretation.
7. **Failed raw residual**: legitimate cancellation artefact, not a hidden constitutive discrepancy (P2).
8. **Negative controls**: P1-1 and P1-2; production-path mutation list in §6.

## 6. Prioritised work plan (necessary vs optional)

**Necessary (before closing BLK-02A-26/27 at the formulation level)**
- W1 (P0→close). Write the LAP lemma M2–M3 with H1–H3 proved (Neumann-series bound for the scaled inverse, Jordan perturbation order). *Acceptance:* a short proof in the report, checked by a numerical bound ‖S(k)⁻¹‖ ≤ K on k/k_ref ∈ [1e-10, 1e-1] for all four pairs plus two new piezo|piezo pairs.
- W2 (P1). Replace the negative controls by constitutive/radiation mutations measured on independent observables: flip p3, β, e33 sign; wrong g placement (with large ωτ); swap outgoing/incoming on one side; drop the D_z/k row; delete the φ row; multiply ε_eff by 2; omit the incident term. *Acceptance:* each mutation changes the compared field or conditioning by > 1e-3 (or fails a count assertion) relative to baseline, and the harness fails if a mutation is undetected.
- W3 (P1). External physics anchor: closed-form LS isotropic thermoelastic dispersion and normal-incidence reflection (no piezo), plus a published piezo benchmark. *Acceptance:* agreement to 1e-12 with the independent operator and to the stated literature tolerance.
- W4 (P1). Frequency sweep, including ωτ_q ~ 0.1, 1, to test H1 and the selector; add flux/group-velocity assertion inside `group()`.
- W5 (P1). Conductor-contact derivation (or remove "metal" language). *Acceptance:* either a derived BC set with D0 argument re-proved, or all manuscript text says "non-piezoelectric dielectric".
- W6 (P1). Natural-scale error reporting (P1-6) and piezo|piezo pairs (P1-5).

**Optional**
- 100-digit and z-profile reruns; tightened acceptance windows; term-wise backward error for the characteristic cross-check; ship the Jordan oracle.

**Production path (only after W1–W5):** keep the 6×6 mechanical/thermal exact-axis solve, treat φ constant as post-processing gauge; implement ACPT02 as "mechanical/thermal fields converge with slope 1, φ jump = 0 at finite k, D_z/k bounded". Do not integrate the experimental 7×7.

## 7. Blocker states
BLK-02A-26 **OPEN** (production code still defective; formulation question largely answered). BLK-02A-27 **OPEN, resolvable** (the extra "condition" is the gauge row; mode count 3+3 stands). SCT02 **UNRESOLVED**. G4 **OPEN**, G5 **OPEN**. No energy results validated; no R/T/A claims.

STOP — review complete.
