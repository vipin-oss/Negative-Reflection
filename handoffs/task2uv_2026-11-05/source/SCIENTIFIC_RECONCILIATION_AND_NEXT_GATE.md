# Phase 2A Task 2U — Referee Reconciliation and Next Scientific Gate

**Status:** Reviewer report assessed against available independent source scripts; no production edits, GitHub push, blocker closures, or energy claims.

**Review source:** CLAUDE_TASK2U_FINITE_K_REFEREE_REPORT.md (5 Nov 2026). **Baseline:** `8713aea631d58bbeff608012ae4f20bfd2e9bdfd`.

## A. Independent confirmations from actual source

1. `negative_controls.py`: `shear_detected=all((7!=B.cols) for _ in [4])` is a constant comparison. With an 8-column matrix it always returns true, regardless of physical shear behavior. It is not a mutation test.
2. `negative_controls.py`: `reversed_R=[m for m in m2 if m['left_out']]` and `all(not m['right_out'] for m in reversed_R)` only reassert mutually exclusive boolean labels created by `group`, not independent radiation/flux admissibility.
3. Wrong-RHS tests solve the unaltered square boundary matrix with an intentionally wrong RHS, then inspect against the old RHS. Detecting a difference is a necessary sanity check, but does not validate the constitutive model, modal selection, or interface physics.
4. `independent_scattering_limit.py`: finite-k scattering and exact-k=0 reference use the same `operator()` from `independent_finite_k_state_space.py`; convergence therefore checks consistency of the limit, not the external correctness of the constitutive equations.
5. `independent_scattering_limit.py` uses `C=fL[2]-fR[2]` followed by `fR[2]+=C` at k=0. Subsequent `phi_jump == 0` is imposed by construction, not a separate predictive test.
6. The available material CSV labels PZT-5H `beta1`, `beta3`, `p3` as **ASSUMED**, and its tau_q as a stated/derived working choice; reported quantitative outcomes must carry this provenance.

## B. Clarifications/corrections to the referee's wording

- **Angular terminology:** with `k=omega sin(theta)/c` and theta measured from the interface normal, `k -> 0+` is the **normal-incidence limit approached from nonzero obliquity**, **not grazing incidence**. Grazing incidence is theta approaching 90 degrees.
- **Strictly identify scope of PASS:** an independent 4x4 characteristic derivation checking the eight eigensolutions is strong support for algebraic correctness **under the stated PDEs**, not validation of measured material parameters or thermodynamic constitutive admissibility across all regimes.
- **Mutation-insensitive convergence is not itself an unexpected failure:** a mutated constitutive law can have the same formal O(k) convergence to its own exact-axis limit. The problem is claiming that convergence validates the unmutated physical law. Add independent absolute observables/closed-form physical benchmarks.
- **Uniform inverse norm:** raw `||S(k)^-1||` is unit-, eigenvector-normalization-, and row-scaling-dependent. Specify a **fixed physically nondimensionalized and normalized** matrix, track smallest singular value, inverse-norm bounds, and no near-misses. This is a finite-range numerical test, not a proof over all k; pair with Neumann-series argument.
- **Al6061:** the present electrical model is a non-piezoelectric epsilon_0 dielectric surrogate; no metal/electrode-contact physical claims are authorized. Mechanical properties of Al6061 may still be used, but the electrical interface assumption must be explicitly named.

## C. Next milestone: Scientific Gate 2U-V (before production editing)

**Deliverable 1: Short mathematical lemma.** Starting from the 8-state operator, show (1) analyticity of six simple outgoing mechanical/thermal branches at k=0 at fixed omega (within documented admissible frequency interval); (2) phi-normalized electrostatic branches lambda_±=±s_e k+O(k^2), mechanical subblock O(k) and D_z=O(k); (3) the scaled-interface matrix after dividing its D_z row by k tends to a block-triangular matrix `S0=[[M6,0],[C,E0]]`, `det(E0)=-(eps_eff,L+eps_eff,R) != 0`; (4) for `||S0^-1 (S(k)-S0)|| < 1`, the Neumann bound `||S(k)^-1|| <= ||S0^-1||/(1-||S0^-1(S(k)-S0)||)`. Specify physical nondimensionalization, normalization and incident radiation assumptions. Do not claim a global-in-frequency theorem.

**Deliverable 2: Real, adversarial negative tests.** Mutation injection must edit the active operator/assembler/selection, re-solve, compare against *independently anchored* physical observables, and fail loudly if undetected. Include: (a) p3 sign or magnitude; (b) beta sign/magnitude; (c) e33 sign, with consistent/otherwise explicitly inconsistent constitutive convention; (d) Lord-Shulman g misplacement at omega*tau~0.1 and ~1; (e) outgoing/incoming swap; (f) remove shear equation in actual assembly; (g) omit incident wave; (h) phi/Dz row mutations; (i) eps_eff electrostatic mismatch. Do not require an arbitrary >1e-3 field change from every mutation: some channels are physically insensitive; declare sensitivity criteria and counterexample cases a priori.

**Deliverable 3: External physics anchors.** Derive a self-contained isotropic one-dimensional Lord–Shulman thermoelastic characteristic equation (no piezoelectric coupling) from the PDEs and compare with the 8-state operator; independently check classical elastic normal-incidence amplitude when thermoelastic/electric coupling is removed. Benchmark LS relaxation factor with artificially enlarged tau AND several frequencies; treat artificial tau as an algorithm stress test, never a measured property. A published piezoelectric benchmark may be added only with source-verified equations and convention mapping.

**Deliverable 4: Stress the material/branch diversity.** Test two genuinely dissimilar piezoelectric materials (e.g. BaTiO3 if provenance adequate), at least one mixed pair with nontrivial electrostatic amplitude. Sweep frequency and k over sufficiently wide, physically interpretable ranges, record root separation and radiation classification by an independent flux/group-velocity check. No automatic extrapolation past branch crossings or non-passivity.

**Deliverable 5: Independent experiment manifest.** Record exact code SHA or hashes, constants, precision, frequencies, k ratio, raw per-row residuals, physical component scales, branch choices, mutant deltas, and failure cases. Flag assumed parameters clearly.

## D. Acceptance criteria before production repair

- Correct normal-incidence approach terminology.
- H1–H3 lemma written and tested with rigorous conditional assumptions.
- Every mutation truly acts on the relevant code path and yields a grounded detection condition; ineffective controls fail the test harness.
- At least one *external analytic* physics benchmark independent of the bulk operator passes with declared tolerance.
- Demonstrate full k>0 ideal-dielectric contact continuity and non-trivial mixed piezoelectric cases.
- No unsupported conductor, energy-conservation or absorption claim.
- Referee's open production blocker BLK-02A-26 remains open until the actual `assembly.py` path is repaired and ACPT02 redefined/validated; BLK-02A-27 remains open until the mathematical hypothesis and mode admissibility proof are accepted. SCT02/G4/G5 remain open.

## E. Division of labor (limited Claude tokens)

- GPT-6: mathematical lemma, independent test design, local numerical runs, artifact assembly and adversarial review.
- Claude: once a narrowly scoped code patch is ready, integrate into a dedicated GitHub branch, run repo CI and push the commit; do not spend limited tokens redoing broad exploratory work.
- Arena: fallback for extensive slow test matrices or alternate implementation only if necessary; avoid overlapping edits.

**Gate decision now:** scientific finite-k milestone **CONDITIONAL**, production **NOT READY**. Proceed to Gate 2U-V; hold production promotion and push.
