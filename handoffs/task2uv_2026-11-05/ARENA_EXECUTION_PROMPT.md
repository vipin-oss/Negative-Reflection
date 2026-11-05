# ARENA EXECUTION BRIEF — TASK 2U-V VALIDATION & TARGETED TASK 2U REPAIR

You are the *implementation and independent testing agent*, NOT the scientific approver. GPT-6 supervises and audits actual commits. Preserve your own Arena chat-specific GitHub branch; don't switch to, reset or push the existing old Arena branch. The user deliberately reserves Claude tokens; do not route routine work to Claude.

## 0. Repository and artifact recovery

Repository: `vipin-oss/Negative-Reflection`
- Baseline: `8713aea631d58bbeff608012ae4f20bfd2e9bdfd`
- Old Arena branch: `arena/a3ede4fe-negative-reflection`
- Handoff branch: `handoff/task2uv-20261105`
- Documentation: `handoffs/task2uv_2026-11-05/ARENA_START_HERE.md`

Before edits, report your own branch, current head, old branch head, handoff commit SHA and ancestry. Import/copy handoff files without overwriting newer repository changes. Read the scientific report, conditional lemma and Claude referee report; inspect source (`foundation_checks.py`, `code/scattering/{assembly,roots,radiation,modal,state}.py`) and claim/blocker registries.

## 1. Scientific objective

Establish a production-usable, correctly scaled and radiation-admissible exact-normal-incidence scattering formulation for ideal **dielectric contact**, not a metal-electrode boundary condition. At exact k=0, formulate the 6×6 mechanical/thermal outgoing-amplitude system (P/thermal and SV decouple); reconstruct a relative potential integration constant by [phi]=0, modulo one global gauge. This is a *proposed* production route and may only be adopted after independent audit of all assumptions. Do not copy the experimental 7×7 as an 8-condition boundary solver.

## 2. Reproduce then attack the handoff evidence

Run scripts in `code/` at their documented precisions. Confirm the 1D independent LS analytic dispersion, elastic impedance benchmark, finite-k 8×8 physical-row mutations and root separation. Compare numeric output to shipped `evidence/` JSON and document the actual executed commands, platform, precision and tolerances. Reject any tautological detector; a mutation must execute a changed physical model through the real assembly and compare a meaningful original physical observable or residual. Test both P and SV and at least two **genuinely dissimilar piezoelectric** materials with provenance limitations explicit.

Do not claim a 40/70-digit check proves material constants or Maxwell physical laws. Record ASSUMED p3, beta, cE ambiguities. Preserve original failures and corrective explanations.

## 3. Proof and independent benchmark requirements

Write a short explicit theorem under H1/H2/H3. For scaled nonzero-k contact matrix S(k), exhibit block limit S0=[[M6,0],[C,E0]], E0=[[1,-1],[-eps_e,L,-eps_e,R]], det E0=-(eps_e,L+eps_e,R), derive the Neumann-series inverse bound and a checkable frequency-dependent domain of validity. This is a **k→0 normal-incidence** limit (angle measured from normal); do not call it grazing incidence or a universal limiting-absorption theorem.

Derive the constitutive-reduced 1D LS quartic with C*, beta*, R* from the actual model; validate g=1-i omega tau_q with synthetic omega tau tests and independent PDE assembly. Explicitly verify outgoing wave classification via flux/group velocity/decay where applicable; include near-degenerate frequency cases. An exact uniform-D0 affine branch may be mathematically allowed at k=0; exclude it only under the specified bounded-amplitude nonzero-k continuation and correctly stated physical assumptions.

## 4. Targeted production changes (only when mathematically justified)

Repair the known production `assembly.py` rule that deletes phi based on incident or basis emptiness. For k=0, use a properly derived 6×6 physical subsystem plus separately reconstructed relative potential gauge; for k>0 keep all ideal dielectric contact conditions including [phi]=[Dz]=0. Do not drop an equation just to make a square matrix. Preserve physical-unit state assembly and explicit row/column equilibration. Validate all four Al6061/PZT-5H pairings (with Al labelled dielectric surrogate) and dissimilar piezo pairs. Implement P and SV incidence with correct incoming RHS and outgoing signs.

The production acceptance test ACPT02 must measure gauge-invariant or gauge-aligned physical observables and k-slopes; not compare arbitrary raw potential amplitudes. Separate [phi] contact residual from global gauge freedom. Add real production-path mutations for phi, Dz, shear, RHS, radiation and LS g.

## 5. Acceptance and stop gates

- Correct physical mode count and complete H1-H3 proof for tested domain, no unlicensed global assertion.
- Finite-k physical 8×8 vs k=0 reduced 6×6 mechanical/thermal convergence with measured per-component orders and precision stability.
- Interface contact residuals explicitly reported by field; no unscaled mixed-SI Euclidean norm.
- Independent analytic LS and elastic benchmark pass under an explicit prespecified tolerance.
- Negative controls detect actual defects; no always-true checks.
- Existing tests and mutations re-run. CI reports actual statuses. Blockers BLK-02A-26/27 close **only if** production and physics evidence genuinely clear them.
- SCT02, G4/G5 and prior-art/material provenance gates remain open; no total energy balance claim, optimization, finite interphase or paper promotion.

Commit and push to **your current Arena-assigned branch**. Return changed paths, exact 40-character SHA, CI run link, reproducible commands, exception cases, all measured residuals and a concise acceptance verdict. **STOP**. Do not spend time on broad literature or rewrite the full manuscript.
