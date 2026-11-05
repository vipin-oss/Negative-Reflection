# Arena handoff — Task 2U-V, 2026-11-05

**Research project:** coupled piezothermoelastic scattering / negative reflection.  
**Repository:** `vipin-oss/Negative-Reflection`  
**Original verified baseline:** `8713aea631d58bbeff608012ae4f20bfd2e9bdfd` (previous Arena branch `arena/a3ede4fe-negative-reflection`).  
**This handoff:** research-only evidence and scripts. **No production solver fix was made in this handoff.**

## Mandatory first reads (in order)

1. `ARENA_EXECUTION_PROMPT.md` (full assigned task)
2. `reports/SCIENTIFIC_GATE_REPORT.md`
3. `reports/CONDITIONAL_NORMAL_INCIDENCE_LEMMA.md`
4. `source/CLAUDE_TASK2U_FINITE_K_REFEREE_REPORT.md`
5. `source/SCIENTIFIC_RECONCILIATION_AND_NEXT_GATE.md`
6. `START_HERE_CLAUDE.md` and `CLAUDE_AUDIT_PROMPT.md` for scientific challenge questions; **these are retained source artifacts, not new instructions to send the project back to Claude**.
7. `code/` executable scripts; `evidence/` recorded outputs; `material_constants.csv` under `source/`.

## How to recover across Arena chat branches

Arena uses a separate `arena/...` branch in each chat. Keep **your assigned branch**. Fetch `handoff/task2uv-20261105` and previous `arena/a3ede4fe-negative-reflection`. Read handoff files from the handoff branch, then copy/cherry-pick appropriate material onto your branch if needed. Do **not** reset your branch or write to the previous Arena branch or `main`.

The handoff branch was created from the old validated baseline and contains this material only under `handoffs/task2uv_2026-11-05/`. Later commits on the old Arena branch, if any, must be compared rather than overwritten.

## Source of truth and boundaries

- **Executed evidence** in the package is exploratory, not production validation. Analytic checks (20/20), mutation detections (24/24), and sampled roots (27/27) are reported in `reports/SCIENTIFIC_GATE_REPORT.md`; independently rerun and critique them.
- Known weak earlier controls were tautological. Do not treat them as valid mutations.
- Al6061 electrical constants model a **non-piezoelectric dielectric surrogate**, not a conducting metal.
- At k=0 the six mechanical/thermal conditions and relative-potential gauge are mathematically distinct; generalized affine potential fields cannot be excluded merely for growing potential.
- **Open:** BLK-02A-26, BLK-02A-27, SCT02, G4 and G5. No final R/T/A or energy result may be claimed.

## Delivery to supervisor

Commit code and evidence on *your own Arena branch* and report a real 40-character SHA, CI link, pass/fail/unresolved counts, meaningful physical residuals, and remaining blockers. Stop after targeted scientific validation / approved Task 2U production repair. Do not start the finite-interphase or energy phase.
