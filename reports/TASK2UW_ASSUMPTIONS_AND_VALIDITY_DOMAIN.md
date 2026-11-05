# Task 2U-W — Assumptions and Validity Domain

## Assumptions (all stated, none hidden)

1. **Linear piezo-thermo-elastic bulk** with Lord-Shulman relaxation
   (foundation_checks material tables; source tau values). The synthetic
   omega*tau = 0.1 and 1 cases in
   `validation/task2uw_modal_completeness_results.json` are **labelled
   synthetic**: algorithmic stress tests of the certificate code, never
   measured relaxation properties.
2. **Al6061 = dielectric surrogate** (eps = eps0), never treated as a metal
   or electrode.
3. **k -> 0 means normal incidence**; small-k expansions use k/kr with
   kr = omega/c_ref of the reference (left) medium. Grazing incidence is a
   different limit and is not covered.
4. **Radiation selection** at k > 0 by the independent classifier's
   phase/decay rule, cross-checked against the production finite-k solve
   (ACPT07). No energy-flux-based radiation condition is used or claimed.
5. **High-precision arithmetic** (dps 40 and 60); residuals below ~1e-20 at
   dps 40 are machine zero. 40/70-digit checks prove implementation
   correctness, not new material constants or Maxwell laws.
6. **Block-triangular scaled limit** uses one fixed diagonal
   row/column equilibration of the reference matrix; the same factors are
   applied to every S(k), so the Neumann series lives in a single fixed
   basis.

## Validity domain

- Materials: Al6061, PZT-5H, BaTiO3 (source tables); pairs Al|Al, PZT|PZT,
  Al|PZT, PZT|Al, PZT|BaTiO3; P and SV incidence; 0.1/1/10 MHz.
- k/kr in [1e-8, 1e-2] for the uniform inverse bound; q_neumann -> 0 as
  k -> 0 and the measured ||S(k)^{-1}||_2 is bounded over the whole range.
- The k = 0 statements (Jordan block, 6 nonzero roots, [[phi]] = 0) are exact
  at k = 0 and verified to machine precision.

## Explicitly out of scope

- No energy-conservation, R/T/A, or matching-layer claims.
- No production code modification (no demonstrated defect required it).
- BLK-02A-26 remains closed and untouched; SCT02, G0, G1, G4, G5,
  LS-STORAGE remain open and are unaffected by this task.

## Open questions

- Whether the supervisor wants the C-block (mech-column Dz/k limit) reported
  entrywise; current deliverables report the block determinants and norms.
- Extension of the uniform inverse bound beyond k/kr = 1e-2 would need a
  coarser-grained basis; not required for normal incidence.
