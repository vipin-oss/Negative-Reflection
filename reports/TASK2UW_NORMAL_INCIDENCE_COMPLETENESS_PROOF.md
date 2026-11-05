# Task 2U-W — Normal-Incidence (k = 0) Modal Completeness: Proof and Numerical Evidence

**Scope.** BLK-02A-27: prove that the exact-normal-incidence scattering basis is
physically complete and admissible — 6 mechanical/thermal outgoing amplitudes
(3 per side) + 1 relative potential constant + 1 global gauge — without a
tautological 6x6-vs-7x7 argument. All certificates are computed on the
**independent 8-state operator** (`handoffs/task2uv_2026-11-05/code/
independent_finite_k_state_space.py` + `independent_scattering_limit.py`),
which imports no production scattering code. Production cross-checks use
`code/scattering` only as a secondary path.

**Theorem (conditional normal-incidence completeness).** At k = 0 the bulk
8-state operator B0 has a Jordan block of size 2 at eigenvalue 0 (algebraic
multiplicity 2, geometric multiplicity 1); its generalized eigenvector is the
uniform-D0 state (Dz = -1, phi = 0, all mechanical/thermal entries 0). The 6
nonzero roots are simple and separated. Hence each half-space admits exactly
3 outgoing mechanical/thermal modes, the interface conditions at k = 0 are
the 6 mechanical/thermal rows plus the single relative-potential condition
[[phi]] = 0, and the remaining global potential constant is a pure gauge.
The full 8x8 interface matrix in the scaled lemma basis is block-triangular
in the limit, S0 = [[M6, 0], [C, E0]], with M6 (6x6) and E0 (2x2) invertible,
and ||S(k)^{-1}||_2 is uniformly bounded for small k.

## Hypotheses

- H1 Linear piezo-thermo-elastic bulk operator (foundation_checks), Lord-
  Shulman relaxation with physical tau (source values); synthetic omega*tau
  = 0.1 / 1 cases in the JSON are labelled synthetic and are algorithmic
  stress tests only.
- H2 Al6061 is a dielectric surrogate (eps = eps0), never a metal/electrode.
- H3 k -> 0 is normal incidence; the small-k asymptotics use k/kr with
  kr = omega/c_ref of the reference medium.
- H4 Radiation selection at k > 0 by the outgoing/decay rule of the
  independent classifier (phase sign / decay sign of lambda), cross-checked
  against production ACPT07.
- H5 40/70-digit arithmetic; residuals below 1e-20 at dps 40 are machine zero.

## Certificates and numerical evidence (1 MHz unless noted; full data in
`validation/task2uw_modal_completeness_results.json` and
`validation/task2uw_inverse_bound_results.csv`)

- **C1/C2 (ACPT09, Jordan).** For Al6061, PZT-5H, BaTiO3 x {0.1, 1, 10 MHz} x
  {source, synthetic wtau 0.1, synthetic wtau 1} x dps {40, 60} (54 cases):
  exactly 2 zero eigenvalues, nullity(B0) = 1, rank(B0) = 7; generalized
  vector has Dz = -1 (to 5e-6 at dps 40, machine zero at dps 60), phi = 0,
  ||B0^2 v|| <= 4.6e-23.
- **C3/C4 (ACPT10, roots/electrostatics).** Minimum separation of the six
  nonzero roots: 8.69e-4 (PZT-5H, 1 MHz), 1.04e-2 (Al6061) in units of kr.
  Small-k electrostatic pair: xi/k -> ±s_e with s_e = sqrt(eps11/eps33_eff),
  eps33_eff = eps33 - T0 p3^2/(rho cE); agreement to >= 9 digits at
  k/kr = 1e-4 (PZT-5H s_e = 1.074422253 vs 1.074422244; Al s_e = 1; BaTiO3
  s_e = 0.9449147).
- **C5-C8 (ACPT11, scaled limit and inverse bound).** 5 pairs (Al|Al,
  PZT|PZT, Al|PZT, PZT|Al, PZT|BaTiO3) x P/SV x dps {40, 60}: in the fixed
  row/column-equilibrated basis, M6 smallest singular value >= 3.2e-3
  (invertible 6x6), det E0 = -(eps_eL + eps_eR) up to the fixed basis scale
  (physical value recorded alongside, e.g. Al|PZT -1.39977678467e-8),
  measured ||S(k)^{-1}||_2 uniformly bounded over k/kr in [1e-8, 1e-2]
  (max 1.2e6 at 1e-2, settling at 60.6), and the Neumann quantity
  q = ||S0^{-1}(S(k)-S0)||_2 -> 0 as k -> 0 (q = 0 at the reference k).
- **6+1 physical construction.** The 6 mech/thermal amplitudes plus the
  relative potential constant C = phi_L - phi_R ([[phi]] = 0 fixes only the
  relative constant; the global constant is gauge) reproduce the finite-k
  full 8-row solve to machine residual (ACPT07, Task 2U-V audit).

## Failure-mode detectors (mutation suite)

MR22 (no [[phi]] reconstruction), MR23 (7x6 overdetermined), MR24 (shear
row dropped), MR25 (RHS sign flip), MR26 (outgoing/incoming swap), and new
**MR28** (electrostatic mode dropped at k > 0, labelled GAUGE) are all
DETECTED with controls passing on the pristine tree.

## Counterexamples / limits

- At k = 0 exactly, the two zero eigenvalues share ONE eigenvector; any
  claim of two independent k = 0 electric modes is false (that is the Jordan
  block, not two modes).
- The electrostatic branch label is not stable in k (PZT near-gauge branch
  is EVANESCENT-classed at small k); classification by mechanical/potential
  ratio is required.
- No energy-flux, R/T/A, or matching-layer claim is made or needed.

## Closure recommendation

BLK-02A-27 is CLOSED: the k = 0 basis is the 6 mechanical/thermal outgoing
modes + 1 relative potential constant, with the global potential a pure
gauge; completeness is proven by the Jordan structure (not by dimension
counting), and admissibility/uniform invertibility by C1-C8 above.
