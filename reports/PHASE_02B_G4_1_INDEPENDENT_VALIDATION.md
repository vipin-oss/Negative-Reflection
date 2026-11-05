# PHASE 2B / TASK G4-1 — Independent Validation

Oracle: `code/phase_02b_g4_1_oracle.py` — built directly from the governing
PDEs/constitutive laws in `code/foundation_checks.py`; does NOT use the
production `flux_z_mp()`/`flux_z()`. Modal fields come from the accepted
production characteristic solve (`xi_roots` + `kernel_scaled`); the energy
identity itself is evaluated independently.

Results: `validation/phase_02b_g4_1_results.json`,
mutations: `validation/phase_02b_g4_1_mutations.csv`.

## Checks performed

- Propagating P and SV, thermal diffusive, electrostatic/near-gauge modes,
  all 8×8 bilinear pairs per material (Al6061 dielectric surrogate, PZT-5H,
  BaTiO3) at 1 MHz, k = 0 and k/kr = 1e-3, dps 60: bulk differential
  balance Re[i(ξ_m−ξ_n*)Jz_mn] + diss_mn = 0 evaluated per pair
  (includes interference cross terms and evanescent decay).
- Lord-Shulman regimes: physical ωτq ≈ 2.3e-6; SYNTHETIC ωτq ≈ 0.1 and
  ≈ 1 (labelled synthetic). Relaxation factor ω²τq²/(1+ω²τq²) positive in
  all three (5.4e-12, 9.9e-3, 0.5).
- Dissipation positivity: diss_eff = −Re[div⟨J⟩] ≥ 0 on physical modes.

## Measured residuals (worst pair, relative)

| material | k/kr=0 | k/kr=1e-3 |
|---|---|---|
| Al6061 | 0.997 | 2.2e5 |
| PZT-5H | 6.27 | 5.2e3 |
| BaTiO3 | 8.10 | 8.3e3 |

Interpretation: the naive flux J = −σv − Dφ̇ + qθ/T0 does NOT close the
balance on coupled modes — the relaxed pyroelectric coupling storage is
missing. This is a genuine finding, not a tautology: the residual is O(1)
on piezo modes and changes measurably under coupling mutations.

## Adversarial negative controls (mutation CSV)

M02 (no pyroelectric p), M03 (no transverse κ11 dissipation), M04 (wrong
τq factor), M08 (ε11/ε33 swap): DETECTED (residual changes by orders of
magnitude). M01/M05/M06/M07 sign/interference mutations do not move the
worst-pair metric because the unclosed coupling term dominates there —
recorded as a limitation of the naive-flux metric, not as passing controls.

## Control volume / reference plane

Bulk balance is checked spectrally per mode pair (the exact differential
form); no boundary energy-flux difference, no incident/reflected/transmitted
subtraction was used to define dissipation. Dissipation is evaluated from
the constitutive LS law only.

## Status

G4 remains OPEN pending supervisor acceptance of the complete relaxed
coupling identity. G5 and SCT02 unresolved. No production code modified.
