# Task 2U-W — Final Analytic-Limit Certificate Check

Narrowly scoped verification of the corrected certificate. No production
scattering code touched; no energy accounting; BLK-02A-27 status unchanged.

## 1. Exact S0 construction and source locations

All in `code/task2uw_certificates.py`, function `scaled_S` (k = 0 branch,
lines ~230–260) and `cert_inverse_bound` (lines ~280–330):

- **M6 block (rows 0–5, mech columns):** analytic — k = 0 nonzero
  eigenvectors of the raw independent operator (`_raw_modes(p, 0, w)`),
  mechanical-max normalized, gauge-fixed by `_fix_sign` (dominant entry
  positive-real). No finite-k data.
- **Top-right block (mech rows × electrostatic columns):** analytic zero
  (potential-normalized electrostatic eigenvector is mechanically clamped,
  mech entries O(k) → 0).
- **E0 block (phi, Dz/k rows × electrostatic columns):** analytic —
  phi = 1, Dz/k = ∓eps_e with eps_e = sqrt(eps11·eps33_eff),
  eps33_eff = eps33 − T0·p3²/(rho·cE) (`eps_eff`, line 55).
- **C block (Dz/k row × mech columns):** NOW FULLY ANALYTIC (this check):
  C_j = −i·e15·ux_j, derived below. The previous version used two-point
  Richardson extrapolation of the on-shell Dz/k from k/kr = 1e-7, 1e-9;
  that path is removed. The corrected `main()` deliverables
  (`validation/task2uw_modal_completeness_results.json`,
  `validation/task2uw_inverse_bound_results.csv`) are regenerated with the
  analytic C-block and give identical q(k) (C ~ 1e-13, numerically
  negligible).
- **Equilibration:** one fixed two-sided diagonal scaling of S0
  (`rsc`, `csc` in `cert_inverse_bound`), applied unchanged to every S(k).

## 2. Independent C-block derivation

Plane-wave convention **exp(i k x − i ω t)**, state y′ = A y, so
d/dz → i·xi (operator header, `independent_finite_k_state_space.py`
line 4). Gauss law ∇·D = 0 gives

    i k Dx + i xi Dz = 0   =>   Dz/k = −Dx/xi        (NO factor of i)

consistent with the operator's own Dz row `Dzp = -I*k*Dx` (line 49):
i·xi·Dz = −i·k·Dx. With Dx = e15·(ux′ + i k uz) − i k eps11·phi
(line 46), at k = 0: Dx = e15·i·xi·ux, hence

    C_j = lim Dz/k = −Dx/xi = −i·e15·ux_j   (mech-normalized branch j).

**Entrywise comparison** (analytic vs the removed Richardson values),
`validation/task2uw_analytic_limit_check.json`:

| material | dps 40 | dps 60 | dps 100 |
|---|---|---|---|
| Al6061 | 1.3e-44 | 8.0e-65 | 1.7e-104 |
| PZT-5H | 7.9e-30 | 7.9e-30 | 7.9e-30 |
| BaTiO3 | 5.2e-30 | 5.2e-30 | 5.2e-30 |

C magnitudes are ~1e-13 (Al ~ 0 since e15 = 0), so agreement to ~1e-30 is
at solver-residual level: the analytic C-block is confirmed.

## 3. Approximate-reference error bound

No longer needed: the analytic C-block replaces the extrapolation, so S0 is
exact (to working precision) and the Neumann interval stands on the exact
reference: q(k) = ‖S0⁻¹(S(k)−S0)‖₂ vanishes linearly (Al|PZT P: 222.5 →
2.2e-6 over k/kr = 1e-2 → 1e-10), **q < 1 for k/kr ≲ 1e-5**, where
‖S(k)⁻¹‖ ≤ ‖S0⁻¹‖/(1−q) (bound 2.08 → 1.617, measured ≤ bound). Rows with
q ≥ 1 remain `bound_holds = False`.

## 4. Factor-of-i verification

Numerically on the modal eigenvectors (dps 60, k/kr = 1e-6): the state's
own Dz/k equals **−Dx/xi** exactly (PZT-5H P-branch: both
6.04881e-14 + 1.95531e-13j); −i·Dx/xi does not match. The previous
certificate code used −i·Dx/xi — a spurious factor of i, now corrected;
it was numerically harmless because the C-block is ~1e-13, but it was
wrong and is fixed.

## 5. Verdict

The mathematical completeness theorem **survives** this check: Jordan
block at k = 0 (alg mult 2, geom mult 1, generalized vector normalized to
Dz = −1 exactly), 6 separated nonzero roots, electrostatic pair ±s_e·k,
exact analytic S0 with invertible M6 (σ_min ≥ 3.2e-3) and E0
(det = −(eps_eL+eps_eR) up to fixed scale), Neumann q < 1 on k/kr ≲ 1e-5.
ACPT09–11 PASS (`validation/tests/test_20_scattering_task2uw.py`).
BLK-02A-27 closure evidence stands; status not altered for CI.
