# Conditional normal-incidence limit lemma — Task 2U-V

**Status:** mathematical proof under stated finite-dimensional hypotheses. This is not an unconditional passivity/radiation theorem, and is not a global-in-frequency assertion.

## Setting

Time-harmonic fields have convention `exp(i k x - i omega t)`; `k=omega sin(theta)/c` with theta measured from the normal. Hence **k→0 is normal incidence**, not grazing incidence.

Let the 8 physical interface fields be `(u_x,u_z,phi,theta,sigma_xz,sigma_zz,D_z,q_z)`. Assume ideal electrically bonded *dielectric* contact, `[phi]=[D_z]=0`, absent free charge/dipole layer, and welded mechanical and thermal contact. Al6061's current numerical electric constants describe an epsilon0 dielectric *surrogate*, NOT a conductor/electrode.

## Exact-axis reduction

At k=0, Gauss gives `D_z'=0`; consider the solution selected as a bounded-amplitude limit of oblique-incidence scattering. For its nonzero-spatial-eigenvalue longitudinal/thermal modes, `D_z=0` individually. From

```
D_z = e33 u_z' - eps33 phi' + p3 theta
sigma_zz = C33 u_z' + e33 phi' - beta3 theta,
```

one obtains in the `D_z=0` sector

```
phi' = (e33/eps33) u_z' + (p3/eps33) theta,
C_* = C33 + e33^2/eps33,
beta_* = beta3 - e33 p3/eps33,
R_* = rho c_E - T0 p3^2/eps33.
```

Then, with spatial mode `exp(lambda z)`, the direct PDE system is

```
(C_* lambda^2 + rho omega^2) U - beta_* lambda Theta = 0,
(i omega T0 beta_* lambda) U + (kappa33 lambda^2 /g + i omega R_*) Theta = 0,
g=1-i omega tau_q.
```

Thus the separate closed-form quartic is

```
(C_* lambda^2 + rho omega^2)
  (kappa33 lambda^2/g + i omega R_*)
  + i omega T0 beta_*^2 lambda^2 = 0.
```

At k=0 the shear wave decouples. Six outgoing amplitudes across two halfspaces (three per side: two P/thermal, one SV) satisfy six displacement/traction/thermal contact equations. Potential integration constants are determined only up to a single global gauge; imposing `[phi]=0` fixes their **relative** offset, but it is not an additional independent mechanical-scattering equation.

## Coalescing electrostatic roots

For each material with positive effective clamped-adiabatic permittivity,

```
eps33_eff = eps33 - T0 p3^2/(rho c_E),
s_e = sqrt(eps11/eps33_eff),
eps_e = sqrt(eps11 eps33_eff).
```

The electrostatic pair has `lambda_± = ± s_e k + o(k)` under the specified fixed-nonzero-omega hypotheses. A potential-normalized state yields `phi=1`, `D_z=O(k)` and nonpotential entries O(k) (possibly smaller), provided branch derivatives are analytic. This is independently supported by the preceding finite-k milestone, not re-derived in this exact-axis analytic code.

The Jordan generalized solution may have `phi` affine in z while other state variables remain finite. Excluding it merely because phi is unbounded is **not a valid argument**. Instead it is excluded **conditionally as a selected k→0⁺ oblique-incidence limit** under bounded inverse assumptions below, rather than as a universal exact-axis solution prohibition.

## Scaled-interface inverse

Order unknown columns `[six outgoing mechanical/thermal; two electrostatic]` and equations `[six mechanical/thermal contact; potential contact; normal displacement contact divided by k]`. Use *fixed* physical nondimensional row and amplitude scales and potential-normalized electrostatic modes. The interface matrix takes the form

```
S(k) → S0 = [ M6   0  ]
              [ C    E0 ],
E0 = [ 1         -1      ]
     [ -eps_e,L  -eps_e,R ],
det E0 = -(eps_e,L + eps_e,R).
```

**Hypotheses**:

- H1: the six mechanical/thermal outgoing roots at fixed frequency are simple, retain their selected radiation branches as k→0⁺, and the `6×6` contact matrix M6 is nonsingular.
- H2: the electrostatic pair is potential-normalized, approaches the above finite block, and the nonpotential entries vanish at least O(k). The thermal/electric constitutive model is nonsingular, including eps33_eff>0.
- H3: S(k) is continuous in the *fixed scaled basis*, its source vector remains bounded, and a radiation/causal branch selection is consistently prescribed for k>0.

Then, if H1–H3 hold, `det S0=det M6 det E0≠0`; for sufficiently small k,

```
||S0^{-1}(S(k)-S0)|| < 1,
S(k)^{-1} = [I + S0^{-1}(S(k)-S0)]^{-1} S0^{-1},
||S(k)^{-1}|| ≤ ||S0^{-1}||/[1-||S0^{-1}(S(k)-S0)||].
```

Therefore electrostatic modal amplitudes remain O(1); the total D_z is O(k) and the six physical mechanical/thermal amplitudes converge to the exact-axis `6×6` solution. The remaining global potential offset is gauge at k=0. A uniform nonzero D0 sector would need amplitude O(1/k), contrary to this bounded-inverse-selected family.

**What this does NOT prove:** H1–H3 at arbitrary frequency/material; an operator-theoretic limiting-absorption theorem under omega→omega+i0; conductor electrical contact; energy positivity/dissipation; correctness of the material constants; or that the previously implemented finite-k radiation classifier can never switch at a degeneracy.

## Executable evidence and limits

`task2uv_analytic_oracle.py` (40/70 digits) checks the LS quartic against an independently assembled **one-dimensional** first-order matrix, and the classical elastic impedance reflection relation. `task2uv_root_gap_certificate.py` checks the **sampled** distinctness of exact-axis P/thermal and SV roots in 27 (material/frequency/tau) combinations. Previous finite-k milestones supply nonpotential convergence and electrostatic block numerical evidence. None is a substitute for a global proof or actual production-path tests.
