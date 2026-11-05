# START HERE — Task 2U-V Claude Scientific Referee Gate

Purpose: adversarial review of new **independent 1D Lord–Shulman dispersion benchmark, actual constitutive perturbations, finite-k physical row mutations**, and conditional normal-incidence block-inverse lemma. User deliberately holds GitHub push and production edit until review.

Read in order:
1. `reports/SCIENTIFIC_GATE_REPORT.md`
2. `reports/CONDITIONAL_NORMAL_INCIDENCE_LEMMA.md`
3. `source/CLAUDE_TASK2U_FINITE_K_REFEREE_REPORT.md` (your previous report)
4. `source/SCIENTIFIC_RECONCILIATION_AND_NEXT_GATE.md`
5. `code/task2uv_analytic_oracle.py`, `code/task2uv_finitek_mutations.py`, `code/task2uv_root_gap_certificate.py`
6. `evidence/*.json`

The two `code/independent_*.py` files are the prior experimental 8-state implementation copied unchanged into this **review harness**. The latest mutation code imports them for finite-k modes. The 1D LS analytic oracle imports *neither* of them.

Reproduce (Python + mpmath installed, from any working directory):

```
python code/task2uv_analytic_oracle.py --digits 40 70 --out evidence/analytic_oracle_recheck.json
python code/task2uv_finitek_mutations.py --digits 40 70 --out evidence/finitek_mutation_recheck.json
python code/task2uv_root_gap_certificate.py
```

Do **not** edit source repository, commit, or push. Examine actual mutations, natural physical units, independent derivations, and any subtle hidden assumptions. Return separate verdicts for analytic LS, 1D mutations, finite-k mutations, conditional lemma and production readiness. Record any failed reproduction with exact script and input.

**Cautions:** Al6061 is epsilon0 dielectric surrogate. Some beta/p/tau values are assumed/artificial. `k→0` is *normal incidence approached obliquely*. A 7×7 algebraic contact identity is not an independent radiation proof. No G4/G5 energy validation is claimed.
