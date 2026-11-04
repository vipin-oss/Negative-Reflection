# Evidence system

This directory is the machine-readable record of what the project claims, what
supports each claim, and what would happen downstream if the support changed. It
exists because Phase 1 contained a genuine methodological failure: documents said
"verified" against tests that did not execute the asserted calculation, and one
document stated a symmetry transformation that is false.

The system is designed around one rule:

> A claim is `VERIFIED` only if a named piece of code, run on a named fixture,
> computes the asserted quantity and produces an output consistent with an
> independent reference. A PASS label is a headline, not a proof.

## Files

| file | what it is | who writes it |
|---|---|---|
| `claim_registry.csv` | 44 research claims with status, exact evidence location and failure mode | human, audited |
| `claim_evidence_map.json` | the machine-checkable binding: claim -> function -> required tokens | human, checked by CI |
| `equation_registry.csv` | 41 equations with assumptions, dimensions, derivation status, oracle, dependencies | human, audited |
| `dependency_graph.json` | nodes and edges from source -> constant -> equation -> code -> claim -> experiment -> figure -> manuscript | human, checked by `--check` |
| `impact_analysis.py` | the tool that answers "if X changes, what breaks" | code |
| `phase_01_claim_audit.md` | the U-11 audit, findings and actions | human |
| `approved_baselines.json` | append-only baseline history, including rejected ones | human |
| `manuscript_evidence_map.json` | stable ids for future manuscript items | human |
| `input_checksums.json` | SHA-256 of every frozen scientific input | `validation/freeze_inputs.py` |

## Status vocabulary

Claims and equations:

- `VERIFIED_STRONG` - executed, independent oracle, tolerance justified
- `VERIFIED_WEAK` - executed, but the oracle is circular or the scope is narrower than the claim
- `UNVERIFIED` - no executed or inspectable basis, though possibly plausible
- `UNRESOLVED` - an execution exists but does not close
- `INVALID` - the statement is false. Retained on purpose.

`INVALID` rows are never deleted. A correction that erases the record of the
failure is worse than the failure.

## The two questions this system answers

**"Is this claim supported?"**
`python3 validation/run_fast.py` runs `test_BIND01`, which resolves every
`VERIFIED` claim to a function and checks that the tokens carrying the asserted
quantity are still present. If the calculation is deleted, the claim is
downgraded automatically and the run fails.

**"If I change this equation, what breaks?"**
`python3 evidence/impact_analysis.py --changed EQ-GOV-08` walks the graph and
lists every dependent node with its current status. Changes that propagate to
nothing numerical - documentation, for example - correctly report zero impact.

## Adding a claim

1. Add the row to `claim_registry.csv` with a new `CLM-<AREA>-<nn>` id.
2. Add a binding to `claim_evidence_map.json` naming the test and at least one
   required token.
3. If there is no test, the status must be `UNVERIFIED`. That is acceptable and
   common. What is not acceptable is `VERIFIED` without a binding.

## Known defects recorded here

These are real and are not hidden by the infrastructure:

- `CLM-C1-04` claims a symbolic null space that `check_C1_structure` no longer
  computes (U-11-1). Status `UNVERIFIED`. The calculation is now done properly in
  `validation/tests/test_02_equations.py::test_EIG03_null_space_is_real`.
- `check_F1_passivity` samples a matrix that is positive semi-definite by
  construction, so the test cannot fail (U-11-2). Replaced by `test_IFC01`, which
  includes a negative control that must be detected.
- `check_I2_interface_rank` returns the literal string `"PASS"` regardless of the
  rank it computes. Retained as evidence; the real measurement is `test_IFC02`.
- `governing_equations.tex` still prints the withdrawn `tau_q` causality floor
  (U-11-4). The document and the code disagree.
