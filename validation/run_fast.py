#!/usr/bin/env python3
"""LEVEL 1 -- FAST. Runs on every commit / every push.

Scope, deliberately narrow so it stays under a minute:
  * imports and syntax of every code module
  * schema conformance of both registries and the material table
  * evidence-artifact presence
  * small symbolic identities (constitutive, parity)
  * one small analytical benchmark (elastic reference via subprocess)
  * nonzero execution of every case
  * registry consistency (no VERIFIED claim without an existing test)
  * targeted dependency-chain integrity

What it does NOT do: energy balances, convergence, mutation tests, figure
regeneration. Those are Level 2/3.

    python3 validation/run_fast.py
    python3 validation/run_fast.py --promotion     # also enforce blocker gate

Exit codes: 0 all pass, 1 failure, 2 unresolved/unknown, 3 nothing executed.
"""
import sys
from _runner import run

if __name__ == "__main__":
    sys.exit(run("1", "FAST", markers=None,
                 description="every-commit gate: syntax, schema, units, symbolic "
                             "identities, one analytical benchmark, evidence "
                             "consistency, dependency chain"))
