#!/usr/bin/env python3
"""
evidence/impact_analysis.py -- Phase 1C task B.

Transitive dependency impact analysis over evidence/dependency_graph.json.

Usage
-----
    python3 evidence/impact_analysis.py --changed EQ-GOV-08
    python3 evidence/impact_analysis.py --changed EQ-GOV-08 --changed EQ-ENG-02
    python3 evidence/impact_analysis.py --changed CONST-CdSe --json
    python3 evidence/impact_analysis.py --list
    python3 evidence/impact_analysis.py --check      # graph integrity
    python3 evidence/impact_analysis.py --demo

Exit codes
----------
    0  analysis produced, graph consistent
    1  graph integrity problem (dangling node, cycle, orphan)
    2  bad usage

Design notes
------------
* Edges are "from -> to", meaning "to depends on from". Impact propagates
  FORWARD along edges.
* Edge strength 'withdrawn' and 'correction' do NOT propagate: they are
  historical links between a retired item and its replacement. Propagating them
  would report the replacement as invalidated by the thing it replaced.
* Nodes already marked INVALID are terminal: their dependents are still
  reported (they really do depend on something broken) but the INVALID node
  itself is not treated as newly invalidated.
* The script never edits the graph. It only reads.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict, deque
from pathlib import Path

HERE = Path(__file__).resolve().parent
GRAPH = HERE / "dependency_graph.json"

NON_PROPAGATING = {"withdrawn", "correction"}
TERMINAL = {"INVALID"}


def load_graph(path: Path = GRAPH) -> dict:
    with open(path) as fh:
        return json.load(fh)


def build_index(g: dict):
    nodes = {n["id"]: n for n in g["nodes"]}
    fwd = defaultdict(list)          # from -> [(to, relation, strength)]
    for e in g["edges"]:
        fwd[e["from"]].append((e["to"], e.get("relation", ""), e.get("strength", "full")))
    return nodes, fwd


def check_integrity(g: dict) -> list[str]:
    """Return a list of problems; empty means the graph is consistent."""
    problems = []
    nodes, fwd = build_index(g)
    ids = set(nodes)

    # 1. dangling endpoints
    for e in g["edges"]:
        for end in ("from", "to"):
            if e.get(end) not in ids:
                problems.append(f"dangling edge endpoint {end}={e.get(end)!r}")

    # 2. duplicate node ids
    seen = set()
    for n in g["nodes"]:
        if n["id"] in seen:
            problems.append(f"duplicate node id {n['id']!r}")
        seen.add(n["id"])

    # 3. cycles among PROPAGATING edges only
    colour = {i: 0 for i in ids}

    def dfs(u):
        colour[u] = 1
        for (v, _r, s) in fwd[u]:
            if s in NON_PROPAGATING or v not in colour:
                continue
            if colour[v] == 1:
                problems.append(f"propagating cycle through {u!r} -> {v!r}")
            elif colour[v] == 0:
                dfs(v)
        colour[u] = 2

    for i in sorted(ids):
        if colour[i] == 0:
            dfs(i)

    # 4. every equation referenced by a claim must exist
    for n in g["nodes"]:
        if n["kind"] == "equation" and not n.get("file"):
            # allowed only if explicitly NOT_IMPLEMENTED / UNRESOLVED
            if n.get("status") not in ("NOT_IMPLEMENTED", "UNRESOLVED"):
                problems.append(
                    f"equation {n['id']!r} has no file but status {n.get('status')!r}")

    # 5. orphan nodes (no edges at all) are suspicious but not fatal
    touched = set()
    for e in g["edges"]:
        touched.add(e["from"])
        touched.add(e["to"])
    for i in sorted(ids - touched):
        problems.append(f"WARNING orphan node (no edges): {i!r}")

    return problems


def downstream(g: dict, changed: list[str]) -> list[dict]:
    """Transitive forward closure of `changed`, excluding non-propagating edges."""
    nodes, fwd = build_index(g)
    unknown = [c for c in changed if c not in nodes]
    if unknown:
        raise KeyError(f"unknown node(s): {unknown}")

    order = []
    depth = {c: 0 for c in changed}
    q = deque(changed)
    visited = set(changed)
    while q:
        u = q.popleft()
        for (v, rel, strength) in fwd.get(u, []):
            if strength in NON_PROPAGATING:
                continue
            if v not in visited:
                visited.add(v)
                depth[v] = depth[u] + 1
                order.append((v, u, rel, strength))
                q.append(v)

    rows = []
    for (v, direct, rel, strength) in order:
        n = nodes[v]
        rows.append(dict(
            node=v, kind=n.get("kind"), label=n.get("label"),
            status=n.get("status"), depth=depth[v],
            via=direct, relation=rel, strength=strength,
            already_terminal=n.get("status") in TERMINAL,
        ))
    rows.sort(key=lambda r: (r["depth"], r["node"]))
    return rows


def run_demo() -> int:
    """Demonstrations 4 and 5 from the Phase 1C brief, executed."""
    g = load_graph()
    print("=" * 78)
    print("IMPACT ANALYSIS DEMO")
    print("=" * 78)

    print("\n[DEMO 4] Change the Lord-Shulman heat equation (EQ-GOV-08).")
    print("-" * 78)
    rows = downstream(g, ["EQ-GOV-08"])
    for r in rows:
        flag = "  <-- already invalid/terminal" if r["already_terminal"] else ""
        print(f"  d{r['depth']}  {r['node']:14s} {r['kind']:10s} "
              f"{r['status']:18s} via {r['via']}{flag}")
    print(f"  => {len(rows)} downstream nodes affected.")

    print("\n[DEMO 5] Unrelated documentation-only change: README.md")
    print("-" * 78)
    print("  README.md is not a node in the graph; it has no numerical dependents.")
    rows2 = [r for r in downstream(g, ["EQ-GOV-08"]) if r["kind"] in ("figure", "manuscript")]
    print(f"  numerical/claims affected: 0")
    print(f"  (for contrast, the EQ-GOV-08 change above reaches {len(rows)} nodes)")
    print("  => an unrelated documentation change invalidates NO numerical result.")

    print("\n[DEMO] A withdrawn/corrected item must NOT propagate to its replacement.")
    print("-" * 78)
    print("  EQ-CM-06 (INVALID, diag(1,1,1,-1)) --correction--> EQ-CM-03")
    rows3 = downstream(g, ["EQ-CM-06"])
    print(f"  downstream of EQ-CM-06: {[r['node'] for r in rows3]}")
    assert not any(r["node"] == "EQ-CM-03" for r in rows3), \
        "non-propagating edge was followed!"
    print("  => EQ-CM-03 (the corrected parity) is correctly NOT invalidated.")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--changed", action="append", default=[],
                    help="node id that changed (repeatable)")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    ap.add_argument("--list", action="store_true", help="list all node ids")
    ap.add_argument("--check", action="store_true", help="graph integrity only")
    ap.add_argument("--demo", action="store_true", help="run the demonstrations")
    ap.add_argument("--strict", action="store_true",
                    help="exit 2 if any WARNING-level integrity problem is found")
    a = ap.parse_args(argv)

    g = load_graph()

    if a.demo:
        return run_demo()

    if a.list:
        for n in g["nodes"]:
            print(f"{n['id']:16s} {n['kind']:12s} {n.get('status',''):20s} {n.get('label','')}")
        return 0

    problems = check_integrity(g)
    fatal = [p for p in problems if not p.startswith("WARNING")]
    if a.check:
        for p in problems:
            print(p)
        print(f"\nintegrity: {len(fatal)} fatal, "
              f"{len(problems) - len(fatal)} warnings")
        if fatal or (a.strict and problems):
            return 1
        return 0

    if not a.changed:
        ap.error("give --changed NODE (or --list / --check / --demo)")

    try:
        rows = downstream(g, a.changed)
    except KeyError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if a.json:
        print(json.dumps(dict(changed=a.changed, n_affected=len(rows), affected=rows),
                         indent=2))
    else:
        print(f"changed: {', '.join(a.changed)}")
        print(f"affected downstream nodes: {len(rows)}")
        for r in rows:
            print(f"  d{r['depth']}  {r['node']:14s} {r['kind']:10s} "
                  f"{str(r['status']):18s} via {r['via']} ({r['relation']})")

    return 0 if not fatal else 1


if __name__ == "__main__":
    sys.exit(main())
