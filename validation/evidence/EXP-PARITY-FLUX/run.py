import json, sys
sys.path.insert(0, "code")
import sympy as sp
import foundation_checks as fc

xi = sp.Symbol("xi")
M = fc.symbolic_M()
Mn = M.subs(xi, -xi)
out = {"survivors": [], "phase1_S_failing_entries": []}
for bits in range(16):
    s = [1 if (bits >> i) & 1 == 0 else -1 for i in range(4)]
    S = sp.diag(*s)
    if sp.simplify(Mn - S * M * S) == sp.zeros(4, 4):
        out["survivors"].append(s)
Sa = sp.diag(1, 1, 1, -1)
out["phase1_S_failing_entries"] = [
    [i, j] for i in range(4) for j in range(4)
    if sp.simplify((Mn - Sa * M * Sa)[i, j]) != 0]
print(json.dumps(out, indent=2))
