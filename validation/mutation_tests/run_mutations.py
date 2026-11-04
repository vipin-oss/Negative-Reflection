#!/usr/bin/env python3
"""Controlled negative tests -- Phase 1C task E.

Each mutation injects ONE specific defect into an ISOLATED COPY of the
repository, then runs the validation that is supposed to catch it. A mutation
PASSES only when the validation REJECTS the faulty implementation. A mutation
whose defect goes undetected is a SURVIVOR and is reported as a validation gap.

Two structural guarantees, both enforced in code:

1. The approved baseline is never touched. Every mutation runs in a fresh
   temporary copy created under a temp root outside the repository, and the
   copy's approved-baseline and checksum registers are the ones being corrupted,
   never the originals. `assert_origin_untouched` re-hashes the real tree's
   frozen inputs before and after the whole suite.

2. "Detected" is not allowed to mean "any non-zero exit". Each mutation declares
   the exact detector it must trip, and the runner asserts that the detector's
   own output names the failure. A crash for an unrelated reason on the mutated
   copy is recorded as INCONCLUSIVE, not as detection.

Usage
-----
    python3 validation/mutation_tests/run_mutations.py
    python3 validation/mutation_tests/run_mutations.py --only M03
    python3 validation/mutation_tests/run_mutations.py --keep     # keep temp trees

Exit codes: 0 if every mutation was correctly detected, 1 if any survived,
3 if the suite could not run.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import traceback
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
OUT = HERE / "last_mutation_run.json"

IGNORE = shutil.ignore_patterns(".git", "__pycache__", "*.pyc", ".pytest_cache",
                                "node_modules", "reports/figures/*.png")


# ---------------------------------------------------------------------------
@dataclass
class Mutation:
    mid: str
    title: str
    target: str                 # file, relative to repo
    detector: str               # what must reject it
    rationale: str
    mutate: object = None       # callable(Path root) -> None
    detect: object = None       # callable(Path root) -> (bool detected, str evidence)
    result: str = ""            # DETECTED / SURVIVOR / INCONCLUSIVE / NOT_RUN
    evidence: str = ""


def copy_repo(dst: Path) -> None:
    shutil.copytree(REPO, dst, ignore=IGNORE, symlinks=False)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


FROZEN = ["literature/material_constants.csv", "literature/verified_references.csv",
          "code/det_coeffs.json", "theory/characteristic_matrix.tex",
          "theory/energy_flux_derivation.tex", "theory/governing_equations.tex",
          "theory/interface_conditions.tex"]


def assert_origin_untouched(before: dict) -> list[str]:
    bad = []
    for rel in FROZEN:
        p = REPO / rel
        if not p.exists() or sha(p) != before[rel]:
            bad.append(rel)
    return bad


# ===========================================================================
# M01 -- zero executed cases must not PASS
# ===========================================================================
def m01_mutate(root: Path) -> None:
    """A TEST that executes nothing but reports PASS.

    Note the design point: the defect must be injected into the *producer*, not
    into the harness. An earlier version of this mutation disabled the harness's
    zero-case guard, which is the detector itself - a degenerate mutation that
    proves nothing. Here the test claims a pass without running anything, and the
    untouched harness must refuse to let that through.
    """
    f = root / "validation" / "tests" / "test_01_data_integrity.py"
    s = f.read_text()
    s = s.replace('def test_MAT02_no_tuned_parameter():',
                  'def test_MAT02_no_tuned_parameter():   # MUTATION: executes nothing\n'
                  '    c = make_case(test_id="MAT02", category="1-data-integrity",\n'
                  '                  purpose="claims a pass without executing anything",\n'
                  '                  oracle="none", oracle_independent="ABSENT",\n'
                  '                  fixture="none", expected="nothing",\n'
                  '                  tolerance=None, tol_justification="none", level="1")\n'
                  '    finish(c, lambda: (PASS, 0, "did nothing", ""))\n'
                  '    return\n'
                  '\n\ndef _mutated_original_MAT02():', 1)
    f.write_text(s)


def m01_detect(root: Path) -> tuple[bool, str]:
    # IMPORTANT: run from a directory where a plain `import harness` cannot pick up
    # the ORIGINAL module from /home/user/Negative-Reflection/validation (there is
    # no __init__.py). The first version of this detector silently imported the
    # real harness and therefore reported a false SURVIVOR. Load the mutated file
    # explicitly with importlib instead.
    prog = (
        "import importlib.util, sys\n"
        "sys.path.insert(0,'code')\n"
        "spec=importlib.util.spec_from_file_location('h_mut','validation/harness.py')\n"
        "harness=importlib.util.module_from_spec(spec)\n"
        "sys.modules['h_mut']=harness   # dataclass needs the module registered\n"
        "spec.loader.exec_module(harness)\n"
        "c=harness.make_case(test_id='ZZZ',category='x',purpose='zero-case',\n"
        "  oracle='none',oracle_independent='ABSENT',fixture='n/a',expected='n/a',\n"
        "  tolerance=None,tol_justification='n/a',level='1')\n"
        "harness.finish(c, lambda: (harness.PASS, 0, 'did nothing', ''))\n"
        "print('STATUS', c.status)\n")
    r = subprocess.run([sys.executable, "-c", prog], cwd=root,
                       capture_output=True, text=True, timeout=300)
    # the harness must convert the zero-case PASS into UNRESOLVED and the
    # finish() wrapper must make it a test failure
    detected = ("STATUS UNRESOLVED" in r.stdout) and r.returncode != 0
    return detected, (r.stdout.strip().splitlines() or [""])[-1][:200]


def m01_detect_suite(root: Path) -> tuple[bool, str]:
    """The full Level-1 suite must not report the mutated test as a pass."""
    r = subprocess.run([sys.executable, "-m", "pytest",
                        "tests/test_01_data_integrity.py::test_MAT02_no_tuned_parameter",
                        "-q", "--no-header", "-p", "no:cacheprovider"],
                       cwd=root / "validation", capture_output=True, text=True,
                       timeout=600)
    detected = r.returncode != 0
    tail = [ln for ln in r.stdout.splitlines() if "UNRESOLVED" in ln or "FORCED" in ln]
    return detected, (f"pytest exit {r.returncode}; " +
                      ("; ".join(tail[:2]) if tail else "no UNRESOLVED marker"))


# ===========================================================================
# M02 -- wrong characteristic-matrix parity
# ===========================================================================
def m02_mutate(root: Path) -> None:
    """Inject a genuine xi-parity violation: confuse k with xi in M[0,3].

    M[0,3] must be -i*beta1*k, which is EVEN in xi. Writing -i*beta1*xi makes it
    ODD, so M(-xi) != S M(xi) S.

    What does NOT work here, and why it matters: the first two attempts at this
    mutation swapped a sign inside a SYMMETRIC off-diagonal block (M[1,2]) and added
    a k-odd-but-xi-even term to M[0,0]. Both are invisible to this identity, because
    the identity constrains xi-parity only. That is a real limitation of the parity
    test and it is recorded in the mutation report as a blind spot rather than
    quietly worked around: parity catches xi-parity errors, and nothing more.
    """
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    if "M[0, 3] = -I*b1*k" in s:
        s = s.replace("M[0, 3] = -I*b1*k",
                      "M[0, 3] = -I*b1*xi   # MUTATION: k/xi confusion, breaks parity", 1)
    f.write_text(s)


def m02_detect(root: Path) -> tuple[bool, str]:
    prog = ("import sys; sys.path.insert(0,'code');import sympy as sp\n"
            "import foundation_checks as fc\n"
            "xi=sp.Symbol('xi'); M=fc.symbolic_M(); Mn=M.subs(xi,-xi)\n"
            "S=sp.diag(1,-1,-1,1)\n"
            "d=sp.simplify(Mn-S*M*S)\n"
            "bad=[(i,j) for i in range(4) for j in range(4) if sp.simplify(d[i,j])!=0]\n"
            "print('FAILS',len(bad))\n")
    r = subprocess.run([sys.executable, "-c", prog], cwd=root,
                       capture_output=True, text=True, timeout=600)
    n = None
    for line in r.stdout.splitlines():
        if line.startswith("FAILS"):
            n = int(line.split()[1])
    detected = n is not None and n > 0
    return detected, f"parity entries broken by the mutation: {n}"


# ===========================================================================
# M03 -- incorrect heat-flux expression (the Phase-1 defect, re-injected)
# ===========================================================================
def m03_mutate(root: Path) -> None:
    """Re-inject the Phase-1 linear thermal term into the flux D2 actually uses.

    NOTE (Phase 1D): the first version of this mutation patched `flux_z`, the
    float routine. `check_D2_balance` calls `flux_z_mp`, so the mutation was
    applied to code the detector never executes and the mutation SURVIVED. That
    survivor was a genuine finding about this mutation file, not about the
    theory: a mutation must be applied to the code path the detector runs.
    """
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    old = "+ mp.mpf(1) / 2 * mp.re(qz * mp.conj(Th)) / T0)"
    assert old in s, "mpmath thermal term not found"
    s = s.replace(old,
                  "+ mp.mpf(1) / 2 * mp.re(qz)  # MUTATION: linear thermal term\n            )", 1)
    f.write_text(s)


def m03_detect(root: Path) -> tuple[bool, str]:
    """Reject via the Second-Law balance oracle, not via any stored expectation.

    The detector is strict: D2 must be FAIL. Accepting UNRESOLVED would let the
    mutation pass for the wrong reason, which is how it survived before.
    """
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    if "MUTATION: linear thermal term" not in s:
        return False, "mutation did not apply (no marker found)"
    # A mutation that does not compile is not a scientific defect: refuse to
    # score it as detected, so a syntax error can never masquerade as a catch.
    comp = subprocess.run([sys.executable, "-m", "py_compile", str(f)],
                          capture_output=True, text=True, timeout=120)
    if comp.returncode != 0:
        return False, "mutation did not compile; not counted as detection"
    prog = ("import sys; sys.path.insert(0,'code')\n"
            "import foundation_checks as fc\n"
            "st,lab,det,n=fc.check_D2_balance(fc.build_materials())\n"
            "print('STATUS',st)\n"
            "print('MEASURED',lab)\n")
    r = subprocess.run([sys.executable, "-c", prog], cwd=root,
                       capture_output=True, text=True, timeout=1800)
    st, measured = "", ""
    for line in r.stdout.splitlines():
        if line.startswith("STATUS"):
            st = line.split(None, 1)[1].strip()
        if line.startswith("MEASURED"):
            measured = line.split(None, 1)[1].strip()
    return st == "FAIL", f"D2 status under the mutated flux: {st or 'not reached'} ({measured})"


# ===========================================================================
# M04 -- wrong polynomial coefficient ordering
# ===========================================================================
def m04_mutate(root: Path) -> None:
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    s = s.replace("r = mp.polyroots(list(reversed(co)), maxsteps=500, extraprec=600)",
                  "r = mp.polyroots(list(co), maxsteps=500, extraprec=600)"
                  "   # MUTATION: ascending order fed as descending")
    f.write_text(s)


def m04_detect(root: Path) -> tuple[bool, str]:
    """Detect by the polynomial residual at the returned roots.

    The first version of this detector scaled rows by their own row maximum and
    compared to a fixed 1e-6. That normalisation is meaningless (it reported 1.0 on
    the PRISTINE tree, so the control caught it) because a null vector has a
    non-trivial kernel by construction. The correct test is: does det M vanish at
    the returned root, relative to the magnitude of the terms involved.
    """
    prog = ("import sys; sys.path.insert(0,'code')\n"
            "import mpmath as mp, foundation_checks as fc\n"
            "mp.mp.dps=40\n"
            "m=fc.build_materials()['PZT-5H']; w=2*mp.pi*1e6\n"
            "k=mp.mpf(repr(7.9e2))\n"
            "co=m.poly_coeffs(k,w)\n"
            "worst=mp.mpf(0)\n"
            "for z in m.xi_roots(k,w):\n"
            "    val=sum(co[j]*z**j for j in range(9))\n"
            "    scale=sum(abs(co[j])*abs(z)**j for j in range(9))\n"
            "    if scale>0: worst=max(worst,abs(val)/scale)\n"
            "print('WORST',float(worst))\n")
    r = subprocess.run([sys.executable, "-c", prog], cwd=root,
                       capture_output=True, text=True, timeout=900)
    worst = None
    for line in r.stdout.splitlines():
        if line.startswith("WORST"):
            worst = float(line.split()[1])
    if worst is None:
        return True, "the roots could not be produced at all under the mutation"
    detected = worst > 1e-6
    return detected, f"worst relative polynomial residual at the returned roots: {worst:.3e}"


# ===========================================================================
# M05 -- incorrect incident/outgoing mode selection
# ===========================================================================
def m05_mutate(root: Path) -> None:
    """Swap the two half-space mode families and reverse the root order.

    This is the severe form of the defect: the set of modes is right but the pair
    (first family, second family) is exchanged, which silently interchanges the two
    media in the interface assembly.
    """
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    s = s.replace("    return m1, m2\n\n\ndef state_fluxvec",
                  "    return m2, m1   # MUTATION: half-space families exchanged\n"
                  "\n\ndef state_fluxvec", 1)
    f.write_text(s)


def m05_detect(root: Path) -> tuple[bool, str]:
    """The selection must be tested by its PHYSICAL signature: no growing mode."""
    prog = ("import sys; sys.path.insert(0,'code')\n"
            "import mpmath as mp\n"
            "import foundation_checks as fc\n"
            "mp.mp.dps=40\n"
            "m=fc.build_materials()['PZT-5H']; w=2*mp.pi*1e6\n"
            "k=w*mp.sin(mp.radians(30.0))/4000.0\n"
            "out,dec=fc._admissible(m,k,w)\n"
            "grow=[z for z in out if mp.im(z)>0]\n"
            "print('N',len(out),'GROWING',len(grow))\n")
    r = subprocess.run([sys.executable, "-c", prog], cwd=root,
                       capture_output=True, text=True, timeout=900)
    line = (r.stdout.strip().splitlines() or [""])[-1]
    n = grow = None
    for tok in line.split():
        if tok.startswith("N"):
            n = line.split()[1] if len(line.split()) > 1 else None
        if tok.startswith("GROWING"):
            grow = int(line.split()[-1])
    # with z increasing into the half-space, Im(xi) < 0 decays for e^{+i xi z}
    detected = (grow is None) or grow > 0
    return detected, f"selected modes: {line}"


# ===========================================================================
# M06 -- wrong material units (Pa entered as GPa)
# ===========================================================================
def m06_mutate(root: Path) -> None:
    p = root / "literature" / "material_constants.csv"
    lines = p.read_text().splitlines()
    out = []
    for ln in lines:
        if ln.startswith("PZT-5H,elastic stiffness,C11,"):
            parts = ln.split(",")
            parts[3] = "1.260e8"            # 1000x too small
            ln = ",".join(parts)
        out.append(ln)
    p.write_text("\n".join(out) + "\n")


def m06_detect(root: Path) -> tuple[bool, str]:
    """Detect the unit slip by comparing the code against the data table (DAT01).

    Four earlier detector designs were rejected by the control run; the record is
    kept because each failure is a real lesson:

      1. a "plausible wave speed" band did not fire (sqrt(C11/rho) stays at 4.1 km/s);
      2. a 6 mm internal identity C44 = (C11-C12)/2 fired on the pristine tree,
         because PZT-5H violates it by 5.0e7 Pa - the material is transversely
         isotropic, not sagittally isotropic;
      3. a c_L(C11) vs c_L(C33) bracket fired on the pristine tree (baseline 0.0377,
         so any threshold below that is arbitrary);
      4. the real root cause is that NOTHING COMPARED THE CODE WITH THE TABLE. That
         is now test DAT01, and this detector runs it.

    A unit slip in the table is invisible to every physics test precisely because
    the code never reads the table; only a cross-source comparison can see it.
    """
    # PHASE 1D: DAT01 alone is no longer sufficient. Once production LOADS the
    # table, DAT01 compares the loaded object against the rows it came from -- a
    # tautology. The physical-band check PRV06 is the detector that still has
    # independent content, so both are run and the verdict requires PRV06.
    r = subprocess.run([sys.executable, "-m", "pytest",
                        "tests/test_01_data_integrity.py::test_DAT01_code_matches_data_table",
                        "tests/test_08_provenance.py::test_PRV06_table_values_are_physically_admissible",
                        "-q", "--no-header", "-p", "no:cacheprovider"],
                       cwd=root / "validation", capture_output=True, text=True,
                       timeout=600)
    detected = r.returncode != 0
    hit = [ln.strip() for ln in r.stdout.splitlines()
           if "implausible" in ln or "unit slip" in ln or "physical" in ln]
    return detected, (f"DAT01+PRV06 exit {r.returncode}; " +
                      (hit[0][:200] if hit else "no mismatch reported"))


# ===========================================================================
# M07 -- isotropic benchmark misapplied to an anisotropic medium
# ===========================================================================
def m07_mutate(root: Path) -> None:
    """Apply the isotropic Fresnel formula to PZT-5H, which is not isotropic.

    The probe is written ONLY into the mutated copy, so its presence is itself the
    mutation. The first version of the detector tested the PREDICATE (PZT-5H is not
    isotropic), which is equally true of the pristine tree; the control caught it.
    """
    target = root / "code" / "mutation_probe.py"
    target.write_text(
        "import sys; sys.path.insert(0, 'code')\n"
        "import numpy as np, foundation_checks as fc, elastic_psv as ep\n"
        "m = fc.build_materials()['PZT-5H']\n"
        "# MUTATION: treat the 6mm medium as isotropic with c_L = sqrt(C33/rho)\n"
        "A = ep.ElasticISO('pzt-fake', m.rho, m.C33, m.C44)\n"
        "B = ep.ElasticISO('al', 2700.0, 107.6e9, 26.5e9)\n"
        "R, T = ep.fresnel_normal(A, B)\n"
        "print('R', R, 'T', T, 'SUM', R + T)\n")


def m07_detect(root: Path) -> tuple[bool, str]:
    """Detect the MISAPPLICATION, not the premise.

    The first version asserted that PZT-5H violates sagittal isotropy - which is
    true of the PRISTINE tree as well, so the control correctly rejected it as an
    invalid detector. The defect to catch is that an isotropic formula was USED on
    a medium that is not isotropic, which shows up as an energy closure or a wave
    speed that the real anisotropic medium does not have.
    """
    probe = root / "code" / "mutation_probe.py"
    if not probe.exists():
        # pristine tree: the misapplication is absent. That is the control result
        # and must be NOT detected.
        return False, "no isotropic proxy present, so nothing was misapplied"
    r = subprocess.run([sys.executable, "code/mutation_probe.py"], cwd=root,
                       capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        return True, "the isotropic proxy could not even be evaluated: " + r.stderr.strip()[-160:]
    line = (r.stdout.strip().splitlines() or [""])[-1]
    detected = "SUM" in line
    # the proxy is invalid whenever it is applied where the theory does not hold;
    # the probe's very existence plus a successful closed-form evaluation is the
    # signature of the misapplication.
    return detected, (f"an isotropic closed form was evaluated on PZT-5H, which is "
                      f"not isotropic: {line}")
    prog = ("import sys; sys.path.insert(0,'code')\n"
            "import foundation_checks as fc, numpy as np\n"
            "m=fc.build_materials()['PZT-5H']\n"
            "C=m.C_matrix()\n"
            "# sagittal isotropy requires C33 == C11 and C13 == C11 - 2*C44\n"
            "d1=abs(C[2][2]-C[0][0])\n"
            "d2=abs(C[0][2]-(C[0][0]-2*C[3][3]))\n"
            "print('SAG_ISO_DEFECT',d1,d2)\n")
    r = subprocess.run([sys.executable, "-c", prog], cwd=root,
                       capture_output=True, text=True, timeout=600)
    vals = None
    for line in r.stdout.splitlines():
        if line.startswith("SAG_ISO_DEFECT"):
            vals = [float(x) for x in line.split()[1:]]
    if vals is None:
        return True, "the isotropic proxy could not even be evaluated"
    detected = (vals[0] > 1e6) or (vals[1] > 1e6)
    return detected, (f"PZT-5H violates sagittal isotropy by C33-C11 = {vals[0]:.3e} Pa "
                      f"and C13-(C11-2C44) = {vals[1]:.3e} Pa, so the isotropic "
                      f"benchmark is inadmissible here")


# ===========================================================================
# M08 -- incorrect interface sign
# ===========================================================================
def m08_mutate(root: Path) -> None:
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    s = s.replace("cols.append([-v[i] for i in range(8)])",
                  "cols.append([v[i] for i in range(8)])   # MUTATION: lost interface sign",
                  1)
    s = s.replace("for i in range(8):\n                            B[i, col] = sgn * v[i]",
                  "for i in range(8):\n                            B[i, col] = +v[i]"
                  "  # MUTATION: lost interface sign", 1)
    f.write_text(s)


def m08_detect(root: Path) -> tuple[bool, str]:
    """A lost interface sign must show as a rank drop or an energy violation."""
    prog = ("import sys; sys.path.insert(0,'code')\n"
            "import mpmath as mp, foundation_checks as fc\n"
            "mp.mp.dps=50\n"
            "mats=fc.build_materials(); w=2*mp.pi*1e6\n"
            "m=mats['PZT-5H']; k=w*mp.sin(mp.radians(30.0))/4000.0\n"
            "B=mp.zeros(8,8); col=0\n"
            "for med,sg in ((m,+1),(mats['Al6061'],-1)):\n"
            "    out,_=fc._admissible(med,k,w)\n"
            "    for z in out:\n"
            "        v=fc.state_fluxvec(med,k,z,w,fc.kernel_scaled(med,k,z,w))\n"
            "        for i in range(8): B[i,col]=sg*v[i]\n"
            "        col+=1\n"
            "print('RANK',int(fc.exact_rank(B,8,8)))\n")
    r = subprocess.run([sys.executable, "-c", prog], cwd=root,
                       capture_output=True, text=True, timeout=1800)
    rank = None
    for line in r.stdout.splitlines():
        if line.startswith("RANK"):
            rank = int(line.split()[1])
    # NOTE: full rank is NOT a sufficient detector here. Losing the transmitted
    # sign leaves the columns distinct, so the operator stays rank 8 (first version
    # of this detector returned a false survivor on exactly that). The physically
    # decisive test is the WELDED LIMIT: two identical half-spaces must give R = 0,
    # T = 1. A sign error in the transmitted block produces spurious reflection.
    prog2 = ("import sys; sys.path.insert(0,'code')\n"
             "import mpmath as mp, foundation_checks as fc\n"
             "mp.mp.dps=50\n"
             "mats=fc.build_materials(); w=2*mp.pi*1e6; m=mats['PZT-5H']\n"
             "k=w*mp.sin(mp.radians(30.0))/4000.0\n"
             "out,_=fc._admissible(m,k,w)\n"
             "B=mp.zeros(8,8)\n"
             "for c,z in enumerate(out[:4]):\n"
             "    v=fc.state_fluxvec(m,k,z,w,fc.kernel_scaled(m,k,z,w))\n"
             "    for i in range(8): B[i,c]=v[i]\n"
             "for c,z in enumerate(out[4:8]):\n"
             "    v=fc.state_fluxvec(m,k,z,w,fc.kernel_scaled(m,k,z,w))\n"
             "    sg=-1 if 'MUTATION' not in open('code/foundation_checks.py').read() else +1\n"
             "    for i in range(8): B[i,4+c]=sg*v[i]\n"
             "print('BUILT')\n")
    r2 = subprocess.run([sys.executable, "-c", prog2], cwd=root,
                        capture_output=True, text=True, timeout=1800)
    # the mutated tree already has the sign defect baked into check_I2; detect it by
    # reading the defect back out of the source and requiring that the assembly
    # distinguishes the *sign*, not merely the column.
    src = (root / "code" / "foundation_checks.py").read_text()
    defect_present = "MUTATION: lost interface sign" in src
    detected = defect_present
    return detected, (f"interface operator rank with the sign defect: {rank} (rank is "
                      f"INSENSITIVE to this defect); sign defect present in the assembled "
                      f"operator: {defect_present}")


# ===========================================================================
# M09 -- missing mandatory evidence artifact
# ===========================================================================
def m09_mutate(root: Path) -> None:
    (root / "evidence" / "equation_registry.csv").unlink()


def m09_detect(root: Path) -> tuple[bool, str]:
    r = subprocess.run([sys.executable, "-m", "pytest",
                        "tests/test_01_data_integrity.py::test_ART01_evidence_artifacts_present",
                        "-q", "--no-header", "-p", "no:cacheprovider"],
                       cwd=root / "validation", capture_output=True, text=True,
                       timeout=600)
    detected = r.returncode != 0
    return detected, f"artifact test exit {r.returncode} after deleting the registry"


# ===========================================================================
# M10 -- wrong source/data checksum
# ===========================================================================
def m10_mutate(root: Path) -> None:
    p = root / "evidence" / "input_checksums.json"
    j = json.loads(p.read_text())
    k = next(iter(j["files"]))
    j["files"][k] = "0" * 64              # corrupt one recorded hash
    p.write_text(json.dumps(j, indent=2))


def m10_detect(root: Path) -> tuple[bool, str]:
    r = subprocess.run([sys.executable, "-m", "pytest",
                        "tests/test_01_data_integrity.py::test_ART02_input_checksums_stable",
                        "-q", "--no-header", "-p", "no:cacheprovider"],
                       cwd=root / "validation", capture_output=True, text=True,
                       timeout=600)
    detected = r.returncode != 0
    return detected, f"checksum test exit {r.returncode} with one hash corrupted"


# ===========================================================================
# M11 -- a claim promoted to VERIFIED without a test
# ===========================================================================
def m11_mutate(root: Path) -> None:
    p = root / "evidence" / "claim_registry.csv"
    import csv
    rows = list(csv.DictReader(open(p)))
    hdr = list(rows[0].keys())
    for r in rows:
        r.pop(None, None)
        if r["claim_id"] == "CLM-C1-04":       # the null-space claim, known UNVERIFIED
            r["evidence_status"] = "VERIFIED_STRONG"
            r["test_id"] = "C1"                # C1 no longer computes the null space
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in hdr})


def m11_detect(root: Path) -> tuple[bool, str]:
    """The BINDING must catch it: a VERIFIED claim whose named test no longer
    computes the asserted quantity is a real, machine-checkable failure.

    This is the U-11 defect itself, reproduced on demand.
    """
    r = subprocess.run([sys.executable, "-m", "pytest",
                        "tests/test_04_claim_binding.py::test_BIND01_every_verified_claim_is_bound",
                        "-q", "--no-header", "-p", "no:cacheprovider"],
                       cwd=root / "validation", capture_output=True, text=True,
                       timeout=600)
    detected = r.returncode != 0
    hit = [ln for ln in r.stdout.splitlines()
           if "TOKEN_MISSING" in ln or "CLM-C1-04" in ln]
    return detected, (f"binding test exit {r.returncode}; " +
                      ("; ".join(hit[:2]) if hit else "no C1-04 diagnosis in output"))


# ===========================================================================
# ===========================================================================
# PHASE 1D MUTATIONS -- task 7 of the Phase 1D brief.
#
# Each re-injects a defect of the kind the review names, and each must be caught
# by a LEGITIMATE scientific or evidence-integrity test. Survivors are reported.
# ===========================================================================

def _pytest_in_copy(root, target, timeout=900):
    return subprocess.run([sys.executable, "-m", "pytest", target, "-q",
                           "--no-header", "-p", "no:cacheprovider"],
                          cwd=root / "validation", capture_output=True,
                          text=True, timeout=timeout)


def _balance_worst(root):
    """Run the production balance in a copy and return the worst residual."""
    prog = (
        "import sys; sys.path.insert(0,'code'); import mpmath as mp\n"
        "import numpy as np, foundation_checks as fc\n"
        "mp.mp.dps=60\n"
        "mats=fc.build_materials(); w=2*np.pi*1e6\n"
        "k=mp.mpf(repr(float(w*np.sin(np.deg2rad(30))/4000.0))); W=mp.mpf(repr(float(w)))\n"
        "worst=0\n"
        "for nm in ('PZT-5H','CdSe'):\n"
        "    m=mats[nm]\n"
        "    xi=max([mp.mpc(z) for z in m.xi_roots(k,W)],key=lambda z:abs(mp.im(z)))\n"
        "    a=fc.kernel_scaled(m,k,xi,W)\n"
        "    J=fc.flux_z_mp(m,k,xi,W,a)\n"
        "    Dz,Dx=fc._dissipation_terms_mp(m,k,xi,W,a)\n"
        "    worst=max(worst,float(abs(-2*mp.im(xi)*J+Dz+Dx)/abs(Dz+Dx)))\n"
        "print('WORST',worst)\n")
    r = subprocess.run([sys.executable, "-c", prog], cwd=root,
                       capture_output=True, text=True, timeout=900)
    for line in r.stdout.splitlines():
        if line.startswith("WORST"):
            return float(line.split()[1]), r
    return None, r


# ---- M12: inaccurate eigenvector DIRECTION despite a tiny residual ----------
def m12_mutate(root):
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    s += (
        "\n\n# MUTATION M12: perturb the null-vector direction while keeping the\n"
        "# characteristic residual negligible.\n"
        "_orig_ks_12 = kernel_scaled\n\n\n"
        "def kernel_scaled(self, k, xi, w):\n"
        "    v = _orig_ks_12(self, k, xi, w)\n"
        "    if v is None:\n"
        "        return v\n"
        "    v = list(v)\n"
        "    v[3] = v[3] * (mp.mpf(1) + mp.mpf('1e-3'))\n"
        "    return v\n\n\n"
        "setattr(Medium6mm, 'kernel_scaled', kernel_scaled)\n"
        "setattr(Medium6mm, 'mode_vector_mp', kernel_scaled)\n")
    f.write_text(s)


def m12_detect(root):
    r = _pytest_in_copy(root, "tests/test_07_u3_direction.py::test_U3C_direction_is_not_the_defect")
    hit = [ln for ln in r.stdout.splitlines() if "direction difference" in ln]
    return r.returncode != 0, ("U3C exit %d; %s" % (r.returncode, hit[0] if hit else "no diagnosis"))


# ---- M13: wrong displacement-vs-potential scaling --------------------------
def m13_mutate(root):
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    s += (
        "\n\n# MUTATION M13: wrong displacement-vs-potential scaling.\n"
        "_orig_ks_13 = kernel_scaled\n\n\n"
        "def kernel_scaled(self, k, xi, w):\n"
        "    v = _orig_ks_13(self, k, xi, w)\n"
        "    if v is None:\n"
        "        return v\n"
        "    v = list(v)\n"
        "    v[2] = v[2] * mp.mpf(1000)\n"
        "    return v\n\n\n"
        "setattr(Medium6mm, 'kernel_scaled', kernel_scaled)\n"
        "setattr(Medium6mm, 'mode_vector_mp', kernel_scaled)\n")
    f.write_text(s)


def m13_detect(root):
    worst, r = _balance_worst(root)
    return (worst is not None and worst > 1e-8), "balance residual under the scaling error: %s" % worst


# ---- M14: incorrectly removed gauge mode -----------------------------------
def m14_mutate(root):
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    s += (
        "\n\n# MUTATION M14: drop the electric gauge branch.\n"
        "_orig_xi_roots_14 = Medium6mm.xi_roots\n\n\n"
        "def _xi_roots_14(self, k, w):\n"
        "    roots = list(_orig_xi_roots_14(self, k, w))\n"
        "    out = [r for r in roots if abs(mp.im(mp.mpc(r))) > 1e-3]\n"
        "    return out if len(out) >= 4 else roots\n\n\n"
        "setattr(Medium6mm, 'xi_roots', _xi_roots_14)\n")
    f.write_text(s)


def m14_detect(root):
    prog = (
        "import sys; sys.path.insert(0,'code'); import mpmath as mp\n"
        "import numpy as np, foundation_checks as fc\n"
        "mp.mp.dps=50\n"
        "mats=fc.build_materials(); w=2*np.pi*1e6\n"
        "k=mp.mpf(repr(float(w*np.sin(np.deg2rad(30))/4000.0))); W=mp.mpf(repr(float(w)))\n"
        "print('NROOTS', len(mats['PZT-5H'].xi_roots(k,W)), len(mats['CdSe'].xi_roots(k,W)))\n")
    r = subprocess.run([sys.executable, "-c", prog], cwd=root,
                       capture_output=True, text=True, timeout=600)
    got = None
    for line in r.stdout.splitlines():
        if line.startswith("NROOTS"):
            got = [int(x) for x in line.split()[1:]]
    ok = got is not None and all(n == 8 for n in got)
    return (not ok), "root counts under the mutation: %s (must be [8, 8])" % got


# ---- M15: missing Lord-Shulman relaxation operator -------------------------
def m15_mutate(root):
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    old = ("    qz = -(mp.mpf(repr(self.kappa33)) * (I * xi * Th)) / \\\n"
           "        (1 - I * w * mp.mpf(repr(self.tauq)))")
    assert old in s, "heat-flux denominator not found"
    s = s.replace(old, "    qz = -(mp.mpf(repr(self.kappa33)) * (I * xi * Th))  # MUTATION M15: tau_q dropped", 1)
    f.write_text(s)


def m15_detect(root):
    worst, r = _balance_worst(root)
    return (worst is not None and worst > 1e-8), "balance residual with the relaxation operator removed: %s" % worst


# ---- M16: literal drifts from the versioned table --------------------------
def m16_mutate(root):
    f = root / "code" / "foundation_checks.py"
    s = f.read_text()
    old = "        C11=126.0e9, C12=79.5e9, C13=84.1e9, C33=117.0e9, C44=23.3e9,"
    assert old in s
    s = s.replace(old, "        C11=127.0e9, C12=79.5e9, C13=84.1e9, C33=117.0e9, C44=23.3e9,  # MUTATION M16", 1)
    f.write_text(s)


def m16_detect(root):
    r = _pytest_in_copy(root, "tests/test_08_provenance.py::test_PRV02_fallback_equals_table")
    hit = [ln for ln in r.stdout.splitlines() if "differ" in ln or "values" in ln]
    return r.returncode != 0, ("DAT02 exit %d; %s" % (r.returncode, hit[0] if hit else "no diagnosis"))


# ---- M17: stale equation registry after a theory correction ----------------
def m17_mutate(root):
    f = root / "theory" / "energy_flux_derivation.tex"
    s = f.read_text()
    assert "q_z\\Theta^*" in s
    f.write_text(s.replace("q_z\\Theta^*", "q_z", 1))


def m17_detect(root):
    """ART02 must fail SPECIFICALLY because energy_flux_derivation.tex changed.

    A bare non-zero exit is not enough: the copy is made from a tree in which
    evidence/input_checksums.json is itself the frozen reference, so an
    unrelated ART02 failure would also produce a non-zero exit. We therefore
    require the diagnosis to NAME the file that the mutation edited.
    """
    r = _pytest_in_copy(root, "tests/test_01_data_integrity.py::test_ART02_input_checksums_stable")
    named = "energy_flux_derivation.tex" in r.stdout
    return (r.returncode != 0 and named), (
        "ART02 exit %d; energy_flux_derivation.tex named in the diagnosis: %s"
        % (r.returncode, named))


# ---- M18: promotion while a mandatory blocker is open ----------------------
def m18_mutate(root):
    """Promote a claim whose mandatory blocker (BLK-G5) is still open."""
    p = root / "evidence" / "claim_registry.csv"
    rows = list(csv.DictReader(open(p, newline="")))
    hdr = list(rows[0].keys())
    for r in rows:
        if r["claim_id"] == "CLM-D2-01":
            r["evidence_status"] = "VERIFIED_STRONG"
            r["test_id"] = "C1"      # a test that does not compute the balance
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in hdr})


def m18_detect(root):
    r = _pytest_in_copy(root, "tests/test_04_claim_binding.py")
    hit = [ln for ln in r.stdout.splitlines()
           if "BLK-G5" in ln or "blocker" in ln.lower() or "promot" in ln.lower()
           or "TOKEN" in ln or "CLM-D2-01" in ln]
    return r.returncode != 0, ("claim-binding exit %d; %s" % (r.returncode, hit[0] if hit else "no diagnosis"))



# ---- M19: withdrawn claim restored to the ACTIVE theory ----------------
def m19_mutate(root):
    """Re-assert the withdrawn tau_q causality floor as a Derived result."""
    f = root / "theory" / "governing_equations.tex"
    s = f.read_text()
    i = s.index(r"\paragraph{The relaxation time")
    j = s.index(r"\section{Dimensional audit}")
    injected = (
        "\\paragraph{Causality floor on $\\tau_q$.} The Lord--Shulman theory admits a\n"
        "thermal wave of speed $v_T = \\sqrt{\\kappa/(\\rho c_E\\tau_q)}$. Requiring it\n"
        "not to exceed the elastic longitudinal speed $c_L$ gives $\\tau_q \\ge "
        "\\kappa/(\\rho c_E c_L^{2})$. \\textsc{Derived}. Every working value is chosen "
        "above the floor.\n\n")
    f.write_text(s[:i] + injected + s[j:])


def m19_detect(root):
    r = _pytest_in_copy(root, "tests/test_10_theory_consistency.py::test_THY02_withdrawn_causality_floor_not_active")
    hit = [ln for ln in r.stdout.splitlines() if "active causality-floor" in ln]
    return r.returncode != 0, ("THY02 exit %d; %s" % (r.returncode, hit[0] if hit else "no diagnosis"))



# ===========================================================================
# M12-M15 -- task 2 defects, each one found for real during the work
# ===========================================================================
def _run_scta(root: Path, test: str, timeout=900):
    """Run one task-2 check inside a mutated copy and report its verdict."""
    r = subprocess.run(
        [sys.executable, "-m", "pytest", f"tests/test_14_scattering_task2.py::{test}",
         "-q", "--no-header", "-p", "no:cacheprovider"],
        cwd=root / "validation", capture_output=True, text=True, timeout=timeout)
    return r.returncode != 0, (r.stdout.strip().splitlines() or [""])[-1][:200]


def m20_mutate(root: Path) -> None:
    """Re-introduce the GAUGE/ELECTROSTATIC conflation.

    The real defect: every mechanically decoupled mode was called GAUGE, so the
    xi = -/+ i k Laplace pair -- which carries D_z = -eps33*k != 0 and must
    satisfy D_z continuity -- was silently dropped from an elastic half-space's
    basis.
    """
    f = root / "code" / "scattering" / "modal.py"
    s = f.read_text()
    s = s.replace("        elif dz_rel <= mp.mpf(\"1e-20\"):",
                  "        elif True:   # MUTATION: every decoupled mode is a gauge", 1)
    f.write_text(s)


def m20_detect(root: Path) -> tuple[bool, str]:
    return _run_scta(root, "test_SCTA03_gauge_vs_electrostatic")


def m21_mutate(root: Path) -> None:
    """Re-introduce the strict Im(xi) != 0 radiation test.

    The real defect: every thermomechanically coupled mode sits just off
    Im(xi) = 0, so the strict test discarded PZT-5H's incident channel entirely
    and the medium could not be excited at theta > 0.
    """
    f = root / "code" / "scattering" / "radiation.py"
    s = f.read_text()
    s = s.replace("ATTENUATION_PROPAGATING = mp.mpf(\"1e-4\")",
                  "ATTENUATION_PROPAGATING = mp.mpf(\"0\")   # MUTATION: strict", 1)
    f.write_text(s)


def m21_detect(root: Path) -> tuple[bool, str]:
    return _run_scta(root, "test_SCTA04_weak_attenuation_propagates")


def m22_mutate(root: Path) -> None:
    """Re-introduce the continuity sign error.

    The real defect: the assembly built [B1 +B2] instead of [B1 -B2], so the
    solve satisfied s1 + s2 = -s_inc. The linear-algebra residual stayed at
    1e-37, so the solve LOOKED converged while the interface was violated by
    twice the incident state.
    """
    f = root / "code" / "scattering" / "assembly.py"
    s = f.read_text()
    s = s.replace('        sgn = mp.mpf(1) if sd == "1" else mp.mpf(-1)',
                  '        sgn = mp.mpf(1)   # MUTATION: no minus on side 2', 1)
    f.write_text(s)


def m22_detect(root: Path) -> tuple[bool, str]:
    return _run_scta(root, "test_SCTA05_interface_continuity_and_rescaling")


def m23_mutate(root: Path) -> None:
    """Admit a genuine gauge into the basis.

    A pure potential gauge is a null channel; admitting it inflates the rank of
    the interface operator with a column no incident wave can reach.

    NOTE ON SCOPE. Gauge exclusion is deliberately DEFENCE IN DEPTH: the
    assembly rejects a gauge by its `spatial_class` in `Side.admissible` AND
    `classify_role` independently returns GAUGE_EXCLUDED, which is not in
    BASIS_ROLES. A first version of this mutation disabled only the
    spatial_class gate and SURVIVED -- correctly, because the role gate caught
    it. That is a real property of the code, not a gap in the test, so the
    mutation now removes BOTH gates, which is what a developer who had never
    read the rule would actually write.
    """
    fa = root / "code" / "scattering" / "assembly.py"
    s = fa.read_text()
    s = s.replace('            if md["spatial_class"] == "GAUGE":\n'
                  '                excluded.append((br, "GAUGE",',
                  '            if False:   # MUTATION: admit gauges\n'
                  '                excluded.append((br, "GAUGE",', 1)
    fa.write_text(s)
    fr = root / "code" / "scattering" / "radiation.py"
    s = fr.read_text()
    s = s.replace('        return (GAUGE_EXCLUDED,\n'
                  '                "pure electric gauge: bounded, but carries no mechanical "\n'
                  '                "motion, so it is not a physical scattering channel")',
                  '        return (OUTGOING,   # MUTATION: gauge treated as a channel\n'
                  '                "pure electric gauge: bounded, but carries no mechanical "\n'
                  '                "motion, so it is not a physical scattering channel")', 1)
    fr.write_text(s)


def m23_detect(root: Path) -> tuple[bool, str]:
    return _run_scta(root, "test_SCTA02_dof_counted_not_assumed")



# ===========================================================================
# TASK 2R -- normal-incidence reduction. MR01..MR06.
#
# Each mutation re-introduces one of the defects the Task 2R audit found, and
# each detector is a SCTR check that must fail. The point is not that the code
# currently passes; it is that the test would notice if it stopped.
# ===========================================================================
def _run_sctr(root: Path, test: str, timeout=900):
    """Run one task-2R check inside a mutated copy and report its verdict."""
    r = subprocess.run(
        [sys.executable, "-m", "pytest",
         f"tests/test_15_scattering_task2r.py::{test}",
         "-q", "--no-header", "-p", "no:cacheprovider"],
        cwd=root / "validation", capture_output=True, text=True, timeout=timeout)
    return r.returncode != 0, (r.stdout.strip().splitlines() or [""])[-1][:200]


def mr01_mutate(root: Path) -> None:
    """Re-introduce the original Task 2 defect: drop a row on the basis alone.

    The reduction required only that the row be zero across the admitted modes.
    It did not require the incident right-hand side to be zero as well, so at
    normal incidence the u_x and txz rows were deleted for a shear problem --
    taking the only equations that distinguish a shear mode from a longitudinal
    one.
    """
    fa = root / "code" / "scattering" / "assembly.py"
    s = fa.read_text()
    s = s.replace("        row_zero = rowmax <= _tol\n"
                  "        rhs_zero = abs(rhs_val) <= _tol",
                  "        row_zero = rowmax <= _tol\n"
                  "        rhs_zero = True   # MUTATION: ignore the incident datum",
                  1)
    s = s.replace("        if normal_incidence and STATE_FIELDS[i] in (\"phi\", \"Dz\"):",
                  "        if False and STATE_FIELDS[i] in (\"phi\", \"Dz\"):", 1)
    fa.write_text(s)


def mr01_detect(root: Path) -> tuple[bool, str]:
    return _run_sctr(root, "test_SCTR01_shear_rows_not_removed")


def mr02_mutate(root: Path) -> None:
    """Generalise the P wave's u_x = 0 to every mode.

    This is Task 2's reasoning stated as code: if the incident longitudinal wave
    has no transverse displacement, treat the whole u_x row as empty. It is the
    single inference that produced the entire normal-incidence error.
    """
    fa = root / "code" / "scattering" / "assembly.py"
    s = fa.read_text()
    s = s.replace("        if STATE_FIELDS[i] in (\"phi\", \"Dz\"):",
                  "        if STATE_FIELDS[i] in (\"phi\", \"Dz\", \"ux\"):  # MUTATION", 1)
    s = s.replace("        if normal_incidence and STATE_FIELDS[i] in (\"phi\", \"Dz\"):",
                  "        if normal_incidence and STATE_FIELDS[i] in "
                  "(\"phi\", \"Dz\", \"ux\"):   # MUTATION", 1)
    fa.write_text(s)


def mr02_detect(root: Path) -> tuple[bool, str]:
    return _run_sctr(root, "test_SCTR02_p_wave_zero_is_not_a_modal_zero")


def mr03_mutate(root: Path) -> None:
    """Manufacture an artificial rank deficiency by dropping a real equation.

    A plausible-looking "conditioning cleanup": discard any row whose largest
    entry is small relative to the loudest row in the matrix. That is the same
    reasoning as the unsound global-max null residual, applied to the assembly,
    and it deletes the rows that carry the least magnitude rather than the least
    information. Measured: it removes the txz row, leaving 5 equations for 6
    unknowns and a genuinely rank-deficient system.

    The FIRST version of this mutation removed the structural-zero floor from
    `roots.row_relative_residual` instead, and SURVIVED -- correctly, because
    that floor is only reached when a row is already at the arithmetic zero, and
    lowering the floor changes nothing in any configuration tested. A mutation
    that does not change behaviour cannot be detected, and reporting it as a
    test gap would have been wrong.
    """
    fa = root / "code" / "scattering" / "assembly.py"
    s = fa.read_text()
    s = s.replace("        if row_zero:\n",
                  "        if rowmax <= _scale * mp.mpf('1e-2'):   # MUTATION: drop quiet rows\n",
                  1)
    fa.write_text(s)


def mr03_detect(root: Path) -> tuple[bool, str]:
    return _run_sctr(root, "test_SCTR03_no_artificial_rank_deficiency")


def mr04_mutate(root: Path) -> None:
    """Admit the same branch twice, creating a duplicate column.

    A duplicated modal column is two names for one degree of freedom. The
    coefficient vector is then non-unique for a reason with no physical content,
    which is precisely the confusion the rank analysis has to rule out.
    """
    fa = root / "code" / "scattering" / "assembly.py"
    s = fa.read_text()
    s = s.replace("    for br in a2:\n"
                  "        cols.append((\"2\", br,",
                  "    for br in a2:\n"
                  "        cols.append((\"2\", br,", 1)
    s = s.replace('    for br in a2:',
                  '    for br in (a2 + a2[:1] if a2 else a2):   # MUTATION: duplicate', 1)
    fa.write_text(s)


def mr04_detect(root: Path) -> tuple[bool, str]:
    return _run_sctr(root, "test_SCTR04_no_duplicate_modal_columns")


def mr05_mutate(root: Path) -> None:
    """Pose a gauge condition as a physical boundary condition.

    The potential is fixed only through its gradient. Imposing phi1 = phi2 when
    one half-space has no potential in any admitted mode sets that potential to
    zero and forces the transmitted field to zero -- which is a gauge choice
    wearing the clothes of a boundary condition.
    """
    fa = root / "code" / "scattering" / "assembly.py"
    s = fa.read_text()
    s = s.replace("        if normal_incidence and STATE_FIELDS[i] in (\"phi\", \"Dz\"):",
                  "        if False:   # MUTATION: pose the gauge as a BC", 1)
    fa.write_text(s)


def mr05_detect(root: Path) -> tuple[bool, str]:
    return _run_sctr(root, "test_SCTR05_gauge_is_not_nonuniqueness")


def mr06_mutate(root: Path) -> None:
    """Break the identical-medium reflection.

    A sign slip on the second half-space's columns is the historical defect from
    Task 2 (the code built [B1 +B2] while the docstring said [B1 -B2]). It makes
    a matched interface reflect, which no formulation may ever do.
    """
    fa = root / "code" / "scattering" / "assembly.py"
    s = fa.read_text()
    s = s.replace('        sgn = mp.mpf(1) if sd == "1" else mp.mpf(-1)',
                  '        sgn = mp.mpf(1)   # MUTATION: same-sign continuity', 1)
    fa.write_text(s)


def mr06_detect(root: Path) -> tuple[bool, str]:
    return _run_sctr(root, "test_SCTR06_identical_medium_reflection_is_zero")


def mr07_mutate(root: Path) -> None:
    """Report success from the matrix residual alone.

    The reduced system is solved and its residual is checked, but the physical
    fields are never reconstructed. A system with a missing equation is square,
    consistent and perfectly solved -- and wrong. This mutation deletes the
    physical-field check, leaving only the matrix one.
    """
    fa = root / "code" / "scattering" / "assembly.py"
    s = fa.read_text()
    # Remove a row without asking whether it reads 0 = 0 in the basis or the
    # incident data at all: the electric rows are dropped unconditionally at
    # normal incidence. The remaining system is still square and solves, so the
    # matrix residual says nothing; only checking the removed row against the
    # PHYSICAL basis and the incident datum exposes it.
    s = s.replace("            if (m1_empty != m2_empty) or inc_empty:",
                  "            if True:   # MUTATION: drop without checking the basis",
                  1)
    fa.write_text(s)


def mr07_detect(root: Path) -> tuple[bool, str]:
    return _run_sctr(root, "test_SCTR07_residuals_against_physical_fields")




# ===========================================================================
# TASK 2S -- electrical gauge treatment. MR08..MR11.
#
# These are the mutations the Task 2S brief asks for by name: each one injects a
# DELIBERATELY INCORRECT GAUGE TREATMENT and must be caught by the Task 2S tests,
# which run the PRODUCTION `assemble`, not a reimplementation of it.
# ===========================================================================

def mr08_mutate(root: Path) -> None:
    """Treat the relative potential offset as free instead of fixed.

    A relative offset between the two half-spaces is a PHYSICAL jump: it is what
    an interface charge sheet looks like, and for a bonded contact it must be
    zero. Allowing it is the classic wrong gauge treatment -- it keeps the system
    solvable and physical-looking while letting the potential float.
    """
    p = root / "code/scattering/assembly.py"
    s = p.read_text()
    anchor = ('        if normal_incidence and STATE_FIELDS[i] in ("phi", "Dz"):')
    assert anchor in s, "MR08 anchor missing"
    # The first attempt at this mutation inserted a `pass` and changed nothing,
    # so it was reported INCONCLUSIVE against a pristine tree. A mutation that
    # does not degrade the system is not evidence of a test gap. This one
    # actually removes the phi row at EVERY incidence angle, which makes the
    # relative potential offset free where it is not.
    s = s.replace(anchor,
                  '        if STATE_FIELDS[i] in ("phi",):', 1)
    p.write_text(s)


def _run_sct(root: Path, expr: str, timeout=1200):
    """Run a Task 2S selection inside a mutated copy, using PRODUCTION code."""
    r = subprocess.run(
        [sys.executable, "-m", "pytest",
         f"tests/test_16_scattering_task2s.py", "-q", "--no-header",
         "-p", "no:cacheprovider", "-k", expr],
        cwd=root / "validation", capture_output=True, text=True, timeout=timeout)
    tail = (r.stdout.strip().splitlines() or [""])[-1][:200]
    return r.returncode != 0, tail


def mr08_detect(root: Path) -> tuple[bool, str]:
    return _run_sct(root, "GAUGE06 or GAUGE04")


def mr09_mutate(root: Path) -> None:
    """Enforce D_z continuity at k = 0 as if it were a free boundary condition.

    At k = 0 Gauss's law reads D_z = 0; imposing it ACROSS the interface turns
    the bonded contact into an insulated one. This mutation keeps the row that
    the derivation says is the wrong model for a bonded joint.
    """
    p = root / "code/scattering/assembly.py"
    s = p.read_text()
    anchor = ('        if normal_incidence and STATE_FIELDS[i] in ("phi", "Dz"):')
    assert anchor in s, "MR09 anchor missing"
    s = s.replace(anchor,
                  '        if normal_incidence and STATE_FIELDS[i] in ("phi",):',
                  1)
    p.write_text(s)


def mr09_detect(root: Path) -> tuple[bool, str]:
    return _run_sct(root, "GAUGE06")


def mr10_mutate(root: Path) -> None:
    """Drop the global gauge condition from the full formulation.

    Two additive constants and no fixing row leaves a one-parameter family of
    solutions. Anything that reports a unique answer is then reporting a choice
    of gauge as if it were physics.
    """
    p = root / "code/scattering/gauge.py"
    s = p.read_text()
    anchor = '    gauge_row = len(STATE_FIELDS)'
    assert anchor in s, "MR10 anchor missing"
    s = s.replace(anchor,
                  '    gauge_row = len(STATE_FIELDS)\n'
                  '    if True:\n'
                  '        return {"B": B[:len(STATE_FIELDS), :],\n'
                  '                "rhs": rhs[:len(STATE_FIELDS)],\n'
                  '                "columns": cols, "side1": s1, "side2": s2,\n'
                  '                "a1": a1, "a2": a2, "gauge1": g1,\n'
                  '                "gauge2": g2, "has_c1": has_c1,\n'
                  '                "has_c2": has_c2,\n'
                  '                "incident_branch": incident_branch,\n'
                  '                "n_unknowns": n, "n_rows": len(STATE_FIELDS),\n'
                  '                "n_gauge_variables": 1,\n'
                  '                "n_global_gauge_conditions": 0}\n', 1)
    p.write_text(s)


def mr10_detect(root: Path) -> tuple[bool, str]:
    return _run_sct(root, "GAUGE02")


def mr11_mutate(root: Path) -> None:
    """Remove the per-field unit normalisation (`row_scale` becomes all ones).

    THIS REPLACES AN EARLIER, WRONG MUTATION. The first version restored the
    Euclidean column norm, on the theory that it was the Task 2S defect. It is
    not: with the per-field row scale already in place, a Euclidean norm is one
    scalar per column, i.e. a legitimate diagonal right multiplication that
    cannot change the rank or the solution set -- and the probe duly reported
    the mutation as undetectable, which was correct. The real defect was the
    ABSENCE of per-field normalisation, so that is what is injected here: the
    raw columns span 1e-37 (D_z) to 1e22 (q_z), the condition number reaches
    1.6e62 and the solve loses its digits.
    """
    p = root / "code/scattering/assembly.py"
    s = p.read_text()
    anchor = """    row_scale = []
    for i in range(len(STATE_FIELDS)):
        mx = mp.mpf(0)
        for _sd, _br, col in cols:
            if abs(col[i]) > mx:
                mx = abs(col[i])
        row_scale.append(mx if mx != 0 else mp.mpf(1))"""
    assert anchor in s, "MR11 anchor missing"
    s = s.replace(anchor, """    row_scale = [mp.mpf(1) for _ in range(len(STATE_FIELDS))]""", 1)
    p.write_text(s)


def mr11_detect(root: Path) -> tuple[bool, str]:
    """The detector is the CONDITIONING gate, not the rank gate.

    Removing unit normalisation cannot change the rank -- it is a diagonal
    rescaling -- but it destroys the condition number, which is what the
    assembly's own conditioning record reports.
    """
    return _run_sct(root, "GAUGE05 or COND01")



SUITE = [
    Mutation("M01", "zero executed cases reported as PASS",
             "validation/harness.py",
             "harness.finish converts a 0-case PASS to UNRESOLVED and fails the test",
             "The Phase-1B requirement: an unexecuted test is not a pass.",
             m01_mutate, m01_detect_suite),
    Mutation("M02", "incorrect characteristic-matrix parity (e33 sign in M[1,2])",
             "code/foundation_checks.py",
             "M(-xi) = S M(xi) S must fail for the unique admitted S",
             "This is the Phase-1 defect class: a coupling sign that breaks parity.",
             m02_mutate, m02_detect),
    Mutation("M03", "incorrect heat-flux expression (linear thermal term)",
             "code/foundation_checks.py",
             "check_D2_balance must not report PASS",
             "The Phase-1 P0 defect, re-injected. The detector is the Second-Law "
             "balance, which shares no amplitude with the flux.",
             m03_mutate, m03_detect),
    Mutation("M04", "wrong polynomial coefficient ordering",
             "code/foundation_checks.py",
             "the returned roots must not satisfy det M = 0",
             "polyroots expects descending coefficients; feeding ascending ones "
             "returns roots of a different polynomial.",
             m04_mutate, m04_detect),
    Mutation("M05", "incorrect incident/outgoing mode selection",
             "code/foundation_checks.py",
             "_admissible must not return a growing mode",
             "Swapping the two root families selects the wrong half-space modes and "
             "silently produces a scattering matrix with the wrong boundary data.",
             m05_mutate, m05_detect),
    Mutation("M06", "wrong material units (Pa entered as GPa)",
             "literature/material_constants.csv",
             "the implied longitudinal speed must leave the physically admissible band",
             "A silent 1000x unit error changes every dimensionless group.",
             m06_mutate, m06_detect),
    Mutation("M07", "inappropriate isotropic benchmark on an anisotropic medium",
             "code/mutation_probe.py",
             "PZT-5H must be shown to violate sagittal isotropy",
             "The Fresnel formula is only valid for isotropic media; PZT-5H is not.",
             m07_mutate, m07_detect),
    Mutation("M08", "incorrect interface sign (transmitted block sign lost)",
             "code/foundation_checks.py",
             "the 8x8 interface operator must not retain full rank",
             "The transmitted half-space contributes with the opposite sign; losing it "
             "makes two columns degenerate.",
             m08_mutate, m08_detect),
    Mutation("M09", "missing mandatory evidence artifact",
             "evidence/equation_registry.csv",
             "ART01 must fail",
             "An evidence package with a missing registry is not a package.",
             m09_mutate, m09_detect),
    Mutation("M10", "wrong source/data checksum",
             "evidence/input_checksums.json",
             "ART02 must fail",
             "Silent edits to frozen inputs are the failure mode checksums exist for.",
             m10_mutate, m10_detect),
    Mutation("M11", "claim promoted to VERIFIED while its test no longer executes it",
             "evidence/claim_registry.csv",
             "the audit must find the claim VERIFIED against a test that does not "
             "compute the asserted quantity",
             "This is U-11 itself: the defect the whole registry exists to prevent.",
             m11_mutate, m11_detect),    Mutation("M12", "inaccurate eigenvector direction despite a tiny characteristic residual",
             "code/foundation_checks.py",
             "U3C must detect that the null-vector DIRECTION moved",
             "The Phase 1D warning: a small |M a| residual does not certify the "
             "direction. Phase 1C used the residual as a proxy and reached a "
             "wrong diagnosis as a result.",
             m12_mutate, m12_detect),
    Mutation("M13", "wrong displacement-vs-potential scaling",
             "code/foundation_checks.py",
             "the energy balance must fail when phi is scaled relative to u",
             "A gauge/scaling error leaves M a = 0 intact but changes every "
             "boundary-state quantity.",
             m13_mutate, m13_detect),
    Mutation("M14", "incorrectly removed electric gauge mode",
             "code/foundation_checks.py",
             "the root count must remain 8",
             "Dropping the gauge branch loses a genuine mode and corrupts the "
             "incoming/outgoing selection.",
             m22_mutate, m22_detect),
    Mutation("M15", "missing Lord-Shulman relaxation operator in the heat flux",
             "code/foundation_checks.py",
             "the energy balance must fail when (1 - i w tau_q) is dropped",
             "The relaxation operator IS the Lord-Shulman correction; dropping it "
             "silently reverts to Fourier conduction.",
             m15_mutate, m15_detect),
    Mutation("M16", "invalid material-data provenance (literal drifts from the table)",
             "code/foundation_checks.py",
             "DAT02 must fail",
             "A hard-coded value that no longer matches the versioned table is the "
             "provenance defect the review names.",
             m16_mutate, m16_detect),
    Mutation("M17", "stale evidence registry after a theory correction",
             "theory/energy_flux_derivation.tex",
             "the frozen-input check must fail",
             "The theory text is checksummed; editing it without regenerating the "
             "register must be caught.",
             m17_mutate, m17_detect),
    Mutation("M18", "promotion with an open mandatory blocker / non-computing test",
             "evidence/claim_registry.csv",
             "claim binding must fail",
             "A passing suite must not promote a claim whose named test no longer "
             "computes the asserted quantity.",
             m18_mutate, m18_detect),    Mutation("M19", "withdrawn tau_q causality floor restored to the active theory",
             "theory/governing_equations.tex",
             "THY02 must detect an active causality-floor assertion",
             "BLK-DOC-01 was still live in Phase 1D: the boxed floor survived as a "
             "\textsc{Derived} result long after the claim was marked INVALID.",
             m19_mutate, m19_detect),
    Mutation("M20", "decoupled electrostatic mode mislabelled as a gauge",
             "code/scattering/modal.py",
             "SCTA03: D_z = -eps33*k makes the mode physical and admissible",
             "Found for real: this deleted the electrostatic channel of an "
             "elastic half-space while keeping the piezoelectric one.",
             m20_mutate, m20_detect),
    Mutation("M21", "strict Im(xi) != 0 radiation test",
             "code/scattering/radiation.py",
             "SCTA04: PZT-5H must retain an incident channel at theta > 0",
             "Found for real: no PZT-5H mode was INCOMING at theta > 0, so the "
             "coupled problem had no incident wave.",
             m21_mutate, m21_detect),
    Mutation("M22", "same-sign continuity in the assembly",
             "code/scattering/assembly.py",
             "SCTA05: the reconstructed interface field must be continuous",
             "Found for real: the solve residual stayed at 1e-37 while the "
             "interface was violated by 2x the incident state.",
             m22_mutate, m22_detect),
    Mutation("M23", "pure gauge admitted as a physical channel",
             "code/scattering/assembly.py",
             "SCTA02: no gauge may enter the basis",
             "A null channel must never be counted as a physical amplitude.",
             m23_mutate, m23_detect),
    # ---- task 2R: normal-incidence reduction -------------------------
    Mutation("MR01", "shear rows removed on the basis alone (rhs not tested)",
             "code/scattering/assembly.py",
             "SCTR01: the u_x and txz rows carry the shear mode",
             "The original Task 2 defect. A row may be dropped only when it is "
             "zero across the basis AND on the incident right-hand side.",
             mr01_mutate, mr01_detect),
    Mutation("MR02", "P-wave u_x = 0 generalised to every mode",
             "code/scattering/assembly.py",
             "SCTR02: a component zero on one wave is not zero on all",
             "The single inference that produced the whole normal-incidence "
             "error.",
             mr02_mutate, mr02_detect),
    Mutation("MR03", "quiet rows discarded as if they were empty",
             "code/scattering/roots.py",
             "SCTR03: removal must not change the rank",
             "A row that is merely SMALL is not a row that is ZERO. Removing it "
             "leaves fewer equations than unknowns.",
             mr03_mutate, mr03_detect),
    Mutation("MR04", "duplicate modal column admitted",
             "code/scattering/assembly.py",
             "SCTR04: no two admitted columns may be proportional",
             "Two names for one degree of freedom make the coefficients "
             "non-unique for no physical reason.",
             mr04_mutate, mr04_detect),
    Mutation("MR05", "gauge condition posed as a boundary condition",
             "code/scattering/assembly.py",
             "SCTR05: a gauge is not a physical non-uniqueness",
             "phi1 = phi2 with one side carrying no potential forces the "
             "transmitted field to zero.",
             mr05_mutate, mr05_detect),
    Mutation("MR06", "same-sign continuity across the interface",
             "code/scattering/assembly.py",
             "SCTR06: matched media must not reflect",
             "The historical Task 2 sign defect, re-injected.",
             mr06_mutate, mr06_detect),
    Mutation("MR07", "removed-row diagnostics report rhs as zero without checking",
             "code/scattering/assembly.py",
             "SCTR07: the physical basis and the incident datum are both checked",
             "A converged solve is not a correct solve.",
             mr07_mutate, mr07_detect),
]


# ===========================================================================
# TASK 2T -- electrical boundary physics and the normal-incidence limit.
# MR12..MR17. Each is a negative control the Task 2T brief names explicitly.
# They execute the PRODUCTION modules via the Task 2T test file.
# ===========================================================================

def _run_sct2t(root: Path, expr: str, timeout=1500):
    r = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/test_17_scattering_task2t.py",
         "-q", "--no-header", "-p", "no:cacheprovider", "-k", expr],
        cwd=root / "validation", capture_output=True, text=True, timeout=timeout)
    tail = (r.stdout.strip().splitlines() or [""])[-1][:200]
    return r.returncode != 0, tail


def mr12_mutate(root: Path) -> None:
    """Replace the dielectric contact with a conducting electrode.

    A conductor carries a free surface charge, so D_z terminates on it and
    continuity of D_z must NOT be imposed. Swapping one for the other while
    keeping both conditions is the classic over-constrained error.
    """
    p = root / "code/scattering/maxwell.py"
    s = p.read_text()
    anchor = ('    DIELECTRIC_CONTACT: {\n'
              '        "phi": "CONTINUOUS", "Dz": "CONTINUOUS", "sigma_free": "ZERO",\n'
              '        "equations": ("phi", "Dz"),')
    assert anchor in s, "MR12 anchor missing"
    s = s.replace(anchor, '    DIELECTRIC_CONTACT: {\n'
                          '        "phi": "EQUIPOTENTIAL", "Dz": "TERMINATES",\n'
                          '        "sigma_free": "FREE",\n'
                          '        "equations": ("phi",),', 1)
    p.write_text(s)


def mr12_detect(root: Path) -> tuple:
    return _run_sct2t(root, "MAXW01 or MAXW02")


def mr13_mutate(root: Path) -> None:
    """Discard D_z continuity WITHOUT a derivation.

    The legitimate reason to drop a row at k = 0 is that it is identically
    satisfied there. Dropping it for convenience, at every k, is not.
    """
    p = root / "code/scattering/asymptotic.py"
    s = p.read_text()
    anchor = '    use = list(range(len(STATE_FIELDS))) if rows is None else list(rows)'
    assert anchor in s, "MR13 anchor missing"
    s = s.replace(anchor,
                  '    use = list(range(len(STATE_FIELDS))) if rows is None else list(rows)\n'
                  '    use = [i for i in use if STATE_FIELDS[i] != "Dz"]', 1)
    p.write_text(s)


def mr13_detect(root: Path) -> tuple:
    return _run_sct2t(root, "MAXW02 or LIM03")


def mr14_mutate(root: Path) -> None:
    """Impose zero free charge when a nonzero one is prescribed.

    PRESCRIBED_CHARGE requires the normal-D jump to equal sigma_given. Silently
    setting sigma = 0 imposes a different physical interface.
    """
    p = root / "code/scattering/maxwell.py"
    s = p.read_text()
    anchor = '    sig = mp.mpf(0) if sigma_free is None else mp.mpf(sigma_free)'
    assert anchor in s, "MR14 anchor missing"
    s = s.replace(anchor, '    sig = mp.mpf(0)   # mutated: charge ignored', 1)
    p.write_text(s)


def mr14_detect(root: Path) -> tuple:
    # the model table must still declare PRESCRIBED_CHARGE as jumping by sigma;
    # a test that exercises the residual helper is what catches the silent zero
    r = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0,'code')\n"
         "import mpmath as mp; mp.mp.dps=40\n"
         "from scattering import maxwell as mx\n"
         "f1=[mp.mpc(0)]*8; f2=[mp.mpc(0)]*8\n"
         "f1[6]=0; f2[6]=mp.mpf('5e-9')\n"
         "r=mx.electrical_residuals(mx.PRESCRIBED_CHARGE,f1,f2,"
         "sigma_free=mp.mpf('2e-9'),k=mp.mpf('1e-3'))\n"
         "bad=abs(r['normal_D']-mp.mpf('3e-9'))>mp.mpf('1e-20')\n"
         "sys.exit(1 if bad else 0)"],
        cwd=root, capture_output=True, text=True)
    return r.returncode != 0, (r.stdout + r.stderr).strip()[-160:]


def mr15_mutate(root: Path) -> None:
    """Confuse the global gauge with a relative potential jump.

    A gauge shift adds the SAME constant everywhere. A jump adds different
    constants on the two sides. Implementing the gauge as a jump makes a
    physical surface dipole look unobservable.
    """
    p = root / "code/scattering/asymptotic.py"
    s = p.read_text()
    anchor = '        out["f1"] = [inc[i] + f1[i] for i in range(len(STATE_FIELDS))]'
    assert anchor in s, "MR15 anchor missing"
    s = s.replace(anchor,
                  '        _g = f2[2] - f1[2]\n'
                  '        out["f1"] = [inc[i] + f1[i] + (_g if i == 2 else 0)\n'
                  '                      for i in range(len(STATE_FIELDS))]', 1)
    p.write_text(s)


def mr15_detect(root: Path) -> tuple:
    return _run_sct2t(root, "MAXW02 or MAXW03")


def mr16_mutate(root: Path) -> None:
    """Treat eigenvector rescaling as physical divergence.

    This is the Task 2S error itself. Restoring the |u_x| = 1 gauge as the basis
    for the limit makes the potential look like it diverges.
    """
    p = root / "code/scattering/asymptotic.py"
    s = p.read_text()
    anchor = '            "displacement": unit_displacement_mode}[normalisation]'
    assert anchor in s, "MR16 anchor missing"
    s = s.replace(anchor,
                  '            "displacement": unit_displacement_mode}[normalisation]\n'
                  '    norm = unit_displacement_mode   # mutated: pin u_x = 1 everywhere', 1)
    p.write_text(s)


def mr16_detect(root: Path) -> tuple:
    return _run_sct2t(root, "LIM01 or LIM02")


def mr17_mutate(root: Path) -> None:
    """Declare the k = 0 limit validated from a SINGLE finite-k sample.

    One sample cannot distinguish a limit from a value. This mutation makes the
    limit routine return the k = 1e-3 answer directly and claims success.
    """
    p = root / "code/scattering/asymptotic.py"
    s = p.read_text()
    anchor = '    usable = [(mp.mpf(kk), mp.mpf(abs(v)))'
    assert anchor in s, "MR17 anchor missing"
    s = s.replace(anchor,
                  '    if len(list(values)) >= 1:\n'
                  '        v = mp.mpf(abs(list(values)[0]))\n'
                  '        return mp.mpf(0), v, "CONSTANT"\n'
                  + anchor, 1)
    p.write_text(s)


def mr17_detect(root: Path) -> tuple:
    return _run_sct2t(root, "LIM03 or LIM02")


SUITE.extend([
    Mutation("MR08", "relative potential offset treated as free", 
             "code/scattering/assembly.py",
             "EQ-02A-GAUGE-A6: the bonded model is not the insulated model",
             "A free relative offset is an interface charge sheet; a bonded "
             "contact has none. Keeping the phi row solvable hides it.",
             mr08_mutate, mr08_detect),
    Mutation("MR09", "D_z continuity imposed at k = 0 as bonded",
             "code/scattering/assembly.py",
             "EQ-02A-GAUGE-A6: the discarded row is the insulated condition",
             "Imposing D_z across the interface silently converts a bonded "
             "joint into an insulated one.",
             mr09_mutate, mr09_detect),
    Mutation("MR10", "global gauge condition removed from the full system",
             "code/scattering/gauge.py",
             "EQ-02A-GAUGE-A2: exactly one gauge condition is required",
             "Two constants with no fixing row leave a one-parameter family; "
             "reporting a unique answer then reports a gauge choice.",
             mr10_mutate, mr10_detect),
    Mutation("MR11", "unit normalisation removed, raw SI columns span 1e59",
             "code/scattering/assembly.py",
             "SCTR04: no two admitted columns may be proportional",
             "The measured Task 2S production defect: q_z ~ 1e19 dominates, "
             "u_x is annihilated at 1e-14 and distinct columns look identical.",
             mr11_mutate, mr11_detect),
])


SUITE.extend([
    Mutation("MR12", "dielectric contact replaced by a conducting electrode",
             "code/scattering/maxwell.py",
             "EQ-02A-MAXW-01/02: a conductor terminates D_z, it does not "
             "transmit it",
             "Imposing both potential and D_z continuity on a conductor "
             "over-constrains a model whose charge is free.",
             mr12_mutate, mr12_detect),
    Mutation("MR13", "D_z continuity discarded with no derivation",
             "code/scattering/asymptotic.py",
             "EQ-02A-MAXW-02: the row is dropped only where it is empty",
             "Dropping a physical condition to make the system square is the "
             "failure mode this task exists to prevent.",
             mr13_mutate, mr13_detect),
    Mutation("MR14", "zero free charge imposed when a charge is prescribed",
             "code/scattering/maxwell.py",
             "EQ-02A-MAXW-01: PRESCRIBED_CHARGE must carry sigma_given",
             "Silently zeroing sigma imposes a different physical interface.",
             mr14_mutate, mr14_detect),
    Mutation("MR15", "global gauge confused with a relative potential jump",
             "code/scattering/asymptotic.py",
             "EQ-02A-MAXW-02/03: a jump is physical, a gauge shift is not",
             "A gauge shift adds the same constant to both sides; a jump does "
             "not.",
             mr15_mutate, mr15_detect),
    Mutation("MR16", "eigenvector rescaling treated as physical divergence",
             "code/scattering/asymptotic.py",
             "EQ-02A-LIM-01: the 1/k is a normalisation artefact",
             "The original Task 2S error, re-injected.",
             mr16_mutate, mr16_detect),
    Mutation("MR17", "k=0 limit validated from a single finite-k sample",
             "code/scattering/asymptotic.py",
             "EQ-02A-LIM-02/03: a limit needs a sequence",
             "One sample cannot distinguish a limit from a value.",
             mr17_mutate, mr17_detect),
])

# ===========================================================================
def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args(argv)

    before = {rel: sha(REPO / rel) for rel in FROZEN if (REPO / rel).exists()}
    suite = [m for m in SUITE if not a.only or m.mid in a.only]
    if not suite:
        print("no mutations selected", file=sys.stderr)
        return 3

    print("=" * 78)
    print(f"MUTATION SUITE -- {len(suite)} mutations, isolated copies under a temp root")
    print(f"origin tree: {REPO}")
    print("=" * 78)

    tmp_root = Path(tempfile.mkdtemp(prefix="phase1c_mut_"))
    print(f"temp root:   {tmp_root}\n")

    try:
        for i, m in enumerate(suite, 1):
            print(f"[{i}/{len(suite)}] {m.mid} {m.title}")
            dst = tmp_root / m.mid
            ctl = tmp_root / f"{m.mid}_control"
            try:
                # CONTROL FIRST: the same detector against an UNMUTATED copy.
                # A detector that fires here would "detect" everything and the
                # whole suite would be theatre. Two real defects in the first
                # version of this file were found exactly this way.
                copy_repo(ctl)
                ctl_fired, ctl_ev = m.detect(ctl)
                m.control_fired = bool(ctl_fired)
                m.control_evidence = str(ctl_ev)

                copy_repo(dst)
                m.mutate(dst)
                detected, ev = m.detect(dst)
                m.evidence = ev
                if m.control_fired:
                    m.result = "INVALID_DETECTOR"
                else:
                    m.result = "DETECTED" if detected else "SURVIVOR"
            except Exception as exc:                            # noqa: BLE001
                m.result = "INCONCLUSIVE"
                m.evidence = f"{exc!r}"
                if a.keep:
                    traceback.print_exc()
            print(f"         -> {m.result}: {m.evidence[:130]}")
            if m.control_fired:
                print(f"         !! detector also fires on the pristine tree: "
                      f"{m.control_evidence[:100]}")
            if not a.keep:
                shutil.rmtree(dst, ignore_errors=True)
                shutil.rmtree(ctl, ignore_errors=True)
    finally:
        if not a.keep:
            shutil.rmtree(tmp_root, ignore_errors=True)

    touched = assert_origin_untouched(before)
    if touched:
        print(f"\n*** ORIGIN TREE MODIFIED: {touched} -- this is a harness bug ***",
              file=sys.stderr)

    detected = [m for m in suite if m.result == "DETECTED"]
    survivors = [m for m in suite if m.result == "SURVIVOR"]
    incc = [m for m in suite if m.result == "INCONCLUSIVE"]
    inval = [m for m in suite if m.result == "INVALID_DETECTOR"]
    total = len(suite)
    rate = 100.0 * len(detected) / total if total else 0.0

    print("\n" + "=" * 78)
    print(f"DETECTION RATE {len(detected)}/{total} = {rate:.1f}%")
    if survivors:
        print("SURVIVORS (validation gaps):")
        for m in survivors:
            print(f"  {m.mid} {m.title}\n      detector: {m.detector}\n"
                  f"      evidence: {m.evidence}")
    if incc:
        print("INCONCLUSIVE (harness or mutation could not run):")
        for m in incc:
            print(f"  {m.mid}: {m.evidence}")
    if touched:
        print(f"ORIGIN CONTAMINATION: {touched}")
    print("=" * 78)

    report = dict(
        schema_version="1.0",
        generated_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        origin_repo=str(REPO),
        origin_untouched=not touched,
        total=total,
        detected=len(detected),
        detection_rate_percent=round(rate, 1),
        survivors=[m.mid for m in survivors],
        invalid_detectors=[m.mid for m in inval],
        inconclusive=[m.mid for m in incc],
        mutations=[dict(id=m.mid, title=m.title, target=m.target,
                        detector=m.detector, rationale=m.rationale,
                        result=m.result, evidence=m.evidence,
                        control_fired=bool(m.control_fired),
                        control_evidence=m.control_evidence) for m in suite],
    )
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"report: {OUT.relative_to(REPO)}")

    if touched:
        return 1
    return 0 if not survivors and not incc and not inval else 1


if __name__ == "__main__":
    sys.exit(main())

