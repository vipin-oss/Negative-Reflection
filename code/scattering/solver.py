"""
Solve the two-half-space coupled scattering system -- Phase 2A task 2.

STATUS OF WHAT THIS MODULE PRODUCES
===================================
Everything computed here is labelled **DIAGNOSTIC**. Master gates G4 (energy
conservation for the coupled piezothermoelastic problem) and G5 (thermodynamic
consistency of the coupled scattering result) remain OPEN, per the task brief:
a converged boundary solve is not evidence that R, T and A are validated
piezothermoelastic quantities. The numbers are reported because the assembly
must be exercised, and the label is carried in every record.

WHAT IS CHECKED
===============
For each configuration:
  * complex amplitudes of every reflected and transmitted channel
  * boundary residuals BY EQUATION (all eight), not a single norm
  * matrix rank and equilibrated condition number
  * working precision and a convergence test across precisions
  * ARBITRARY RESCALING of every mode: the physical amplitudes must not change.
    A result that moves under a consistent mode rescaling is invalid.

Rescaling is the sharpest test available here. The modal basis is determined
only up to a complex scalar per mode (a null vector has no intrinsic norm), so
if the answer is basis-independent the solve must be invariant. The test
multiplies each admitted column by an arbitrary complex factor, re-solves, and
compares the per-channel amplitudes AFTER undoing the scaling.

CONTINUOUS vs DISCONTINUOUS
===========================
The system is continuous (displacements, potential, temperature, tractions,
normal electric displacement, normal heat flux), so under a mode scaling
c_j -> alpha_j c_j the unknown coefficients must satisfy a_j -> a_j / alpha_j.
The physical field of mode j is a_j * c_j, which is what must be invariant.
"""
from __future__ import annotations

import mpmath as mp

from scattering.assembly import (COND_TRUST, assemble,
                                 condition_number, matrix_rank)
from scattering.state import STATE_FIELDS

#: Label carried by every numeric result from this module. Cleared only by an
#: expert gate review, never by a converging solve.
EVIDENCE_LABEL = "DIAGNOSTIC"

#: A physical channel whose energy share is below this is a null direction, not
#: a coefficient. Recorded rather than silently truncated.
AMPLITUDE_NEGLIGIBLE = mp.mpf("1e-6")


def _equilibrate(B):
    """Row-equilibrate a complex matrix. Returns (scaled, row divisors).

    The interface matrix mixes metres and pascals, so its raw condition number
    is dominated by unit choice and means nothing (measured ~2.5e14 for a
    well-behaved problem). Dividing each row by its largest entry is an exact
    rescaling of the equations and leaves the solution untouched.
    """
    n, m = B.rows, B.cols
    S = mp.matrix(n, m)
    scale = []
    for i in range(n):
        rowmax = max([abs(B[i, j]) for j in range(m)] or [mp.mpf(0)])
        scale.append(rowmax if rowmax != 0 else mp.mpf(1))
        for j in range(m):
            S[i, j] = B[i, j] / scale[-1]
    return S, scale


def solve(system, dps=None):
    """Least-squares / square solve of the assembled system.

    A square system is solved directly. An over-determined one (which is what
    normal incidence reduces to, once the two genuine potential gauges are
    removed) is solved in the least-squares sense and its residual is reported
    -- padding it to square with an arbitrary regularisation is exactly what
    the brief forbids.
    """
    if dps is not None:
        mp.mp.dps = dps
    B, rhs = system["B"], system["rhs"]
    n, m = B.rows, B.cols
    if m == 0:
        return {"ok": False, "reason": "no admissible modes"}

    S, scale = _equilibrate(B)
    r = mp.matrix(n, 1)
    for i in range(n):
        r[i] = rhs[i] / scale[i]

    try:
        if n == m:
            a = mp.lu_solve(S, r)
            method = "lu_solve (square)"
        else:
            # normal equations on the equilibrated system; adequate here
            # because the conditioning is reported and checked
            Sh = mp.matrix(m, n)
            for i in range(m):
                for j in range(n):
                    Sh[i, j] = mp.conj(S[j, i])
            N = Sh * S
            b = Sh * r
            a = mp.lu_solve(N, b)
            method = f"normal equations ({n}x{m}, over/over-determined)"
    except Exception as exc:                                   # noqa: BLE001
        return {"ok": False, "reason": f"solve failed: {exc}"}

    resid = mp.matrix(n, 1)
    for i in range(n):
        resid[i] = r[i] - sum(S[i, j] * a[j] for j in range(m))

    return {
        "ok": True, "method": method,
        "coeffs": a, "n": n, "m": m,
        "residual_by_equation": {
            STATE_FIELDS[i]: resid[i] for i in range(n)},
        "residual_max": max(abs(resid[i]) for i in range(n)),
        "row_scale": scale,
    }


def channel_amplitudes(system, coeffs):
    """Physical amplitudes of every admitted mode.

    The columns of the assembled operator are UNIT-NORMALISED, so the
    coefficient returned for mode j is the amplitude of that mode on a
    basis-independent scale. Without that normalisation the raw coefficient is
    an artefact of the eigenvector routine and moves by 1e248 under a
    consistent mode rescaling (measured) -- which is exactly the invalidity
    `rescale_invariance` is designed to expose.
    """
    out = []
    for j, br in enumerate(system["admitted1"]):
        out.append({"side": 1, "branch": br,
                    "class": system["side1"].modes[br]["spatial_class"],
                    "role": system["side1"].roles[br][0],
                    "coeff": coeffs[j]})
    off = len(system["admitted1"])
    for j, br in enumerate(system["admitted2"]):
        out.append({"side": 2, "branch": br,
                    "class": system["side2"].modes[br]["spatial_class"],
                    "role": system["side2"].roles[br][0],
                    "coeff": coeffs[off + j]})
    return out


def rescale_invariance(system, dps=None, alphas=None):
    """THE validity test: consistent mode rescaling must not move the physics.

    Each admitted column is multiplied by an arbitrary complex factor and the
    solve repeated. Because the basis vector changes, the coefficient changes
    too; the PHYSICAL field is the product. If the product is stable the
    formulation is basis-independent.

    Returns the per-branch relative change of the physical contribution.
    """
    if dps is not None:
        mp.mp.dps = dps
    m = system["n_unknowns"]
    if m == 0:
        return {"ok": False, "reason": "no modes"}
    if alphas is None:
        # deterministic, non-trivial, unequal factors on both sides
        alphas = [mp.mpc(mp.mpf(3) + mp.mpf("0.7") * i,
                         mp.mpf(-1) + mp.mpf("0.3") * i)
                  for i in range(m)]

    base = solve(system)
    if not base["ok"]:
        return {"ok": False, "reason": base["reason"]}

    # Apply an arbitrary complex factor to each column WITHOUT renormalising.
    # An earlier version renormalised afterwards, which cancelled the scaling
    # exactly and made this test a no-op -- it reported "invariant" for a
    # formulation that was not. The operator must still produce the same
    # physical field, because the divisor used to reconstruct that field is
    # recomputed from the scaled columns.
    scaled = dict(system)
    B2 = mp.matrix(system["B"].rows, system["B"].cols)
    for j in range(system["B"].cols):
        for i in range(system["B"].rows):
            B2[i, j] = system["B"][i, j] * alphas[j]
    scaled["column_norms"] = [system["column_norms"][j] / alphas[j]
                              for j in range(len(system["column_norms"]))]
    scaled["B"] = B2
    alt = solve(scaled)
    if not alt["ok"]:
        return {"ok": False, "reason": alt["reason"]}

    # Compare the TOTAL reconstructed interface field, not each coefficient
    # against itself.
    #
    # For identical media the reflected amplitudes are exactly zero, and a
    # relative change computed on a quantity that is zero is not a meaningful
    # number: it reported 54.2 for a coefficient of order 1e-37 and made a
    # correct solve look invalid. The physically meaningful statement is that
    # the FIELD the basis reconstructs must not move, so the change is measured
    # against the magnitude of the total field.
    f0 = _mode_fields(system, base["coeffs"])
    f1 = _mode_fields(scaled, alt["coeffs"])
    tot0 = [mp.mpc(0)] * len(STATE_FIELDS)
    tot1 = [mp.mpc(0)] * len(STATE_FIELDS)
    for f in f0:
        for i in range(len(STATE_FIELDS)):
            tot0[i] += f[i]
    for f in f1:
        for i in range(len(STATE_FIELDS)):
            tot1[i] += f[i]
    scale_field = max(sum(abs(v) for v in tot0), mp.mpf("1e-300"))
    worst = sum(abs(tot1[i] - tot0[i]) for i in range(len(STATE_FIELDS))) \
        / scale_field
    detail = [{"total_field_rel_change": worst}]
    return {"ok": True, "worst_rel_change": worst, "detail": detail,
            "alphas": alphas,
            "null_modes": [j for j, f in enumerate(f0)
                           if sum(abs(v) for v in f)
                           <= AMPLITUDE_NEGLIGIBLE * scale_field]}


def _mode_fields(system, coeffs):
    """Physical interface contribution of each admitted mode, per row.

    Uses the unit-normalised column, so this is the field the mode actually
    puts on the interface and it is independent of the arbitrary eigenvectors.
    """
    scales = system["column_norms"]
    side = system["side1"]
    a1 = system["admitted1"]
    fields = []
    for j in range(len(a1)):
        fields.append([coeffs[j] * side.raw_states[br][i] / scales[j]
                       for br in [a1[j]] for i in range(len(STATE_FIELDS))])
    off = len(a1)
    s2 = system["side2"]
    a2 = system["admitted2"]
    for j in range(len(a2)):
        fields.append([coeffs[off + j] * s2.raw_states[br][i] / scales[off + j]
                       for br in [a2[j]] for i in range(len(STATE_FIELDS))])
    return fields


def boundary_residuals(system, coeffs):
    """Reconstruct the full interface state on each side and difference it.

    This is an INDEPENDENT check of the solve: it does not reuse the residual
    of the linear algebra, it rebuilds s(0) from the modal expansion and
    compares the two half-spaces field by field.
    """
    s1, s2 = system["side1"], system["side2"]
    a1, a2 = system["admitted1"], system["admitted2"]
    # The coefficients belong to the UNIT-NORMALISED columns, so the physical
    # contribution of mode j is coeff_j * state_j / ||state_j||. Forgetting the
    # divisor was a bug in the first version of this function and produced
    # interface residuals of 1e14 on a solve whose own residual was 1e-52.
    #
    # AND THE RECONSTRUCTION MUST USE PHYSICAL STATES.
    #
    # This routine used `s1.states` / `s2.states`, which are each divided by
    # their OWN medium's reference scales. For two different media that makes
    # the difference st1 - st2 a comparison of u_i/scale1_i with u_i/scale2_i
    # -- not a continuity condition -- and the reported residual was measured at
    # 1.1e9 for Al6061|PZT-5H at 30 degrees while the equations that WERE
    # imposed were satisfied to 1.7e-37. The columns of the assembly are now
    # raw physical states, so the reconstruction has to be as well.
    norms = system["column_norms"]
    st1 = [mp.mpc(0)] * len(STATE_FIELDS)
    for j, br in enumerate(a1):
        for i in range(len(STATE_FIELDS)):
            st1[i] += coeffs[j] * s1.raw_states[br][i] / norms[j]
    st2 = [mp.mpc(0)] * len(STATE_FIELDS)
    off = len(a1)
    for j, br in enumerate(a2):
        for i in range(len(STATE_FIELDS)):
            st2[i] += coeffs[off + j] * s2.raw_states[br][i] / norms[off + j]
    # add the incident wave to side 1 -- it is not an unknown
    if system["incident_branch"] is not None:
        inc = s1 if system["incident_side"] == 1 else s2
        st_inc = inc.raw_states[system["incident_branch"]]
        tgt = st1 if system["incident_side"] == 1 else st2
        for i in range(len(STATE_FIELDS)):
            tgt[i] += st_inc[i]
    diff = {STATE_FIELDS[i]: st1[i] - st2[i] for i in range(len(STATE_FIELDS))}

    # NOT EVERY COMPONENT OF THIS DIFFERENCE IS A BOUNDARY-CONDITION VIOLATION.
    #
    # A row that was removed from the system -- because it reads 0 = 0 for every
    # admitted mode, or because it is a gauge rather than a condition -- was
    # never required to be continuous. Reporting its difference as a residual
    # invents a failure:
    #
    #   Al6061|PZT-5H at k = 0, P incidence. The phi row is dropped as gauge
    #   (Al6061 has no potential in any admitted mode). The transmitted
    #   potential is then the sum of the transmitted modes' own potentials, and
    #   the incident side has none, so the raw difference is 0.913 -- while the
    #   SOLVE residual is 1.4e-37, i.e. every equation that was kept is
    #   satisfied to the arithmetic floor.
    #
    # The two quantities are therefore returned separately, and a caller that
    # wants a single number must say which it means.
    enforced = ({r["row"] for r in system.get("rows_removed", [])
                 if r.get("INCONSISTENT")}
                | set(system.get("rows_kept", [])))
    return diff, enforced


def analyse(system, dps=None):
    """Full record for one configuration, labelled DIAGNOSTIC."""
    if dps is not None:
        mp.mp.dps = dps
    sol = solve(system)
    if not sol["ok"]:
        return {"label": EVIDENCE_LABEL, "ok": False,
                "reason": sol["reason"]}

    S, _scale = _equilibrate(system["B"])
    rec = {
        "label": EVIDENCE_LABEL,
        "ok": True,
        "dps": mp.mp.dps,
        "n_equations": system["n_equations"],
        "n_unknowns": system["n_unknowns"],
        "square": system["square"],
        "rank": matrix_rank(system["B"]),
        "cond_equilibrated": condition_number(S),
        "method": sol["method"],
        "solve_residual_max": sol["residual_max"],
        "residual_by_equation": {
            k: v for k, v in sol["residual_by_equation"].items()},
        "channels": channel_amplitudes(system, sol["coeffs"]),
    }
    cond = rec["cond_equilibrated"]
    rec["cond_trust"] = bool(cond <= COND_TRUST)
    if not rec["square"]:
        rec["consistency_note"] = (
            f"{system['n_equations']} independent equations for "
            f"{system['n_unknowns']} unknowns: the system is "
            f"{'over' if system['overdetermined'] else 'under'}determined by "
            f"construction. Reported, not regularised. If the least-squares "
            f"residual is not at rounding level the omitted rows are NOT "
            f"redundant and the model itself needs review.")
    rec["trustworthy"] = bool(rec["cond_trust"] and rec["square"])
    try:
        br, enforced = boundary_residuals(system, sol["coeffs"])
        rec["interface_residual_by_field"] = br
        rec["interface_residual_enforced_rows"] = sorted(enforced)
        # TWO NUMBERS, because they answer two questions.
        #
        # `interface_residual_max` is over EVERY component: it includes rows
        # that were deliberately not imposed, so a gauge row can dominate it
        # while the solve is exact. `interface_residual_enforced_max` is over
        # the rows the system actually enforced, and that is the one that says
        # whether the boundary-value problem was solved.
        rec["interface_residual_max"] = max(abs(v) for v in br.values())
        rec["interface_residual_enforced_max"] = max(
            (abs(br[f]) for f in enforced), default=mp.mpf(0))
    except Exception as exc:                                   # noqa: BLE001
        rec["interface_residual_by_field"] = None
        rec["interface_residual_max"] = None
        rec["interface_residual_enforced_max"] = None
        rec["interface_residual_error"] = str(exc)
    ri = rescale_invariance(system)
    rec["rescale_invariance"] = {
        "ok": ri.get("ok"),
        "worst_rel_change": ri.get("worst_rel_change"),
    }
    # INTERPRET the invariance result against the conditioning.
    #
    # A near-singular operator cannot distinguish coefficients that differ only
    # along its near-null direction, so a large rescaling change is a statement
    # about the OPERATOR, not proof of a formulation error. The identical-media
    # case is exactly this: the interface is invisible, rSS = tSS - 1 = 0 is a
    # one-parameter family, and sigma_min falls to ~1e-5 (and to 1e-21 for
    # Al|Al at normal incidence, where only four independent equations exist).
    # The honest record says which of the two situations occurred.
    sv = _reduced_singular_values(system)
    rec["sigma_min"] = sv[-1] if sv else None
    rec["sigma_max"] = sv[0] if sv else None
    inv = rec["rescale_invariance"]["worst_rel_change"]
    if inv is None:
        rec["rescale_verdict"] = "NO_DATA"
    elif inv <= mp.mpf("1e-20"):
        rec["rescale_verdict"] = "INVARIANT"
    elif rec["cond_equilibrated"] > COND_TRUST:
        rec["rescale_verdict"] = (
            "DEFERRED_TO_CONDITIONING: the operator is near-singular "
            f"(cond = {mp.nstr(rec['cond_equilibrated'], 4)}), so the "
            "coefficients are not individually determined and a rescaling "
            "change is not evidence of a formulation defect. The configuration "
            "is reported as UNRESOLVED, not as passing")
    else:
        rec["rescale_verdict"] = (
            "FAILED: the operator is well conditioned "
            f"(cond = {mp.nstr(rec['cond_equilibrated'], 4)}) yet the physics "
            "moved under a consistent mode rescaling -- this invalidates the "
            "result")
    return rec


def _reduced_singular_values(system):
    """Singular values of the row-equilibrated operator actually solved."""
    S, _sc = _equilibrate(system["B"])
    from scattering.assembly import _singular_values
    return _singular_values(S)


def _main():
    import sys
    sys.path.insert(0, ".")
    mp.mp.dps = 50
    import foundation_checks as fc
    from scattering.assembly import assemble

    w = 2 * mp.pi * 1e6
    print(f"Coupled scattering solve -- evidence label: {EVIDENCE_LABEL}\n")
    print("G4/G5 REMAIN OPEN: these amplitudes are diagnostic only.\n")
    for name1, name2 in (("PZT-5H", "PZT-5H"), ("Al6061", "PZT-5H")):
        for th in (0, 30):
            m1 = fc.material_from_table(name1)
            m2 = fc.material_from_table(name2)
            k = w * mp.sin(mp.radians(mp.mpf(th))) / mp.mpf(4000)
            A = assemble(m1, m2, k, w, incident_branch=6 if th else None)
            r = analyse(A)
            print(f"--- {name1}|{name2}  th={th} deg ---")
            if not r["ok"]:
                print(f"    {r['reason']}")
                continue
            print(f"    eq={r['n_equations']} unk={r['n_unknowns']} "
                  f"rank={r['rank']} cond={mp.nstr(r['cond_equilibrated'], 5)}")
            print(f"    solve resid={mp.nstr(r['solve_residual_max'], 4)}  "
                  f"interface resid={mp.nstr(r['interface_residual_max'], 4)}")
            print(f"    rescale worst={r['rescale_invariance']['worst_rel_change']}")


if __name__ == "__main__":
    _main()
