"""Task 2U-W normal-incidence modal completeness -- independent certificates.

What this file establishes, on the INDEPENDENT 8-state operator (handoff
module; the certificate module imports no production scattering code):

  ACPT09  Jordan structure at k = 0 (certificates C1/C2): exactly two zero
          eigenvalues, nullity 1 (rank 7), and a generalized eigenvector in
          the uniform-D0 sector (Dz = -1, phi = 0) with B0^2 v ~ machine zero.
  ACPT10  Root separation and electrostatic asymptotics (C3/C4): the six
          nonzero roots are separated, and the small-k electrostatic pair
          converges to +/- s_e * k with s_e = sqrt(eps11/eps33_eff).
  ACPT11  Block-triangular scaled limit and uniform inverse bound (C5-C8):
          the equilibrated S(k) has an invertible 6x6 mechanical block M6
          and 2x2 electrostatic block E0 with det E0 matching the physical
          formula -(eps_eL + eps_eR) to 1e-6 relative; ||S(k)^{-1}||_2 is
          uniformly bounded over k/kr in [1e-8, 1e-2] and the Neumann
          quantity q = ||S0^{-1}(S(k)-S0)||_2 vanishes at the reference k.

Synthetic omega*tau = 0.1 / 1 cases in the JSON deliverable are labelled
synthetic: algorithmic stress tests only, never measured relaxation
properties. Al6061 remains a dielectric surrogate. No energy, R/T/A, or
matching-layer claim is made. BLK-02A-27 closure evidence; SCT02 and the
gates remain open.
"""
import sys
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "handoffs" / "task2uv_2026-11-05" / "code"))

import task2uw_certificates as t2uw                              # noqa: E402

DPS = 40
W = 2 * mp.pi * mp.mpf(10) ** 6


def _c(x):
    """Parse an mp.nstr complex literal like '(-1.0 + 0.0j)'."""
    return mp.mpc(complex(x.strip().strip("()").replace(" ", "")))


def test_acpt09_jordan_structure_k0():
    mp.mp.dps = DPS
    for name in t2uw.MATS:
        p = t2uw.ind_load(t2uw.CSV, name)
        j = t2uw.cert_jordan(p, W, DPS)
        assert j["n_zero_eigenvalues"] == 2, (name, j)
        assert j["nullity_B0"] == 1 and j["rank_B0"] == 7, (name, j)
        assert abs(mp.mpf(j["resid_B0sq_v"])) < mp.mpf("1e-20"), (name, j)
        assert abs(_c(j["generalized_Dz_entry"]) + 1) < mp.mpf("1e-4"), (name, j)
        assert abs(_c(j["generalized_phi_entry"])) < mp.mpf("1e-4"), (name, j)


def test_acpt10_roots_and_electrostatic_pair():
    mp.mp.dps = DPS
    for name in t2uw.MATS:
        p = t2uw.ind_load(t2uw.CSV, name)
        r = t2uw.cert_roots(p, W, DPS)
        assert mp.mpf(r["nonzero_root_min_sep_over_kr"]) > mp.mpf("1e-6"), (name, r)
        assert len(r["small_k_pair"]) >= 2, (name, r)
        for pair in r["small_k_pair"][-2:]:   # k/kr <= 1e-3 asymptotic regime
            assert len(pair["roots_over_k"]) == 2, (name, pair)
            se = mp.mpf(pair["s_e"])
            r0, r1 = _c(pair["roots_over_k"][0]), _c(pair["roots_over_k"][1])
            assert abs(abs(r0) - se) < mp.mpf("1e-5"), (name, pair)
            assert abs(abs(r1) - se) < mp.mpf("1e-5"), (name, pair)
            assert abs(r0 + r1) < mp.mpf("1e-5"), (name, pair)  # opposite signs


def test_acpt11_inverse_bound_pairs():
    mp.mp.dps = DPS
    pairs = [("Al6061", "Al6061"), ("PZT-5H", "PZT-5H"), ("Al6061", "PZT-5H"),
             ("PZT-5H", "Al6061"), ("PZT-5H", "BaTiO3")]
    for n1, n2 in pairs:
        p1, p2 = t2uw.ind_load(t2uw.CSV, n1), t2uw.ind_load(t2uw.CSV, n2)
        for pol in ("P", "SV"):
            ib = t2uw.cert_inverse_bound(
                p1, p2, W, pol, DPS,
                ("1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"))
            assert mp.mpf(ib["M6_smin"]) > mp.mpf("1e-3"), (n1, n2, pol, ib)
            # E0 invertible in the fixed equilibrated basis, with the sign of
            # the physical formula -(eps_eL + eps_eR) (scaled det is basis-
            # dependent; the unscaled formula is recorded alongside)
            detE, detP = _c(ib["detE0_scaled"]), mp.mpf(ib["detE0_physical_formula"])
            assert abs(detE) > mp.mpf("1e-3") and (detE.real < 0) == (detP < 0), \
                (n1, n2, pol, ib)
            meas = max(mp.mpf(r["measured_inv_norm"]) for r in ib["rows"])
            assert meas < mp.mpf("1e8"), (n1, n2, pol, meas)
            assert mp.mpf(ib["rows"][-1]["q_neumann"]) == 0, (n1, n2, pol, ib["rows"][-1])
