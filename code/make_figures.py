"""make_figures.py -- figures built ONLY from data computed in this project.
Run: python3 code/make_figures.py
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plate_spectral import (Plate, isotropic_tensors, crystal_6mm, bond_rotate,
                            energy_velocity, rl_roots, eigen_full)

os.makedirs("reports/figures", exist_ok=True)
EPS0 = 8.8541878128e-12
E, NU, RHO = 71.7e9, 0.34, 2790.0
mu = E / (2 * (1 + NU)); lam = E * NU / ((1 + NU) * (1 - 2 * NU))
cL, cT = np.sqrt((lam + 2 * mu) / RHO), np.sqrt(mu / RHO)
C_iso, Z0, ED, _ = isotropic_tensors(E, NU, RHO)
C_pz, e_pz, eps_pz = crystal_6mm(126e9, 79.5e9, 84.1e9, 117e9, 23.0e9,
                                17.0, -6.5, 23.3, 1700 * EPS0, 1470 * EPS0)
Rn = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], float)
C_pzn, e_pzn, eps_pzn = bond_rotate(C_pz, e_pz, eps_pz, Rn)

# ---------------------------------------------------------------- FIGURE 1
# The S1/S2b fold and the sign change of the energy flux.
# No branch tracking is used: at every wavenumber ALL modes in the band are
# evaluated.  This is deliberate -- eigenvector-continuity tracking loses the
# branch at the near-degeneracy, which is itself part of the point.
h = 0.75e-3
pl = Plate(h, C_iso, Z0, ED, RHO, N=26)
ks = np.linspace(880.0, 1240.0, 61)
K, F, V = [], [], []
for k in ks:
    f = pl.spectrum(k)
    f = f[(f > 1.80e6) & (f < 1.96e6)]
    for ff in f:
        ve = energy_velocity(pl, k, ff)[0]
        K.append(k); F.append(ff); V.append(ve)
K, F, V = map(np.array, (K, F, V))

fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
ax[0].scatter(K, F / 1e6, c=np.sign(V), cmap="coolwarm", s=12,
              vmin=-1.2, vmax=1.2)
ax[0].axvline(1061.0, color="k", ls=":", lw=1.1)
ax[0].set_xlabel("wavenumber $k$  (1/m)")
ax[0].set_ylabel("frequency $f$  (MHz)")
ax[0].set_title("(a) $S_1$/$S_{2b}$ fold, Duralumin-like 1.5 mm plate\n"
                "(colour: sign of the energy flux)")
ax[0].grid(alpha=.3)
ax[1].scatter(K, V, c=np.sign(V), cmap="coolwarm", s=12, vmin=-1.2, vmax=1.2)
ax[1].axhline(0, color="k", lw=1)
ax[1].axvline(1061.0, color="k", ls=":", lw=1.1)
ax[1].set_xlabel("wavenumber $k$  (1/m)")
ax[1].set_ylabel("$v_e = P/(2\\langle T\\rangle)$  (m/s)")
ax[1].set_title("(b) the flux changes sign exactly at the ZGV\n"
                "$k_{ZGV}$ marked; no branch tracking used")
ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig("reports/figures/fig1_dispersion_ve.pdf", dpi=200)
plt.savefig("reports/figures/fig1_dispersion_ve.png", dpi=130)
plt.close()

# ---------------------------------------------------------------- FIGURE 2
# validation summary
d1 = json.load(open("data/check_01_elastic_validation.json"))
d2 = json.load(open("data/check_02_flux_and_zgv.json"))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
labels = ["cutoffs\nvs exact", "dispersion\nvs exact", "ZGV: spectral\nvs exact theory",
          "ZGV vs published\nDuralumin value"]
vals = [d1["V1_max_rel_err"], d1["V2b_root_mismatch"],
        d1["V4_vs_V3_rel_diff"], abs(d1["V5_rel_dev"])]
ax[0].bar(range(4), vals, color=["#2b6cb0"] * 3 + ["#c05621"])
ax[0].set_yscale("log"); ax[0].set_ylabel("relative difference")
ax[0].set_xticks(range(4)); ax[0].set_xticklabels(labels, fontsize=7.5)
ax[0].set_ylim(2e-7, 4e-2)
ax[0].set_title("(a) solver validation against exact theory")
ax[0].axhline(1e-2, color="r", ls="--", lw=1)
ax[0].text(0.02, 1.25e-2, "2% tolerance", color="r", fontsize=8)
for i, v in enumerate(vals):
    ax[0].text(i, v * 1.5, f"{v:.1e}", ha="center", fontsize=8)
ax[0].text(3, 4e-4, "literature\ncomparison, not\na solver error", ha="center",
           fontsize=6.5, color="#c05621")
ax[0].grid(alpha=.3, axis="y")

tab = d2["C_table"]
th = [r["thickness_mm"] for r in tab]
keys = ["isotropic reference", "PZT elastic, no coupling",
        "PZT coupled, OPEN", "PZT coupled, SHORT"]
for kk, mk in zip(keys, ["o", "s", "^", "v"]):
    ys = [r[kk]["f_MHz"] if r.get(kk) else np.nan for r in tab]
    ax[1].plot(th, ys, mk + "-", label=kk.replace("PZT elastic, no coupling",
                                                  "PZT elastic, coupling OFF"))
ax[1].set_xscale("log"); ax[1].set_yscale("log")
ax[1].set_xlabel("plate thickness (mm)"); ax[1].set_ylabel("$f_{ZGV}$ (MHz)")
ax[1].set_title("(b) coupling shifts $f_{ZGV}$ by 35%; open vs short does not", fontsize=9.5)
ax[1].legend(fontsize=7); ax[1].grid(alpha=.3, which="both")
plt.tight_layout(); plt.savefig("reports/figures/fig2_validation_premise.pdf", dpi=200)
plt.savefig("reports/figures/fig2_validation_premise.png", dpi=130)
plt.close()

# ---------------------------------------------------------------- FIGURE 3
# Which branches are electrically active, and why the S1/S2b one is not.
d2 = json.load(open("data/check_02_flux_and_zgv.json"))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))

# (a) per-branch sensitivity to the electrical surface condition, at fixed k
kk0 = 900.0
po = Plate(h, C_pzn, e_pzn, eps_pzn, 7500.0, N=24, electrical="open")
ps = Plate(h, C_pzn, e_pzn, eps_pzn, 7500.0, N=24, electrical="short")
fo = po.spectrum(kk0); fs_ = ps.spectrum(kk0)
nd = min(len(fo), len(fs_))
rel = np.abs(fs_[:nd] - fo[:nd]) / fo[:nd] * 100
idx = np.arange(nd)
ax[0].bar(idx, np.maximum(rel, 1e-6))
ax[0].set_yscale("log")
ax[0].set_xlabel("branch index at $k = 900$ 1/m  (increasing frequency)")
ax[0].set_ylabel("$|\\Delta f/f|$ , OPEN vs SHORT  (%)")
ax[0].set_title("(a) most branches are strongly electrically active;\n"
                "the $S_1$/$S_{2b}$ branch (arrow) is not")
zf = d2["C_table"][-1]["PZT coupled, OPEN"]["f_MHz"]
j = int(np.argmin(rel))          # the electrically inactive branch
ax[0].annotate("$S_1$/$S_{2b}$ branch", xy=(j, max(rel[j], 1e-6)),
               xytext=(j + 2.5, 2.0e-3),
               arrowprops=dict(arrowstyle="->", lw=1), fontsize=8,
               bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", lw=.6))
ax[0].grid(alpha=.3, axis="y")

# (b) the mode's own structure: normalised displacement and potential
u, phi = eigen_full(po, 808.95, 1.12975e6)
n = po.n
y = np.linspace(-1, 1, n)
u1 = u[:n].real; u2 = u[n:2 * n].real
sc = np.max(np.abs(u2)) + 1e-30
ph = phi.real / (np.max(np.abs(phi.real)) + 1e-30)
ax[1].plot(u1 / sc, y, label="$u_1$ (scaled)", lw=1.6)
ax[1].plot(u2 / sc, y, label="$u_2$ (scaled)", lw=1.6)
ax[1].plot(ph, y, "--", label="$\\varphi$ (normalised)", lw=1.6)
ax[1].axvline(0, color="k", lw=.8)
ax[1].axhline(0, color="k", lw=.8)
ax[1].set_xlabel("normalised field amplitude")
ax[1].set_ylabel("$x_2 / h$")
ax[1].set_title("(b) the $S_1$/$S_{2b}$ ZGV mode:\n"
                "antisymmetric potential, small surface charge")
ax[1].legend(fontsize=8, loc="lower left"); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig("reports/figures/fig3_premise_limits.pdf", dpi=200)
plt.savefig("reports/figures/fig3_premise_limits.png", dpi=130)
plt.close()
print("wrote reports/figures/fig1..fig3")
