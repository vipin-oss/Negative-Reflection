# 03 — Negative Reflection: A Working Taxonomy

**Document status:** Phase 2 deliverable
**Companion data:** `data/references.csv` (60 verified records), `data/excluded_records.csv`

---

## 0. Why a taxonomy is needed at all

The phrase "negative reflection" is used in the literature for **at least four physically distinct
mechanisms**, and it is routinely confused with four other phenomena that are not reflection at all.
A proposal that says "we will find negative reflection in material X" is unfalsifiable unless the
mechanism is named. This document fixes the vocabulary, then states the mathematical test that
distinguishes each case.

The single most important sentence in this document:

> **A negative reflection coefficient is neither necessary nor sufficient for negative reflection.**

- A coefficient can be driven negative by a 180° phase shift at a hard/soft contrast with no
  negative-angle behaviour whatsoever.
- Negative reflection can occur with a reflection coefficient of modulus 1 and any phase, because
  what makes it "negative" is the **direction of the reflected wave vector relative to the incident
  one**, not the sign of a number.

---

## 1. The ten phenomena, distinguished

Notation. Interface normal **n**; incident wave vector **k**ᵢ with tangential part **k**∥;
reflected wave vector **k**ᵣ; phase velocity **v**φ = (ω/|k|²)**k**; group velocity
**v**g = ∇_**k** ω; time-averaged Poynting/power flux **P**.

### 1.1 Ordinary specular reflection

**Definition.** **k**ᵣ,∥ = **k**ᵢ,∥ and the normal component reverses sign; the reflected beam lies on
the *opposite* side of the normal from the incident beam; **P**ᵣ·**n** = −**P**ᵢ·**n**; and the
tangential energy flow reverses sign together with the wave vector.

**Test.** sign[(**P**ᵣ·**t**)(**k**ᵣ·**t**)] = sign[(**P**ᵢ·**t**)(**k**ᵢ·**t**)] — phase and energy
travel *together* on both sides.

**Comment.** This is the baseline. Everything else below is defined by which of these statements
fails.

### 1.2 Negative reflection (NR) — the primary definition used in this project

**Definition.** The reflected beam lies on the **same side of the normal** as the incident beam:

$$\theta_r < 0 \quad \text{when} \quad \theta_i > 0,$$

i.e. **k**ᵣ,∥ and **k**ᵢ,∥ have the *same sign*, while the reflected energy flow has the *opposite*
tangential sign to the incident energy flow:

$$(\mathbf{P}_r\cdot\mathbf{t})(\mathbf{P}_i\cdot\mathbf{t}) < 0 .$$

For an isotropic medium this forces the reflected mode to be **backward** (phase velocity
anti-parallel to energy flux) on the reflected side. It is therefore *not* a property of the
interface alone — it requires (i) an interface, (ii) a backward branch on one side, and
(iii) tangential-wavenumber matching between a forward incident mode and that backward branch.

**Test (three conditions, all necessary):**
1. **Tangential matching:** **k**ᵣ,∥ = **k**ᵢ,∥ (Snell, exact).
2. **Backward character of the reflected mode:** v_g,∥ · k_∥ < 0 for the reflected mode at that (ω, k).
3. **Energy reversal:** (**P**ᵣ·**t**)(**P**ᵢ·**t**) < 0 with |P| > 0.

Condition 3 is *not* implied by condition 2 when the medium is coupled, anisotropic, or lossy, because
the flux is then not simply proportional to the group velocity. This is the technical heart of
`05_Physical_Feasibility_Foundation.md`.

**Established demonstrations (elastic guided waves):** Germano et al. 2012
(`10.1103/PhysRevB.85.012102`); Gérardin et al. 2016 (`10.1121/1.4959024`); Veres et al. 2016
(`10.1103/PhysRevB.93.174304`); Gérardin et al. 2019 (`10.1038/s41598-019-38480-3`).

### 1.3 Negative refraction

**Definition.** At an interface between two media, the *transmitted* beam lies on the same side of
the normal as the incident beam: **k**ₜ,∥ and **k**ᵢ,∥ have the same sign, with the transmitted energy
flow continuing forward.

**Established demonstrations:** Bramhavar et al. 2011 (`10.1103/PhysRevB.83.014106`);
Legrand et al. 2018 (`10.1103/PhysRevB.98.214114`); Philippe et al. 2015 (`10.1038/srep11112`);
Legrand et al. 2021 (`10.1038/s41598-021-03146-6`).

**Why the distinction matters.** Refraction and reflection are *different boundary conditions*
(transmission versus a free edge) and are realised by different experiments. A result about
negative refraction does not transfer to negative reflection, although both rest on the same
forward↔backward conversion.

### 1.4 Negative group velocity

**Definition.** ∂²ω/∂k² changes sign such that v_g,∥ < 0 at k > 0 for a propagating mode:
**backward-wave propagation**. It is a property of a *mode of a waveguide*, not of any interface.

**Why it is not NR.** A waveguide full of backward modes with no interface exhibits no reflection at
all. Conversely, NR requires a specific interface at which conversion is efficient.

**Important corollary, verified numerically in this project** (`code/check_02_flux_and_zgv.py`):
the sign of v_g equals the sign of the flux **only** for a lossless, non-coupled, non-dissipative
mode; and even then the identity is v_e = P/(2⟨T⟩) = dω/dk, i.e. one must compare with dω/dk and
not with df/dk. A 2π bookkeeping slip in this comparison produced an apparent factor-of-2π
disagreement during this work before being corrected. Dispersion curves alone are not evidence.

### 1.5 Backward-wave propagation

Same content as 1.4, stated as a phenomenon rather than a derivative. The canonical elastic
example is the S₂ᵦ branch of a plate above its ZGV frequency, coexisting with the forward S₁
branch. **Sufficient for NR only in combination with an interface and matching.**

### 1.6 Phase conjugation

**Definition.** A mirror that reverses both the propagation direction *and* the phase of an incident
field: it maps a field ψ to ψ*. For a divergent incident wave, the conjugated wave converges back
to the source.

**Relation to NR.** NR with a reflection coefficient of modulus 1 and phase π acts as a *passive*
phase-conjugating mirror for the matched mode pair. Gérardin et al. 2016 make this equivalence
explicit; Gérardin et al. 2019 exploit it to obtain back-focusing in a chaotic cavity.

**Caution.** The equivalence holds only near the ZGV point, where the forward and backward wavenumbers
coincide. Away from it, the mismatch k₁ ≠ k_{2b} gives a *tunable* focal distance, not a conjugated
field. Phase conjugation is therefore a *consequence*, not a definition.

### 1.7 Reflectionless scattering

**Definition.** R → 0 at an interface, typically through impedance matching, a graded transition,
or an exceptional point of a non-Hermitian system.

**Why it must not be confused with NR.** Reflectionless is |R| → 0; NR is about the *direction* of a
strongly reflected field. A perfect NR mirror has |R| ≈ 1.

### 1.8 Negative-valued reflection amplitude

**Definition.** Re(R) < 0 at some frequency/angle. This happens constantly in ordinary acoustics
(a pressure-release surface has R = −1).

**Verdict.** **Not** evidence of NR. This is the most common error in the literature and the one a
referee will test first. A ±1 amplitude sign carries no directional information; only the
**tangential wave-vector relation and the flux direction** do.

### 1.9 Mode-converted reflection

**Definition.** The incident mode leaves the interface in a *different* mode. This is the norm, not the
exception: a single Lamb mode at a free edge scatters into a large family of propagating,
inhomogeneous and evanescent modes.

**Relation to NR.** NR of guided elastic waves is *always* mode-converted reflection (S₁ → S₂ᵦ or the
reverse). But mode conversion alone is ordinary anomalous reflection: if k_r ≠ k_i the reflected
angle simply differs from the incident one and the beam stays on the opposite side of the normal.

### 1.10 Temporal reflection

**Definition.** Reflection from a *time* interface — the medium's properties are switched in time at
fixed position. Momentum is conserved, frequency is not; the reflected wave has a shifted frequency.

**Why it must be separated.** It is sometimes called "time reflection" and appears in recent elastic
metamaterial work (e.g. `10.1038/s41467-025-64530-8`). It has *no* wavenumber matching and *no*
negative-angle meaning. Conflating it with spatial NR would be a category error.

---

## 2. The decisive physical picture: why NR needs a fold, not a sign

The mechanism established by the elastic guided-wave literature is the following, and it is worth
stating precisely because every viable proposal in this project rests on it.

1. A plate supports an infinite family of modes ω_n(k). At the zero-group-velocity point of a
   branch, ∂ω/∂k = 0 at k ≠ 0.
2. Two branches of the *same symmetry* with close cutoff frequencies repel each other. Above the
   repulsion, a **forward** branch (S₁) and a **backward** branch (S₂ᵦ) coexist over a frequency
   band. At the ZGV frequency their wavenumbers coincide: k₁ = k_{2b}.
3. Because k₁ = k_{2b} there, the two modes have nearly identical stress–displacement fields and
   differ essentially only in the sign of their energy flux.
4. At a free edge, the traction-free condition can then be satisfied by a combination of the two
   *alone*, giving a reflection coefficient of modulus 1 and phase ≈ π for that pair: ρ ≈ −1.
5. Because the tangential wavenumbers are equal and the Poynting vectors are opposite, the reflected
   beam retraces the incident path — retroreflection, back-focusing, passive phase conjugation.

**The two conditions that a new material must satisfy** are therefore:

- **(C1) Existence of a forward/backward pair with coincident wavenumber** — a ZGV fold of the
  appropriate symmetry. This is a property of the *dispersion relation*, so it can be checked before
  any scattering calculation.
- **(C2) Near-degeneracy of the stress–displacement fields of that pair**, so that the boundary
  condition can be met by the pair alone. This is what makes ρ ≈ −1 rather than merely ≠ 0.

Condition (C1) is easy to check numerically and is checked in this project for piezoelectric plates
(`reports/figures/fig1_dispersion_ve.png`, `fig2_validation_premise.png`).
Condition (C2) is *not* checked here — it requires the edge calculation, which is deliberately out
of scope of the feasibility audit. It is the first thing the proposed research must do.

---

## 3. The four mechanisms that can produce NR, and their literature footprint

| Mechanism | How NR arises | Representative verified studies | Available in a plate? |
|---|---|---|---|
| **M1. Forward–backward mode conversion at a boundary** | Incident forward mode converts to a backward mode of the same wavenumber at a free edge | Germano 2012; Gérardin 2016; Veres 2016; Gérardin 2019 | **Yes** — the adopted mechanism |
| **M2. Forward–backward conversion at an interface between two waveguides** | Crossing of dispersion curves across a thickness step or material change | Bramhavar 2011; Legrand 2018; Philippe 2015; Legrand 2021 | Yes, but that is *refraction*, not reflection |
| **M3. Diffraction-order / phase-gradient metasurface** | A phase gradient adds a tangential wave-vector kick; the −1 order appears on the incident side | Song 2017 (`10.1038/s41598-017-14387-9`); Liu 2018 (`10.1038/lsa.2018.8`); Han 2025 (`10.1063/5.0283804`) | Requires a patterned surface, not a homogeneous plate |
| **M4. Anisotropy / hyperbolic dispersion** | The equifrequency contour is non-convex, so the group velocity points into the "wrong" quadrant | Álvarez-Pérez 2022 (`10.1126/sciadv.abp8486`); spin-wave work `10.1021/acsami.5c21700` | Only with strong anisotropy, not in a 6mm plate along a symmetry axis |

**Consequence for problem design.** Only M1 is compatible with the "homogeneous (or graded) coupled
plate with a straight free edge" configuration that matches Dr. Gupta's analytical strengths.
M3 requires fabrication; M4 requires a specific anisotropic cut whose parameters would have to be
justified. Any proposal must therefore be explicit that it is an **M1** study, and must be judged on
whether the coupled material *changes M1's efficiency*, not merely whether a backward branch exists.

---

## 4. Classification of the verified corpus

60 records in `data/references.csv`. Classified by the scheme required in the brief:

| Class (brief's scheme) | Count | Notes |
|---|---|---|
| **Direct negative-reflection study** | **8** | 4 M1 (free-edge forward/backward conversion) + 4 M2 (interface conversion → negative *refraction*). All DOI-verified. |
| **Related forward/backward or ZGV study** | **5** | ZGV in isotropic, anisotropic and graded plates; backward modes under prestress; backward waves in piezoelectric plates. |
| **Related scattering / mode-conversion study** | **7** | Oblique-incidence edge scattering, anomalous reflection in metamaterials, non-Hermitian media, temporal reflection, nonreciprocal/Willis layers. |
| **Relevant mathematical or material extension** | **6** | MEE/nonlocal plates, MEE imperfect interfaces, gradient electro-mechanics/flexoelectric Lamb waves, flexomagnetic reflection, size-dependent piezoelectric review. |
| **Dr. Gupta's own verified work** | **28** | Profile material for Phase 1. |
| **NR in other wave systems** (context only) | **6** | Metasurface, metagrating, anisotropy, spin-wave, magnetically tuned, plus reviews. Included **only** to show the breadth of the term; these are **not** evidence about elastic plates. |
| **TOTAL** | **60** | Tabulated in `reports/02_Verified_Literature.xlsx`, sheet "Verified Literature". |

**Explicit warning, as required by the brief:** the six "other wave systems" papers are **not**
treated as demonstrations of negative reflection of elastic guided waves. Their mechanisms (M3, M4)
are different, and the fact that "negative reflection" appears in their titles is the exact
terminological hazard this taxonomy exists to control.

---

## 5. Honest statement of what the corpus does and does not cover

**Covered, with verified primary sources:**
- M1 mechanism in homogeneous isotropic plates, up to full experimental demonstration, including
  chaotic and disordered geometries.
- M2 mechanism at thickness steps and material interfaces.
- Exact definitions and existence conditions for ZGV points in isotropic and anisotropic plates.
- Backward acoustic waves in *piezoelectric* plates, with ZGV points, established as an
  *experimental* fact (YX-LiNbO₃, `10.3390/s23020648`).

**Not covered by any verified source found:**
- A modal scattering matrix at a free edge of a **coupled** (piezoelectric / piezomagnetic /
  magneto-electro-elastic) plate, i.e. an M1 study in a coupled medium.
- The effect of a **strain-gradient or nonlocal internal length** on the edge reflection
  coefficients of a plate.
- The effect of **imperfect (spring-type) edge or interface conditions** on forward↔backward
  conversion efficiency.

These three absences are the raw material for Phases 4–5. They are *absences of evidence*, which is
weaker than *evidence of absence*; they are graded accordingly in `06_Candidate_Gap_Matrix.md` and
`04_Research_Gap_Evidence.md`.

---

## 6. Search strategy and how to reproduce it

The literature search was run as follows; the queries are recorded so the committee can repeat them.

**Anchor.** Gérardin et al. 2019, `10.1038/s41598-019-38480-3`, read at abstract + results + methods
level (open access), plus its arXiv version `1710.06150`.

**Backward references followed:** Pagneux 2006 (edge resonance); Prada, Clorennec & Royer 2008 (ZGV
definition); Germano 2012 (first NR observation); Bramhavar 2011 and Legrand 2018 (conversion at an
interface); Meitzler-type backward-wave work cited therein.

**Forward citations followed:** Legrand 2021 (devices); Veres 2016 (broad-angle); the anisotropic
ZGV work; the acoustoelastic backward-wave work; the recent spin-wave and THz NR papers that cite
the elastic corpus.

**Query families actually executed** (each with variant phrasings):
1. `negative reflection` + {elastic guided, Lamb, plate, free edge, backward}
2. `zero group velocity` + {Lamb, plate, piezoelectric, anisotropic, orthotropic}
3. `backward wave` + {piezoelectric plate, composite laminate, Lamb, S2b}
4. `forward backward mode conversion` + {edge, thickness step, interface}
5. `negative reflection` + {piezoelectric, piezomagnetic, magneto-electro-elastic, functionally graded}
6. `negative reflection` + {metasurface, gradient, phononic crystal, water wave, spin wave, polariton}
7. `nonlocal strain gradient` + {Lamb wave, piezoelectric plate, flexoelectric}
8. `imperfect interface` + {magneto-electro-elastic, membrane-type, spring model}
9. `exceptional point` + {elastic, Lamb, PT-symmetric, non-Hermitian}
10. `Willis coupling`, `time-modulated` / `temporal reflection` + {elastic, plate}

**Yield and limits.** Query families 1–4 and 7 returned dense, highly relevant results. Queries 5
and 8 returned extensive *material* literature but **no combination with negative reflection** —
which is the gap, but only weakly evidenced: absence in a search index is not absence in the
literature. Queries 6 and 9 returned large literatures on *other* mechanisms (M3, M4, non-Hermitian
scattering) that must be kept terminologically separate.

**Search was limited to what this sandbox can reach:** OpenAlex (full API), publisher landing pages,
arXiv, and a general web index. Scopus and Web of Science were **not** accessible, so no claim here
rests on them. Crossref's REST API was blocked; OpenAlex DOI resolution was used and cross-checked
against publisher pages for the core references.

---

## 7. Summary judgement

The negative-reflection literature for **elastic guided waves** is mature, physically well
understood, and rests on a single clear mechanism (M1: forward–backward conversion at a boundary
where the two branches have equal wavenumber). Its mature status is good news for feasibility and bad
news for novelty: "do M1 in material X" is not, by itself, a research contribution.

What is genuinely open is the *coupling* question — whether, and through which term, a
multi-field constitutive law **changes the efficiency or the conditions of M1**, and specifically
whether the higher-order boundary conditions that gradient theories require can be satisfied by the
mode pair at all. That is a mechanism question, it is answerable analytically, and it is the one the
proposals in Phase 5 are built around.
