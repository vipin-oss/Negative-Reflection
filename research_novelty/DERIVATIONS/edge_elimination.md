# DERIVATIONS/edge_elimination.md

Arena edge Lagrangian (time-harmonic, r,Q ∝ e^{i(kx x − ωt)}):

  L = ½ m ṙ² + ½ L Q̇² − ½ Kc (r−u_edge)² − ½ A |∇_edge r|² − (Q−αr)²/(2C)

Electrical EOM:  −Lω² Q + (Q − αr)/C = 0  ⇒  Q = −αr/(CLω² − 1).

Mechanical edge EOM:
  −mω² r + Kc(r−u) − A kx² r + α(Q − αr)/C = 0.

Elimination (sympy, verified) gives

  r/u = −Kc(CLω²−1) / (ACL kx²ω² − A kx² − C Kc L ω² + CLmω⁴ + Kc + Lα²ω² − mω²),

  K_eff(ω,kx)/Kc = −2A kx²(CLω²−1) / (same denominator).

Interpretation (derived result): K_eff is a rational function of ω² plus a
single A kx² gradient term — i.e. a standard local dynamic boundary impedance
(frequency-dependent stiffness) of the Ingard–Myers/dynamic-BC class. No new
governing law; the electrical DOFs collapse to a lumped resonator impedance.
Double-counting check: if the plate itself is piezoelectric, the (Q−αr)²/2C
term duplicates the plate's own e-coupling — a further consistency objection.
