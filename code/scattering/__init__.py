"""Phase 2A scattering package: modal classification and the null-vector route.

`modal` answers "what kind of mode is this?" with four INDEPENDENT fields rather
than one label, and `nullspace` supplies a null vector whose result does not
depend on the working precision. The second module exists because the first
cannot be trusted on top of a route that intermittently selects the wrong
branch.
"""
from .modal import (  # noqa: F401
    ATTENUATION_EVANESCENT,
    ATTENUATION_WEAK,
    ENERGY_STATUSES,
    GAUGE_MECHANICAL,
    SPATIAL_CLASSES,
    THERMAL_DISPERSION_TOL,
    classify_mode,
    screen_halfspaces,
    verify_convention,
)
from .nullspace import null_vector, null_vector_cramer, null_vector_gram  # noqa: F401

__all__ = [
    "ATTENUATION_EVANESCENT", "ATTENUATION_WEAK", "ENERGY_STATUSES",
    "GAUGE_MECHANICAL", "SPATIAL_CLASSES", "THERMAL_DISPERSION_TOL",
    "classify_mode", "screen_halfspaces", "verify_convention",
    "null_vector", "null_vector_cramer", "null_vector_gram",
]
