from enum import Enum

__all__ = ["DTIBasedMeasure"]


class DTIBasedMeasure(str, Enum):
    """Possible DTI measures."""

    FRACTIONAL_ANISOTROPY = "FA"
    MEAN_DIFFUSIVITY = "MD"
    AXIAL_DIFFUSIVITY = "AD"
    RADIAL_DIFFUSIVITY = "RD"
