"""This module contains utilities for PET data handling."""
from enum import Enum

__all__ = [
    "Tracer",
    "SUVRReferenceRegion",
    "ReconstructionMethod",
    "get_pet_tracer_from_filename",
]


class Tracer(str, Enum):
    """BIDS label for PET tracers.

    Follows the convention proposed in the PET section of the BIDS specification.

    See: https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/09-positron-emission-tomography.html
    """

    PIB = "11CPIB"
    AV1451 = "18FAV1451"
    AV45 = "18FAV45"
    FBB = "18FFBB"
    FDG = "18FFDG"
    FMM = "18FFMM"


class SUVRReferenceRegion(str, Enum):
    PONS = "pons"
    CEREBELLUM_PONS = "cerebellumPons"
    PONS2 = "pons2"
    CEREBELLUM_PONS2 = "cerebellumPons2"


class ReconstructionMethod(str, Enum):
    """BIDS label for PET reconstruction methods.

    Follows the convention proposed in the PET section of the BIDS specification.

    See: https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/09-positron-emission-tomography.html#pet-recording-data

    For ADNI specific reconstruction methods, see:

    https://adni.loni.usc.edu/data-samples/adni-data/neuroimaging/pet/
    """

    # Reconstruction methods defined in the BIDS specifications
    STATIC = "nacstat"
    DYNAMIC = "nacdyn"
    STATIC_ATTENUATION_CORRECTION = "acstat"
    DYNAMIC_ATTENUATION_CORRECTION = "acdyn"

    # ADNI specific reconstruction methods
    CO_REGISTERED_DYNAMIC = "coregdyn"  # Corresponds to ADNI processing steps 1
    CO_REGISTERED_AVERAGED = "coregavg"  # Corresponds to ADNI processing steps 2
    CO_REGISTERED_STANDARDIZED = "coregstd"  # Corresponds to ADNI processing steps 3
    COREGISTERED_ISOTROPIC = "coregiso"  # Corresponds to ADNI processing steps 4


def get_pet_tracer_from_filename(filename: str) -> Tracer:
    """Return the PET tracer from the provided filename.

    Parameters
    ----------
    filename : str
        The filename from which to extract the PET tracer.

    Returns
    -------
    tracer : Tracer
        The PET tracer.

    Raises
    ------
    ValueError
        If no tracer found in the filename.
    """
    import re

    tracer = None
    for entity in ("trc", "acq"):
        m = re.search(rf"({entity}-[a-zA-Z0-9]+)", filename)
        if m:
            tracer = m.group(1)[4:].upper()
    if tracer is None:
        raise ValueError(
            f"Could not extract the PET tracer from the following file name {filename}."
        )

    return Tracer(tracer)
