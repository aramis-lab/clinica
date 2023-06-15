"""This module contains utilities for PET data handling."""
import os
import typing as ty
from enum import Enum

import pandas as pd


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


class ReconstructionMethod(str, Enum):
    """BIDS label for PET reconstruction methods.

    Follows the convention proposed in the PET section of the BIDS specification.

    See: https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/09-positron-emission-tomography.html#pet-recording-data

    For ADNI specific reconstruction methods, see:

    https://adni.loni.usc.edu/methods/pet-analysis-method/pet-analysis/
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


def read_psf_information(
    pvc_psf_tsv: os.PathLike,
    subject_ids: ty.List[str],
    session_ids: ty.List[str],
    pet_tracer: str,
) -> ty.List[ty.List[int]]:
    """Read PSF information from TSV file.

    Parameters
    ----------
    pvc_psf_tsv : PathLike
        Path to the TSV file containing the following columns:
        'participant_id', 'session_id', 'acq_label', 'psf_x',
        'psf_y', and 'psf_z'

    subject_ids : List[str]
        List of participant IDs.
        ex: ['sub-CLNC01', 'sub-CLNC01']

        .. warning::
            Must have the same length as `session_ids`.

    session_ids : List[str]
        List of session IDs.
        ex: ['ses-M000', 'ses-M018']

        .. warning::
            Must have the same length as `subject_ids`.

    pet_tracer : str
        Tracer we want to select in the 'acq_label' column.
        Other tracers will not be read in this function

    Returns
    -------
    psf : List[List[int]]
        The PSF information following [subject_ids, session_ids] order.

    Examples
    --------
    Example of pvc_psf_tsv:

    participant_id    session_id     acq_label     psf_x    psf_y    psf_z
    sub-CLNC01        ses-M000        FDG           8        9        10
    sub-CLNC01        ses-M018        FDG           8        9        10
    sub-CLNC01        ses-M000        AV45          7        6        5
    sub-CLNC02        ses-M000        FDG           8        9        10
    sub-CLNC03        ses-M000        FDG           8        9        10
    """
    valid_columns = {
        "participant_id",
        "session_id",
        "acq_label",
        "psf_x",
        "psf_y",
        "psf_z",
    }
    psf_df = pd.read_csv(pvc_psf_tsv, sep="\t")
    diff = valid_columns.symmetric_difference(set(psf_df.columns))
    if len(diff) > 0:
        raise IOError(
            f"The file {pvc_psf_tsv} must contain the following columns (separated by tabulations):\n"
            f"participant_id, session_id, acq_label, psf_x, psf_y, psf_z\n"
            f"{str(list(psf_df.columns))}\n"
            f"Pay attention to the spaces (there should be none)."
        )
    psf = []
    for sub, ses in zip(subject_ids, session_ids):
        result = psf_df.query(
            f"participant_id == '{sub}' and session_id == '{ses}' and acq_label == '{pet_tracer}'"
        )
        if len(result) == 0:
            raise RuntimeError(
                f"Subject {sub} with session {ses} and tracer {pet_tracer} "
                f"that you want to proceed was not found in the TSV file containing "
                f"PSF specifications ({pvc_psf_tsv})."
            )
        if len(result) > 1:
            raise RuntimeError(
                f"Subject {sub} with session {ses} and tracer {pet_tracer} "
                f"that you want to proceed was found multiple times "
                f"in the TSV file containing PSF specifications ({pvc_psf_tsv})."
            )
        psf.append(result[["psf_x", "psf_y", "psf_z"]].values.tolist()[0])

    return psf


LIST_SUVR_REFERENCE_REGIONS = ["pons", "cerebellumPons", "pons2", "cerebellumPons2"]


def get_suvr_mask(suvr_reference_region: str) -> os.PathLike:
    """Returns the path to the SUVR mask from SUVR reference region label.

    Parameters
    ----------
    suvr_reference_region : str
        The label of the SUVR reference region.
        Supported labels are: 'pons', 'cerebellumPons', 'pons2', and 'cerebellumPons2'

    Returns
    -------
    PathLike :
        The path to the SUVR mask.
    """
    from pathlib import Path

    current_dir = Path(os.path.realpath(__file__))
    masks_dir = current_dir.parent.parent / "resources" / "masks"

    suvr_reference_region_labels_to_filenames = {
        "pons": "region-pons_eroded-6mm_mask.nii.gz",
        "cerebellumPons": "region-cerebellumPons_eroded-6mm_mask.nii.gz",
        "pons2": "region-pons_remove-extrabrain_eroded-2it_mask.nii.gz",
        "cerebellumPons2": "region-cerebellumPons_remove-extrabrain_eroded-3it_mask.nii.gz",
    }
    if suvr_reference_region not in suvr_reference_region_labels_to_filenames:
        raise ValueError(
            f"SUVR reference region label {suvr_reference_region} is not supported. "
            f"Supported values are : {list(suvr_reference_region_labels_to_filenames.keys())}."
        )
    return masks_dir / suvr_reference_region_labels_to_filenames[suvr_reference_region]
