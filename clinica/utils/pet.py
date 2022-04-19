"""This module contains utilities for PET data handling."""
from enum import Enum


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


def read_psf_information(pvc_psf_tsv, subject_ids, session_ids, pet_tracer):
    """Read PSF information from TSV file.

    Args:
        pvc_psf_tsv: TSV file containing participant_id, session_id, acq_label,
            psf_x, psf_y & psf_z columns
        subject_ids: list of participant IDs (e.g. ['sub-CLNC01', 'sub-CLNC01'])
        session_ids: list of session IDs (e.g. ['ses-M00', 'ses-M18'])
        pet_tracer: Tracer we want to select in acq_label column. Other tracers
            will not be read in this function

    Example of pvc_psf_tsv:

    participant_id    session_id     acq_label     psf_x    psf_y    psf_z
    sub-CLNC01        ses-M00        FDG           8        9        10
    sub-CLNC01        ses-M18        FDG           8        9        10
    sub-CLNC01        ses-M00        AV45          7        6        5
    sub-CLNC02        ses-M00        FDG           8        9        10
    sub-CLNC03        ses-M00        FDG           8        9        10

    Returns:
        PSF information following [subject_ids, session_ids] order
    """
    import os

    from pandas.io.parsers import read_csv

    if not os.path.isfile(pvc_psf_tsv):
        raise FileNotFoundError(f"Could not find the psf_tsv file {pvc_psf_tsv}")
    try:
        psf_df = read_csv(pvc_psf_tsv, sep="\t")
    except (IOError, UnicodeDecodeError):
        raise RuntimeError("An error while reading {pvc_psf_tsv} happened")

    if any(
        elem
        not in ["participant_id", "session_id", "acq_label", "psf_x", "psf_y", "psf_z"]
        for elem in list(psf_df.columns)
    ):
        raise IOError(
            f"The file {pvc_psf_tsv} must contain the following columns (separated by tabulations):\n"
            f"participant_id, session_id, acq_label, psf_x, psf_y, psf_z\n"
            f"{str(list(psf_df.columns))}\n"
            f"Pay attention to the spaces (there should be none)."
        )

    subjects_psf = list(psf_df.participant_id)
    sessions_psf = list(psf_df.session_id)
    pet_tracer_psf = list(psf_df.acq_label)
    idx_reordered = []
    for i, sub in enumerate(subject_ids):
        current_ses = session_ids[i]
        idx_sub = [
            j
            for j in range(len(subjects_psf))
            if (sub == subjects_psf[j])
            and (current_ses == sessions_psf[j])
            and (pet_tracer == pet_tracer_psf[j])
        ]
        if len(idx_sub) == 0:
            raise RuntimeError(
                f"Subject {sub} with session {current_ses} and tracer {pet_tracer} "
                f"that you want to proceed was not found in the TSV file containing "
                f"PSF specifications ({pvc_psf_tsv})."
            )
        if len(idx_sub) > 1:
            raise RuntimeError(
                f"Subject {sub} with session {current_ses} and tracer {pet_tracer} "
                f"that you want to proceed was found multiple times "
                f"in the TSV file containing PSF specifications ({pvc_psf_tsv})."
            )
        idx_reordered.append(idx_sub[0])

    psf_x = list(psf_df.psf_x)
    psf_y = list(psf_df.psf_y)
    psf_z = list(psf_df.psf_z)
    iterables_psf = [[psf_x[i], psf_y[i], psf_z[i]] for i in idx_reordered]
    return iterables_psf


LIST_SUVR_REFERENCE_REGIONS = ["pons", "cerebellumPons", "pons2", "cerebellumPons2"]


def get_suvr_mask(suvr_reference_region):
    """Get path of the SUVR mask from SUVR reference region label.

    Args:
        suvr_reference_region: Label of the SUVR reference region

    Returns:
        Path of the SUVR mask
    """
    import os

    suvr_reference_region_to_suvr = {
        "pons": os.path.join(
            os.path.split(os.path.realpath(__file__))[0],
            "..",
            "resources",
            "masks",
            "region-pons_eroded-6mm_mask.nii.gz",
        ),
        "cerebellumPons": os.path.join(
            os.path.split(os.path.realpath(__file__))[0],
            "..",
            "resources",
            "masks",
            "region-cerebellumPons_eroded-6mm_mask.nii.gz",
        ),
        "pons2": os.path.join(
            os.path.split(os.path.realpath(__file__))[0],
            "..",
            "resources",
            "masks",
            "region-pons_remove-extrabrain_eroded-2it_mask.nii.gz",
        ),
        "cerebellumPons2": os.path.join(
            os.path.split(os.path.realpath(__file__))[0],
            "..",
            "resources",
            "masks",
            "region-cerebellumPons_remove-extrabrain_eroded-3it_mask.nii.gz",
        ),
    }
    return suvr_reference_region_to_suvr[suvr_reference_region]
