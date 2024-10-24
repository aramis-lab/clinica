from pathlib import Path
from typing import List, Union

import pandas as pd

from clinica.utils.pet import SUVRReferenceRegion, Tracer

__all__ = [
    "get_suvr_mask",
    "read_psf_information",
]


def get_suvr_mask(region: Union[str, SUVRReferenceRegion]) -> Path:
    """Returns the path to the SUVR mask from SUVR reference region label.

    Parameters
    ----------
    region : str or SUVRReferenceRegion
        The label of the SUVR reference region.
        Supported labels are: 'pons', 'cerebellumPons', 'pons2', and 'cerebellumPons2'

    Returns
    -------
    Path :
        The path to the SUVR mask.
    """
    masks_dir = Path(__file__).resolve().parents[2] / "resources" / "masks"

    return masks_dir / _get_suvr_reference_region_labels_filename(
        SUVRReferenceRegion(region)
    )


def get_mni_mask() -> Path:
    """ """
    return Path(__file__).resolve().parents[2] / "resources" / "masks" / ""


def _get_suvr_reference_region_labels_filename(region: SUVRReferenceRegion) -> str:
    if region == SUVRReferenceRegion.PONS:
        return "region-pons_eroded-6mm_mask.nii.gz"
    if region == SUVRReferenceRegion.CEREBELLUM_PONS:
        return "region-cerebellumPons_eroded-6mm_mask.nii.gz"
    if region == SUVRReferenceRegion.PONS2:
        return "region-pons_remove-extrabrain_eroded-2it_mask.nii.gz"
    if region == SUVRReferenceRegion.CEREBELLUM_PONS2:
        return "region-cerebellumPons_remove-extrabrain_eroded-3it_mask.nii.gz"


def read_psf_information(
    pvc_psf_tsv_path: Path,
    subject_ids: List[str],
    session_ids: List[str],
    pet_tracer: Union[str, Tracer],
) -> List[List[int]]:
    """Read PSF information from TSV file.

    Parameters
    ----------
    pvc_psf_tsv_path : Path
        The path to the TSV file containing the following columns:
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

    pet_tracer : str or Tracer
        The tracer we want to select in the 'acq_label' column.
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
    pet_tracer = Tracer(pet_tracer)
    df = _read_psf_dataframe(pvc_psf_tsv_path)
    psf = []
    for subject, session in zip(subject_ids, session_ids):
        result = df.query(
            f"participant_id == '{subject}' and session_id == '{session}' and acq_label == '{pet_tracer.value}'"
        )
        if len(result) == 0:
            raise RuntimeError(
                f"Subject {subject} with session {session} and tracer {pet_tracer.value} "
                f"that you want to proceed was not found in the TSV file containing "
                f"PSF specifications ({pvc_psf_tsv_path})."
            )
        if len(result) > 1:
            raise RuntimeError(
                f"Subject {subject} with session {session} and tracer {pet_tracer.value} "
                f"that you want to proceed was found multiple times "
                f"in the TSV file containing PSF specifications ({pvc_psf_tsv_path})."
            )
        psf.append(result[["psf_x", "psf_y", "psf_z"]].values.tolist()[0])

    return psf


def _read_psf_dataframe(filename: Path) -> pd.DataFrame:
    valid_columns = {
        "participant_id",
        "session_id",
        "acq_label",
        "psf_x",
        "psf_y",
        "psf_z",
    }
    df = pd.read_csv(filename, sep="\t")
    column_mismatch = valid_columns.symmetric_difference(set(df.columns))
    if column_mismatch:
        raise IOError(
            f"The file {filename} must contain the following columns (separated by tabulations):\n"
            f"{valid_columns}\n"
            f"The provided TSV file contains the following columns instead:{list(df.columns)}\n"
            f"Pay attention to the spaces (there should be none)."
        )
    return df
