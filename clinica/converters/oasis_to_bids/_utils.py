import os
from pathlib import Path
from typing import Iterable, Optional, Union

import nibabel as nb
import numpy as np
import pandas as pd

from clinica.converters.study_models import (
    OASISBIDSSubjectID,
    StudyName,
    bids_id_factory,
)

__all__ = [
    "create_participants_df",
    "get_subjects_list",
    "create_sessions_df",
    "write_sessions_tsv",
    "write_scans_tsv",
    "get_first_image",
    "mapping_diagnosis",
    "get_image_with_good_orientation",
]


def create_participants_df(
    clinical_specifications_folder: Path,
    clinical_data_dir: Path,
    bids_ids: list[str],
    delete_non_bids_info: bool = True,
) -> pd.DataFrame:
    """Create the file participants.tsv.

    Parameters
    ----------
    clinical_specifications_folder : Path
        The path to the clinical file.

    clinical_data_dir : Path
        The path to the directory where the clinical data are stored.

    bids_ids : list of str
        The list of bids ids.

    delete_non_bids_info : bool, optional
        If True delete all the rows of the subjects that are not available in the BIDS dataset.
        Default=True.

    Returns
    -------
    pd.DataFrame :
        A pandas dataframe that contains the participants data.
    """
    from clinica.converters.study_models import bids_id_factory

    study_name = StudyName.OASIS
    location = f"{study_name.value} location"

    participants_specs = pd.read_csv(
        clinical_specifications_folder / "participant.tsv", sep="\t"
    )[["BIDS CLINICA", study_name, location]].dropna()
    excel_location = participants_specs[location].unique()[0]

    participant_df = pd.DataFrame()
    for _, row in participants_specs.iterrows():
        if row[location] != excel_location:
            excel_location = row[location]
        file_to_read = pd.read_excel(clinical_data_dir / excel_location, sheet_name=0)
        participant_df = pd.concat(
            [
                participant_df,
                file_to_read[[row[study_name]]].rename(
                    {row[study_name]: row["BIDS CLINICA"]}, axis=1
                ),
            ],
            axis=1,
        )

    # OASIS provides several MRI for the same session
    participant_df = participant_df[
        ~participant_df.alternative_id_1.str.endswith("_MR2")
    ]
    # Adding participant_id column with BIDS ids
    participant_df["participant_id"] = participant_df["alternative_id_1"].apply(
        lambda x: bids_id_factory(study_name).from_original_study_id(x)
    )

    # Delete all the rows of the subjects that are not available in the BIDS dataset
    if delete_non_bids_info:
        participant_df = participant_df.set_index("participant_id", drop=False).loc[
            bids_ids
        ]
        # todo : what happens if subject not in clinical data ?
    participant_df.reset_index(inplace=True, drop=True)
    return participant_df.fillna("n/a")


def _get_subjects_list_from_file(subjects_list_path: Path) -> list[OASISBIDSSubjectID]:
    """Gets the list of subjects folders names from the subjects list file.

    Parameters
    ----------
    subjects_list_path : Path
        The path to the subjects list file.

    Returns
    -------
    list[OASISBIDSSubjectID] :
        List of subjects folders names.
    """
    return [
        OASISBIDSSubjectID(OASISBIDSSubjectID.from_original_study_id(subject))
        for subject in subjects_list_path.read_text().splitlines()
    ]


def _get_subjects_list_from_data(source_dir: Path) -> list[OASISBIDSSubjectID]:
    """Gets the list of subjects folders names from the input dataset folder.

    Parameters
    ----------
    source_dir : Path
        The path to the input dataset folder.

    Returns
    -------
    list[OASISBIDSSubjectID] :
        List of subjects folders names.
    """
    return [
        OASISBIDSSubjectID(OASISBIDSSubjectID.from_original_study_id(folder.name))
        for folder in source_dir.iterdir()
        if not folder.name.startswith(".")
    ]


def _filter_oasis_subjects(
    source_dir: Path, subjects_list: list[OASISBIDSSubjectID]
) -> list[Path]:
    """Filters and builds the paths to the subjects folders.

    Parameters
    ----------
    source_dir : Path
        The path to the input dataset folder.

    subjects_list : list[OASISBIDSSubjectID]
        The list of subjects folders names.

    Returns
    -------
    list[Path] :
        List of paths to the subjects folders.
    """
    return list(
        filter(
            lambda path: path.is_dir(),
            [source_dir / subj.to_original_study_id() for subj in subjects_list],
        )
    )


def get_subjects_list(
    source_dir: Path,
    subjs_list_path: Optional[Path] = None,
) -> list[Path]:
    """Gets the list of paths to the subjects folders from the raw dataset.

    Parameters
    ----------
    source_dir : Path
        The path to the input dataset folder.

    subjs_list_path : Optional[Path]
        The path to the subjects list file.

    Returns
    -------
    list[Path] :
        List of paths to the subjects folders.
    """

    if subjs_list_path:
        return _filter_oasis_subjects(
            source_dir, _get_subjects_list_from_file(subjs_list_path)
        )

    return _filter_oasis_subjects(source_dir, _get_subjects_list_from_data(source_dir))


def _convert_cdr_to_diagnosis(cdr: Union[int, str]) -> str:
    if cdr == 0:
        return "CN"
    elif isinstance(cdr, int) and cdr > 0:
        return "AD"
    else:
        return "n/a"


def create_sessions_df(
    clinical_data_dir: Path,
    clinical_specifications_folder: Path,
    bids_ids: Iterable[str],
) -> pd.DataFrame:
    """Extract the information regarding sessions M000 and store them in a dataframe.

    Parameters
    ----------
    clinical_data_dir : Path
        The path to the input folder.

    clinical_specifications_folder : Path
        The path to the clinical file folder.

    bids_ids : list of str
        The list of bids ids which are in the BIDS directory.

    Returns
    -------
    pd.Dataframe :
        Session df.
    """

    study = StudyName.OASIS.value
    location = f"{study} location"
    spec = pd.read_csv(clinical_specifications_folder / "sessions.tsv", sep="\t")[
        [study, location, "BIDS CLINICA"]
    ].dropna()

    sessions_df = pd.DataFrame()
    if len(spec[location].unique()) == 1:
        loc = spec[location].unique()[0]
    else:
        raise ValueError(
            f"OASIS1 metadata is supposed to be contained in only 1 file, {len(spec[location].unique())} were detected : {spec[location].unique()}"
        )

    file = pd.read_excel(clinical_data_dir / loc)
    file["BIDS ID"] = file.ID.apply(
        lambda x: bids_id_factory(StudyName.OASIS).from_original_study_id(x)
    )
    file.set_index("BIDS ID", drop=True, inplace=True)

    for _, row in spec[spec[location] == loc].iterrows():
        sessions_df[row["BIDS CLINICA"]] = file[row[[study]]]

    missing_subjects = set(bids_ids) - set(sessions_df.index)
    for ms in missing_subjects:
        sessions_df.loc[ms] = ["n/a" for _ in sessions_df.columns]

    sessions_df = sessions_df.loc[bids_ids]

    sessions_df["diagnosis"] = sessions_df["diagnosis"].apply(
        lambda x: _convert_cdr_to_diagnosis(x)
    )

    sessions_df.insert(loc=0, column="session_id", value="ses-M000")

    return sessions_df


def write_sessions_tsv(bids_dir: Path, sessions_df: pd.DataFrame) -> None:
    """Writes the content of the function `clinica.iotools.bids_utils.create_sessions_df`
    in several TSV files following the BIDS specification.

    Parameters
    ----------
    bids_dir : Path
        The path to the BIDS directory.

    sessions_df : DataFrame
        Contains sessions metadata.

        .. note::
            This is the output of the function
            `clinica.iotools.bids_utils.create_sessions_df`.

    See also
    --------
    create_sessions_df
    """
    for subject, data in sessions_df.iterrows():
        session_path = bids_dir / subject
        data.to_frame().T.to_csv(
            session_path / f"{subject}_sessions.tsv",
            sep="\t",
            encoding="utf8",
            index=False,
        )


def write_scans_tsv(bids_dir: Path) -> None:
    """
    Write the scans.tsv file at the root of baseline sessions (ses-M000).

    Parameters
    ----------
    bids_dir : Path to the BIDS output
    """
    for subject_path in bids_dir.rglob("sub-*"):
        if subject_path.is_dir():
            to_write = pd.DataFrame(
                {
                    "filename": [
                        f"{path.parent.name}/{path.name}"
                        for path in subject_path.rglob("*ses-M000*.nii.gz")
                    ]
                }
            )

            to_write.to_csv(
                subject_path / "ses-M000" / f"{subject_path.name}_ses-M000_scans.tsv",
                sep="\t",
                index=False,
            )


def get_first_image(input_folder: Path) -> Path:
    """Get the first .img file in the folder given as parameter.
    Throw an exception if no file is found.

    Parameters
    ----------
    input_folder : Path
        The path to the input folder.

    Returns
    -------
    Path :
        The path to the first image found in the provided folder.
    """
    from clinica.utils.stream import log_and_raise

    try:
        img_file_path = next(input_folder.glob("*.img"))
        return img_file_path

    except StopIteration:
        log_and_raise(
            f"No file ending in .img found in {input_folder}.",
            FileNotFoundError,
        )


def mapping_diagnosis(diagnosis_bl: float) -> str:
    """Mapping diagnosis label to corresponding number.

    Parameters
    ----------
    diagnosis_bl : float
        Diagnosis number.

    Returns
    -------
    str :
        Diagnosis label if matching a number, Not Available (n/a) otherwise.
    """

    if diagnosis_bl == 0.0:
        return "CN"

    elif diagnosis_bl in (0.5, 1.0, 1.5, 2.0):
        return "AD"

    return "n/a"


def get_image_with_good_orientation(image_path: Path) -> nb.Nifti1Image:
    """Convert an image from its path to Nifti with the correct orientation.

    Parameters
    ----------
    image_path : Path
        The path to the input image.

    Returns
    -------
    nb.Nifti1Image :
        The converted Nifti image with the correct orientation.
    """
    # First, convert to Nifti so that we can extract the s_form with NiBabel
    # (NiBabel creates an 'Spm2AnalyzeImage' object that does not contain 'get_sform' method
    img_with_wrong_orientation_analyze = nb.load(image_path)

    # OASIS-1 images have the same header but sform is incorrect
    # To solve this issue, we use header from images converted with FreeSurfer
    # to generate a 'clean hard-coded' header
    # affine:
    # [[   0.    0.   -1.   80.]
    #  [   1.    0.    0. -128.]
    #  [   0.    1.    0. -128.]
    #  [   0.    0.    0.    1.]]
    # fmt: off
    affine = np.array(
        [
            0, 0, -1, 80,
            1, 0, 0, -128,
            0, 1, 0, -128,
            0, 0, 0, 1
        ]
    ).reshape(4, 4)
    # fmt: on
    s_form = affine.astype(np.int16)

    hdr = nb.Nifti1Header()
    hdr.set_data_shape((256, 256, 160))
    hdr.set_data_dtype(np.int16)
    hdr["bitpix"] = 16
    hdr.set_sform(s_form, code="scanner")
    hdr.set_qform(s_form, code="scanner")
    hdr["extents"] = 16384
    hdr["xyzt_units"] = 10

    img_with_good_orientation_nifti = nb.Nifti1Image(
        np.round(img_with_wrong_orientation_analyze.get_fdata(dtype="float32")).astype(
            np.int16
        ),
        s_form,
        header=hdr,
    )
    # Header correction to obtain dim0 = 3
    return nb.funcs.four_to_three(img_with_good_orientation_nifti)[0]
