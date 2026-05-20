"""Module for converting FDG PET of ADNI."""

from functools import partial
from os import PathLike
from pathlib import Path
from typing import Iterable, List, Set, Tuple

import pandas as pd

__all__ = [
    "convert_fdg_pet",
]

from clinica.converters.adni_to_bids.modality_converters._pet_utils import (
    ADNITracer,
)

from .._modality import ADNIModalityConverter, ADNIPETPreprocessingStep


def convert_fdg_pet(
    source_dir: Path,
    csv_dir: Path,
    destination_dir: Path,
    conversion_dir: Path,
    subjects: Iterable[str],
    force_new_extraction: bool = False,
    n_procs: int = 1,
    pet_preprocessing_step: int = 4,
):
    """Convert FDG PET images of ADNI into BIDS format.

    Parameters
    ----------
    source_dir : PathLike
        Path to the ADNI directory.

    csv_dir : PathLike
        Path to the clinical data directory.

    destination_dir : PathLike
        Path to the destination BIDS directory.

    conversion_dir : PathLike
        Path to the TSV files including the paths to original images.

    pet_preprocessing_step : int, optional
        ADNI processing step. Default = 4

    subjects : List of str, optional
        List of subjects.

    force_new_extraction : bool
        If True, pre-existing images in the BIDS directory
        will be erased and extracted again.

    n_procs : int, optional
        The requested number of processes.
        If specified, it should be between 1 and the number of available CPUs.
        Default=1.
    """
    from clinica.utils.stream import cprint

    from .._utils import paths_to_bids

    cprint(
        "Calculating paths of FDG PET images. "
        f"Output will be stored in {conversion_dir}."
    )

    pet_preprocessing_step = ADNIPETPreprocessingStep.from_step_value(
        pet_preprocessing_step
    )

    images = _compute_fdg_pet_paths(
        source_dir, csv_dir, subjects, conversion_dir, pet_preprocessing_step
    )

    cprint("Paths of FDG PET images found. Exporting images into BIDS ...")
    modality = ADNIModalityConverter.PET_FDG
    paths_to_bids(
        images,
        destination_dir,
        modality,
        force_new_extraction=force_new_extraction,
        n_procs=n_procs,
        pet_preprocessing_step=pet_preprocessing_step,
    )
    cprint(msg=f"{modality.value} conversion done.", lvl="debug")


def _compute_fdg_pet_paths(
    source_dir: Path,
    csv_dir: Path,
    subjects: Iterable[str],
    conversion_dir: Path,
    pet_preprocessing_step: ADNIPETPreprocessingStep = ADNIPETPreprocessingStep.STEP4_8MM,
) -> pd.DataFrame:
    """Compute the paths to the FDG PET images and store them in a TSV file.

    Parameters
    ----------
    source_dir : PathLike
        Path to the ADNI directory.

    csv_dir : PathLike
        Path to the clinical data directory. It must contain the following
        CSV files:
            - PETQC.csv
            - PETC3.csv
            - PET_META_LIST.csv

    subjects : list of str
        List of subjects.

    conversion_dir : PathLike
        Path to the TSV files including the paths to original images.

    pet_preprocessing_step : PreprocessingStep
        ADNI processing step, is an int between 0 and 5.

    Returns
    -------
    images : pd.DataFrame
        DataFrame with all the paths to the PET images that will be converted into BIDS.
    """
    from clinica.utils.pet import Tracer

    from ._image_path_utils import find_image_path

    tracer = Tracer.FDG.value
    pet_fdg_df = _get_pet_fdg_df(csv_dir, subjects, pet_preprocessing_step)
    images = find_image_path(pet_fdg_df, source_dir, tracer)
    images.to_csv(
        conversion_dir / f"{tracer}_pet_paths.tsv",
        sep="\t",
        index=False,
    )

    return images


def _get_pet_fdg_df(
    csv_dir: Path,
    subjects: Iterable[str],
    pet_preprocessing_step: ADNIPETPreprocessingStep,
) -> pd.DataFrame:
    """Build a DataFrame for the PET FDG images for the provided list of subjects."""
    dfs = []
    for subject in subjects:
        dfs.extend(
            _get_images_pet_for_subject(
                subject,
                _get_csv_data(Path(csv_dir)),
                pet_preprocessing_step,
            )
        )
    if len(dfs) == 0:
        return pd.DataFrame(columns=_get_pet_fdg_columns())
    pet_fdg_df = pd.concat(dfs, ignore_index=True)
    return _remove_known_conversion_errors(pet_fdg_df)


def _get_pet_fdg_columns() -> List[str]:
    return [
        "Phase",
        "Subject_ID",
        "VISCODE",
        "Visit",
        "Sequence",
        "Scan_Date",
        "Study_ID",
        "Series_ID",
        "Image_ID",
        "Original",
    ]


def _get_csv_data(csv_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load needed data from .csv files in csv_dir folder."""
    return (
        _get_pet_qc_df(csv_dir),
        _get_qc_adni_3_df(csv_dir),
        _get_meta_list_df(csv_dir),
    )


def _load_df_with_column_check(
    csv_dir: Path, filename: str, required_columns: Set[str]
) -> pd.DataFrame:
    """Load the requested CSV file in a dataframe and check that the requested columns are present."""
    from clinica.converters._utils import load_clinical_csv

    df = load_clinical_csv(csv_dir, filename)
    if not required_columns.issubset(set(df.columns)):
        raise ValueError(
            f"Missing column(s) from {filename} file."
            f"Required columns for this file are {required_columns}."
        )
    return df


_get_pet_qc_df = partial(
    _load_df_with_column_check,
    filename="PETQC",
    required_columns={"PASS", "RID"},
)
_get_qc_adni_3_df = partial(
    _load_df_with_column_check,
    filename="PETC3",
    required_columns={"SCANQLTY", "RID", "SCANDATE"},
)
_get_meta_list_df = partial(
    _load_df_with_column_check,
    filename="PET_META_LIST",
    required_columns={"Subject"},
)


def _get_images_pet_for_subject(
    subject: str,
    csv_data: Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
    pet_preprocessing_step: ADNIPETPreprocessingStep = ADNIPETPreprocessingStep.STEP2,
) -> List[pd.DataFrame]:
    """Filter the PET images' QC dataframes for the given subject."""
    from clinica.converters.adni_to_bids.modality_converters._pet_utils import (
        define_pet_processing_step_with_tracer,
    )

    from ._pet_utils import get_images_pet

    pet_qc_df, pet_qc_adni_3_df, pet_meta_list_df = csv_data
    subject_pet_metadata = pet_meta_list_df[pet_meta_list_df["Subject"] == subject]

    if subject_pet_metadata.empty:
        return []

    subject_pet_qc = _build_pet_qc_all_studies_for_subject(
        subject, pet_qc_df, pet_qc_adni_3_df
    )

    return get_images_pet(
        subject,
        subject_pet_qc,
        subject_pet_metadata,
        _get_pet_fdg_columns(),
        "FDG-PET",
        [
            define_pet_processing_step_with_tracer(
                ADNITracer.FDG, pet_preprocessing_step
            )
        ],
    )


def _build_pet_qc_all_studies_for_subject(
    subject: str,
    qc_df: pd.DataFrame,
    qc_adni_3_df: pd.DataFrame,
) -> pd.DataFrame:
    """Filter QC for FDG PET images and for the given subject, for studies ADNI 1, GO, 2, and 3.
    Merge everything in a single dataframe.
    """
    subject_rid = _convert_subject_to_rid(subject)
    qc_1go2 = qc_df[(qc_df.PASS == 1) & (qc_df.RID == subject_rid)]
    qc_3 = qc_adni_3_df[
        (qc_adni_3_df.SCANQLTY == 1) & (qc_adni_3_df.RID == subject_rid)
    ]
    qc_3.insert(0, "EXAMDATE", qc_3.SCANDATE.to_list())
    return pd.concat([qc_1go2, qc_3], axis=0, ignore_index=True, sort=False)


def _convert_subject_to_rid(subject: str) -> int:
    """Get the QC RID from the subject string identifier.
    Examples
    --------
    >>> _convert_subject_to_rid("123_S_4567")
    4567
    """

    from re import search

    match = search(r"\d{3}_S_(\d{4})", subject)

    if match:
        return int(match.group(1))
    else:
        raise ValueError(
            f"Cannot convert the subject '{subject}' identifier into a RID "
            "for PET QC filtering. The expected format for the subject is XXX_S_XXXX."
        )


def _remove_known_conversion_errors(df: pd.DataFrame) -> pd.DataFrame:
    """Remove known exceptions from images to convert."""
    if df.empty:
        return df
    error_ind = df.index[df.apply(lambda x: _is_visit_a_conversion_error(x), axis=1)]
    df.drop(error_ind, inplace=True)
    return df


def _is_visit_a_conversion_error(row: pd.Series) -> bool:
    """Return whether a given visit row is among the known exceptions or not."""
    conversion_errors = [  # NONAME.nii
        ("031_S_0294", "bl"),
        ("037_S_1421", "m36"),
        ("037_S_1078", "m36"),
        # Empty folders
        ("941_S_1195", "m48"),
        ("005_S_0223", "m12"),
    ]
    return (row.Subject_ID, row.VISCODE) in conversion_errors
