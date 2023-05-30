# coding: utf-8

"""Module for converting FDG PET of ADNI."""

from enum import Enum
from functools import partial
from os import PathLike
from typing import Optional, Union


class ADNIPreprocessingStep(Enum):
    """ADNI preprocessing steps."""

    STEP0 = "ADNI Brain PET: Raw FDG"
    STEP1 = "Co-registered Dynamic"
    STEP2 = "Co-registered, Averaged"
    STEP3 = "Coreg, Avg, Standardized Image and Voxel Size"
    STEP4 = "Coreg, Avg, Std Img and Vox Siz, Uniform Resolution"
    STEP5 = "Coreg, Avg, Std Img and Vox Siz, Uniform 6mm Res"

    @classmethod
    def from_step_value(cls, step_value: Union[int, str]):
        """Accept step specification in raw integer (0, 1, ..., 5) and
        string forms ("1", "3", "step2", "step_4"....).
        """
        error_msg = (
            f"Step value {step_value} is not a valid ADNI preprocessing step value."
            f"Valid values are {list(ADNIPreprocessingStep)}."
        )
        if isinstance(step_value, str) and step_value.startswith("step"):
            step_value = step_value.lstrip("step").lstrip("_")
        try:
            step_value = int(step_value)
        except Exception:
            raise ValueError(error_msg)
        if 0 <= step_value <= 5:
            return cls[f"STEP{step_value}"]
        raise ValueError(error_msg)


def _convert_adni_fdg_pet(
    source_dir: PathLike,
    csv_dir: PathLike,
    destination_dir: PathLike,
    conversion_dir: PathLike,
    preprocessing_step: ADNIPreprocessingStep,
    subjects: Optional[list] = None,
    mod_to_update: bool = False,
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

    preprocessing_step : ADNIPreprocessingStep
        ADNI processing step.

    subjects : List, optional
        List of subjects.

    mod_to_update : bool
        If True, pre-existing images in the BIDS directory
        will be erased and extracted again.
    """
    from pathlib import Path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import paths_to_bids
    from clinica.utils.stream import cprint

    if subjects is None:
        adni_merge = pd.read_csv(
            Path(csv_dir) / "ADNIMERGE.csv", sep=",", low_memory=False
        )
        subjects = list(adni_merge.PTID.unique())
    cprint(
        "Calculating paths of FDG PET images. "
        f"Output will be stored in {conversion_dir}."
    )
    images = _compute_fdg_pet_paths(
        source_dir, csv_dir, subjects, conversion_dir, preprocessing_step
    )
    cprint("Paths of FDG PET images found. Exporting images into BIDS ...")
    modality = _get_modality_from_adni_preprocessing_step(preprocessing_step)
    paths_to_bids(images, destination_dir, modality, mod_to_update=mod_to_update)
    cprint(msg="FDG PET conversion done.", lvl="debug")


def _get_modality_from_adni_preprocessing_step(step: ADNIPreprocessingStep) -> str:
    if step == ADNIPreprocessingStep.STEP2:
        return "fdg"
    if step == ADNIPreprocessingStep.STEP4:
        return "fdg_uniform"
    raise ValueError(
        f"The ADNI preprocessing step {step} is not (yet) supported by the converter."
        f"The converter only supports {ADNIPreprocessingStep.STEP2} and "
        f"{ADNIPreprocessingStep.STEP4} for now."
    )


convert_adni_fdg_pet = partial(
    _convert_adni_fdg_pet, preprocessing_step=ADNIPreprocessingStep.STEP2
)
convert_adni_fdg_pet_uniform = partial(
    _convert_adni_fdg_pet, preprocessing_step=ADNIPreprocessingStep.STEP4
)


def _compute_fdg_pet_paths(
    source_dir: PathLike,
    csv_dir: PathLike,
    subjects: list,
    conversion_dir: PathLike,
    preprocessing_step: ADNIPreprocessingStep,
):
    """Compute the paths to the FDG PET images and store them in a TSV file.

    Parameters
    ----------
    source_dir : PathLike
        Path to the ADNI directory.

    csv_dir : PathLike
        Path to the clinical data directory.

    subjects : list
        List of subjects.

    conversion_dir : PathLike
        Path to the TSV files including the paths to original images.

    preprocessing_step : PreprocessingStep
        ADNI processing step, is an int between 0 and 5.

    Returns
    -------
    images: a dataframe with all the paths to the PET images that will be converted into BIDS
    """
    from pathlib import Path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        find_image_path,
        get_images_pet,
    )
    from clinica.utils.pet import Tracer

    pet_fdg_col = [
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
    csv_dir = Path(csv_dir)
    pet_fdg_df = pd.DataFrame(columns=pet_fdg_col)
    pet_fdg_dfs_list = []

    # Loading needed .csv files
    petqc = pd.read_csv(csv_dir / "PETQC.csv", sep=",", low_memory=False)
    petqc3 = pd.read_csv(csv_dir / "PETC3.csv", sep=",", low_memory=False)
    pet_meta_list = pd.read_csv(
        csv_dir / "PET_META_LIST.csv", sep=",", low_memory=False
    )

    for subject in subjects:
        # PET images metadata for subject
        subject_pet_meta = pet_meta_list[pet_meta_list["Subject"] == subject]

        if subject_pet_meta.empty:
            continue

        # QC for FDG PET images for ADNI 1, GO and 2
        pet_qc_1go2_subj = petqc[(petqc.PASS == 1) & (petqc.RID == int(subject[-4:]))]

        # QC for FDG PET images for ADNI 3
        pet_qc3_subj = petqc3[
            (petqc3.SCANQLTY == 1) & (petqc3.RID == int(subject[-4:]))
        ]
        pet_qc3_subj.insert(0, "EXAMDATE", pet_qc3_subj.SCANDATE.to_list())

        # Concatenating visits in both QC files
        pet_qc_subj = pd.concat(
            [pet_qc_1go2_subj, pet_qc3_subj], axis=0, ignore_index=True, sort=False
        )

        subj_dfs_list = get_images_pet(
            subject,
            pet_qc_subj,
            subject_pet_meta,
            pet_fdg_col,
            "FDG-PET",
            [preprocessing_step.value],
        )
        if subj_dfs_list:
            pet_fdg_dfs_list += subj_dfs_list

    if pet_fdg_dfs_list:
        # Concatenating dataframes into one
        pet_fdg_df = pd.concat(pet_fdg_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [  # NONAME.nii
        ("031_S_0294", "bl"),
        ("037_S_1421", "m36"),
        ("037_S_1078", "m36"),
        # Empty folders
        ("941_S_1195", "m48"),
        ("005_S_0223", "m12"),
    ]

    # Removing known exceptions from images to convert
    if not pet_fdg_df.empty:
        error_ind = pet_fdg_df.index[
            pet_fdg_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1
            )
        ]
        pet_fdg_df.drop(error_ind, inplace=True)

    images = find_image_path(pet_fdg_df, source_dir, "FDG", "I", "Image_ID")
    images.to_csv(
        Path(conversion_dir) / f"{Tracer.FDG}_pet_paths.tsv", sep="\t", index=False
    )

    return images
