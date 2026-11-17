"""Module for converting DWI of ADNI."""

from pathlib import Path
from typing import Any, Iterable, Optional

import pandas as pd

__all__ = ["convert_dwi"]


def convert_dwi(
    source_dir: Path,
    csv_dir: Path,
    destination_dir: Path,
    conversion_dir: Path,
    subjects: Iterable[str],
    force_new_extraction: bool = False,
    n_procs: int = 1,
):
    """Convert DW images of ADNI into BIDS format.

    Parameters
    ----------
    source_dir : Path
        The path to the ADNI directory.

    csv_dir : Path
        The path to the clinical data directory.

    destination_dir : Path
        The path to the destination BIDS directory.

    conversion_dir : Path
        The path to the TSV files including the paths to original images.

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

    from .._utils import ADNIModalityConverter, paths_to_bids

    cprint(
        f"Calculating paths of DWI images. Output will be stored in {conversion_dir}.",
        lvl="info",
    )
    images = _compute_dwi_paths(source_dir, csv_dir, subjects, conversion_dir)
    cprint("Paths of DWI images found. Exporting images into BIDS ...", lvl="info")
    paths_to_bids(
        images,
        destination_dir,
        ADNIModalityConverter.DWI,
        force_new_extraction=force_new_extraction,
        n_procs=n_procs,
    )
    cprint(msg="DWI conversion done.", lvl="debug")


def _compute_dwi_paths(
    source_dir: Path,
    csv_dir: Path,
    subjects: Iterable[str],
    conversion_dir: Path,
) -> pd.DataFrame:
    """Compute paths to DW images to convert to BIDS.

    Parameters
    ----------
    source_dir : Path
        The path to the ADNI directory.

    csv_dir : Path
        The path to the clinical data directory.

    subjects : list of str
        The subjects list.

    conversion_dir : Path
        The path to the TSV files including the paths to original images.

    Returns
    -------
    images : pd.DataFrame
        A pandas dataframe that contains the path to all the DWI images to convert.
    """
    from .._utils import load_clinical_csv
    from ._image_path_utils import find_image_path
    from ._mri_utils import load_all_images_metadata

    dwi_df = _initialize_dwi_df()
    dwi_dfs_list = []
    # Loading needed .csv files
    mayo_mri_qc = load_clinical_csv(csv_dir, "MAYOADIRL_MRI_IMAGEQC_05_07_15")
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == "DTI"]
    all_images = load_all_images_metadata(csv_dir)

    # Selecting only DTI images that are not Multiband, processed or enhanced images
    all_images = all_images[
        all_images.modality.str.contains("DTI", case=False, na=False)
    ]
    all_images = all_images[
        all_images.image_description.map(
            lambda x: not any(
                subs in x for subs in ("MB", "ADC", "FA", "TRACEW", "Enhanced", "Reg")
            )
        )
    ]
    for subject in subjects:
        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and
        # sort the rows/visits by examination date
        all_images_subj = all_images[all_images.subject_id == subject]
        all_images_subj = all_images_subj.sort_values("image_date")

        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subject[-4:])]

        for _, row in all_images_subj.iterrows():
            timepoint = row["image_visit"]
            visit_all_images = all_images_subj[all_images_subj.image_visit == timepoint]
            axial = _select_dwi_image(
                subject, timepoint, visit_all_images, mayo_mri_qc_subj
            )
            if axial is not None:
                row_to_append = pd.DataFrame(axial, index=["i"])
                dwi_dfs_list.append(row_to_append)
    if dwi_dfs_list:
        dwi_df = pd.concat(dwi_dfs_list, ignore_index=True)

    # Removing known exceptions from images to convert
    if not dwi_df.empty:
        error_ind = dwi_df.index[
            dwi_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in _get_known_conversion_errors()),
                axis=1,
            )
        ]
        dwi_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(dwi_df, source_dir, "DWI")
    images.to_csv(conversion_dir / "dwi_paths.tsv", sep="\t", index=False)

    return images


def _initialize_dwi_df() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "Subject_ID",
            "VISCODE",
            "Sequence",
            "Scan_Date",
            "Study_ID",
            "Series_ID",
            "Image_ID",
            "Field_Strength",
        ]
    )


def _get_known_conversion_errors() -> Iterable[tuple[str, str]]:
    return [
        ("029_S_2395", "m60"),
        ("029_S_0824", "m108"),
        ("029_S_0914", "m108"),
        ("027_S_2219", "m36"),
        ("129_S_2332", "m12"),
        ("029_S_4384", "m48"),
        ("029_S_4385", "m48"),
        ("029_S_4585", "m48"),
        ("016_S_4591", "m24"),
        ("094_S_4630", "m06"),
        ("094_S_4649", "m06"),
        ("029_S_5219", "m24"),
        ("094_S_2238", "m48"),
        ("129_S_4287", "bl"),
        ("007_S_4611", "m03"),
        ("016_S_4638", "bl"),
        ("027_S_5118", "bl"),
        ("098_S_4018", "bl"),
        ("098_S_4003", "m12"),
        ("016_S_4584", "m24"),
        ("016_S_5007", "m12"),
        ("129_S_2347", "m06"),
        ("129_S_4220", "bl"),
        ("007_S_2058", "m12"),
        ("016_S_2007", "m06"),
        ("020_S_6358", "bl"),
        ("114_S_6039", "m12"),
        ("114_S_6057", "bl"),
        ("153_S_6274", "bl"),
        ("006_S_4485", "m84"),
        ("153_S_6237", "bl"),
        ("153_S_6336", "bl"),
        ("153_S_6450", "bl"),
        ("003_S_4441", "m12"),
        # Eq_1 Images
        ("006_S_6252", "m12"),
        ("010_S_0419", "m156"),
        ("127_S_0259", "m156"),
        ("129_S_6830", "bl"),
        ("135_S_6104", "m24"),
        ("021_S_5237", "m84"),
        ("027_S_2245", "m120"),
        ("099_S_6396", "m24"),
        ("126_S_4507", "m96"),
        ("126_S_4891", "m96"),
        ("126_S_4896", "m96"),
        ("126_S_6559", "m24"),
        ("126_S_6724", "m12"),
        ("127_S_6203", "m24"),
        ("127_S_6330", "m24"),
        ("127_S_6512", "m24"),
        ("127_S_6549", "m24"),
        ("129_S_6459", "m24"),
        ("129_S_6763", "m12"),
        ("129_S_6784", "m12"),
        ("129_S_6830", "m12"),
        ("135_S_6446", "m24"),
        ("301_S_6224", "m36"),
        # Missing .bval and .bvec
        ("020_S_5203", "m72"),
        ("020_S_6185", "m24"),
        ("020_S_6513", "m12"),
        ("021_S_0178", "m156"),
        ("153_S_6336", "m12"),
        ("153_S_6755", "bl"),
        ("020_S_6227", "m24"),
        ("020_S_6449", "m24"),
        ("020_S_6513", "m24"),
        # Volume mismatch between .nii and .bvec / .bval
        ("006_S_6610", "bl"),
        ("006_S_6682", "bl"),
        ("006_S_6696", "bl"),
        ("006_S_6770", "bl"),
        ("029_S_6289", "bl"),
        ("130_S_6043", "bl"),
        ("130_S_6329", "bl"),
        ("027_S_6183", "m24"),
        ("123_S_6891", "bl"),
        # Several output images
        ("029_S_2395", "m72"),
    ]


def _select_dwi_image(
    subject_id: str,
    timepoint: str,
    visit_all_images: list,
    mri_qc_subj: pd.DataFrame,
) -> Optional[dict[str, Any]]:
    """Select a DWI image based on QC information.

    One image among those in the input list is chosen according to QC
    and then corresponding metadata is extracted to a dictionary.

    Parameters
    ----------
    subject_id : str
        The subject identifier.

    timepoint : str
        The visit code.

    visit_all_images : list
        The list of images metadata.

    mri_qc_subj : pd.Dataframe
        DataFrame containing list of QC of scans for the subject.

    Returns
    -------
    dict or None:
        Contains image metadata.
        Returns None if no image was found.
    """
    from clinica.utils.filemanip import replace_special_characters_with_symbol

    from ._qc_utils import select_image_qc

    selected_image_id = select_image_qc(
        list(visit_all_images.image_id),
        mri_qc_subj,
    )
    if selected_image_id is None:
        return None
    selected_scan = visit_all_images[
        visit_all_images.image_id == selected_image_id
    ].iloc[0]

    return {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Sequence": replace_special_characters_with_symbol(
            selected_scan.image_description, symbol="_"
        ),
        "Scan_Date": selected_scan["image_date"],
        "Study_ID": str(int(selected_scan.study_id)),
        "Series_ID": str(int(selected_scan.series_id)),
        "Image_ID": str(int(selected_scan.image_id)),
        "Field_Strength": selected_scan.MagneticFieldStrength,
    }
