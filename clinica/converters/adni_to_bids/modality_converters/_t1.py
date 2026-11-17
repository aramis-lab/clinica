"""Module for converting T1 of ADNI."""

from pathlib import Path
from typing import Any, Iterable, Optional

import pandas as pd

__all__ = ["convert_t1"]


def convert_t1(
    source_dir: Path,
    csv_dir: Path,
    destination_dir: Path,
    conversion_dir: Path,
    subjects: Iterable[str],
    force_new_extraction: bool = False,
    n_procs: int = 1,
):
    """Convert T1 MR images of ADNI into BIDS format.

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
        (
            f"Calculating paths of {ADNIModalityConverter.T1.value} images. "
            f"Output will be stored in {conversion_dir}."
        ),
        lvl="info",
    )
    images = _compute_t1_paths(source_dir, csv_dir, subjects, conversion_dir)
    cprint(
        f"Paths of {ADNIModalityConverter.T1.value} images found. Exporting images into BIDS ...",
        lvl="info",
    )
    paths_to_bids(
        images,
        destination_dir,
        ADNIModalityConverter.T1,
        force_new_extraction=force_new_extraction,
        n_procs=n_procs,
    )
    cprint(msg=f"{ADNIModalityConverter.T1.value} conversion done.", lvl="debug")


def _compute_t1_paths(
    source_dir: Path,
    csv_dir: Path,
    subjects: Iterable[str],
    conversion_dir: Path,
) -> pd.DataFrame:
    """Compute the paths to T1 MR images and store them in a TSV file.

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
        A dataframe with all the paths to the T1 MR images that will be converted into BIDS.
    """
    from clinica.utils.stream import cprint

    from .._utils import load_clinical_csv
    from ._image_path_utils import find_image_path
    from ._mri_utils import load_all_images_metadata

    t1_dfs_list = []
    t1_df = _initialize_t1_df()
    all_images = load_all_images_metadata(csv_dir)
    mri_quality = load_clinical_csv(csv_dir, "MRIQUALITY")
    mayo_mri_qc = load_clinical_csv(csv_dir, "MAYOADIRL_MRI_IMAGEQC_05_07_15")

    # Keep only T1 scans
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == "T1"]

    # We will convert the images for each subject in the subject list
    for subject in subjects:
        # Filter ADNIMERGE, MPRAGE METADATA and QC for only one subject and sort the rows/visits by examination date

        all_images_subj = all_images[
            (all_images.subject_id == subject) & (all_images.modality == "MRI")
        ]
        all_images_subj = all_images_subj.sort_values("image_date")

        mri_quality_subj = mri_quality[mri_quality.RID == int(subject[-4:])]
        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subject[-4:])]

        for _, row in all_images_subj.iterrows():
            cohort = row["phase"]
            timepoint = row["image_visit"]
            selected_scan = None

            if cohort in ("ADNI1", "ADNIGO", "ADNI2"):
                selected_scan = _select_preferred_scan_for_subject_in_adni1go2(
                    subject_id=subject,
                    timepoint=timepoint,
                    all_images_subj=all_images_subj,
                    mri_quality_subj=mri_quality_subj,
                    mayo_mri_qc_subj=mayo_mri_qc_subj,
                    preferred_field_strength=1.5 if cohort == "ADNI1" else 3.0,
                )
            elif cohort == "ADNI3":
                selected_scan = _select_preferred_scan_for_subject_in_adni3(
                    subject_id=subject,
                    timepoint=timepoint,
                    all_images_subj=all_images_subj,
                    mayo_mri_qc_subj=mayo_mri_qc_subj,
                )
            elif cohort == "ADNI4":
                selected_scan = _select_preferred_scan_for_subject_in_adni4(
                    subject_id=subject,
                    timepoint=timepoint,
                    all_images_subj=all_images_subj,
                )
            else:
                cprint(
                    f"Subject {subject} visit {timepoint} belongs to an unknown cohort: {cohort}"
                )
            if selected_scan is not None:
                t1_dfs_list.append(pd.DataFrame(selected_scan, index=["i"]))
                if "Subject_ID" not in selected_scan:
                    print("SELECTED SCAN: \n", selected_scan)
                    print("SUBJECT: ", subject, "VISCODE: ", timepoint)
    if t1_dfs_list:
        t1_df = pd.concat(t1_dfs_list, ignore_index=True)

    # Removing known exceptions from images to convert
    if not t1_df.empty:
        error_indices = t1_df.index[
            t1_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in _get_known_conversion_errors()),
                axis=1,
            )
        ]
        t1_df.drop(error_indices, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(t1_df, source_dir, modality="T1")
    images.to_csv(conversion_dir / "t1_paths.tsv", sep="\t", index=False)

    return images


def _initialize_t1_df() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "Subject_ID",
            "VISCODE",
            "Visit",
            "Sequence",
            "Scan_Date",
            "Study_ID",
            "Field_Strength",
            "Series_ID",
            "Image_ID",
            "Original",
        ]
    )


def _get_known_conversion_errors() -> Iterable[tuple[str, str]]:
    return [  # Eq_1
        ("031_S_0830", "m48"),
        ("100_S_0995", "m18"),
        ("031_S_0867", "m48"),
        ("100_S_0892", "m18"),
        ("003_S_6264", "m12"),
        ("011_S_4105", "m72"),
        # Empty folders
        # ('029_S_0845', 'm24'),
        # ('094_S_1267', 'm24'),
        # ('029_S_0843', 'm24'),
        # ('027_S_0307', 'm48'),
        # ('057_S_1269', 'm24'),
        # ('036_S_4899', 'm03'),
        # ('033_S_1016', 'm120'),
        # ('130_S_4984', 'm12'),
        # ('027_S_4802', 'm06'),
        # ('131_S_0409', 'bl'),
        # ('082_S_4224', 'm24'),
        # ('006_S_4960', 'bl'),
        # ('006_S_4960', 'm03'),
        # ('006_S_4960', 'm06'),
        # ('006_S_4960', 'm12'),
        # ('006_S_4960', 'm24'),
        # ('006_S_4960', 'm36'),
        # ('006_S_4960', 'm72'),
        # ('022_S_5004', 'bl'),
        # ('022_S_5004', 'm03'),
        # Several images: T1wa ...
        ("006_S_4485", "m84"),
        ("029_S_2395", "m72"),
        ("114_S_6039", "bl"),
        ("016_S_4952", "m48"),
    ]


def _select_preferred_scan_for_subject_in_adni1go2(
    subject_id: str,
    timepoint: str,
    all_images_subj: pd.DataFrame,
    mri_quality_subj: pd.DataFrame,
    mayo_mri_qc_subj: pd.DataFrame,
    preferred_field_strength: float = 3.0,
) -> Optional[dict[str, Any]]:
    """Select the preferred scan for a subject in a visit, given the subject belongs to ADNI 1, Go or 2 cohorts.

    Parameters
    ----------
    subject_id : str
        A string containing subject ID.

    timepoint : str
        A string of visit code in months.

    all_images_subj : pd.DataFrame
        A DataFrame of all user selected images metadata corresponding to the subject.

    mri_quality_subj : pd.DataFrame
        A DatFrame of MR image quality of images corresponding to the subject.

    mayo_mri_qc_subj : pd.DatFrame
        A DataFrame of MAYO Clinic MR image quality of images corresponding to the subject.

    preferred_field_strength : float, optional
        Field strength that is preferred in case there are several image acquisitions.

    Returns
    -------
    dict or None:
        Dictionary containing selected scan information.
        Returns None if no scan found.
    """
    from clinica.utils.filemanip import replace_special_characters_with_symbol

    # filter out images that do not pass QC
    all_images_subj = all_images_subj[
        all_images_subj.apply(
            lambda x: _check_qc(x, subject_id, timepoint, mri_quality_subj), axis=1
        )
    ]
    filtered_all_images = _get_last_processed_scan_for_visit(all_images_subj, timepoint)

    # If no N3 processed image found (it means there are no processed images at all), get best original image
    if filtered_all_images.empty:
        return _select_best_original_image(
            subject_id,
            timepoint,
            all_images_subj,
            mayo_mri_qc_subj,
            preferred_field_strength,
        )

    # If there are images with different magnetic field strength, prefer 1.5T images for ADNI1, 3.0T otherwise
    if len(filtered_all_images.MagneticFieldStrength.unique()) > 1:
        filtered_all_images = filtered_all_images[
            filtered_all_images.MagneticFieldStrength == preferred_field_strength
        ]

    # Sort by Series ID in case there are several images, so we keep the one acquired first
    filtered_all_images = filtered_all_images.sort_values("series_id")
    scan = filtered_all_images.iloc[0]

    n3 = scan.image_description.find("N3")
    # Sequence ends in 'N3' or in 'N3m'
    sequence = scan.image_description[
        : n3 + 2 + int(scan.image_description[n3 + 2] == "m")
    ]
    sequence = replace_special_characters_with_symbol(sequence, symbol="_")

    return {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Sequence": sequence,
        "Scan_Date": scan.image_date,
        "Study_ID": str(scan.study_id),
        "Series_ID": str(scan.series_id),
        "Image_ID": str(scan.image_id),
        "Field_Strength": scan.MagneticFieldStrength,
        "Original": False,
    }


def _get_last_processed_scan_for_visit(
    all_images_subj: pd.DataFrame,
    timepoint: str,
    unwanted_series_id: Optional[tuple[str]] = None,
) -> pd.DataFrame:
    """Filter content of Dataframe containing MPRAGE metadata to obtain last processed image for the specified visit.

    Parameters
    ----------
    all_images_subj : pd.DataFrame
        A DataFrame of all user selected images metadata corresponding to the subject.

    timepoint : str
        The viscode.

    unwanted_series_id : tuple or str, optional
        Series ID of scans not passing QC.

    Returns
    -------
    pd.Dataframe
    """
    unwanted_series_id = unwanted_series_id or ()
    all_images_subj = all_images_subj[
        ~all_images_subj.series_id.isin(unwanted_series_id)
    ]
    # Get the preferred scan (image series that has been Scaled)
    filtered_all_images = all_images_subj[
        (all_images_subj["image_type"] == "Pre-processed")
        & (all_images_subj.image_visit == timepoint)
        & all_images_subj.image_description.str.endswith("Scaled", na=False)
    ]
    # If no preferred image found, get N3 processed image (N3m)
    if filtered_all_images.empty:
        filtered_all_images = all_images_subj[
            (all_images_subj["image_type"] == "Pre-processed")
            & (all_images_subj.image_visit == timepoint)
            & all_images_subj.image_description.str.endswith("N3m", na=False)
        ]

    return filtered_all_images


def _select_preferred_scan_for_subject_in_adni3(
    subject_id: str,
    timepoint: str,
    all_images_subj: pd.DataFrame,
    mayo_mri_qc_subj: pd.DataFrame,
) -> Optional[dict[str, Any]]:
    """Select the preferred scan for a subject in a visit, given the subject belongs to ADNI 3 cohort.

    Parameters
    ----------
    subject_id : str
        A string containing subject ID.

    timepoint : str
        A string of visit code in months.

    all_images_subj : pd.DataFrame
        A DataFrame of all user selected images metadata corresponding to the subject.

    mayo_mri_qc_subj : pd.DataFrame
        A DatFrame of MAYO Clinic MR image quality of images corresponding to the subject.

    Returns
    -------
    dict or None :
        Dictionary containing selected scan information.
        Returns None if no scan found.
    """
    from clinica.utils.filemanip import replace_special_characters_with_symbol
    from clinica.utils.stream import cprint

    filtered_scan = all_images_subj[
        (all_images_subj["image_type"] == "Original")
        & (all_images_subj.image_visit == timepoint)
        & all_images_subj.image_description.str.contains("Accel", case=False, na=False)
        & ~all_images_subj.image_description.str.lower().str.endswith("_nd", na=False)
    ]
    if filtered_scan.empty:
        cprint(f"NO MPRAGE Meta for ADNI3: {subject_id} for visit {timepoint}")
        return None
    scan = _select_scan_from_qc(
        filtered_scan, mayo_mri_qc_subj, preferred_field_strength=3.0
    )
    sequence = replace_special_characters_with_symbol(
        scan.image_description, symbol="_"
    )

    return {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Sequence": sequence,
        "Scan_Date": scan.image_date,
        "Study_ID": str(scan.study_id),
        "Series_ID": str(scan.series_id),
        "Image_ID": str(scan.image_id),
        "Field_Strength": scan.MagneticFieldStrength,
        "Original": True,
    }


def _select_preferred_scan_for_subject_in_adni4(
    subject_id: str,
    timepoint: str,
    all_images_subj: pd.DataFrame,
) -> Optional[dict[str, Any]]:
    """Select the preferred scan for a subject in a visit, given the subject belongs to ADNI 3 cohort.

    Parameters
    ----------
    subject_id : str
        A string containing subject ID.

    timepoint : str
        A string of visit code in months.

    all_images_subj : pd.DataFrame
        A DataFrame of all user selected images metadata corresponding to the subject.

    Returns
    -------
    dict or None :
        Dictionary containing selected scan information.
        Returns None if no scan found.
    """
    from clinica.utils.filemanip import replace_special_characters_with_symbol
    from clinica.utils.stream import cprint

    filtered_scan = all_images_subj[
        (all_images_subj["image_type"] == "Original")
        & (all_images_subj.image_visit == timepoint)
        & all_images_subj.image_description.str.contains("Accel", case=False, na=False)
        & ~all_images_subj.image_description.str.lower().str.endswith("_nd", na=False)
    ]
    if filtered_scan.empty:
        cprint(f"NO MPRAGE Meta for ADNI4: {subject_id} for visit {timepoint}")
        return None
    scan = _select_scan_no_qc(filtered_scan)
    sequence = replace_special_characters_with_symbol(
        scan.image_description, symbol="_"
    )

    return {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Sequence": sequence,
        "Scan_Date": scan.image_date,
        "Study_ID": str(scan.study_id),
        "Series_ID": str(scan.series_id),
        "Image_ID": str(scan.image_id),
        "Field_Strength": scan.MagneticFieldStrength,
        "Original": True,
    }


def _select_best_original_image(
    subject_id: str,
    timepoint: str,
    all_images_subj: pd.DataFrame,
    mayo_mri_qc_subj: pd.DataFrame,
    preferred_field_strength: float = 3.0,
) -> Optional[dict[str, Any]]:
    """Select the preferred scan for a subject in a visit, given the scan is not preprocessed.

    Parameters
    ----------
    subject_id : str
        A string containing subject ID.

    timepoint : str
        A string of visit code in months.

    all_images_subj : pd.DataFrame
        A DataFrame of all user selected images metadata corresponding to the subject.

    mayo_mri_qc_subj : pd.DataFrame
        A DatFrame of MAYO Clinic MR image quality of images corresponding to the subject.

    preferred_field_strength : float, optional
        Field strength that is preferred in case there are several image acquisitions.

    Returns
    -------
    dict or None :
        Dictionary containing selected scan information.
        Returns None if no scan found.
    """
    from clinica.utils.filemanip import replace_special_characters_with_symbol
    from clinica.utils.stream import cprint

    all_images_subj_orig = all_images_subj[all_images_subj["image_type"] == "Original"]
    cond_mprage = (
        (all_images_subj_orig.image_visit == timepoint)
        & all_images_subj_orig.image_description.str.contains(
            "MP-RAGE|MPR|MPRAGE", case=False, na=False
        )
        & ~all_images_subj_orig.image_description.str.contains(
            "2", na=False
        )  # what is the meaning behind 2?
    )
    cond_spgr = (
        (all_images_subj_orig.image_visit == timepoint)
        & all_images_subj_orig.image_description.str.contains(
            "SPGR", case=False, na=False
        )  # what is the meaning behind spgr?
        & ~all_images_subj_orig.image_description.str.contains(
            "ACC", case=False, na=False
        )  # what is the meaning behind acc?
    )
    filtered_scan = all_images_subj_orig[cond_mprage | cond_spgr]

    if filtered_scan.empty:
        cprint(f"NO MPRAGE Meta: {subject_id} for visit {timepoint}")
        return None
    scan = _select_scan_from_qc(
        filtered_scan, mayo_mri_qc_subj, preferred_field_strength
    )
    sequence = replace_special_characters_with_symbol(
        scan.image_description, symbol="_"
    )

    return {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Sequence": scan.image_description,
        "Scan_Date": scan.image_date,
        "Study_ID": str(scan.study_id),
        "Series_ID": str(scan.series_id),
        "Image_ID": str(scan.image_id),
        "Field_Strength": scan.MagneticFieldStrength,
        "Original": True,
    }


def _select_scan_from_qc(
    scans_meta: pd.DataFrame,
    mayo_mri_qc_subj: pd.DataFrame,
    preferred_field_strength: float,
) -> pd.DataFrame:
    """Select a scan from a list of scans taking into account available QC.

    Parameters
    ----------
    scans_meta : pd.DataFrame
        A DataFrame containing the metadata for images to choose among.

    mayo_mri_qc_subj : pd.DataFrame
        A DatFrame of MAYO Clinic MR image quality of images corresponding to the subject.

    preferred_field_strength : pd.DataFrame
        Field strength that is preferred in case there are several image acquisitions.

    Returns
    -------
    pd.DataFrame :
        A DataFrame row containing selected scan.
    """
    import numpy as np

    multiple_mag_strength = len(scans_meta.MagneticFieldStrength.unique()) > 1

    # Select preferred_field_strength images
    if multiple_mag_strength:
        # Save for later the scans with the non preferred magnetic field strength
        not_preferred_scan = scans_meta[
            scans_meta.MagneticFieldStrength != preferred_field_strength
        ]
        # Filtering to keep only the scans with the preferred magnetic field strength
        scans_meta = scans_meta[
            scans_meta.MagneticFieldStrength == preferred_field_strength
        ]
    else:
        not_preferred_scan = None

    if scans_meta.MagneticFieldStrength.unique()[0] == 3.0:
        id_list = scans_meta.image_id.unique()
        image_ids = ["I" + str(image_id) for image_id in id_list]
        int_ids = [int(image_id) for image_id in id_list]
        images_qc = mayo_mri_qc_subj[mayo_mri_qc_subj.loni_image.isin(image_ids)]

        if not images_qc.empty:
            selected_image = None
            # Check if there is only one selected series image.
            if np.sum(images_qc.series_selected) == 1:
                selected_image = (
                    images_qc[images_qc.series_selected == 1].iloc[0].loni_image[1:]
                )
            # Otherwise, select the one with the best QC if available.
            else:
                images_not_rejected = images_qc[images_qc.series_quality < 4]

                if images_not_rejected.empty:
                    # There are no images that passed the qc,
                    # so we'll try to see if there are other images without qc.
                    # Otherwise, return None.
                    qc_ids = set(
                        [int(qc_id[1:]) for qc_id in images_qc.loni_image.unique()]
                    )
                    no_qc_ids = list(set(int_ids) - qc_ids)

                    # If none of images passed qc the scan is None, otherwise:
                    if no_qc_ids:
                        no_qc_scans_meta = scans_meta[
                            scans_meta.image_id.isin(no_qc_ids)
                        ]
                        return _select_scan_no_qc(no_qc_scans_meta)
                else:
                    # We select the image with the best (lower) QC.
                    # If no positive QC available we choose image with -1 (not performed)
                    series_quality = [
                        q if q > 0 else 4
                        for q in list(images_not_rejected.series_quality)
                    ]
                    best_q = np.amin(series_quality)
                    if best_q == 4:
                        best_q = -1
                    images_best_qc = images_not_rejected[
                        images_not_rejected.series_quality == best_q
                    ]
                    if len(images_best_qc) == 1:
                        selected_image = images_best_qc.iloc[0].loni_image[1:]
                    else:
                        best_ids = [
                            int(x[1:]) for x in images_best_qc.loni_image.unique()
                        ]
                        best_qc_meta = scans_meta[scans_meta.image_id.isin(best_ids)]
                        return _select_scan_no_qc(best_qc_meta)

            # If we did not find an image passing QC for the preferred magnetic field strength,
            # then we will choose between the other available images with other magnetic field strength
            if not selected_image and multiple_mag_strength:
                scans_meta = not_preferred_scan
            else:
                scan = scans_meta[scans_meta.image_id == int(selected_image)].iloc[0]
                return scan

    # 1.5T (no available QC)
    return _select_scan_no_qc(scans_meta)


def _select_scan_no_qc(scans_meta: pd.DataFrame) -> pd.DataFrame:
    """Select a scan from available scans in case there is not an available QC.

    Parameters
    ----------
    scans_meta : pd.DataFrame
        A DataFrame containing the metadata for images to choose among.

    Returns
    -------
    pd.DataFrame :
        DataFrame row containing selected scan.
    """
    # We choose the first scan (not containing repeat in name)
    selected_scan = scans_meta[
        ~scans_meta.image_description.str.contains("REPEAT|MPR-R", case=False, na=False)
    ]
    if selected_scan.empty:
        selected_scan = scans_meta
    scan = selected_scan.iloc[0]
    return scan


def _check_qc(
    scan: pd.DataFrame, subject_id: str, timepoint: str, mri_quality_subj: pd.DataFrame
) -> bool:
    """Check if a scan has passed check quality, if it was performed.

    Parameters
    ----------
    scan : pd.DataFrame
        A DataFrame row containing scan information.

    subject_id : str
        A string containing subject ID.

    timepoint : str
        A string containing VISCODE.

    mri_quality_subj : pd.DataFrame
        A DatFrame of MR image quality of images corresponding to the subject.

    Returns
    -------
    bool :
        Boolean, True if image passed QC or if there is no available QC, False otherwise.
    """
    from clinica.utils.stream import cprint

    # Check if QC exists for image series
    qc = mri_quality_subj[mri_quality_subj.LONIUID == "S" + str(scan.series_id)]

    # If QC exists and failed we keep the other scan (in case 2 scans were performed)
    if not qc.empty and qc.iloc[0].PASS != 1:
        cprint("QC found but NOT passed", lvl="info")
        cprint(
            msg=(
                f"Subject {subject_id} for visit {timepoint} "
                f"- Series: {str(scan.series_id)} - Study: {str(scan.study_id)}"
            ),
            lvl="info",
        )
        return False
    return True
