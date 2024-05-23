"""Module for converting field maps of ADNI."""
import json
import os
import re
import shutil
from os import PathLike
from pathlib import Path
from typing import List, Optional

import pandas as pd

from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv
from clinica.utils.stream import cprint


def convert_adni_fmap(
    source_dir: PathLike,
    csv_dir: PathLike,
    destination_dir: PathLike,
    conversion_dir: PathLike,
    subjects: Optional[List[str]] = None,
    mod_to_update: bool = False,
    n_procs: Optional[int] = 1,
):
    """Convert field map images of ADNI into BIDS format.

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

    subjects : List of str, optional
        List of subjects.

    mod_to_update : bool
        If True, pre-existing images in the BIDS directory
        will be erased and extracted again.
    """

    from os import path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import paths_to_bids
    from clinica.utils.stream import cprint

    csv_dir = Path(csv_dir)
    source_dir = Path(source_dir)
    conversion_dir = Path(conversion_dir)

    if not subjects:
        adni_merge = load_clinical_csv(csv_dir, "ADNIMERGE")
        subjects = list(adni_merge.PTID.unique())

    cprint(
        f"Calculating paths of fMRI field maps (FMAPs). Output will be stored in {conversion_dir}.",
        lvl="debug",
    )

    images = compute_fmap_path(source_dir, csv_dir, subjects, conversion_dir)

    cprint("Paths of field maps found. Exporting images into BIDS ...")

    paths_to_bids(images, destination_dir, "fmap", mod_to_update=mod_to_update)
    reorganize_fmaps(Path(destination_dir))

    cprint(msg="Field maps conversion done.", lvl="debug")


def compute_fmap_path(
    source_dir: Path, csv_dir: Path, subjs_list: list[str], conversion_dir: Path
) -> pd.DataFrame:
    """Compute the paths to fMR images.

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        subjs_list: subjects list
        conversion_dir: path to the TSV files including the paths to original images

    Returns:
        pandas Dataframe containing the path for each fmri
    """
    from os import path
    from pathlib import Path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        find_image_path,
        visits_to_timepoints,
    )

    fmap_col = [
        "Subject_ID",
        "VISCODE",
        "Visit",
        "Sequence",
        "Scan_Date",
        "Study_ID",
        "Field_Strength",
        "Series_ID",
        "Image_ID",
    ]
    fmap_df = pd.DataFrame(columns=fmap_col)
    fmap_dfs_list = []
    adni_merge = load_clinical_csv(csv_dir, "ADNIMERGE")

    mayo_mri_qc = load_clinical_csv(csv_dir, "MAYOADIRL_MRI_IMAGEQC_12_08_15")
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == "fMRI"]
    mayo_mri_qc.columns = [x.upper() for x in mayo_mri_qc.columns]

    mayo_mri_qc3 = load_clinical_csv(csv_dir, "MAYOADIRL_MRI_QUALITY_ADNI3")
    mayo_mri_qc3 = mayo_mri_qc3[mayo_mri_qc3.SERIES_TYPE == "EPB"]

    # Concatenating visits in both QC files
    mayo_mri_qc = pd.concat(
        [mayo_mri_qc, mayo_mri_qc3], axis=0, ignore_index=True, sort=False
    )

    mri_list = load_clinical_csv(csv_dir, "MRILIST")

    # Selecting only fMRI images that are not Multiband
    mri_list = mri_list[
        mri_list.SEQUENCE.str.contains("apping")
    ]  # 'apping' includes all field map scans, but not others

    for subj in subjs_list:
        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        adnimerge_subj = adnimerge_subj.sort_values("EXAMDATE")

        mri_list_subj = mri_list[mri_list.SUBJECT == subj]
        mri_list_subj = mri_list_subj.sort_values("SCANDATE")

        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subj[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints(subj, mri_list_subj, adnimerge_subj, "FMAP")

        for visit_info, visit_str in visits.items():
            timepoint = visit_info[0]

            visit_mri_list = mri_list_subj[mri_list_subj.VISIT == visit_str]

            image = fmap_image(
                subj, timepoint, visits[visit_info], visit_mri_list, mayo_mri_qc_subj
            )

            if image is not None:
                row_to_append = pd.DataFrame(
                    image,
                    index=[
                        "i",
                    ],
                )
                fmap_dfs_list.append(row_to_append)

        if fmap_dfs_list:
            fmap_df = pd.concat(fmap_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [
        ("006_S_4485", "m84"),
        ("123_S_4127", "m96"),
        # Eq_1
        ("094_S_4503", "m24"),
        ("009_S_4388", "m72"),
        ("036_S_6088", "bl"),
        ("036_S_6134", "bl"),
        ("016_S_6802", "bl"),
        ("016_S_6816", "bl"),
        ("126_S_4891", "m84"),
        ("177_S_6448", "m24"),
        ("023_S_4115", "m126"),
        # Multiple images
        ("029_S_2395", "m72"),
    ]

    # Removing known exceptions from images to convert
    if not fmap_df.empty:
        error_ind = fmap_df.index[
            fmap_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1
            )
        ]
        fmap_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(fmap_df, source_dir, "FMAP", "S", "Series_ID")
    images.to_csv(conversion_dir / "fmap_paths.tsv", sep="\t", index=False)

    return images


def fmap_image(
    subject_id: str,
    timepoint: str,
    visit_str: str,
    visit_mri_list: list[pd.DataFrame],
    mri_qc_subj: pd.DataFrame,
) -> Optional[dict]:
    """
    One image among those in the input list is chosen according to QC
    and then corresponding metadata is extracted to a dictionary.

    Args:
        subject_id: Subject identifier
        timepoint: Visit code
        visit_str: Visit name
        visit_mri_list: List of images metadata
        mri_qc_subj: Dataframe containing list of QC of scans for the subject

    Returns: dictionary - contains image metadata
             None - no image was able to be selected from the MRI list based on IDs provided
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        replace_sequence_chars,
        select_image_qc,
    )

    mri_qc_subj.columns = [x.lower() for x in mri_qc_subj.columns]
    sel_image = select_image_qc(list(visit_mri_list.IMAGEUID), mri_qc_subj)

    if not sel_image:
        return None

    sel_scan = visit_mri_list[visit_mri_list.IMAGEUID == sel_image].iloc[0]

    image_dict = {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Visit": visit_str,
        "Sequence": replace_sequence_chars(sel_scan.SEQUENCE),
        "Scan_Date": sel_scan["SCANDATE"],
        "Study_ID": str(int(sel_scan.STUDYID)),
        "Series_ID": str(int(sel_scan.SERIESID)),
        "Image_ID": str(int(sel_scan.IMAGEUID)),
        "Field_Strength": sel_scan.MAGSTRENGTH,
    }

    return image_dict


def renaming_fmap_extensions_case1(old_extension: str) -> str:
    # Assuming case 1 with 1 or 2 magnitudes have the same outputs from dcm2nix
    if old_extension == "e1":
        new_extension = "magnitude1"
    elif old_extension == "e2":
        new_extension = "magnitude2"
    elif old_extension == "e2_ph" or old_extension == "e1_ph":
        new_extension = "phasediff"
    return new_extension


def case1_1phase_2mag(fmap_path: Path):
    "Performs the checks and renaming associated to BIDS specifications, case 1 for fieldmaps"
    files = [f for f in os.listdir(fmap_path) if not f.startswith(".")]
    js = [f for f in files if f.endswith("ph.json")][0]

    # Checking for Echotime keys in phase json
    with open(fmap_path / js, "r") as file:
        json_data = json.load(file)
    if "EchoTime1" not in json_data or "EchoTime2" not in json_data:
        cprint(f'Invalid file {js}, missing "EchoTime" key.', lvl="warning")
        unrecognized_fmap_case(fmap_path)
        return
    # Renaming
    pattern = r"(.*_)fmap_(.*?)(\..*)"
    for previous_filename in files:
        rgx = re.search(pattern, previous_filename)
        cut, extension, type = rgx.group(1), rgx.group(2), rgx.group(3)
        new_name = cut + renaming_fmap_extensions_case1(extension) + type
        os.rename(fmap_path / previous_filename, fmap_path / new_name)


def renaming_fmap_extensions_case2(old_extension: str) -> str:
    if old_extension == "e1":
        new_extension = "magnitude1"
    elif old_extension == "e2":
        new_extension = "magnitude2"
    elif old_extension == "e1_ph":
        new_extension = "phase1"
    elif old_extension == "e2_ph":
        new_extension = "phase2"
    return new_extension


def case2_2phase_2mag(fmap_path: Path):
    "Performs the checks and renaming for BIDS spec case 2 for fieldmaps"
    files = [f for f in os.listdir(fmap_path) if not f.startswith(".")]
    check_json = [f for f in files if "ph.json" in f]
    # Checking for Echotime keys in phase jsons
    for js in check_json:
        with open(fmap_path / js, "r") as file:
            json_data = json.load(file)
        if "EchoTime" not in json_data:
            cprint(f'Invalid file {js}, missing "EchoTime" key.', lvl="warning")
            unrecognized_fmap_case(fmap_path)
            return
    # Renaming
    pattern = r"(.*_)fmap_(.*?)(\..*)"
    for previous_filename in files:
        rgx = re.search(pattern, previous_filename)
        cut, extension, type = rgx.group(1), rgx.group(2), rgx.group(3)
        new_name = cut + renaming_fmap_extensions_case2(extension) + type
        os.rename(fmap_path / previous_filename, fmap_path / new_name)


def renaming_fmap_extensions_case3(old_extension: str) -> str:
    if old_extension.endswith("ph"):
        new_extension = "fieldmap"
    else:
        new_extension = "magnitude"
    return new_extension


def direct_fieldmap(fmap_path: Path):
    """Performs the checks and renaming for BIDS spec case 3 for fieldmaps"""

    files = [f for f in os.listdir(fmap_path) if not f.startswith(".")]
    # Assuming the extension of the file ends with _ph
    check_json = [f for f in files if "ph.json" in f]
    js = [f for f in files if f.endswith("ph.json")][0]

    # Checking for Unit key in fmap json
    with open(fmap_path / js, "r") as file:
        json_data = json.load(file)
    if "Units" not in json_data:
        print(
            f'Invalid file {js} for Direct Fieldmapping, missing "Units" key.'
            f"Does not correspond to BIDS Case 3."
        )
        case1_1phase_2mag(fmap_path)
        return
    # Renaming
    pattern = r"(.*_)fmap_(.*?)(\..*)"
    for previous_filename in files:
        rgx = re.search(pattern, previous_filename)
        cut, extension, type = rgx.group(1), rgx.group(2), rgx.group(3)
        new_name = cut + renaming_fmap_extensions_case3(extension) + type
        os.rename(fmap_path / previous_filename, fmap_path / new_name)


def unrecognized_fmap_case(fmap_path: Path):
    """Deletes fmap directory"""
    to_delete = os.listdir(fmap_path)
    cprint(
        f"The following files were found in {fmap_path} : {to_delete}."
        f"They will be deleted as they are not usable.",
        lvl="info",
    )
    for file in to_delete:
        os.remove(fmap_path / file)


def check_case_fmap(fmap_path: Path) -> str or None:
    sub = fmap_path.parent.parent.name
    ses = fmap_path.parent.name

    extensions = [
        re.search(r"fmap_(.*).json", f).group(1)
        for f in os.listdir(fmap_path)
        if re.search(r"fmap_(.*).json", f)
    ]

    files = [f for f in os.listdir(fmap_path) if not f.startswith(".")]
    nb_files = len(files)

    if nb_files == 0:
        cprint(f"Folder for {sub}, {ses} is empty", lvl="warning")
        return None
    elif nb_files == 4:
        # todo : verify it can handle both case 3 and 1
        cprint(
            f"BIDS Case 1 or Case 3: expecting 1 magnitude and 1 phase or fieldmap files for {sub}, {ses}.",
            lvl="info",
        )
        return "case3"
    elif nb_files == 6 and set(extensions) == set(["e1", "e2", "e2_ph"]):
        cprint(
            f"BIDS Case 1 : expecting 1 phase and 2 magnitude files for {sub}, {ses}.",
            lvl="info",
        )
        return "case1"
    elif nb_files == 8 and set(extensions) == set(["e1", "e2", "e1_ph", "e2_ph"]):
        cprint(
            f"BIDS Case 2 : expecting 2 phase and 2 magnitude files for {sub}, {ses}.",
            lvl="info",
        )
        return "case2"
    else:
        cprint(f"No BIDS case was recognized for {sub}, {ses}.", lvl="warning")
        return "unrecognized"


def reorganize_fmaps(bids_path: Path):
    """
    Performs the renaming and suppression of fmap files according to BIDS specifications.

    Parameters
    ----------
    bids_path : Path to the BIDS dataset

    """

    bids_case = {
        "case1": case1_1phase_2mag,
        "case2": case2_2phase_2mag,
        "case3": direct_fieldmap,
        "unrecognized": unrecognized_fmap_case,
    }

    lst_subjects = [d for d in os.listdir(bids_path) if "sub-ADNI" in d]

    for subject in lst_subjects:
        lst_ses = [
            d
            for d in os.listdir(bids_path / subject)
            if "ses" in d and (bids_path / subject / d).is_dir()
        ]
        for ses in lst_ses:
            fmap_path = bids_path / subject / ses / "fmap"
            case = check_case_fmap(fmap_path)
            if case:
                bids_case[case](fmap_path)


def bids_guess(bids_path: Path):
    # todo : WIP - add where necessary to check

    bids_path = Path(bids_path)

    for file_path in bids_path.rglob(pattern=r"*.json"):
        str_file_path = str(file_path)
        if "fmap" in file_path.name:
            with open(file_path, "r") as f:
                file_json = json.load(f)
            if "BidsGuess" in file_json:
                bids_guess = file_json["BidsGuess"][-1].split("_")[-1]

                cut = re.search(r"\S*_fmap", str_file_path).group(0)
                os.rename(src=str_file_path, dst=cut + "_" + bids_guess + ".json")
                os.rename(
                    src=str_file_path.removesuffix(".json") + ".nii.gz",
                    dst=cut + "_" + bids_guess + ".nii.gz",
                )
            else:
                cprint(
                    msg=f"The FMAP file {file_path.name} could not be renamed "
                    f"since dcm2nix did not find any BIDS correspondence",
                    lvl="warning",
                )
                # todo : delete files
