"""Module for converting field maps of ADNI."""
from os import PathLike
from pathlib import Path
from typing import List, Optional

import pandas as pd

from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv


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
    # rename_fmaps(destination_dir)

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

    # csv_dir = "/Users/alice.joubert/clinicaQC/data/ADNI/ADNI_clinical_data"
    #
    # cond = mri_list.copy()
    # grouped = cond.groupby(["SUBJECT", "SCANDATE"])
    # filtered = grouped.filter(lambda x:x['IMAGEUID'].count() == 6)

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


def rename_fmaps(destination_dir: Path):
    import json
    import os

    dir_name = destination_dir

    for root, dirs, _ in os.walk(destination_dir):
        for dir_name in dirs:
            if dir_name.startswith("sub-ADNI"):
                for subroot, subdir, filenames in os.walk(os.path.join(root, dir_name)):
                    for filename in filenames:
                        fmap_filename = os.path.join(subroot, filename)

                        print("Dealing with file :", Path(fmap_filename).name)

                        cut, extension = fmap_filename.split("_fmap")

                        if extension.startswith("."):
                            print("Not handled yet")
                            # todo : means what format ??
                            return

                        elif extension.startswith("_"):
                            extension = extension.removeprefix("_")

                        # todo : json and nii handled differently ?
                        # todo : need to check if files exist ?? where ?
                        if extension == "e1.nii.gz":
                            print("Renaming in magnitude1")
                            os.rename(fmap_filename, cut + "_magnitude1.nii.gz")
                        if extension == "e2.nii.gz":
                            os.rename(fmap_filename, cut + "_magnitude2.nii.gz")
                            print("Renaming in magnitude2")
                        if extension == "e1.json" or extension == "e2.json":
                            print("Deleting")
                            os.remove(fmap_filename)
                        if extension == "e1_ph.json" or extension == "e2_ph.json":
                            print("Check json")
                            check_json(fmap_filename, filename, subroot)
                        if extension == "e1_ph.nii.gz":
                            if os.path.exists(cut + "_phase1.json"):
                                print("Renaming in phase1")
                                os.rename(fmap_filename, cut + "_phase1.nii.gz")
                            elif os.path.exists(cut + "_fmap_e1_ph.json"):
                                json_filepath = cut + "_fmap_e1_ph.json"
                                with open(json_filepath, "r") as file:
                                    fmap_data = json.load(file)
                                    if (
                                        "EchoTime1" in fmap_data
                                        and "EchoTime2" in fmap_data
                                    ):
                                        print("Renaming in phasediff")
                                        os.rename(
                                            fmap_filename, cut + "_phasediff.nii.gz"
                                        )
                                        os.rename(
                                            json_filepath, cut + "_phasediff.json"
                                        )
                                    else:
                                        print("Renaming in phase1")
                                        os.rename(fmap_filename, cut + "_phase1.nii.gz")
                                        os.rename(json_filepath, cut + "_phase1.json")
                            elif os.path.exists(cut + "_phasediff.json"):
                                json_filepath = cut + "_phasediff.json"
                                with open(json_filepath, "r") as file:
                                    fmap_data = json.load(file)
                                    if (
                                        "EchoTime1" in fmap_data
                                        and "EchoTime2" in fmap_data
                                    ):
                                        print("Renaming in phasediff")
                                        os.rename(
                                            fmap_filename, cut + "_phasediff.nii.gz"
                                        )
                        if extension == "e2_ph.nii.gz":
                            if os.path.exists(cut + "_phase2.json"):
                                print("Renaming in phase2")
                                os.rename(fmap_filename, cut + "_phase2.nii.gz")
                            elif os.path.exists(cut + "_fmap_e2_ph.json"):
                                json_filepath = cut + "_fmap_e2_ph.json"
                                with open(json_filepath, "r") as file:
                                    fmap_data = json.load(file)
                                    if (
                                        "EchoTime1" in fmap_data
                                        and "EchoTime2" in fmap_data
                                    ):
                                        print("Renaming in phasediff")
                                        os.rename(
                                            fmap_filename, cut + "_phasediff.nii.gz"
                                        )
                                        os.rename(
                                            json_filepath, cut + "_phasediff.json"
                                        )
                                    else:
                                        print("Renaming in phase2")
                                        os.rename(fmap_filename, cut + "_phase2.nii.gz")
                                        os.rename(json_filepath, cut + "_phase2.json")
                            elif os.path.exists(cut + "_phasediff.json"):
                                json_filepath = cut + "_phasediff.json"
                                with open(json_filepath, "r") as file:
                                    fmap_data = json.load(file)
                                    if (
                                        "EchoTime1" in fmap_data
                                        and "EchoTime2" in fmap_data
                                    ):
                                        print("Renaming in phasediff")
                                        os.rename(
                                            fmap_filename, cut + "_phasediff.nii.gz"
                                        )


def check_json(json_filename, filename, subroot):
    """

    Parameters
    ----------
    json_filename : json
    filename
    subroot

    Returns
    -------

    """
    import json
    import os

    nifti_file = json_filename.rsplit(".", 1)[0] + ".nii.gz"
    cut, extension = json_filename.split("_fmap")

    # todo : attention !! CUT will work only if was not renamed bc renaming removes 'fmap'

    if extension.startswith("."):
        print("Not handled yet")
        # todo : means what format ??
        return
    elif extension.startswith("_"):
        extension = extension.removeprefix("_")

    if os.path.exists(json_filename):
        with open(json_filename, "r") as file:
            fmap_data = json.load(file)
            if "EchoTime1" in fmap_data and "EchoTime2" in fmap_data:
                os.rename(json_filename, cut + "_phasediff.json")
                if os.path.exists(nifti_file):
                    os.rename(nifti_file, cut + "_phasediff.nii.gz")
            elif extension == "e1_ph.json":
                os.rename(json_filename, cut + "_phase1.json")
                if os.path.exists(nifti_file):
                    os.rename(nifti_file, cut + "_phase1.nii.gz")
            elif extension == "e2_ph.json":
                os.rename(json_filename, cut + "_phase2.json")
                if os.path.exists(nifti_file):
                    os.rename(nifti_file, cut + "_phase2.nii.gz")

    else:
        print("The following file was already modified : ", Path(json_filename).name)


def bids_guess(bids_path: Path):
    import json
    import os
    import re
    from pathlib import Path

    from clinica.utils.stream import cprint

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
                # todo : adapt this bc breaks code later if suffix not usual


def reorganize_fmap_files(bids_path: Path):
    # todo : should suppress sessions that contain only one nifti
    # should verify that fmap got sessions. then what, abort process?
    # also should pay attention to what's written in scans.tsv and sessions.tsv
    return
