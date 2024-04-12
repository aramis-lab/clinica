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
        f"Calculating paths of fMRI field maps (FMAPs). Output will be stored in {conversion_dir}."
    )

    images = compute_fmap_path(source_dir, csv_dir, subjects, conversion_dir)
    
    cprint("Paths of field maps found. Exporting images into BIDS ...")

    paths_to_bids(images, destination_dir, "fmap", mod_to_update=mod_to_update)
    rename_fmaps(destination_dir)

    cprint(msg="Field maps conversion done.", lvl="debug")

def compute_fmap_path(source_dir: Path, csv_dir: Path, subjs_list: list[str], conversion_dir: Path) -> pd.DataFrame:
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

    import pandas as pd
    from pathlib import Path

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
    import os
    import json

    dir_name = destination_dir

    for root, dirs, _ in os.walk(destination_dir):
        for dir_name in dirs:
            if dir_name.startswith('sub-ADNI'):
                for subroot, subdir, filenames in os.walk(os.path.join(root, dir_name)):
                    for filename in filenames:
                        fmap_filename = os.path.join(subroot, filename)
                        if filename[-11:] == "fmap.nii.gz" and os.path.exists(fmap_filename):
                            mag_filename = os.path.join(subroot, filename[:-11] + "magnitude1.nii.gz")
                            os.rename(fmap_filename, mag_filename)
                        if filename[-14:] == "fmap_e2.nii.gz" and os.path.exists(fmap_filename):
                            mag_filename = os.path.join(subroot, filename[:-14] + "magnitude2.nii.gz")
                            os.rename(fmap_filename, mag_filename)
                        if filename[-9:] == "fmap.json" and os.path.exists(fmap_filename):
                            os.remove(fmap_filename)
                        if filename[-12:] == "fmap_e2.json" and os.path.exists(fmap_filename):
                            os.remove(fmap_filename)
                        if filename[-12:] == "fmap_ph.json" and os.path.exists(fmap_filename):
                            check_json(fmap_filename, filename, subroot)
                        if filename[-15:] == "fmap_e2_ph.json" and os.path.exists(fmap_filename):
                            check_json(fmap_filename, filename, subroot)               
                        if filename[-14:] == "fmap_ph.nii.gz":
                            if os.path.exists(fmap_filename[:-14] + "phase1.json"):
                                os.rename(fmap_filename, fmap_filename[:-14] + "phase1.nii.gz")
                            elif os.path.exists(fmap_filename[:-14] + "fmap_ph.json"):
                                json_filepath = fmap_filename[:-14] + "fmap_ph.json"
                                with open(json_filepath, 'r') as file:
                                    fmap_data = json.load(file)
                                    if "EchoTime1" in fmap_data and "EchoTime2" in fmap_data:
                                        os.rename(fmap_filename, fmap_filename[:-14] + "phasediff.nii.gz")
                                        os.rename(json_filepath, json_filepath[:-12] + "phasediff.json")
                                    else:
                                        os.rename(fmap_filename, fmap_filename[:-14] + "phase1.nii.gz")
                                        os.rename(json_filepath, json_filepath[:-12] + "phase1.json")
                            elif os.path.exists(fmap_filename[:-17] + "phasediff.json"):
                                json_filepath = fmap_filename[:-17] + "phasediff.json"
                                with open(json_filepath, 'r') as file:
                                    fmap_data = json.load(file)
                                    if "EchoTime1" in fmap_data and "EchoTime2" in fmap_data:
                                        os.rename(fmap_filename, fmap_filename[:-17] + "phasediff.nii.gz")
                        if filename[-17:] == "fmap_e2_ph.nii.gz":
                            if os.path.exists(fmap_filename[:-17] + "phase2.json"):
                                os.rename(fmap_filename, fmap_filename[:-17] + "phase2.nii.gz")
                            elif os.path.exists(fmap_filename[:-17] + "fmap_e2_ph.json"):
                                json_filepath = fmap_filename[:-17] + "fmap_e2_ph.json"
                                with open(json_filepath, 'r') as file:
                                    fmap_data = json.load(file)
                                    if "EchoTime1" in fmap_data and "EchoTime2" in fmap_data:
                                        os.rename(fmap_filename, fmap_filename[:-17] + "phasediff.nii.gz")
                                        os.rename(json_filepath, json_filepath[:-15] + "phasediff.json")
                                    else:
                                        os.rename(fmap_filename, fmap_filename[:-17] + "phase2.nii.gz")
                                        os.rename(json_filepath, json_filepath[:-15] + "phase2.json")
                            elif os.path.exists(fmap_filename[:-17] + "phasediff.json"):
                                json_filepath = fmap_filename[:-17] + "phasediff.json"
                                with open(json_filepath, 'r') as file:
                                    fmap_data = json.load(file)
                                    if "EchoTime1" in fmap_data and "EchoTime2" in fmap_data:
                                        os.rename(fmap_filename, fmap_filename[:-17] + "phasediff.nii.gz")

def check_json(fmap_filename, filename, subroot):
    import os
    import json

    nifti_file = fmap_filename[:-5] + "nii.gz"
    with open(fmap_filename, 'r') as file:
        fmap_data = json.load(file)
        if "EchoTime1" in fmap_data and "EchoTime2" in fmap_data:
            phase_json = os.path.join(subroot, filename[:-15] + "phasediff.json")
            os.rename(fmap_filename, phase_json)
            if os.path.exists(nifti_file):
                phase_nifti = os.path.join(subroot, filename[:-15] + "_phasediff.nii.gz")
                os.rename(nifti_file, phase_nifti)
        else:
            if filename[:-12] == "fmap_ph.json":
                phase_json = os.path.join(subroot, filename[:-12] + "phase1.json")
                os.rename(fmap_filename, phase_json)
                if os.path.exists(nifti_file):
                    phase_nifti = os.path.join(subroot, filename[:-12] + "phase1.nii.gz")
                    os.rename(nifti_file, phase_nifti)
            if filename[:-15] == "fmap_e2_ph.json":
                phase_json = os.path.join(subroot, filename[:-15] + "phase2.json")
                os.rename(fmap_filename, phase_json)
                if os.path.exists(nifti_file):
                    phase_nifti = os.path.join(subroot, filename[:-15] + "phase2.nii.gz")
                    os.rename(nifti_file, phase_nifti)            
            