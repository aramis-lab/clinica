"""Module for converting field maps of ADNI."""
import json
import os
import re
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import Callable, Iterable, Optional

import pandas as pd

from clinica.utils.stream import cprint

from .._utils import load_clinical_csv

__all__ = ["convert_fmap"]


def convert_fmap(
    source_dir: PathLike,
    csv_dir: PathLike,
    destination_dir: PathLike,
    conversion_dir: PathLike,
    subjects: Iterable[str],
    force_new_extraction: bool = False,
    n_procs: int = 1,
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

    force_new_extraction : bool
        If True, pre-existing images in the BIDS directory
        will be erased and extracted again.
    """
    from clinica.utils.stream import cprint

    from .._utils import ADNIModalityConverter, paths_to_bids

    csv_dir = Path(csv_dir)
    source_dir = Path(source_dir)
    conversion_dir = Path(conversion_dir)

    cprint(
        f"Calculating paths of fMRI field maps (FMAPs). Output will be stored in {conversion_dir}.",
        lvl="debug",
    )

    images = compute_fmap_path(source_dir, csv_dir, subjects, conversion_dir)

    cprint("Paths of field maps found. Exporting images into BIDS ...")

    paths_to_bids(
        images,
        destination_dir,
        ADNIModalityConverter.FMAP,
        force_new_extraction=force_new_extraction,
    )
    reorganize_fmaps(destination_dir)

    cprint(msg="Field maps conversion done.", lvl="debug")


def compute_fmap_path(
    source_dir: Path, csv_dir: Path, subjs_list: Iterable[str], conversion_dir: Path
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

    from ._image_path_utils import find_image_path
    from ._visits_utils import visits_to_timepoints

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
        # Multiple images
        ("029_S_2395", "m72"),
        # Real/Imaginary
        ("002_S_1261", "m60"),
        ("002_S_1261", "m72"),
        ("002_S_1261", "m84"),
        ("002_S_1261", "m96"),
        ("006_S_4485", "bl"),
        ("006_S_4485", "m03"),
        ("006_S_4485", "m06"),
        ("006_S_4485", "m12"),
        ("006_S_4485", "m24"),
        ("006_S_4485", "m48"),
        # Unrecognized BIDSCase
        ("006_S_4485", "m78"),
        ("009_S_4388", "m03"),
        ("009_S_4388", "m06"),
        ("009_S_4388", "m12"),
        ("009_S_4388", "m24"),
        ("009_S_4388", "m48"),
        ("023_S_4115", "bl"),
        ("023_S_4115", "m03"),
        ("023_S_4115", "m06"),
        ("023_S_4115", "m12"),
        ("023_S_4115", "m24"),
        ("023_S_4115", "m48"),
        ("123_S_4127", "bl"),
        ("123_S_4127", "m12"),
        ("123_S_4127", "m24"),
        ("123_S_4127", "m36"),
        # Missing EchoTime Keys
        ("006_S_4485", "m90"),
        ("036_S_6088", "m12"),
        ("123_S_4127", "m84"),
        # Missing DICOMS
        ("023_S_4115", "m126"),
        ("177_S_6448", "m24"),
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
    images = find_image_path(fmap_df, source_dir, "FMAP")
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
    from clinica.utils.filemanip import replace_special_characters_with_symbol

    from ._qc_utils import select_image_qc

    mri_qc_subj.columns = [x.lower() for x in mri_qc_subj.columns]
    sel_image = select_image_qc(list(visit_mri_list.IMAGEUID), mri_qc_subj)

    if not sel_image:
        return None

    sel_scan = visit_mri_list[visit_mri_list.IMAGEUID == sel_image].iloc[0]

    image_dict = {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Visit": visit_str,
        "Sequence": replace_special_characters_with_symbol(
            sel_scan.SEQUENCE, symbol="_"
        ),
        "Scan_Date": sel_scan["SCANDATE"],
        "Study_ID": str(int(sel_scan.STUDYID)),
        "Series_ID": str(int(sel_scan.SERIESID)),
        "Image_ID": str(int(sel_scan.IMAGEUID)),
        "Field_Strength": sel_scan.MAGSTRENGTH,
    }

    return image_dict


class BIDSFMAPCase(str, Enum):
    ALREADY_RENAMED = "already_renamed"
    EMPTY_FOLDER = "empty_folder"
    ONE_PHASE_TWO_MAGNITUDES = "one_phase_two_magnitudes"
    TWO_PHASES_TWO_MAGNITUDES = "two_phases_two_magnitudes"
    DIRECT_FIELDMAPS = "direct_fieldmaps"
    NOT_SUPPORTED = "not_supported"


def phase_magnitude_renamer(old_extension: str, case: BIDSFMAPCase) -> str:
    # Assuming case 1 with 1 or 2 magnitudes have the same outputs from dcm2nix
    if (
        case != BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES
        and case != BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES
    ):
        raise NotImplementedError(
            f"No renaming should be performed for case {case.value}"
        )

    magnitude_match = re.match(pattern="^e([1-2])$", string=old_extension)
    if magnitude_match:
        return f"magnitude{magnitude_match.group(1)}"

    phase_match = re.match(pattern="^e([1-2])_ph$", string=old_extension)
    if phase_match:
        if case == BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES:
            return "phasediff"
        elif case == BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES:
            return f"phase{phase_match.group(1)}"

    raise ValueError(f"Extension {old_extension} not taken in charge.")


def rename_files(fmap_path: Path, case: BIDSFMAPCase):
    filenames = [f.name for f in fmap_path.iterdir() if not f.name.startswith(".")]
    pattern = r"(.*_)fmap_(.*?)(\..*)$"
    for previous_filename in filenames:
        rgx = re.search(pattern, previous_filename)
        if rgx:
            cut, extension, type = rgx.group(1), rgx.group(2), rgx.group(3)
            new_name = f"{cut}{phase_magnitude_renamer(extension, case)}{type}"
            os.rename(fmap_path / previous_filename, fmap_path / new_name)
        else:
            raise ValueError(f"Invalid file {previous_filename} was found.")


def fmap_case_handler_factory(case: BIDSFMAPCase) -> Callable:
    if case == BIDSFMAPCase.ALREADY_RENAMED:
        return already_renamed_handler
    if case == BIDSFMAPCase.EMPTY_FOLDER:
        return empty_folder_handler
    if case == BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES:
        return one_phase_two_magnitudes_handler
    if case == BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES:
        return two_phases_two_magnitudes_handler
    if case == BIDSFMAPCase.DIRECT_FIELDMAPS:
        return direct_fieldmaps_handler
    if case == BIDSFMAPCase.NOT_SUPPORTED:
        return not_supported_handler


def get_json_file_matching_pattern(fmap_path: Path, pattern: str) -> Path:
    json_file = [
        f
        for f in fmap_path.iterdir()
        if not f.name.startswith(".") and f.name.endswith(pattern)
    ]
    if len(json_file) != 1:
        msg = f"Expected only a single JSON file ending in '{pattern}' in {fmap_path}, but got:"
        msg += "\n".join([f.name for f in json_file])
        raise ValueError(msg)
    return json_file[0]


def check_json_contains_keys(json_file: Path, keys: Iterable[str]) -> bool:
    with open(json_file, "r") as file:
        json_data = json.load(file)
    missing_keys = [key for key in keys if key not in json_data]
    if missing_keys:
        cprint(
            f"Invalid file {json_file}, missing the following keys: {missing_keys}.",
            lvl="warning",
        )
        return False
    return True


def already_renamed_handler(fmap_path: Path):
    cprint(
        f"Files for subject {fmap_path.parent.parent.name},"
        f"session {fmap_path.parent.name} were already renamed.",
        lvl="info",
    )


def empty_folder_handler(fmap_path: Path):
    cprint(
        f"Folder for subject {fmap_path.parent.parent.name},"
        f"session {fmap_path.parent.name} is empty",
        lvl="warning",
    )


def one_phase_two_magnitudes_handler(fmap_path: Path):
    """Performs the checks and renaming associated to BIDS specifications, case 1 for fieldmaps"""
    cprint(
        f"BIDS Case 1 : expecting 1 phasediff and at least 1 magnitude files for subject"
        f"{fmap_path.parent.parent.name}, session {fmap_path.parent.name}.",
        lvl="info",
    )
    json_file = get_json_file_matching_pattern(fmap_path, pattern="ph.json")
    if check_json_contains_keys(json_file, ("EchoTime1", "EchoTime2")):
        rename_files(fmap_path, BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES)
    else:
        not_supported_handler(fmap_path)


def two_phases_two_magnitudes_handler(fmap_path: Path):
    """Performs the checks and renaming for BIDS spec case 2 for fieldmaps"""
    cprint(
        f"BIDS Case 2 : expecting 2 phase and 2 magnitude files for subject {fmap_path.parent.parent.name},"
        f"session {fmap_path.parent.name}.",
        lvl="info",
    )
    json_files = [
        f
        for f in fmap_path.iterdir()
        if not f.name.startswith(".") and f.name.endswith("ph.json")
    ]
    if all([check_json_contains_keys(f, ("EchoTime",)) for f in json_files]):
        rename_files(fmap_path, BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES)
    else:
        not_supported_handler(fmap_path)


def direct_fieldmaps_handler(fmap_path: Path):
    cprint(
        f"Assuming BIDS Case 3 for subject {fmap_path.parent.parent.name},"
        f"session {fmap_path.parent.name}. Not supported yet.",
        lvl="info",
    )
    # Assuming filename ends with "fmap" for the fieldmap
    json_file = get_json_file_matching_pattern(fmap_path, pattern="fmap.json")
    check_json_contains_keys(json_file, ("Units",))
    not_supported_handler(fmap_path)


def not_supported_handler(fmap_path: Path):
    to_delete = os.listdir(fmap_path)
    cprint(
        f"No currently supported BIDS case was recognized for subject {fmap_path.parent.parent.name},"
        f"session {fmap_path.parent.name}."
        f"The following files will be deleted : {to_delete}.",
        lvl="warning",
    )
    for file in to_delete:
        os.remove(fmap_path / file)


def infer_case_fmap(fmap_path: Path) -> BIDSFMAPCase:
    files = [f.name for f in fmap_path.iterdir() if not f.name.startswith(".")]
    nb_files = len(files)

    extensions = set(
        re.search(r"fmap(.*).json", f).group(1)
        for f in files
        if re.search(r"fmap(.*).json", f)
    )

    if nb_files == 0:
        return BIDSFMAPCase.EMPTY_FOLDER
    elif all([re.search(r"magnitude|phase", f) for f in files]):
        return BIDSFMAPCase.ALREADY_RENAMED
    elif nb_files == 4:
        if extensions == {""}:
            return BIDSFMAPCase.DIRECT_FIELDMAPS
        else:
            return BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES
    elif nb_files == 6 and extensions == {"_e1", "_e2", "_e2_ph"}:
        return BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES
    elif nb_files == 8 and extensions == {"_e1", "_e2", "_e1_ph", "_e2_ph"}:
        return BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES
    else:
        return BIDSFMAPCase.NOT_SUPPORTED


def reorganize_fmaps(bids_path: Path):
    for subject in (
        folder for folder in bids_path.iterdir() if folder.name.startswith("sub-ADNI")
    ):
        for session in (
            folder
            for folder in subject.iterdir()
            if folder.name.startswith("ses-") and folder.is_dir()
        ):
            fmap_path = bids_path / subject / session / "fmap"
            if fmap_path.exists():
                fmap_case_handler_factory(infer_case_fmap(fmap_path))(fmap_path)
