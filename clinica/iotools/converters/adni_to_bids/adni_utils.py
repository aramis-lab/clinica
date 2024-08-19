from enum import Enum
from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd

from clinica.utils.pet import ReconstructionMethod, Tracer
from clinica.utils.stream import cprint

__all__ = [
    "get_subjects_list",
    "ADNIModality",
    "ADNIModalityConverter",
    "ADNIStudy",
    "correct_diagnosis_sc_adni3",
    "create_adni_sessions_dict",
    "create_adni_scans_files",
    "paths_to_bids",
    "load_clinical_csv",
]


class ADNIStudy(str, Enum):
    """Possible versions of ADNI studies."""

    ADNI1 = "ADNI1"
    ADNI2 = "ADNI2"
    ADNI3 = "ADNI3"
    ADNIGO = "ADNIGO"


def _define_subjects_list(
    source_dir: Path,
    subjs_list_path: Optional[Path] = None,
) -> List[str]:
    import re

    from clinica.utils.stream import cprint

    if subjs_list_path:
        cprint("Loading a subjects lists provided by the user...")
        return subjs_list_path.read_text().splitlines()

    cprint(f"Using the subjects contained in the ADNI dataset at {source_dir}")
    rgx = re.compile(r"\d{3}_S_\d{4}")
    return list(filter(rgx.fullmatch, [folder.name for folder in source_dir.iterdir()]))


def _check_subjects_list(
    subjs_list: List[str],
    clinical_dir: Path,
) -> List[str]:
    from copy import copy

    from clinica.utils.stream import cprint

    subjs_list_copy = copy(subjs_list)
    adni_merge = load_clinical_csv(clinical_dir, "ADNIMERGE")
    # Check that there are no errors in subjs_list given by the user
    for subj in subjs_list_copy:
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        if len(adnimerge_subj) == 0:
            cprint(
                msg=f"Subject with PTID {subj} does not have corresponding clinical data."
                f"Please check your subjects list or directory.",
                lvl="warning",
            )
            subjs_list.remove(subj)
    del subjs_list_copy

    if not subjs_list:
        cprint(f"Processing an empty list of subjects.", lvl="warning")

    return subjs_list


def get_subjects_list(
    source_dir: Path,
    clinical_dir: Path,
    subjs_list_path: Optional[Path] = None,
) -> List[str]:
    return _check_subjects_list(
        _define_subjects_list(source_dir, subjs_list_path), clinical_dir
    )


class ADNIModality(str, Enum):
    """Possible modalities supported by the ADNI-to-BIDS converter.

    These are the modalities exposed to the user. There is not a
    one-to-one relationship with the modality converters. That is,
    some modalities, like PET_AMYLOID, are associated with multiple
    converters, while others are associated with only one converter.
    For this reason, the ADNIModalityConverter enumeration exists,
    and there is a mapping between a ADNIModality and an iterable of
    ADNIModalityConverter variants.
    """

    T1 = "T1"
    PET_FDG = "PET_FDG"
    PET_AMYLOID = "PET_AMYLOID"
    PET_TAU = "PET_TAU"
    DWI = "DWI"
    FLAIR = "FLAIR"
    FMRI = "fMRI"
    FMAP = "FMAP"


class ADNIModalityConverter(str, Enum):
    """Possible modality converters for ADNI.

    These are not exposed to the user. However, there is a one-to-one
    relationship with the modality converters. That is, each variant
    has a corresponding converter.
    """

    T1 = "T1"
    PET_FDG = "PET_FDG"
    PET_FDG_UNIFORM = "PET_FDG_UNIFORM"
    PET_PIB = "PET_PIB"
    PET_AV45 = "PET_AV45"
    PET_TAU = "PET_TAU"
    DWI = "DWI"
    FLAIR = "FLAIR"
    FMRI = "fMRI"
    FMAP = "FMAP"


def correct_diagnosis_sc_adni3(
    clinical_data_dir: Path, participants_df: pd.DataFrame
) -> pd.DataFrame:
    """Add missing diagnosis at screening in ADNI3 due to changes in DX_bl.

    See https://groups.google.com/g/adni-data/c/whYhKafN_0Q.
    The identification of SMC, EMCI, LMCI is lost in ADNI3.

    Args:
        clinical_data_dir: path to the input folder containing clinical data.
        participants_df: DataFrame containing metadata at the participant level.

    Returns:
        Corrected participants_df.
    """
    diagnosis_dict = {1: "CN", 2: "MCI", 3: "AD"}
    # DXSUM_PDXCONV_ADNIALL has been renamed to DXSUM_PDXCONV
    # this ensures old ADNI downloads still work with recent versions of Clinica
    try:
        dxsum_df = load_clinical_csv(
            clinical_data_dir, "DXSUM_PDXCONV_ADNIALL"
        ).set_index(["PTID", "VISCODE2"])
    except OSError:
        dxsum_df = load_clinical_csv(clinical_data_dir, "DXSUM_PDXCONV").set_index(
            ["PTID", "VISCODE2"]
        )
    missing_sc = participants_df[
        participants_df.original_study == ADNIStudy.ADNI3.value
    ]
    participants_df.set_index("alternative_id_1", drop=True, inplace=True)
    for alternative_id in missing_sc.alternative_id_1.values:
        try:
            diagnosis_sc = diagnosis_dict[
                dxsum_df.loc[(alternative_id, "sc"), "DIAGNOSIS"].values[0]
            ]
            participants_df.loc[alternative_id, "diagnosis_sc"] = diagnosis_sc
        except KeyError:
            cprint(
                msg=f"Unknown screening diagnosis for subject {alternative_id}.",
                lvl="warning",
            )
            participants_df.loc[alternative_id, "diagnosis_sc"] = "n/a"

    participants_df.reset_index(inplace=True, drop=False)
    return participants_df


def _write_adni_sessions_tsv(
    df_subj_sessions: pd.DataFrame, bids_subjs_paths: list[Path]
):
    """Write the result of method create_session_dict into several TSV files.

    Args:
        df_subj_sessions: global dataframe containing clinical sessions data for all subjects
        bids_subjs_paths: a list with the path to all bids subjects
    """

    df_subj_sessions["adas_memory"] = (
        df_subj_sessions["adas_Q1"]
        + df_subj_sessions["adas_Q4"]
        + df_subj_sessions["adas_Q7"]
        + df_subj_sessions["adas_Q8"]
        + df_subj_sessions["adas_Q9"]
    )  # / 45
    df_subj_sessions["adas_language"] = (
        df_subj_sessions["adas_Q2"]
        + df_subj_sessions["adas_Q5"]
        + df_subj_sessions["adas_Q10"]
        + df_subj_sessions["adas_Q11"]
        + df_subj_sessions["adas_Q12"]
    )  # / 25
    df_subj_sessions["adas_praxis"] = (
        df_subj_sessions["adas_Q3"] + df_subj_sessions["adas_Q6"]
    )  # / 10
    df_subj_sessions["adas_concentration"] = df_subj_sessions["adas_Q13"]  # / 5

    df_subj_sessions = df_subj_sessions.fillna("n/a")

    # compute the amyloid and ptau status
    df_subj_sessions["a_stat"] = df_subj_sessions.apply(
        lambda x: _compute_amyloid_status(
            x[["adni_av45", "adni_pib", "adni_abeta"]].apply(
                pd.to_numeric, errors="coerce"
            )
        ),
        axis=1,
    )
    df_subj_sessions["tau_stat"] = df_subj_sessions.apply(
        lambda x: _compute_ptau_status(
            x[["adni_ptau"]].apply(pd.to_numeric, errors="coerce")
        ),
        axis=1,
    )
    for subject_path in bids_subjs_paths:
        subject_path.mkdir(parents=True, exist_ok=True)
        if subject_path.name == "conversion_info":
            continue
        df_tmp = df_subj_sessions[df_subj_sessions["RID"] == subject_path.name[12:]]
        df_tmp.to_csv(
            subject_path / f"{subject_path.name}_sessions.tsv",
            sep="\t",
            index=False,
            encoding="utf-8",
        )


def _compute_amyloid_status(row: pd.DataFrame) -> str:
    if (
        pd.isnull(row["adni_av45"])
        and pd.isnull(row["adni_pib"])
        and isinstance(row["adni_abeta"], float)
    ):
        return "Au"
    if row["adni_av45"] > 1.1 or row["adni_pib"] > 1.5 or row["adni_abeta"] < 192:
        return "A+"
    if (
        (pd.isnull(row["adni_av45"]) or row["adni_av45"] < 1.1)
        and (pd.isnull(row["adni_pib"]) or row["adni_pib"] < 1.5)
        and (isinstance(row["adni_abeta"], float) or row["adni_abeta"] > 192)
    ):
        return "A-"
    return "n/a"


def _compute_ptau_status(row: pd.DataFrame) -> str:
    if pd.isnull(row["adni_ptau"]):
        return "Tu"
    if row["adni_ptau"] > 23:
        return "T+"
    return "T-"


def _filter_subj_bids(
    df_files: pd.DataFrame, location: str, bids_ids: list[str]
) -> pd.DataFrame:
    from clinica.iotools.bids_utils import remove_space_and_symbols

    # Depending on the file that needs to be open, identify and
    # preprocess the column that contains the subjects ids.
    # todo : use id class here ?
    bids_ids = [x[8:] for x in bids_ids if "sub-ADNI" in x]
    if location == "ADNIMERGE.csv":
        df_files["RID"] = df_files["PTID"].apply(
            lambda x: remove_space_and_symbols(x)[4:]
        )
    else:
        df_files["RID"] = df_files["RID"].apply(lambda x: _pad_id(x))
    return df_files.loc[df_files["RID"].isin([x[4:] for x in bids_ids])]


def _update_age(row: pd.DataFrame) -> float:
    """Update age with time passed since bl to current visit"""
    from datetime import datetime

    if row["session_id"] != "ses-M000":
        exam_date = datetime.strptime(row["EXAMDATE"], "%Y-%m-%d")
        exam_date_baseline = datetime.strptime(row["EXAMDATE_bl"], "%Y-%m-%d")
        delta = exam_date - exam_date_baseline
        return round(
            float(row["AGE"]) + (delta.days / 365.25),
            1,
        )
    return row["AGE"]


def _pad_id(ref_id: str) -> str:
    """Add leading zeros to the RID to keep à length of 4"""
    rid = str(ref_id)
    # Fill the rid with the needed number of zero
    if 4 - len(rid) > 0:
        zeros_to_add = 4 - len(rid)
        rid = "0" * zeros_to_add + rid
    return rid


def _compute_session_id(df: pd.DataFrame, csv_filename: str) -> pd.DataFrame:
    """Compute session ID from visit code.

    The CSV file name is used to determine in which column of
    the dataframe the visit code should be found.
    """
    visit_code_column = _get_visit_code_column_name(csv_filename)
    if visit_code_column not in df:
        raise ValueError(
            f"DataFrame does not contain a column named '{visit_code_column}', "
            "which is supposed to encode the visit code. The columns present in "
            f"the DataFrame are: {df.columns}."
        )
    return df.assign(
        session_id=lambda _df: _df[visit_code_column].apply(
            lambda x: _get_session_id_from_visit_code(x)
        )
    )


def _get_visit_code_column_name(csv_filename: str) -> str:
    """Return the name of the column containing the visit code."""
    if _is_a_visit_code_2_type(csv_filename):
        return "VISCODE2"
    if _is_a_time_point_type(csv_filename):
        return "Timepoint"
    return "VISCODE"


def _is_a_visit_code_2_type(csv_filename: str) -> bool:
    """If the csv file is among these files, then the visit code column is 'VISCODE2'."""
    return csv_filename in {
        "ADAS_ADNIGO2.csv",
        "DXSUM_PDXCONV.csv",
        "DXSUM_PDXCONV_ADNIALL.csv",
        "CDR.csv",
        "NEUROBAT.csv",
        "GDSCALE.csv",
        "MODHACH.csv",
        "MOCA.csv",
        "NPIQ.csv",
        "MEDHIST.csv",
        "VITALS.csv",
        "UWNPSYCHSUM_03_07_19.csv",
        "UWNPSYCHSUM_03_26_20.csv",
        "ECOGPT.csv",
        "ECOGSP.csv",
        "FCI.csv",
        "CCI.csv",
        "NPIQ.csv",
        "NPI.csv",
    }


def _is_a_time_point_type(csv_filename: str) -> bool:
    """If the csv file is among these files, then the visit code column is 'Timepoint'."""
    return csv_filename in {
        "BHR_EVERYDAY_COGNITION.csv",
        "BHR_BASELINE_QUESTIONNAIRE.csv",
        "BHR_LONGITUDINAL_QUESTIONNAIRE.csv",
    }


def _get_session_id_from_visit_code(visit_code: str) -> Optional[str]:
    """Checks that the visit code found in the data is supported by
    the converter utility function `viscode_to_session` before using it.
    There are a few special cases that are handled here.
    """
    from clinica.iotools.converter_utils import viscode_to_session

    if _is_visit_code_not_supported(visit_code):
        return None
    if visit_code == "sc":
        return "sc"
    return viscode_to_session(visit_code)


def _is_visit_code_not_supported(visit_code: str) -> bool:
    """Return True is the visit code is not supported by the converter.

    There are a few known values like "f" or "uns1" which are present in
    ADNI data that are not supported for a mapping to a session ID.
    """
    unsupported_values = {"f", "uns1"}
    return pd.isnull(visit_code) or visit_code in unsupported_values


def create_adni_sessions_dict(
    bids_ids: list[str],
    clinical_specifications_folder: Path,
    clinical_data_dir: Path,
    bids_subjs_paths: list[Path],
):
    """Extract all the data required for the sessions files and organize them in a dictionary.

    Args:
        bids_ids: list of subject IDs to process
        clinical_specifications_folder: path to the specifications for converting the clinical data
        clinical_data_dir: path to the clinical data folder
        bids_subjs_paths: a list with the path to all the BIDS subjects
    """
    from clinica.utils.stream import cprint

    df_sessions = pd.read_csv(clinical_specifications_folder / "sessions.tsv", sep="\t")
    files = list(df_sessions["ADNI location"].unique())
    files = [x for x in files if not pd.isnull(x)]
    df_subj_session = pd.DataFrame()
    # write line to get field_bids = sessions['BIDS CLINICA'] without the null values

    # Iterate over the metadata files

    for location in files:
        location = location.split("/")[0]
        try:
            df_file = load_clinical_csv(clinical_data_dir, location.split(".")[0])
        except IOError:
            continue
        df_filtered = _filter_subj_bids(df_file, location, bids_ids).copy()
        if not df_filtered.empty:
            df_filtered = _compute_session_id(df_filtered, location)
            # Filter rows with invalid session IDs.
            df_filtered.dropna(subset="session_id", inplace=True)

            if location == "ADNIMERGE.csv":
                df_filtered["AGE"] = df_filtered.apply(lambda x: _update_age(x), axis=1)
            df_subj_session = _update_sessions_df(
                df_subj_session, df_filtered, df_sessions, location
            )
        else:
            cprint(
                f"Clinical dataframe extracted from {location} is empty after filtering.",
                lvl="debug",
            )
            dict_column_correspondence = dict(
                zip(df_sessions["ADNI"], df_sessions["BIDS CLINICA"])
            )
            df_filtered.rename(columns=dict_column_correspondence, inplace=True)
            df_filtered = df_filtered.loc[
                :, (~df_filtered.columns.isin(df_subj_session.columns))
            ]
            df_subj_session = pd.concat([df_subj_session, df_filtered], axis=1)
    if df_subj_session.empty:
        cprint(
            "Empty dataset detected. Clinical data cannot be extracted.", lvl="warning"
        )
        return
    # Nv/None refer to sessions whose session is undefined. "sc" is the screening session with unreliable (incomplete)
    # data.
    df_subj_session = df_subj_session[
        df_subj_session.session_id.str.contains(r"ses-M\d+")
    ]
    _write_adni_sessions_tsv(df_subj_session, bids_subjs_paths)


def _update_sessions_df(
    df_subj_session: pd.DataFrame,
    df_filtered: pd.DataFrame,
    df_sessions: pd.DataFrame,
    location: str,
) -> pd.DataFrame:
    """Update the sessions dataframe with data of current subject.

    Args:
        df_subj_session: dataframe containing aggregate sessions
        df_filtered: dataframe with current subject sessions data
        df_sessions: dataframe with the the metadata information
        location: the clinica data filename
    """
    df_columns_to_add = df_sessions[
        df_sessions["ADNI location"].str.contains(location, na=False)
    ][["BIDS CLINICA", "ADNI"]]

    df_columns_to_add = df_columns_to_add[
        (df_columns_to_add["ADNI"].isin(df_filtered.columns))
        & (~df_columns_to_add["BIDS CLINICA"].isin(df_subj_session.columns))
    ]
    df_temp = df_filtered[
        ["RID", "session_id"] + list(df_columns_to_add["ADNI"])
    ].copy()
    df_temp.columns = ["RID", "session_id"] + list(df_columns_to_add["BIDS CLINICA"])

    # if error in adni data (duplicate session id), keep only the first row
    df_temp.drop_duplicates(subset=["RID", "session_id"], keep="first", inplace=True)
    try:
        df_temp["diagnosis"] = df_temp["diagnosis"].apply(
            lambda x: _convert_diagnosis_code(x)
        )
    except KeyError:
        pass
    if df_subj_session.empty:
        df_subj_session = df_temp
    else:
        df_subj_session = pd.merge(
            df_temp, df_subj_session, on=["RID", "session_id"], how="outer"
        )
    return df_subj_session


def _convert_diagnosis_code(diagnosis_code: str) -> str:
    """Convert the string field DX contained in ADNIMERGE.csv into a code that identify the diagnosis.

    Args:
        diagnosis_code: a string that represents the diagnosis

    Returns:
        A code that identify a diagnosis
    """
    diagnosis = {"CN": "CN", "MCI": "MCI", "Dementia": "AD"}

    if pd.isna(diagnosis_code):
        return diagnosis_code
    try:
        return diagnosis[diagnosis_code]
    except KeyError:
        raise KeyError(f"Unknown diagnosis code {diagnosis_code}.")


def create_adni_scans_files(conversion_path: Path, bids_subjs_paths: list[Path]):
    """Create scans.tsv files for ADNI.

    Args:
        conversion_path: path to the conversion_info folder
        bids_subjs_paths: list of bids subject paths
    """
    from os import path

    from clinica.iotools.bids_utils import StudyName, bids_id_factory
    from clinica.utils.stream import cprint

    scans_fields_bids = ["filename", "scan_id", "mri_field"]

    conversion_versions = [
        folder.name
        for folder in conversion_path.iterdir()
        if folder.name[0] == "v" and (conversion_path / folder).is_dir()
    ]
    conversion_versions = sorted(
        conversion_versions, key=lambda x: int(x.split("v")[1])
    )
    older_version = conversion_versions[-1]
    converted_dict = dict()
    for tsv_path in (conversion_path / older_version).iterdir():
        modality = tsv_path.name.split("_paths")[0]
        df = pd.read_csv(conversion_path / older_version / tsv_path, sep="\t")
        df.set_index(["Subject_ID", "VISCODE"], inplace=True, drop=True)
        converted_dict[modality] = df

    for bids_subject_path in bids_subjs_paths:
        # Create the file
        bids_id = bids_subject_path.resolve().name
        subject_id = bids_id_factory(StudyName.ADNI)(bids_id).to_original_study_id()
        for session_path in bids_subject_path.glob("ses-*"):
            viscode = _session_label_to_viscode(session_path.name[4::])
            tsv_name = f"{bids_id}_{session_path.name}_scans.tsv"
            (session_path / tsv_name).unlink(missing_ok=True)
            scans_tsv = open(session_path / tsv_name, "a")
            scans_df = pd.DataFrame(columns=scans_fields_bids)
            scans_df.to_csv(scans_tsv, sep="\t", index=False, encoding="utf-8")

            # Extract modalities available for each subject
            for mod in session_path.glob("*"):
                for file in mod.glob("*"):
                    scans_df = pd.DataFrame(index=[0], columns=scans_fields_bids)
                    scans_df["filename"] = path.join(mod.name, file.name)
                    converted_mod = _find_conversion_mod(file.name)
                    conversion_df = converted_dict[converted_mod]
                    try:
                        scan_id = conversion_df.loc[(subject_id, viscode), "Image_ID"]
                        scans_df["scan_id"] = scan_id
                        if "Field_Strength" in conversion_df.columns.values:
                            field_strength = conversion_df.loc[
                                (subject_id, viscode), "Field_Strength"
                            ]
                            scans_df["mri_field"] = field_strength
                    except KeyError:
                        cprint(
                            msg=f"No information found for file {file.name}",
                            lvl="warning",
                        )
                    scans_df = scans_df.fillna("n/a")
                    scans_df.to_csv(
                        scans_tsv, header=False, sep="\t", index=False, encoding="utf-8"
                    )


def _find_conversion_mod(file_name: str) -> str:
    from clinica.utils.pet import Tracer

    suffix = file_name.split("_")[-1].split(".")[0]

    if suffix == "T1w":
        return "t1"
    if suffix == "FLAIR":
        return "flair"
    if suffix == "dwi":
        return "dwi"
    if suffix == "bold":
        return "fmri"
    if suffix in (
        "fmap",
        "ph",
        "e1",
        "e2",
        "phase1",
        "phase2",
        "phasediff",
        "magnitude1",
        "magnitude2",
    ):
        return "fmap"
    if suffix == "pet":
        tracer = file_name.split("trc-")[1].split("_")[0]
        if tracer in (Tracer.AV45, Tracer.FBB):
            return "amyloid_pet"
        return f"{tracer}_pet"
    raise ValueError(f"Conversion modality could not be found for file {file_name}")


def paths_to_bids(
    images: pd.DataFrame,
    bids_dir: Path,
    modality: ADNIModalityConverter,
    mod_to_update: bool = False,
    n_procs: Optional[int] = 1,
) -> List[Path]:
    """Images in the list are converted and copied to directory in BIDS format.

    Parameters
    ----------
    images : pd.DataFrame
        List of images metadata and paths.

    bids_dir : Path
        The path to the output BIDS directory.

    modality : str
        Imaging modality.

    mod_to_update : bool
        If True, pre-existing images in the BIDS directory will be erased
        and extracted again.

    n_procs : int, optional
        The number of processes to use for conversion.
        Default=1.

    Returns
    -------
    output_file_treated : list of paths
        List of path to image files created.
        This list contains None values for files where the
        conversion wasn't successful.
    """
    from functools import partial
    from multiprocessing import Pool

    images_list = list([data for _, data in images.iterrows()])
    create_file_ = partial(
        _create_file,
        modality=modality,
        bids_dir=bids_dir,
        mod_to_update=mod_to_update,
    )
    # If n_procs==1 do not rely on a Process Pool to enable classical debugging
    if n_procs == 1:
        return [create_file_(image) for image in images_list]
    with Pool(processes=n_procs) as pool:
        output_file_treated = pool.map(create_file_, images_list)
    return output_file_treated


def _get_images_with_suffix(
    folder: Path,
    suffixes: Iterable[str],
) -> Iterable[Path]:
    import re

    if suffixes:
        rgx = re.compile(r"|".join(suffixes))
        images = list(filter(rgx.search, [f.name for f in folder.iterdir()]))
        return [folder / image for image in images]
    return []


def _remove_existing_images_if_necessary(
    folder: Path,
    suffixes: Iterable[str],
    mod_to_update: bool = False,
) -> bool:
    """Returns whether it was possible to delete the existing images or not."""
    from clinica.utils.stream import cprint

    images_to_remove = _get_images_with_suffix(folder, suffixes)
    if not images_to_remove:
        return True
    elif mod_to_update:
        cprint(f"Removing old images : {images_to_remove}", lvl="info")
        for path_to_unlink in images_to_remove:
            (path_to_unlink).unlink()
        return True
    cprint(
        f"There already exist images : {images_to_remove}. "
        "The parameter 'mod_to_update' is set to False so that "
        "they cannot be overwritten.",
        lvl="warning",
    )
    return False


def _remove_files_with_unsupported_suffixes(
    folder: Path, suffixes: Iterable[str]
) -> Iterable[Path]:
    from clinica.utils.stream import cprint

    invalid_suffix_images = _get_images_with_suffix(folder, suffixes)
    for image_to_remove in invalid_suffix_images:
        cprint(
            f"Image with bad suffix {', '.join(suffixes)} was generated by dcm2nix. These are not"
            f"supported by Clinica so the image will NOT be converted.",
            lvl="warning",
        )
        cprint(f"Removing image {folder / image_to_remove}", lvl="info")
        (folder / image_to_remove).unlink()
    return [folder / image for image in invalid_suffix_images]


def _create_file(
    image: pd.Series,
    modality: ADNIModalityConverter,
    bids_dir: Path,
    mod_to_update: bool,
) -> Optional[Path]:
    """Creates an image file at the corresponding output folder.

    The image file is the result of image conversion (DICOM to NIFTI)
    and centering, as specified for each modality

    Parameters
    ----------
    image : pd.Series
        The image metadata.

    modality : str
        Imaging modality.

    bids_dir : Path
        The path to the output BIDS directory.

    mod_to_update : bool
        If True, pre-existing images in the BIDS directory will be
        erased and extracted again.

    Returns
    -------
    output_image : Path
        The path to the output image created.
        If the conversion wasn't successful, this function returns None.
    """
    import glob
    import os
    import re
    import shutil

    import numpy as np

    from clinica.cmdline import setup_clinica_logging
    from clinica.iotools.bids_utils import StudyName, bids_id_factory, run_dcm2niix
    from clinica.iotools.converter_utils import viscode_to_session
    from clinica.iotools.utils.data_handling import center_nifti_origin
    from clinica.utils.stream import cprint

    # This function is executed in a multiprocessing context
    # such that we need to re-configure the clinica logger in the child processes.
    # Note that logging messages could easily be lost (for example when logging
    # to a file from two different processes). A better solution would be to
    # implement a logging process consuming logging messages from a multiprocessing.Queue...
    setup_clinica_logging("INFO")

    subject = image.Subject_ID
    viscode = image.VISCODE

    image_tracer = None
    logging_header = f"[{modality.value}]"
    if modality == ADNIModalityConverter.PET_AV45:
        image_tracer = Tracer(image.Tracer)
        logging_header = f"[{modality.value} - {image_tracer.value}]"
    if not image.Path:
        # todo : check for message redundancy
        cprint(
            f"{logging_header} No path specified for {subject} in session {viscode}",
            lvl="info",
        )
        return None
    cprint(
        f"{logging_header} Processing subject {subject} in session {viscode}",
        lvl="info",
    )
    session = viscode_to_session(viscode)
    image_path = Path(image.Path)
    image_id = image.Image_ID
    # If the original image is a DICOM, check if contains two DICOM inside the same folder
    if image.Is_Dicom:
        image_path = _check_two_dcm_folder(image_path, bids_dir, image_id)
    bids_id = bids_id_factory(StudyName.ADNI).from_original_study_id(subject)
    output_path = bids_dir / bids_id / session / _get_output_path(modality)
    output_filename = (
        f"{bids_id}_{session}{_get_output_filename(modality, image_tracer)}"
    )
    output_path.mkdir(parents=True, exist_ok=True)

    # Check if old images exist and if updated mode is set to TRUE remove them
    if not _remove_existing_images_if_necessary(
        output_path,
        (_get_output_filename(modality, image_tracer), "magnitude", "phase"),
        mod_to_update,
    ):
        return None

    file_without_extension = output_path / output_filename
    output_image = file_without_extension.with_suffix(".nii.gz")

    if image.Is_Dicom:
        success = run_dcm2niix(
            input_dir=image_path
            if modality != ADNIModalityConverter.FMAP
            else Path(image_path).parent,
            output_dir=output_path,
            output_fmt=output_filename,
            compress=not _should_be_centered(modality),
            bids_sidecar=_write_json_sidecar(modality),
        )
        if not success:
            cprint(
                f"{logging_header} Error converting image {image_path} for subject {subject} and session {session}",
                lvl="warning",
            )

        # If "_t" - the trigger delay time - exists in dcm2niix output filename, we remove it
        for trigger_time in output_path.glob(f"{output_filename}_t[0-9]*"):
            res = re.search(r"_t\d+\.", str(trigger_time))
            no_trigger_time = str(trigger_time).replace(
                trigger_time[res.start() : res.end()], "."
            )
            os.rename(trigger_time, no_trigger_time)

        # Removing images with unsupported suffixes if generated by dcm2niix
        # todo : maybe later use output to remind user what failed
        _remove_files_with_unsupported_suffixes(
            output_path, ("ADC", "real", "imaginary")
        )

        # Conditions to check if output NIFTI files exists,
        # and, if DWI, if .bvec and .bval files are also present
        nifti_exists = (
            file_without_extension.with_suffix(".nii").is_file()
            or output_image.is_file()
        )
        if modality == ADNIModalityConverter.FMAP:
            nifti_exists = nifti_exists or any(
                [
                    file_without_extension.with_name(
                        f"{file_without_extension.name}_{suffix}.nii.gz"
                    ).is_file()
                    for suffix in {"ph", "e2", "e1", "e2_ph", "e1_ph"}
                ]
            )
        dwi_bvec_and_bval_exist = not (modality == ADNIModalityConverter.DWI) or (
            file_without_extension.with_suffix(".bvec").is_file()
            and file_without_extension.with_suffix(".bval").is_file()
        )

        # Check if conversion worked (output files exist)
        if not nifti_exists or not dwi_bvec_and_bval_exist:
            cprint(
                msg=f"{logging_header} Conversion with dcm2niix failed for subject {subject} and session {session}",
                lvl="warning",
            )
            return np.nan

        # Case when JSON file was expected, but not generated by dcm2niix
        json_exists = file_without_extension.with_suffix(".json").exists()

        if modality == ADNIModalityConverter.FMAP:
            json_exists = json_exists or any(
                [
                    file_without_extension.with_name(
                        f"{file_without_extension.name}_{suffix}.json"
                    ).is_file()
                    for suffix in {"ph", "e2", "e1", "e2_ph", "e1_ph"}
                ]
            )

        if _write_json_sidecar(modality) and not json_exists:
            cprint(
                msg=f"{logging_header} JSON file not generated by dcm2niix for subject {subject} and session {session}",
                lvl="warning",
            )

        if _should_be_centered(modality):
            output_image, error_msg = center_nifti_origin(
                file_without_extension.with_suffix(".nii"), output_image
            )
            if error_msg:
                cprint(msg=error_msg, lvl="error")
                raise ValueError(error_msg)
            file_without_extension.with_suffix(".nii").unlink()

    else:
        if _should_be_centered(modality):
            try:
                output_image, error_msg = center_nifti_origin(image_path, output_image)
            except Exception:
                cprint(
                    msg=f"No output image specified.",
                    lvl="error",
                )
                output_image = ""
                error_msg = False
            if error_msg:
                cprint(
                    msg=(
                        f"For subject {subject} in session {session}, "
                        f"an error occurred whilst recentering Nifti image: {image_path}"
                        f"The error is: {error_msg}"
                    ),
                    lvl="error",
                )
                raise ValueError(error_msg)
        else:
            shutil.copy(image_path, output_image)

    # Check if there is still the folder tmp_dcm_folder and remove it
    _remove_tmp_dmc_folder(bids_dir, image_id)
    return output_image


def _get_output_path(modality: ADNIModalityConverter) -> str:
    if modality == ADNIModalityConverter.T1:
        return "anat"
    if modality == ADNIModalityConverter.DWI:
        return "dwi"
    if modality == ADNIModalityConverter.FLAIR:
        return "anat"
    if modality == ADNIModalityConverter.FMRI:
        return "func"
    if modality == ADNIModalityConverter.FMAP:
        return "fmap"
    if modality in (
        ADNIModalityConverter.PET_FDG,
        ADNIModalityConverter.PET_FDG_UNIFORM,
        ADNIModalityConverter.PET_PIB,
        ADNIModalityConverter.PET_AV45,
        ADNIModalityConverter.PET_TAU,
    ):
        return "pet"


def _get_output_filename(
    modality: ADNIModalityConverter, tracer: Optional[Tracer] = None
) -> str:
    if modality == ADNIModalityConverter.T1:
        return "_T1w"
    if modality == ADNIModalityConverter.DWI:
        return "_dwi"
    if modality == ADNIModalityConverter.FLAIR:
        return "_FLAIR"
    if modality == ADNIModalityConverter.FMRI:
        return "_task-rest_bold"
    if modality == ADNIModalityConverter.FMAP:
        return "_fmap"
    if modality == ADNIModalityConverter.PET_FDG:
        return f"_trc-{Tracer.FDG.value}_rec-{ReconstructionMethod.CO_REGISTERED_AVERAGED.value}_pet"
    if modality == ADNIModalityConverter.PET_FDG_UNIFORM:
        return f"_trc-{Tracer.FDG.value}_rec-{ReconstructionMethod.COREGISTERED_ISOTROPIC.value}_pet"
    if modality == ADNIModalityConverter.PET_PIB:
        return f"_trc-{Tracer.PIB.value}_pet"
    if modality == ADNIModalityConverter.PET_AV45:
        return f"_trc-{tracer.value}_pet"
    if modality == ADNIModalityConverter.PET_TAU:
        return f"_trc-{Tracer.AV1451.value}_pet"


def _should_be_centered(modality: ADNIModalityConverter) -> bool:
    if modality == ADNIModalityConverter.T1:
        return True
    if modality == ADNIModalityConverter.DWI:
        return False
    if modality == ADNIModalityConverter.FLAIR:
        return True
    if modality == ADNIModalityConverter.FMRI:
        return False
    if modality == ADNIModalityConverter.FMAP:
        return False
    if modality in (
        ADNIModalityConverter.PET_FDG,
        ADNIModalityConverter.PET_FDG_UNIFORM,
        ADNIModalityConverter.PET_PIB,
        ADNIModalityConverter.PET_AV45,
        ADNIModalityConverter.PET_TAU,
    ):
        return True


def _write_json_sidecar(modality: ADNIModalityConverter) -> bool:
    if modality == ADNIModalityConverter.T1:
        return False
    if modality == ADNIModalityConverter.DWI:
        return True
    if modality == ADNIModalityConverter.FLAIR:
        return True
    if modality == ADNIModalityConverter.FMRI:
        return True
    if modality == ADNIModalityConverter.FMAP:
        return True
    if modality in (
        ADNIModalityConverter.PET_FDG,
        ADNIModalityConverter.PET_FDG_UNIFORM,
        ADNIModalityConverter.PET_PIB,
        ADNIModalityConverter.PET_AV45,
        ADNIModalityConverter.PET_TAU,
    ):
        return False


def _session_label_to_viscode(session_name: str) -> str:
    """Replace the session name passed as input with the session label 'bl' or 'mXXX'.

    Parameters
    ----------
    session_name: str
        Name of the session (MXXX).

    Returns
    -------
    str:
        'bl' if is the baseline session or the original session name.
    """
    if session_name == "M000":
        return "bl"
    return f"m{(int(session_name[1:])):02d}"


def _check_two_dcm_folder(dicom_path: Path, bids_folder: Path, image_uid: str) -> Path:
    """[summary].

    Check if a folder contains more than one DICOM and if yes, copy the DICOM related to
    image id passed as parameter into a temporary folder called tmp_dicom_folder.

    Parameters
    ----------
    dicom_path : Path
        The path to the DICOM folder.

    bids_folder : Path
        The path to the BIDS folder where the dataset will be stored.

    image_uid : str
        The image id of the fMRI.

    Returns
    -------
    Path :
        The path to the original DICOM folder or the path to a temporary folder called
        tmp_dicom_folder where only the DICOM to convert is copied.
    """
    import shutil
    from shutil import copy

    dest_path = bids_folder / f"tmp_dcm_folder_{image_uid.strip(' ')}"
    # Check if there are dicom files inside the folder not belonging to the desired image
    dicom_list = [f for f in dicom_path.glob("*.dcm")]
    image_list = [f for f in dicom_path.glob(f"*{image_uid}.dcm")]
    if len(dicom_list) != len(image_list):
        # Remove the precedent tmp_dcm_folder if present.
        if dest_path.exists():
            shutil.rmtree(dest_path)
        dest_path.mkdir(parents=True)
        for d in dicom_path.glob(f"*{str(image_uid)}.dcm*"):
            copy(d, dest_path)
        return dest_path
    return dicom_path


def _remove_tmp_dmc_folder(bids_dir: Path, image_id: str):
    """Remove the folder tmp_dmc_folder created by the method check_two_dcm_folder (if existing).

    Args:
        bids_dir: path to the BIDS directory
        image_id:
    """
    from shutil import rmtree

    tmp_dcm_folder_path = bids_dir / f"tmp_dcm_folder_{image_id.strip(' ')}"
    if tmp_dcm_folder_path.exists():
        rmtree(tmp_dcm_folder_path)


def load_clinical_csv(clinical_dir: Path, filename: str) -> pd.DataFrame:
    """Load the clinical csv from ADNI. This function is able to find the csv in the
    different known format available, the old format with just the name, and the new
    format with the name and the date of download.

    Parameters
    ----------
    clinical_dir: str
        Directory containing the csv.

    filename: str
        name of the file without the suffix.

    Returns
    -------
    pd.DataFrame:
        Dataframe corresponding to the filename.
    """
    import re

    pattern = filename + r"(_\d{1,2}[A-Za-z]{3}\d{4})?.csv"
    files_matching_pattern = [
        f for f in clinical_dir.rglob("*.csv") if re.search(pattern, f.name)
    ]
    if len(files_matching_pattern) != 1:
        raise IOError(
            f"Expecting to find exactly one file in folder {clinical_dir} "
            f"matching pattern {pattern}. {len(files_matching_pattern)} "
            f"files were found instead : \n{'- '.join(str(files_matching_pattern))}"
        )
    try:
        return pd.read_csv(files_matching_pattern[0], sep=",", low_memory=False)
    except Exception:
        raise ValueError(
            f"File {str(files_matching_pattern[0])} was found but could not "
            "be loaded as a DataFrame. Please check your data."
        )
