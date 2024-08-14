import json
import re
import shutil
from pathlib import Path
from typing import List, Optional

import nibabel as nib
import pandas as pd
from nilearn.image import concat_imgs

from clinica.iotools.bids_utils import StudyName, bids_id_factory
from clinica.utils.stream import cprint

__all__ = [
    "read_clinical_data",
    "define_participants",
    "write_ixi_subject_data",
    "write_sessions",
    "write_scans",
]


def _get_subjects_list_from_data(data_directory: Path) -> List[str]:
    return list(
        dict.fromkeys(
            re.match(r"IXI\d{3}", path.name).group(0)
            for path in data_directory.rglob(pattern="IXI*.nii.gz")
            if re.match(r"IXI\d{3}", path.name)
        )
    )


def _filter_subjects_list(
    subjects_list: List[str], clinical_data: pd.DataFrame
) -> List[str]:
    return [
        subject
        for subject in subjects_list
        if subject in clinical_data["ixi_id"].values
    ]


def define_participants(
    data_directory: Path,
    clinical_data: pd.DataFrame,
    subjs_list_path: Optional[Path] = None,
) -> List[str]:
    """
    Defines the actual list of participants based on the (provided) subjects list filtered
    using clinical data.

    Parameters
    ----------
    data_directory : Path to the raw data directory.
    clinical_data : Dataframe containing the formatted clinical data of the IXI study.
    subjs_list_path : [Optional] Path to a text file containing a list of specific subjects to extract.

    Returns
    -------
    The list of ids that were either present in the data directory or asked for with a specific text file, provided
    associated clinical data actually exists.

    """
    if subjs_list_path:
        cprint("Loading a subjects list provided by the user...")
        subjects_to_filter = subjs_list_path.read_text().splitlines()
    else:
        subjects_to_filter = _get_subjects_list_from_data(data_directory)
    return _filter_subjects_list(subjects_to_filter, clinical_data)


def read_clinical_data(clinical_data_path: Path) -> pd.DataFrame:
    """
    Reads and formats IXI clinical data.

    Parameters
    ----------
    clinical_data_path : Path to the directory where the clinical data .xls is stored.

    Returns
    -------
    A dataframe containing the clinical data, with some modifications: padded study id and added session id.

    """
    # todo : what keys ?
    clinical_data = pd.read_excel(clinical_data_path / "IXI.xls")
    clinical_data.dropna(inplace=True)
    clinical_data.drop("DATE_AVAILABLE", axis=1, inplace=True)
    clinical_data.rename(lambda x: x.lower(), axis=1, inplace=True)
    clinical_data["ixi_id"] = clinical_data.ixi_id.apply(
        lambda x: "IXI" + "0" * (3 - len(str(x))) + str(x)
    )
    clinical_data["session_id"] = "ses-M000"
    return clinical_data


def _rename_modalities(input_mod: str) -> str:
    if input_mod == "T1":
        return "T1w"
    elif input_mod == "T2":
        return "T2w"
    elif input_mod == "MRA":
        return "angio"
    elif input_mod == "PD":
        return "PDw"
    elif input_mod == "DTI":
        return "dti"
    else:
        raise ValueError(
            f"The modality {input_mod} is not recognized in the IXI dataset."
        )


def _define_magnetic_field(hospital: str) -> str:
    # todo : raise error ?
    if hospital == "Guys" or hospital == "IOP":
        return "1.5"
    if hospital == "HH":
        return "3"


def _get_img_data(data_directory: Path) -> pd.DataFrame:
    """Finds paths for all images that are not DTI data and processes the info contained in their names"""
    df = pd.DataFrame(
        {
            "img_path": [
                path
                for path in data_directory.rglob(pattern="IXI*.nii.gz")
                if re.search(r"IXI\d{3}(-\w*){3}.nii.gz$", str(path))
            ]
        }
    )
    df = (
        df.assign(img_name=lambda df: df.img_path.apply(lambda x: x.name))
        .assign(img_name_no_ext=lambda df: df.img_name.apply(lambda x: x.split(".")[0]))
        .assign(subject=lambda df: df.img_name_no_ext.apply(lambda x: x.split("-")[0]))
        .assign(
            bids_id=lambda df: df.subject.apply(
                lambda x: bids_id_factory(StudyName.IXI).from_original_study_id(x)
            )
        )
        .assign(hospital=lambda df: df.img_name_no_ext.apply(lambda x: x.split("-")[1]))
        .assign(
            modality=lambda df: df.img_name_no_ext.apply(
                lambda x: _rename_modalities(x.split("-")[3])
            )
        )
        .assign(field=lambda df: df.hospital.apply(lambda x: _define_magnetic_field(x)))
        .assign(session="ses-M000")
    )
    return df


def _select_subject_data(data_df: pd.DataFrame, subject: str) -> pd.DataFrame:
    """Select subject images to copy based on IXI id"""
    return data_df[data_df["subject"] == subject]


def _get_bids_filename_from_image_data(img: pd.Series) -> str:
    return f"{img['bids_id']}_{img['session']}_{img['modality']}"


def write_ixi_subject_data(
    bids_dir: Path, participant: str, path_to_dataset: Path
) -> None:
    """
    Writes the data of the IXI subject in the BIDS directory following BIDS specifications.

    Parameters
    ----------
    bids_dir : Path to the output BIDS directory.
    participant : Current converted subject study id (str).
    path_to_dataset : Path to the raw dataset directory.
    """
    _write_subject_no_dti(
        _select_subject_data(_get_img_data(path_to_dataset), participant), bids_dir
    )
    _write_subject_dti_if_exists(bids_dir, participant, path_to_dataset)


def _write_ixi_json_image(writing_path: Path, hospital: str, field: str) -> None:
    """
    Writes a json associated to one IXI image.

    Parameters
    ----------
    writing_path : Path indicating under what name to write the json.
    hospital : identifier for the hospital responsible for the acquisition, determined from the image filename.
    field :  magnetic field used for the acquisition of the image, determined from the hospital.
    """
    with open(writing_path, "w") as f:
        json.dump(
            {
                "InstitutionName": hospital,
                "MagneticFieldStrength (T)": field,
            },
            f,
            indent=4,
        )


def _write_subject_no_dti(subject_df: pd.DataFrame, bids_path: Path) -> None:
    """Copies all subject data but DTI"""
    for _, row in subject_df.iterrows():
        cprint(
            f"Converting modality {row['modality']} for subject {row['subject']}.",
            lvl="debug",
        )
        filename = _get_bids_filename_from_image_data(row)
        data_path = bids_path / row["bids_id"] / row["session"] / "anat"
        data_path.mkdir(parents=True, exist_ok=True)
        shutil.copy2(row["img_path"], f"{data_path}/{filename}.nii.gz")
        _write_ixi_json_image(
            (data_path / filename).with_suffix(".json"), row["hospital"], row["field"]
        )
        # todo : 1 json per image or 1 json per subject ?


def _write_subject_dti_if_exists(
    bids_path: Path, subject: str, data_directory: Path
) -> None:
    """Processes DTI data if found for a subject"""
    if dti_paths := _find_subject_dti_data(data_directory, subject):
        cprint(f"Converting modality DTI for subject {subject}.", lvl="debug")
        dti_to_save = _merge_dti(dti_paths)
        bids_id = bids_id_factory(StudyName.IXI).from_original_study_id(subject)
        data_path = bids_path / bids_id / "ses-M000" / "dwi"
        data_path.mkdir(parents=True, exist_ok=True)
        filename = Path(f"{bids_id}_ses-M000_dwi")
        dti_to_save.to_filename(data_path / filename.with_suffix(".nii.gz"))
        hospital = dti_paths[0].name.split("-")[1]
        _write_ixi_json_image(
            data_path / filename.with_suffix(".json"),
            hospital,
            _define_magnetic_field(hospital),
        )
    else:
        cprint(f"No DTI data was found for IXI subject {subject}.", lvl="warning")


def _find_subject_dti_data(data_directory: Path, subject: str) -> List[Path]:
    pattern = subject + r"(-\w*){4}.nii.gz$"
    return [
        path
        for path in data_directory.rglob(pattern="IXI*.nii.gz")
        if re.search(pattern, str(path))
    ]


def _merge_dti(dti_images: List[Path]):
    return concat_imgs([nib.load(img) for img in dti_images])


def write_scans(bids_dir: Path, participant: str) -> None:
    """
    Write the scans.tsv for the only session (ses-M000) of a IXI subject.

    Parameters
    ----------
    bids_dir : Path to the output BIDS directory.
    participant : Current converted subject study id (str).
    """
    bids_id = bids_id_factory(StudyName.IXI).from_original_study_id(participant)
    to_write = pd.DataFrame(
        {
            "filename": [
                f"{path.parent.name}/{path.name}"
                for path in (bids_dir / bids_id).rglob("*.nii.gz")
            ]
        }
    )
    to_write.to_csv(
        bids_dir / bids_id / f"ses-M000/{bids_id}_ses-M000_scans.tsv", sep="\t"
    )


def write_sessions(
    bids_dir: Path, clinical_data: pd.DataFrame, participant: str
) -> None:
    """
    Writes the sessions.tsv for a IXI subject.

    Parameters
    ----------
    bids_dir : Path to the output BIDS directory.
    clinical_data : Dataframe containing the formatted clinical data of the IXI study.
    participant : Current converted subject study id (str).
    """
    line = clinical_data[clinical_data["ixi_id"] == participant]
    bids_id = bids_id_factory(StudyName.IXI).from_original_study_id(participant)
    line.to_csv(bids_dir / bids_id / f"{bids_id}_sessions.tsv", sep="\t")
