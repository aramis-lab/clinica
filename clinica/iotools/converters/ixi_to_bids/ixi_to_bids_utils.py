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


def get_subjects_list_from_data(data_directory: Path) -> List[str]:
    return list(
        dict.fromkeys(
            re.match(r"IXI\d{3}", path.name).group(0)
            for path in data_directory.rglob(pattern="IXI*.nii.gz")
        )
    )


def filter_subjects_list(
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
    if subjs_list_path:
        cprint("Loading a subjects list provided by the user...")
        subjects_to_filter = subjs_list_path.read_text().splitlines()
    else:
        subjects_to_filter = get_subjects_list_from_data(data_directory)
    return filter_subjects_list(subjects_to_filter, clinical_data)


def read_ixi_clinical_data(clinical_data_path: Path) -> pd.DataFrame:
    clinical_data = pd.read_excel(clinical_data_path / "IXI.xls")
    clinical_data.dropna(inplace=True)
    clinical_data.drop("DATE_AVAILABLE", axis=1, inplace=True)
    clinical_data.rename(lambda x: x.lower(), axis=1, inplace=True)
    clinical_data["ixi_id"] = clinical_data.ixi_id.apply(
        lambda x: "IXI" + "0" * (3 - len(str(x))) + str(x)
    )
    clinical_data["session_id"] = "ses-M000"
    return clinical_data


def rename_ixi_modalities(input_mod: str) -> str:
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


def define_magnetic_field(hospital: str) -> str:
    if hospital == "Guys" or hospital == "IOP":
        return "1.5"
    if hospital == "HH":
        return "3"


def get_img_data_df(data_directory: Path) -> pd.DataFrame:
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
                lambda x: rename_ixi_modalities(x.split("-")[3])
            )
        )
        .assign(field=lambda df: df.hospital.apply(lambda x: define_magnetic_field(x)))
        .assign(session="ses-M000")
    )
    return df


def select_subject_data(data_df: pd.DataFrame, subject: str) -> pd.DataFrame:
    """Select subject images to copy based on IXI id"""
    return data_df[data_df["subject"] == subject]


def bids_filename_from_ixi(img: pd.Series) -> str:
    return f"{img['bids_id']}_{img['session']}_{img['modality']}"


def write_ixi_json_subject(writing_path: str, hospital: str, field: str) -> None:
    with open(writing_path, "w") as f:
        json.dump(
            {
                "InstitutionName": hospital,
                "MagneticFieldStrength (T)": field,
            },
            f,
            indent=4,
        )


def write_subject_no_dti(subject_df: pd.DataFrame, bids_path: Path) -> None:
    for _, row in subject_df.iterrows():
        cprint(
            f"Converting modality {row['modality']} for subject {row['subject']}.",
            lvl="debug",
        )
        filename = bids_filename_from_ixi(row)
        data_path = bids_path / row["bids_id"] / row["session"] / "anat"
        data_path.mkdir(parents=True, exist_ok=True)
        shutil.copy2(row["img_path"], f"{data_path}/{filename}.nii.gz")
        write_ixi_json_subject(
            f"{data_path}/{filename}.json", row["hospital"], row["field"]
        )


def write_subject_dti_if_exists(
    bids_path: Path, subject: str, data_directory: Path
) -> None:
    if dti_paths := find_subject_dti_data(data_directory, subject):
        cprint(f"Converting modality DTI for subject {subject}.", lvl="debug")
        dti_to_save = merge_dti(dti_paths)
        bids_id = bids_id_factory(StudyName.IXI).from_original_study_id(subject)
        data_path = bids_path / bids_id / "ses-M000" / "dwi"
        data_path.mkdir(parents=True, exist_ok=True)
        filename = f"{bids_id}_ses-M000_dwi"
        dti_to_save.to_filename(f"{data_path}/{filename}.nii.gz")
        hospital = dti_paths[0].name.split("-")[1]
        write_ixi_json_subject(
            f"{data_path}/{filename}.json", hospital, define_magnetic_field(hospital)
        )
    else:
        cprint(f"No DTI data was found for IXI subject {subject}.", lvl="warning")


def find_subject_dti_data(data_directory: Path, subject: str) -> List[Path]:
    pattern = subject + r"(-\w*){4}.nii.gz$"
    return [
        path
        for path in data_directory.rglob(pattern="IXI*.nii.gz")
        if re.search(pattern, str(path))
    ]


def merge_dti(dti_images: List[Path]):
    return concat_imgs([nib.load(img) for img in dti_images])


def write_ixi_scans(bids_dir: Path, participant: str):
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
        f"{bids_dir}/{bids_id}/ses-M000/{bids_id}_ses-M000_scans.tsv", sep="\t"
    )


def write_ixi_sessions(bids_dir: Path, clinical_data: pd.DataFrame, participant: str):
    line = clinical_data[clinical_data["ixi_id"] == participant]
    bids_id = bids_id_factory(StudyName.IXI).from_original_study_id(participant)
    line.to_csv(f"{bids_dir}/{bids_id}/{bids_id}_sessions.tsv", sep="\t")
