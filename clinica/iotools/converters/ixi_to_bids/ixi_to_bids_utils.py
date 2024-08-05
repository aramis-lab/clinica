# todo : verify typing compatibility with user input
import json
import re
import shutil
from pathlib import Path
from typing import List

import nibabel as nib
import pandas as pd
from nilearn.image import concat_imgs

from clinica.iotools.bids_utils import StudyName, bids_id_factory


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
        if subject in clinical_data["IXI_ID"].values
    ]


# todo : filter img data based on this


def read_ixi_clinical_data(clinical_data_path: Path) -> pd.DataFrame:
    clinical_data = pd.read_excel(clinical_data_path / "IXI.xls")
    clinical_data.dropna(inplace=True)
    clinical_data["IXI_ID"] = clinical_data.IXI_ID.apply(
        lambda x: "IXI" + "0" * (3 - len(str(x))) + str(x)
    )
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
    elif input_mod == "DWI":
        return "dwi"
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
    # todo : works for all except DTI !!! (I mean we should drop it after this)
    df = pd.DataFrame(
        {"img_path": [path for path in data_directory.rglob(pattern="IXI*.nii.gz")]}
    )
    # todo : regex pattern for 3 groups IXI\d{3}(-\w*){3}.nii.gz vs 4 groups IXI\d{3}(-\w*){4}.nii.gz
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
    # todo : has to work for both nifti and json
    return f"{img['bids_id']}_{img['session']}_{img['modality']}"


def write_subject(subject_df: pd.DataFrame, bids_path: Path) -> None:
    for _, row in subject_df.iterrows():
        filename = bids_filename_from_ixi(row)

        if row["modality"] in ["T1w", "T2w", "PDw", "angio"]:
            data_path = bids_path / row["bids_id"] / row["session"] / "anat"
        elif row["modality"] == "dwi":
            data_path = bids_path / row["bids_id"] / row["session"] / "func"
        else:
            raise ValueError(
                f"Modality {row['modality']} not recognized for IXI dataset."
            )

        data_path.mkdir(parents=True, exist_ok=True)
        shutil.copy2(row["img_path"], f"{data_path}/{filename}.nii.gz")

        with open(f"{data_path}/{filename}.json", "w") as f:
            json.dump(
                {
                    "InstitutionName": row["hospital"],
                    "MagneticFieldStrength (T)": row["field"],
                },
                f,
                indent=4,
            )


def merge_dti(dti_df: pd.DataFrame) -> pd.DataFrame:
    pass


# todo : write only one json per subject ?
# todo : test all
# todo : question DTI - sessions différentes ou images différentes ?
# todo : dataset descr? scans/sessions/...
