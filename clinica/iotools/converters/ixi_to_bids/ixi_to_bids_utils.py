# todo : verify typing compatibility with user input
import re
from pathlib import Path
from typing import List

import pandas as pd

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


def get_data_df(data_directory: Path) -> pd.DataFrame:
    # todo : works for all except DTI !!! (I mean we should drop it after this)
    df = pd.DataFrame(
        {"img_path": [path for path in data_directory.rglob(pattern="IXI*.nii.gz")]}
    )
    # todo : rename modalities inside
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
        .assign(session="ses-M000")
    )
    return df


def select_subject_data(data_df: pd.DataFrame, subject: str) -> pd.DataFrame:
    return data_df[data_df["subject"] == subject]


def bids_filename_from_ixi(img: pd.Series) -> str:
    # todo : adapter for modalities
    return f"{img['bids_id']}_{img['session']}"


def write_subject(subject_df: pd.DataFrame, bids_path: Path) -> None:
    # todo : place in subject folder, anat folder (T1,T2,PD) or dwi (dwi) folder if modality there ; where angio ??
    pass


# todo : function to place in files
# todo : function to create jsons
# todo : test all

# todo : question DTI - sessions différentes ou images différentes ?
