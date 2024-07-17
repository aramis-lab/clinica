# todo : verify typing compatibility with user input
import re
from pathlib import Path
from typing import List

import pandas as pd


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


def get_data_df(data_directory: Path) -> pd.DataFrame:
    # todo : works for all except DTI !!! (I mean we should drop it after this)
    df = pd.DataFrame(
        {"img_path": [path for path in data_directory.rglob(pattern="IXI*.nii.gz")]}
    )
    df = (
        df.assign(img_name=lambda df: df.img_path.apply(lambda x: x.name))
        .assign(img_name_no_ext=lambda df: df.img_name.apply(lambda x: x.split(".")[0]))
        .assign(subject=lambda df: df.img_name_no_ext.apply(lambda x: x.split("-")[0]))
        .assign(hospital=lambda df: df.img_name_no_ext.apply(lambda x: x.split("-")[1]))
        .assign(modality=lambda df: df.img_name_no_ext.apply(lambda x: x.split("-")[3]))
        .assign(session="ses-M000")
    )
    return df


# todo : question DTI - sessions différentes ou images différentes ?
