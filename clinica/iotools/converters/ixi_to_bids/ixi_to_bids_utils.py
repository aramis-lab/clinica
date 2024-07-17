# todo : verify typing compatibility with user input
from pathlib import Path

import pandas as pd


def get_subjects_list_from_data(data_directory: Path) -> List[str]:
    img_list = list(data_directory.rglob(pattern="*.nii.gz"))
    return


def filter_subjects_list(
    subjects_list: List[str], clinical_data: pd.DataFrame
) -> List[str]:
    return [subject for subject in subjects_list if subject in clinical_data["IXI_ID"]]


def read_ixi_clinical_data(clinical_data_path: Path) -> pd.DataFrame:
    clinical_data = pd.read_excel(clinical_data_path / "IXI.xls")
    clinical_data.dropna(inplace=True)
    clinical_data["IXI_ID"] = clinical_data.IXI_ID.apply(
        lambda x: "IXI" + "0" * (3 - len(str(x))) + str(x)
    )
    return clinical_data
