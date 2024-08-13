from pathlib import Path
from typing import List

import pandas as pd
import pytest

from clinica.iotools.converters.ixi_to_bids.ixi_to_bids_utils import (
    filter_subjects_list,
    get_subjects_list_from_data,
    rename_ixi_modalities,
)


def test_get_subjects_list_from_data(tmp_path):
    for filename in ["IXI1", "IXI123", "IXIaaa", "foo"]:
        Path(f"{tmp_path}/{filename}_T1w.nii.gz").touch()
    assert get_subjects_list_from_data(tmp_path) == ["IXI123"]


def test_filter_subjects_list():
    clinical_data = pd.DataFrame({"ixi_id": ["IXI123"]})
    subjects_list = ["IXI123", "IXI000", "IXI001"]
    assert filter_subjects_list(
        subjects_list=subjects_list, clinical_data=clinical_data
    ) == ["IXI123"]


@pytest.mark.parametrize(
    "input, expected",
    [
        ("T1", "T1w"),
        ("T2", "T2w"),
        ("MRA", "angio"),
        ("PD", "PDw"),
        ("DTI", "dti"),
    ],
)
def test_rename_ixi_modalities_success(input, expected):
    assert rename_ixi_modalities(input) == expected


@pytest.mark.parametrize("input", ["t1", "foo", "T1w"])
def test_rename_ixi_modalities_error(input):
    with pytest.raises(
        ValueError, match=f"The modality {input} is not recognized in the IXI dataset."
    ):
        rename_ixi_modalities(input)
