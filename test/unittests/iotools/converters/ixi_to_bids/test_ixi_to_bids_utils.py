from pathlib import Path

import pandas as pd
import pytest

from clinica.iotools.converters.ixi_to_bids.ixi_to_bids_utils import (
    _filter_subjects_list,
    _get_subjects_list_from_data,
    _rename_modalities,
    define_participants,
)


def test_get_subjects_list_from_data(tmp_path):
    for filename in ["IXI1", "IXI123", "IXIaaa", "foo"]:
        Path(f"{tmp_path}/{filename}_T1w.nii.gz").touch()
    assert _get_subjects_list_from_data(tmp_path) == ["IXI123"]


def test_filter_subjects_list():
    clinical_data = pd.DataFrame({"ixi_id": ["IXI123"]})
    subjects_list = ["IXI123", "IXI000", "IXI001"]
    assert _filter_subjects_list(
        subjects_list=subjects_list, clinical_data=clinical_data
    ) == ["IXI123"]


@pytest.mark.parametrize(
    "input_str, expected",
    [
        ("T1", "T1w"),
        ("T2", "T2w"),
        ("MRA", "angio"),
        ("PD", "PDw"),
        ("DTI", "dti"),
    ],
)
def test_rename_ixi_modalities_success(input_str, expected):
    assert _rename_modalities(input_str) == expected


@pytest.mark.parametrize("input_str", ["t1", "foo", "T1w"])
def test_rename_ixi_modalities_error(input_str):
    with pytest.raises(
        ValueError,
        match=f"The modality {input_str} is not recognized in the IXI dataset.",
    ):
        _rename_modalities(input_str)


def test_define_participants_from_file(tmp_path):
    pass
