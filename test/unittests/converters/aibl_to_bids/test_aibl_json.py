import pytest
from pydicom.dataset import DataElement, Dataset


def build_dicom_header():
    header = Dataset()
    header.add_new("Time", "TM", "00:00")
    return header


@pytest.mark.parametrize("keys, expected", [((), None), (("Time",), "00:00")])
def test_get_dcm_value_from_header(keys, expected):
    from clinica.converters.aibl_to_bids.utils.json import _get_dcm_value_from_header

    assert expected == _get_dcm_value_from_header(keys, build_dicom_header())
