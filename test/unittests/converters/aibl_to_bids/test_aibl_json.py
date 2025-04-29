import pytest
from pydicom.dataset import DataElement, Dataset
from pydicom.multival import MultiValue

from clinica.converters.aibl_to_bids.utils.json import (
    _get_dicom_tags_and_defaults_base,
    _get_dicom_tags_and_defaults_pet,
)


def build_dicom_header():
    header = Dataset()
    header.add_new("Time", "TM", "093015")
    header.BeamSequence = [Dataset()]
    header.BeamSequence[0].Manufacturer = "Linac, co."
    return header


@pytest.mark.parametrize(
    "keys, expected",
    [
        ((), None),
        (("Time",), "093015"),
        (("BeamSequence", "Manufacturer"), "Linac, co."),
    ],
)
def test_fetch_dcm_data_from_header(keys, expected):
    from clinica.converters.aibl_to_bids.utils.json import _fetch_dcm_data_from_header

    assert expected == _fetch_dcm_data_from_header(keys, build_dicom_header())


@pytest.mark.parametrize(
    "element, expected",
    [
        (None, None),
        ("093015", "093015"),
        (2, 2),
        (MultiValue(int, [1, 2]), [1, 2]),
    ],
)
def test_check_dcm_value(element, expected):
    from clinica.converters.aibl_to_bids.utils.json import _check_dcm_value

    assert expected == _check_dcm_value(element)


@pytest.mark.parametrize(
    "get_default, length, tag",
    [
        (_get_dicom_tags_and_defaults_pet, 38, "AttenuationCorrection"),
        (_get_dicom_tags_and_defaults_base, 9, "Units"),
    ],
)
def test_get_dicom_tags_and_defaults_pet(get_default, length, tag):
    defaults = get_default()
    assert len(defaults) == length
    assert set(defaults.columns) == {"BIDSname", "DCMtag", "Value"}
    assert tag in defaults["BIDSname"]
