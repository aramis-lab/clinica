from unittest.mock import patch

import pandas as pd
import pytest
from pydicom.dataset import Dataset
from pydicom.multival import MultiValue

from clinica.converters.aibl_to_bids.utils.json import (
    _get_dicom_tags_and_defaults_base,
    _get_dicom_tags_and_defaults_pet,
    _update_metadata_from_image_dicoms,
)


def build_dicom_header():
    header = Dataset()
    header.add_new("Time", "TM", "093015")
    header.add_new("SoftwareVersions", "LO", ["clinica", "0.9.4"])
    header.BeamSequence = [Dataset()]
    header.BeamSequence[0].Manufacturer = "foo"
    return header


@pytest.mark.parametrize(
    "keys, expected",
    [
        ((), None),
        (("Time",), "093015"),
        (("BeamSequence", "Manufacturer"), "foo"),
    ],
)
def test_fetch_dcm_data_from_header(keys, expected):
    from clinica.converters.aibl_to_bids.utils.json import _fetch_dcm_data_from_header

    assert expected == _fetch_dcm_data_from_header(keys, build_dicom_header())


def test_get_dcm_value_from_header_multivalue():
    from clinica.converters.aibl_to_bids.utils.json import _get_dcm_value_from_header

    assert ["clinica", "0.9.4"] == _get_dcm_value_from_header(
        ("SoftwareVersions",), build_dicom_header()
    )


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
        (_get_dicom_tags_and_defaults_pet, 33, "AttenuationCorrection"),
        (_get_dicom_tags_and_defaults_base, 8, "Units"),
    ],
)
def test_get_dicom_tags_and_defaults_pet(get_default, length, tag):
    defaults = get_default()
    assert len(defaults) == length
    assert set(defaults.columns) == {"BIDSname", "DCMtag", "Value"}
    assert tag in defaults["BIDSname"]


@pytest.mark.parametrize(
    "modality, length",
    [
        ("av45", 33),
        ("flute", 33),
        ("pib", 33),
        ("foo", 8),
    ],
)
def test_get_default_for_modality(modality, length):
    from clinica.converters.aibl_to_bids.utils.json import _get_default_for_modality

    assert len(_get_default_for_modality(modality)) == length


@patch("clinica.converters.aibl_to_bids.utils.json.cprint")
def test_update_metadata_from_image_dicoms_error(mock_cprint, tmp_path):
    _update_metadata_from_image_dicoms(pd.DataFrame(), tmp_path)
    mock_cprint.assert_called_once_with(
        msg=f"No DICOM found at {tmp_path}, the image json will be filled with default values",
        lvl="warning",
    )


# todo : success : need to mock header + defaults
