from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from pydicom.dataset import Dataset
from pydicom.multival import MultiValue

from clinica.converters.aibl_to_bids.utils.bids import Modality
from clinica.converters.aibl_to_bids.utils.json import (
    _format_time,
    _get_dicom_tags_and_defaults_base,
    _get_dicom_tags_and_defaults_pet,
    _update_metadata_from_image_dicoms,
)


@pytest.mark.parametrize(
    "input, expected", [("112233", "11:22:33"), ("112233.00", "11:22:33")]
)
def test_format_time_success(input, expected):
    assert _format_time(input) == expected


@pytest.mark.parametrize("input", ["foo", "112233,0000", "12", 112233])
def test_format_time_error(input):
    with pytest.raises((ValueError, AttributeError)):
        _format_time(input)


@pytest.mark.parametrize(
    "reference, event, expected",
    [("01:00:00", "01:00:01", 1), ("01:00:01", "01:00:00", -1)],
)
def test_substract_formatted_times(reference, event, expected):
    from clinica.converters.aibl_to_bids.utils.json import _substract_formatted_times

    assert expected == _substract_formatted_times(reference, event)


@pytest.mark.parametrize(
    "input, expected", [("112233", "11:22:33"), (112233, "n/a"), ("112233,00", "n/a")]
)
def test_format_timezero(input, expected):
    from clinica.converters.aibl_to_bids.utils.json import _format_timezero

    metadata = pd.DataFrame({"BIDSname": ["TimeZero"], "Value": [input]}).set_index(
        "BIDSname", drop=False
    )
    _format_timezero(metadata)
    assert expected == metadata.loc["TimeZero", "Value"]


@pytest.mark.parametrize(
    "input, expected", [("010001", 1), (112233, np.nan), ("112233,00", np.nan)]
)
def test_set_time_relative_to_zero(input, expected):
    from clinica.converters.aibl_to_bids.utils.json import _set_time_relative_to_zero

    metadata = pd.DataFrame(
        {"BIDSname": ["TimeZero", "Test"], "Value": ["01:00:00", input]}
    ).set_index("BIDSname", drop=False)
    _set_time_relative_to_zero(metadata, "Test")
    if not np.isnan(expected):
        assert expected == metadata.loc["Test", "Value"]
    else:
        assert np.isnan(metadata.loc["Test", "Value"])


@pytest.mark.parametrize(
    "input, expected",
    [("START", True), ("ADMIN", True), ("foo", False), ("NONE", False)],
)
def test_check_decay_correction(input, expected):
    from clinica.converters.aibl_to_bids.utils.json import _check_decay_correction

    metadata = pd.DataFrame(
        {"BIDSname": ["ImageDecayCorrected"], "Value": [input]}
    ).set_index("BIDSname", drop=False)
    _check_decay_correction(metadata)
    assert metadata.loc["ImageDecayCorrected", "Value"] == expected


@pytest.mark.parametrize(
    "input, expected",
    [("START", "scan"), ("ADMIN", "injection"), ("foo", None), ("NONE", None)],
)
def test_set_decay_time(input, expected):
    from clinica.converters.aibl_to_bids.utils.json import _set_decay_time

    metadata = pd.DataFrame(
        {
            "BIDSname": ["ImageDecayCorrectionTime", "ScanStart", "InjectionStart"],
            "Value": [input, "scan", "injection"],
        }
    ).set_index("BIDSname", drop=False)
    _set_decay_time(metadata)
    assert metadata.loc["ImageDecayCorrectionTime", "Value"] == expected


@pytest.mark.parametrize(
    "injected, specific, expected",
    [(5, "2", 2.5), (None, "2", "n/a"), (12, "n/a", "n/a")],
)
def test_update_injected_mass(injected, specific, expected):
    from clinica.converters.aibl_to_bids.utils.json import _update_injected_mass

    metadata = pd.DataFrame(
        {
            "BIDSname": [
                "InjectedRadioactivity",
                "SpecificRadioactivity",
                "InjectedMass",
            ],
            "Value": [injected, specific, "n/a"],
        }
    ).set_index("BIDSname", drop=False)
    _update_injected_mass(metadata)

    assert metadata.loc["InjectedMass", "Value"] == expected


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

    assert expected == _fetch_dcm_data_from_header(build_dicom_header(), keys)


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
        (_get_dicom_tags_and_defaults_base, 8, "BodyPart"),
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
        (Modality.AV45, 33),
        (Modality.FLUTE, 33),
        (Modality.PIB, 33),
        ("foo", 8),
    ],
)
def test_get_default_for_modality(modality, length):
    from clinica.converters.aibl_to_bids.utils.json import _get_dicom_tags_and_defaults

    assert len(_get_dicom_tags_and_defaults(modality)) == length


@patch("clinica.converters.aibl_to_bids.utils.json.cprint")
def test_update_metadata_from_image_dicoms_error(mock_cprint, tmp_path):
    _update_metadata_from_image_dicoms(pd.DataFrame(), tmp_path)
    mock_cprint.assert_called_once_with(
        msg=f"No DICOM found at {tmp_path}, the image json will be filled with default values",
        lvl="warning",
    )


def test_update_metadata_from_image_dicoms_success(tmp_path, mocker):
    metadata = pd.DataFrame(
        {"BIDSname": ["BodyPart"], "DCMtag": ["Body"], "Value": ["n/a"]}
    ).set_index("BIDSname", drop=False)

    mocker.patch("clinica.converters.aibl_to_bids.utils.json.next", return_value=[])
    mocker.patch(
        "clinica.converters.aibl_to_bids.utils.json.dcmread",
        return_value=build_dicom_header(),
    )
    mocker.patch(
        "clinica.converters.aibl_to_bids.utils.json._get_dcm_value_from_header",
        return_value=1,
    )
    _update_metadata_from_image_dicoms(metadata, tmp_path)
    assert metadata.loc["BodyPart", "Value"] == 1
