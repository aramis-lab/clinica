"""
This module handles json writing from DICOM for PET modality in AIBL for BIDS 1.10 compliance
"""

from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
from pydicom.dataset import DataElement, FileDataset
from pydicom.multival import MultiValue

from clinica.utils.stream import cprint

# todo : __all__

# todo : can this be used for other than AIBL ?
# todo : write docstrings

# todo : would using a class be beneficial ? I am starting to get a lot of : if modality do something


def write_json(json_path: Path, json_data: pd.DataFrame):
    # todo : why n/a written as n><a ? or np.nan as null ?
    json_data["Value"].to_json(json_path, indent=4)


def _format_time(timestamp: str) -> str:
    # todo : test
    from time import strftime, strptime

    return strftime("%H:%M:%S", strptime(timestamp.split(".")[0], "%H%M%S"))


def _substract_formatted_times(reference: str, action: str) -> float:
    # todo : test
    from datetime import datetime

    return (
        datetime.strptime(action, "%H:%M:%S") - datetime.strptime(reference, "%H:%M:%S")
    ).total_seconds()


def _get_dicom_tags_and_defaults_pet() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["BIDSname", "DCMtag", "Value"],
        data=[
            # Hardware
            ["DeviceSerialNumber", ("DeviceSerialNumber",), "n/a"],
            ["ManufacturersModelName", ("ManufacturerModelName",), "n/a"],
            ["SoftwareVersions", ("SoftwareVersions",), "n/a"],
            ["BodyPart", ("BodyPartExamined",), "n/a"],
            ["MagneticFieldStrength", ("MagneticFieldStrength",), "n/a"],
            ["Units", ("Units",), "n/a"],
            # Institution
            ["InstitutionName", ("InstitutionName",), "n/a"],
            ["InstitutionAddress", ("InstitutionAddress",), "n/a"],
            [
                "InstitutionalDepartmentName",
                ("InstitutionalDepartmentName",),
                "n/a",
            ],
            # Time
            ["TimeZero", ("SeriesTime",), "n/a"],
            ["ScanStart", ("AcquisitionTime",), np.nan],
            [
                "InjectionStart",
                (
                    "RadiopharmaceuticalInformationSequence",
                    "RadiopharmaceuticalStartTime",
                ),
                np.nan,
            ],
            ["FrameTimesStart", ("FrameReferenceTime",), []],
            ["FrameDuration", ("ActualFrameDuration",), []],
            # Radiochemistry
            [
                "TracerName",
                (
                    "RadiopharmaceuticalInformationSequence",
                    "RadionuclideCodeSequence",
                    "MappingResourceName",
                ),
                "n/a",
            ],
            [
                "TracerRadionuclide",
                (
                    "RadiopharmaceuticalInformationSequence",
                    "RadionuclideCodeSequence",
                    "CodeMeaning",
                ),
                "n/a",
            ],
            [
                "InjectedRadioactivity",
                (
                    "RadiopharmaceuticalInformationSequence",
                    "RadionuclideTotalDose",
                ),
                np.nan,
            ],
            ["InjectedRadioactivityUnits", (), "MBq"],
            ["InjectedMass", (), "n/a"],
            ["InjectedMassUnits", (), "n/a"],
            [
                "SpecificRadioactivity",
                (
                    "RadiopharmaceuticalInformationSequence",
                    "RadiopharmaceuticalSpecificActivity",
                ),
                "n/a",
            ],
            ["SpecificRadioactivityUnits", (), "Bq/micromole"],
            [
                "ModeOfAdministration",
                (),
                "n/a",
            ],
            # Reconstruction
            [
                "AcquisitionMode",
                (),
                "n/a",
            ],  # todo : could not find / (0054, 1000) Series Type ?
            ["ImageDecayCorrected", ("DecayCorrection",), "n/a"],
            ["ImageDecayCorrectionTime", ("DecayCorrection",), np.nan],
            ["ReconMethodName", ("ReconstructionMethod",), "n/a"],
            [
                "ReconMethodParameterLabels",
                ("ReconstructionMethod",),
                [""],
            ],  # todo : could not find
            ["ReconMethodParameterUnits", ("ReconstructionMethod",), [""]],
            ["ReconMethodParameterValues", ("ReconstructionMethod",), [""]],
            ["ReconFilterType", ("ConvolutionKernel",), "n/a"],
            [
                "ReconFilterSize",
                ("ConvolutionKernel",),
                np.nan,
            ],  # todo : could not find
            ["AttenuationCorrection", ("AttenuationCorrectionMethod",), "n/a"],
        ],
    ).set_index(keys="BIDSname", drop=False)


def _get_dicom_tags_and_defaults_base() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["BIDSname", "DCMtag", "Value"],
        data=[
            # Hardware
            ["DeviceSerialNumber", ("DeviceSerialNumber",), "n/a"],
            ["ManufacturersModelName", ("ManufacturerModelName",), "n/a"],
            ["SoftwareVersions", ("SoftwareVersions",), "n/a"],
            ["BodyPart", ("BodyPartExamined",), "n/a"],
            ["MagneticFieldStrength", ("MagneticFieldStrength",), "n/a"],
            # Institution
            ["InstitutionName", ("InstitutionName",), "n/a"],
            ["InstitutionAddress", ("InstitutionAddress",), "n/a"],
            [
                "InstitutionalDepartmentName",
                ("InstitutionalDepartmentName",),
                "n/a",
            ],
        ],
    ).set_index(keys="BIDSname", drop=False)


def _check_dcm_value(
    dicom_data: Union[None, str, list, float],
) -> Union[None, str, list, float]:
    # Handles case where result is MultiValue, which is not JSON serializable and would be replaced by a blank space
    if dicom_data:
        if type(dicom_data) == MultiValue:
            return list(dicom_data)
        return dicom_data
    return None


def _fetch_dcm_data_from_header(
    keys_list: tuple[str], dicom_header: Union[FileDataset, DataElement]
) -> Optional[DataElement]:
    """
    Get the value from the dicom header corresponding to a dicom Tag/key list

    Parameters
    ----------
    keys_list :
        Tuple representing where to get the expected information in the dicom header

    dicom_header :
        Dicom header to parse

    Returns
    -------
        The value of the tag obtained from the header

    """
    if not keys_list:
        return None
    header = dicom_header.get(keys_list[0])
    if len(keys_list) == 1:
        return header
    if header:
        return _get_dcm_value_from_header(keys_list[1:], header[0])
    return None


def _get_dcm_value_from_header(
    keys_list: tuple[str], dicom_header: Union[FileDataset, DataElement]
) -> Union[None, str, list, float]:
    return _check_dcm_value(_fetch_dcm_data_from_header(keys_list, dicom_header))


def _update_metadata_from_image_dicoms(metadata: pd.DataFrame, dcm_dir: Path) -> None:
    from pydicom import dcmread

    try:
        dicom_header = dcmread(next(dcm_dir.rglob("*.dcm")))
    except StopIteration:
        cprint(
            msg=f"No DICOM found at {dcm_dir}, the image json will be filled with default values",
            lvl="warning",
        )
    else:
        # todo : prints in debug ?
        for index, data in metadata.iterrows():
            if dcm_value := _get_dcm_value_from_header(data.DCMtag, dicom_header):
                metadata.loc[index, "Value"] = dcm_value


def _format_timezero(dcm_result: pd.DataFrame) -> None:
    try:
        dcm_result.loc["TimeZero", "Value"] = _format_time(
            dcm_result.loc["TimeZero", "Value"]
        )
    except ValueError:
        return


def _set_time_relative_to_zero(dcm_result: pd.DataFrame, field: str) -> None:
    if (zero := dcm_result.loc["TimeZero", "Value"]) != "n/a":
        try:
            dcm_result.loc[field, "Value"] = _format_time(
                dcm_result.loc[field, "Value"]
            )
        except ValueError:
            return
        dcm_result.loc[field, "Value"] = _substract_formatted_times(
            reference=zero, action=dcm_result.loc[field, "Value"]
        )


def _set_time_from_ms_to_seconds(dcm_result: pd.DataFrame, field: str) -> None:
    if start := dcm_result.loc[field, "Value"]:
        dcm_result.loc[field, "Value"] = start / 1000


def _get_admin_injection_time_to_zero(dcm_result: pd.DataFrame) -> None:
    """The administration mode should be defined as bolus if < 10 min and infusion if > 30 min."""
    # todo : confirm ?
    if (injection := dcm_result.loc["InjectionStart", "Value"]) != "n/a":
        if abs(injection) / 60 < 10:
            dcm_result.loc["ModeOfAdministration", "Value"] = "bolus"
        elif abs(injection) / 60 > 30:
            dcm_result.loc["ModeOfAdministration", "Value"] = "infusion"
        else:
            dcm_result.loc["ModeOfAdministration", "Value"] = "bolus-infusion"


def _check_decay_correction(dcm_result: pd.DataFrame) -> None:
    corrected = dcm_result.loc["ImageDecayCorrected", "Value"]
    if corrected == "NONE":
        dcm_result.loc["ImageDecayCorrected", "Value"] = False
    elif corrected:
        dcm_result.loc["ImageDecayCorrected", "Value"] = True


def _set_decay_time(dcm_result: pd.DataFrame) -> None:
    correction = dcm_result.loc["ImageDecayCorrectionTime", "Value"]

    if correction == "START":
        dcm_result.loc["ImageDecayCorrectionTime", "Value"] = dcm_result.loc[
            "ScanStart", "Value"
        ]
    elif correction == "ADMIN":
        dcm_result.loc["ImageDecayCorrectionTime", "Value"] = dcm_result.loc[
            "InjectionStart", "Value"
        ]
    else:
        dcm_result.loc["ImageDecayCorrectionTime", "Value"] = None


def _update_injected_mass(dcm_result: pd.DataFrame) -> None:
    injected = dcm_result.loc["InjectedRadioactivity", "Value"]
    specific = dcm_result.loc["SpecificRadioactivity", "Value"]
    if injected and specific != "n/a":
        dcm_result.loc["InjectedMass", "Value"] = injected / float(specific)
        dcm_result.loc["InjectedMassUnits", "Value"] = "mole"


def _get_default_for_modality(modality: str) -> Optional[pd.DataFrame]:
    # todo : might need to adapt if other converters involved one day
    # todo : more restrictive on modality ? dataclass or something
    if modality in ("av45", "flute", "pib"):
        return _get_dicom_tags_and_defaults_pet()
    return _get_dicom_tags_and_defaults_base()


def _postprocess_for_pet(metadata: pd.DataFrame):
    _format_timezero(metadata)
    # Scan Start in BIDS is the difference between Time Zero and the beginning of image acquisition
    _set_time_relative_to_zero(metadata, field="ScanStart")
    # Injection Start in BIDS is the difference between Time Zero and the beginning of injection
    # RQ : can be negative
    _set_time_relative_to_zero(metadata, field="InjectionStart")
    _set_time_from_ms_to_seconds(metadata, "FrameTimesStart")
    _set_time_from_ms_to_seconds(metadata, "FrameDuration")
    _get_admin_injection_time_to_zero(metadata)

    _set_decay_time(metadata)
    _check_decay_correction(metadata)

    _update_injected_mass(metadata)


def get_json_data(dcm_dir: Path, modality: str) -> pd.DataFrame:
    metadata = _get_default_for_modality(modality)
    _update_metadata_from_image_dicoms(metadata, dcm_dir)
    if modality in ("av45", "flute", "pib"):
        _postprocess_for_pet(metadata)
    return metadata
