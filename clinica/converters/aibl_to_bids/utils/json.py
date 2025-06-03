"""
This module handles json writing from DICOM for PET modality in AIBL for BIDS 1.10 compliance
"""

import json
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np
import pandas as pd
from pydicom import dcmread
from pydicom.dataset import DataElement, FileDataset
from pydicom.multival import MultiValue

from clinica.utils.stream import cprint

from .bids import Modality

__all__ = ["write_json", "get_json_data"]


def write_json(json_path: Path, json_data: pd.DataFrame):
    with open(json_path, "w") as f:
        f.write(
            json.dumps(
                {row["BIDSname"]: row["Value"] for _, row in json_data.iterrows()},
                indent=4,
            )
        )


def _format_time(timestamp: str) -> str:
    """Format the timestamps that comes under the format HHMMSS or HHMMSS.0000 to HH:MM:SS."""
    from time import strftime, strptime

    return strftime("%H:%M:%S", strptime(timestamp.split(".")[0], "%H%M%S"))


def _substract_formatted_times(reference_time: str, event_time: str) -> float:
    from datetime import datetime

    return (
        datetime.strptime(event_time, "%H:%M:%S")
        - datetime.strptime(reference_time, "%H:%M:%S")
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
            ],  # todo : could not find
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
    dicom_data: Optional[Any],
) -> Union[None, str, list, float]:
    """Handles case where result is MultiValue, which is not JSON serializable and would be replaced by a blank space"""
    if dicom_data:
        if type(dicom_data) == MultiValue:
            return list(dicom_data)
        return dicom_data
    return None


def _fetch_dcm_data_from_header(
    dicom_header: Union[FileDataset, DataElement], dicom_tags: Optional[tuple[str, ...]]
) -> Optional[Any]:
    """
    Gets the value from the dicom header corresponding to a dicom Tag/key list

    Parameters
    ----------
    dicom_tags :
        Tuple representing where to get the expected information in the dicom header

    dicom_header :
        Dicom header to parse

    Returns
    -------
        The value of the tag obtained from the header

    """
    if not dicom_tags:
        return None
    header = dicom_header.get(dicom_tags[0])
    if len(dicom_tags) == 1:
        return header
    if header:
        return _fetch_dcm_data_from_header(header[0], dicom_tags[1:])
    return None


def _get_dcm_value_from_header(
    dicom_tags: tuple[str, ...], dicom_header: Union[FileDataset, DataElement]
) -> Union[None, str, list, float]:
    return _check_dcm_value(_fetch_dcm_data_from_header(dicom_header, dicom_tags))


def _update_metadata_from_image_dicoms(metadata: pd.DataFrame, dcm_dir: Path) -> None:
    try:
        dicom_header = dcmread(next(dcm_dir.rglob("*.dcm")))
    except StopIteration:
        cprint(
            msg=f"No DICOM found at {dcm_dir}, the image json will be filled with default values",
            lvl="warning",
        )
    else:
        for index, data in metadata.iterrows():
            if dcm_value := _get_dcm_value_from_header(data.DCMtag, dicom_header):
                cprint(msg=f"Value for {data.DCMtag} is : {dcm_value}", lvl="debug")
                metadata.loc[index, "Value"] = dcm_value
            else:
                cprint(
                    msg=f"Value for {data.DCMtag} could not be retrieved.", lvl="debug"
                )


def _format_timezero(dcm_result: pd.DataFrame) -> None:
    """Formats TimeZero if possible, else uses the default value."""
    try:
        dcm_result.loc["TimeZero", "Value"] = _format_time(
            dcm_result.loc["TimeZero", "Value"]
        )
    except (ValueError, AttributeError) as e:
        dcm_result.loc["TimeZero", "Value"] = "n/a"


def _set_time_relative_to_zero(dcm_result: pd.DataFrame, field: str) -> None:
    """Computes the difference from the time inside the field compared to TimeZero which should be formatted."""
    if (zero := dcm_result.loc["TimeZero", "Value"]) != "n/a":
        try:
            dcm_result.loc[field, "Value"] = _format_time(
                dcm_result.loc[field, "Value"]
            )
        except (ValueError, AttributeError) as e:
            dcm_result.loc[field, "Value"] = np.nan
        else:
            dcm_result.loc[field, "Value"] = _substract_formatted_times(
                reference_time=zero, event_time=dcm_result.loc[field, "Value"]
            )


def _set_time_from_ms_to_seconds(dcm_result: pd.DataFrame, field: str) -> None:
    """Set time to seconds for FrameTimeStart and FrameDuration which are encoded as list(float), given in ms."""
    if start := dcm_result.loc[field, "Value"]:
        dcm_result.loc[field, "Value"] = start / 1000


def _check_decay_correction(dcm_result: pd.DataFrame) -> None:
    dcm_result.loc["ImageDecayCorrected", "Value"] = dcm_result.loc[
        "ImageDecayCorrected", "Value"
    ] in ("START", "ADMIN")


def _set_decay_time(dcm_result: pd.DataFrame) -> None:
    correction = dcm_result.loc["ImageDecayCorrectionTime", "Value"]

    decay_time = None

    if correction == "START":
        decay_time = dcm_result.loc["ScanStart", "Value"]
    elif correction == "ADMIN":
        decay_time = dcm_result.loc["InjectionStart", "Value"]

    dcm_result.loc["ImageDecayCorrectionTime", "Value"] = decay_time


def _update_injected_mass(dcm_result: pd.DataFrame) -> None:
    injected = dcm_result.loc["InjectedRadioactivity", "Value"]
    specific = dcm_result.loc["SpecificRadioactivity", "Value"]
    if injected and specific != "n/a":
        dcm_result.loc["InjectedMass", "Value"] = injected / float(specific)
        dcm_result.loc["InjectedMassUnits", "Value"] = "mole"


def _get_dicom_tags_and_defaults(modality: Modality) -> pd.DataFrame:
    if modality in (Modality.AV45, Modality.FLUTE, Modality.PIB):
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

    _set_decay_time(metadata)
    _check_decay_correction(metadata)

    _update_injected_mass(metadata)


def get_json_data(dcm_dir: Path, modality: Modality) -> pd.DataFrame:
    metadata = _get_dicom_tags_and_defaults(modality)
    _update_metadata_from_image_dicoms(metadata, dcm_dir)
    if modality in (Modality.AV45, Modality.FLUTE, Modality.PIB):
        _postprocess_for_pet(metadata)
    return metadata
