"""
This module handles json writing from DICOM for PET modality in AIBL for BIDS 1.10 compliance
"""

from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
from pydicom.dataset import DataElement, FileDataset
from pydicom.multival import MultiValue
from pydicom.tag import Tag

from clinica.utils.stream import cprint

# todo : __all__

# todo : can this be used for other than AIBL ?
# todo : adapt to modality, rn computed for all though tags are for PET
# todo : write docstrings


# todo : would using a class be beneficial ? I am starting to get a lot of if modality truc...


def write_json(json_path: Path, json_data: pd.DataFrame):
    # todo : date format ?
    # todo : why n/a written as n><a ? or np.nan as null ?
    json_data["Value"].to_json(json_path, indent=4)


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
            # todo : definitions ?? Study Time != Series time != Content Time != Acquisition Time (but dates same)
            ["TimeZero", ("AcquisitionTime",), "n/a"],
            ["ScanStart", (), "n/a"],  # same as acquisition time
            [
                "InjectionStart",
                ("RadiopharmaceuticalStartTime",),
                np.nan,
            ],  # diff between "RadiopharmaceuticalStartTime" and acquisition time
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
            ["InjectedRadioactivity", ("RadionuclideTotalDose",), np.nan],
            ["InjectedRadioactivityUnits", (), "n/a"],
            ["InjectedMass", (), "n/a"],
            ["InjectedMassUnits", (), "n/a"],
            [
                "SpecificRadioactivity",
                ("RadiopharmaceuticalSpecificActivity",),
                "n/a",
            ],
            ["SpecificRadioactivityUnits", (), "n/a"],
            [
                "ModeOfAdministration",
                ("RadiopharmaceuticalStartTime",),
                "n/a",
            ],  # uses "RadiopharmaceuticalStartTime" to decide on admin mode
            # todo :  and (0054, 1002) Counts Source
            ["InfusionRadioactivity", (), np.nan],  # todo : corresp ?
            ["InfusionStart", ("RadiopharmaceuticalStartTime",), np.nan],
            ["InfusionSpeed", (), np.nan],  # todo : corresp ?
            ["InfusionSpeedUnits", (), "n/a"],  # todo : corresp ?
            ["InjectedVolume", ("RadiopharmaceuticalVolume",), np.nan],
            # Reconstruction
            # todo : see(0054, 1000) Series Type ; (0054, 1001) Units ; (0028, 0051) Corrected Image
            # see (0028, 1052) Rescale Intercept ; (0028, 1053) Rescale Slope
            # (0018, 1100) Reconstruction Diameter ; (0018, 1149) Field of View Dimension(s) ; (0018, 1147) Field of View Shape
            ["AcquisitionMode", (), "n/a"],
            ["ImageDecayCorrected", ("DecayCorrection",), "n/a"],
            ["ImageDecayCorrectionTime", ("DecayCorrection",), np.nan],
            ["ReconMethodName", ("ReconstructionMethod",), "n/a"],
            [
                "ReconMethodParameterLabels",
                ("ReconstructionMethod",),
                [""],
            ],  # todo : corresp more precise ?
            ["ReconMethodParameterUnits", ("ReconstructionMethod",), [""]],
            ["ReconMethodParameterValues", ("ReconstructionMethod",), [""]],
            ["ReconFilterType", ("ConvolutionKernel",), "n/a"],
            ["ReconFilterSize", ("ConvolutionKernel",), np.nan],  # todo : corresp ?
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
            ["Units", ("Units",), "n/a"],
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
        _get_dcm_value_from_header(keys_list[1:], header[0])
    return None


def _get_dcm_value_from_header(
    keys_list: tuple[str], dicom_header: Union[FileDataset, DataElement]
) -> Union[None, str, list, float]:
    test = _fetch_dcm_data_from_header(keys_list, dicom_header)
    return _check_dcm_value(test)


def _update_metadata_from_image_dicoms(metadata: pd.DataFrame, dcm_dir: Path) -> None:
    from pydicom import dcmread

    try:
        dicom_header = dcmread(next(dcm_dir.rglob("*.dcm")))
    except IndexError:
        cprint(
            msg=f"No DICOM found at {dcm_dir}, the image json will be filled with default values",
            lvl="warning",
        )
    else:
        # todo : prints in debug ?
        for index, data in metadata.iterrows():
            if dcm_value := _get_dcm_value_from_header(data.DCMtag, dicom_header):
                metadata.loc[index, "Value"] = dcm_value


def _set_scan_start(dcm_result: pd.DataFrame) -> None:
    dcm_result.loc["ScanStart", "Value"] = "0"


def _set_injection_start(dcm_result: pd.DataFrame) -> None:
    if (injection := dcm_result.loc["InjectionStart", "Value"]) and (
        acquisition := dcm_result.loc["TimeZero", "Value"]
    ):
        dcm_result.loc["InjectionStart", "Value"] = injection - float(acquisition)


def _get_admin_mode_from_start_time(dcm_result: pd.DataFrame) -> None:
    if (injection := dcm_result.loc["ModeOfAdministration", "Value"]) and (
        acquisition := dcm_result.loc["TimeZero", "Value"]
    ):
        if injection - float(acquisition):
            dcm_result.loc["ModeOfAdministration", "Value"] = "bolus-infusion"
        dcm_result.loc["ModeOfAdministration", "Value"] = "bolus"


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


def _update_default_units(dcm_result: pd.DataFrame) -> None:
    dcm_result.loc["InjectedRadioactivityUnits", "Value"] = "MBq"
    dcm_result.loc["SpecificRadioactivityUnits", "Value"] = "Bq/micromole"


def _update_injected_mass(dcm_result: pd.DataFrame) -> None:
    injected = dcm_result.loc["InjectedRadioactivity", "Value"]
    specific = dcm_result.loc["SpecificRadioactivity", "Value"]
    if injected and specific:
        dcm_result.loc["InjectedMass", "Value"] = injected / float(specific)
        dcm_result.loc["InjectedMassUnits", "Value"] = "mole"


def _get_default_for_modality(modality: str) -> Optional[pd.DataFrame]:
    # todo : might need to adapt if other converters involved one day
    if modality in ("av45", "flute", "pib"):
        return _get_dicom_tags_and_defaults_pet()
    return _get_dicom_tags_and_defaults_base()


def _postprocess_for_pet(metadata: pd.DataFrame):
    _set_decay_time(metadata)
    _set_scan_start(metadata)

    # _get_admin_mode_from_start_time(df)
    # _update_default_units(df)
    # _update_injected_mass(df)
    # _check_decay_correction(df)
    # _set_injection_start(df)


def get_json_data(dcm_dir: Path, modality: str) -> pd.DataFrame:
    metadata = _get_default_for_modality(modality)
    _update_metadata_from_image_dicoms(metadata, dcm_dir)
    if modality in ("av45", "flute", "pib"):
        _postprocess_for_pet(metadata)
    return metadata
