import json
from pathlib import Path
from turtledemo.penrose import start
from typing import Tuple, Union

import numpy as np
import pandas as pd
from pydicom.dataset import FileDataset
from pydicom.tag import Tag

from clinica.utils.stream import cprint

from .bids_utils import StudyName

# todo : can this be used for other than AIBL ?
# todo : test
# todo : should do specific construction depending on modality ?
# todo : clarify class and use, what it will return...
# Currently done there for AIBL, PET (applied for all though)


def write_json(json_path: Path, json_dict: dict):
    # todo : maybe exists somewhere else in the code base
    with open(json_path, "w") as f:
        f.write(json.dumps(json_dict, indent=4))


class BidsCompliantJson:
    """For BIDS 1.10 compliance ; json from dcm"""

    def __init__(self):
        # todo : finish filling default value and verify dicom correspondence
        # todo : there : optional + required fields for PET
        self.meta_collection = pd.DataFrame(
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
                # #todo : definitions ??
                ["TimeZero", ("AcquisitionTime",), "n/a"],
                ["ScanStart", ("AcquisitionTime",), "n/a"],
                ["InjectionStart", ("RadiopharmaceuticalStartTime",), np.nan],
                ["FrameTimesStart", ("FrameReferenceTime",), []],
                ["FrameDuration", ("ActualFrameDuration",), []],
                # Radiochemistry
                [
                    "TracerName",
                    (
                        "RadiopharmaceuticalInformationSequence",
                        "RadionuclideCodeSequence",
                        "MappingResource",
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
                ["ModeOfAdministration", (), "n/a"],
                ["InfusionRadioactivity", (), np.nan],  # todo : corresp ?
                ["InfusionStart", ("RadiopharmaceuticalStartTime",), np.nan],
                ["InfusionSpeed", (), np.nan],  # todo : corresp ?
                ["InfusionSpeedUnits", (), "n/a"],  # todo : corresp ?
                ["InjectedVolume", ("RadiopharmaceuticalVolume",), np.nan],
                # Reconstruction
                ["AcquisitionMode", (), "n/a"],  # todo : corresp ?
                # todo : DecayCorrection ?
                [
                    "ImageDecayCorrected",
                    ("DecayCorrection",),
                    np.nan,
                ],  # todo : if DecayCorrection 'NONE' then False
                ["ImageDecayCorrectionTime", (), np.nan],
                ["ReconMethodName", ("ReconstructionMethod",), "n/a"],
                [
                    "ReconMethodParameterLabels",
                    ("ReconstructionMethod",),
                    [""],
                ],  # todo : may use .VM on element to check for length
                ["ReconMethodParameterUnits", ("ReconstructionMethod",), [""]],
                ["ReconMethodParameterValues", ("ReconstructionMethod",), [""]],
                ["ReconFilterType", ("ConvolutionKernel",), "n/a"],
                ["ReconFilterSize", ("ConvolutionKernel",), np.nan],
                ["AttenuationCorrection", ("AttenuationCorrectionMethod",), "n/a"],
            ],
        )
        self.meta_collection.set_index(keys="BIDSname", drop=False, inplace=True)

    @staticmethod
    def _get_dcm_tag(
        keys_list: Tuple[str], dicom_header: FileDataset
    ) -> Union[str, int, list, None]:
        # todo : test
        # todo : for now don't call with exception to see how it works
        if not keys_list:
            return None
        if len(keys_list) == 1:
            get_result = dicom_header.get(Tag(keys_list[0]))
            return get_result.value if get_result else None
        dh = dicom_header.copy()
        try:
            for key in keys_list[:-1]:
                dh = dh.get(Tag(key))[0]
            get_result = dh.get(Tag(keys_list[-1]))
        except TypeError:
            get_result = None
        return get_result.value if get_result else None

    def _update_from_dcm(self, dcm_dir: Path) -> None:
        from pydicom import dcmread

        try:
            dcm_path = [f for f in dcm_dir.glob("*.dcm")][0]
            dicom_header = dcmread(dcm_path)
        except IndexError:
            cprint(
                msg=f"No DICOM found at {dcm_dir}, the image json will be filled with default values",
                lvl="warning",
            )
            # todo : kinda weird to be here and have no dicom though, does it ever happen ?
            return

        # todo : use exception ?
        for _, metadata in self.meta_collection.iterrows():
            if dcm_value := self._get_dcm_tag(metadata.DCMtag, dicom_header):
                self.meta_collection.loc[metadata.BIDSname, "Value"] = dcm_value

    def _get_admin_mode_from_study(self) -> None:
        # todo : test
        start_time = self.meta_collection.loc["TimeZero", "Value"]
        injection_time = self.meta_collection.loc["InjectionStart", "Value"]

        # todo : more check on injection_time ?? (np.isnan does not always work, mb write another function)

        if start_time != "n/a" and not bool(np.isnan(injection_time)):
            if injection_time - float(start_time):
                self.meta_collection.loc[
                    "ModeOfAdministration", "Value"
                ] = "bolus-infusion"
            self.meta_collection.loc["ModeOfAdministration", "Value"] = "bolus"

    def _set_decay_tag(self) -> None:
        corrected = self.meta_collection.loc["ImageDecayCorrected", "Value"]
        self.meta_collection.loc["ImageDecayCorrected", "Value"] = (
            not bool(np.isnan(corrected)) and corrected != "NONE"
        )
        # todo : more check on corrected

    def _update_default_units(self) -> None:
        self.meta_collection.loc["InjectedRadioactivityUnits", "Value"] = "MBq"
        self.meta_collection.loc["SpecificRadioactivityUnits", "Value"] = "Bq/micromole"

    def _update_injected_mass(self) -> None:
        # todo :test
        injected = self.meta_collection.loc["InjectedRadioactivity", "Value"]
        specific = self.meta_collection.loc["SpecificRadioactivity", "Value"]

        if not bool(np.isnan(injected)) and specific != "n/a":
            self.meta_collection.loc["InjectedMass", "Value"] = injected / float(
                specific
            )
            self.meta_collection.loc["InjectedMassUnits", "Value"] = "mole"

    def _translate_to_dict(self) -> dict:
        return {
            data.BIDSname: data.Value for _, data in self.meta_collection.iterrows()
        }

    def build_dict(self, dcm_dir: Path) -> dict:
        self._update_from_dcm(dcm_dir)
        self._get_admin_mode_from_study()
        self._update_default_units()
        self._update_injected_mass()
        breakpoint()
        return self._translate_to_dict()
