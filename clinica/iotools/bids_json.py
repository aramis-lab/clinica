import json
from pathlib import Path
from typing import Tuple, Union

import numpy as np
import pandas as pd
from pydicom.dataset import FileDataset
from pydicom.tag import Tag

from clinica.utils.stream import cprint

# todo : can this be used for other than AIBL ?
# todo : test
# todo : clarify class and use, what it will return...
# Currently done there for AIBL, PET (applied for all though)


def write_json(json_path: Path, json_dict: dict):
    # todo : maybe exists somewhere else in the code base
    # todo : date format ?
    with open(json_path, "w") as f:
        f.write(json.dumps(json_dict, indent=4))


class BidsCompliantJson:
    """For BIDS 1.10 compliance ; json from dcm"""

    def __init__(self):
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
                # todo : definitions ??
                ["TimeZero", ("AcquisitionTime",), "n/a"],
                ["ScanStart", (), "n/a"],  # same as acquisition time
                ["InjectionStart", ("RadiopharmaceuticalStartTime",), np.nan],
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
                ["ModeOfAdministration", (), "n/a"],
                ["InfusionRadioactivity", (), np.nan],  # todo : corresp ?
                ["InfusionStart", ("RadiopharmaceuticalStartTime",), np.nan],
                ["InfusionSpeed", (), np.nan],  # todo : corresp ?
                ["InfusionSpeedUnits", (), "n/a"],  # todo : corresp ?
                ["InjectedVolume", ("RadiopharmaceuticalVolume",), np.nan],
                # Reconstruction
                ["AcquisitionMode", (), "n/a"],
                ["ImageDecayCorrected", ("DecayCorrection",), "n/a"],
                ["ImageDecayCorrectionTime", ("DecayCorrection",), np.nan],
                ["ReconMethodName", ("ReconstructionMethod",), "n/a"],
                # todo : may use .VM on element to check for length
                ["ReconMethodParameterLabels", ("ReconstructionMethod",), [""]],
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

    def _from_image(self, dcm_dir: Path) -> dict:
        from pydicom import dcmread

        dict_json = {key: None for key in self.meta_collection.index}

        try:
            dcm_path = [f for f in dcm_dir.glob("*.dcm")][0]
            dicom_header = dcmread(dcm_path)
        except IndexError:
            cprint(
                msg=f"No DICOM found at {dcm_dir}, the image json will be filled with default values",
                lvl="warning",
            )
            # todo : kinda weird to be here and have no dicom though, does it ever happen ?
        else:
            # todo : use exception ?
            for index, metadata in self.meta_collection.iterrows():
                if dcm_value := self._get_dcm_tag(metadata.DCMtag, dicom_header):
                    dict_json.update({index: dcm_value})
        return dict_json

    @staticmethod
    def _set_scan_start(dcm_result: dict) -> dict:
        dcm_result["ScanStart"] = "0"
        return dcm_result

    @staticmethod
    def _set_injection_start(dcm_result: dict) -> dict:
        if (injection := dcm_result["InjectionStart"]) and (
            acquisition := dcm_result["TimeZero"]
        ):
            dcm_result["InjectionStart"] = injection - float(acquisition)
        return dcm_result

    @staticmethod
    def _get_admin_mode_from_start_time(dcm_result: dict) -> dict:
        # todo : test
        if (injection := dcm_result["InjectionStart"]) and (
            acquisition := dcm_result["TimeZero"]
        ):
            if injection - float(acquisition):
                dcm_result["ModeOfAdministration"] = "bolus-infusion"
            dcm_result["ModeOfAdministration"] = "bolus"
        return dcm_result

    @staticmethod
    def _check_decay_correction(dcm_result: dict) -> dict:
        # todo : test
        corrected = dcm_result["ImageDecayCorrected"]
        if corrected == "NONE":
            dcm_result["ImageDecayCorrected"] = False
        elif corrected:
            dcm_result["ImageDecayCorrected"] = True
        return dcm_result

    @staticmethod
    def _set_decay_time(dcm_result: dict) -> dict:
        correction = dcm_result["ImageDecayCorrectionTime"]
        if correction == "START":
            dcm_result["ImageDecayCorrectionTime"] = dcm_result["ScanStart"]
        elif correction == "ADMIN":
            dcm_result["ImageDecayCorrectionTime"] = dcm_result["InjectionStart"]
        return dcm_result

    @staticmethod
    def _update_default_units(dcm_result: dict) -> dict:
        dcm_result["InjectedRadioactivityUnits"] = "MBq"
        dcm_result["SpecificRadioactivityUnits"] = "Bq/micromole"
        return dcm_result

    @staticmethod
    def _update_injected_mass(dcm_result: dict) -> dict:
        # todo :test
        injected = dcm_result["InjectedRadioactivity"]
        specific = dcm_result["SpecificRadioactivity"]
        if injected and specific:
            dcm_result["InjectedMass"] = injected / float(specific)
            dcm_result["InjectedMassUnits"] = "mole"
        return dcm_result

    def _get_default_bids(self, dcm_result: dict) -> dict:
        for key in dcm_result:
            if not dcm_result[key]:
                dcm_result[key] = self.meta_collection.loc[key, "Value"]
        return dcm_result

    def build_dict(self, dcm_dir: Path) -> dict:
        dict_json = self._from_image(dcm_dir)

        dict_json = self._get_admin_mode_from_start_time(dict_json)
        dict_json = self._update_default_units(dict_json)
        dict_json = self._update_injected_mass(dict_json)
        dict_json = self._check_decay_correction(dict_json)
        dict_json = self._set_decay_time(dict_json)
        dict_json = self._set_scan_start(dict_json)
        dict_json = self._set_injection_start(dict_json)

        dict_json = self._get_default_bids(dict_json)
        return dict_json
