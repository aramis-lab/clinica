import json
from pathlib import Path
from types import NoneType
from typing import Iterable, Tuple, Union

import numpy as np
import pandas as pd
from pydicom.dataset import FileDataset
from pydicom.tag import Tag

from clinica.utils.stream import cprint

from .bids_utils import StudyName

# todo : can this be used for other than AIBL ?
# todo : test
# todo : should do specific construction depending on modality ?


class BidsCompliantJson:
    """For BIDS 1.10 compliance ; json from dcm"""

    def __init__(self):
        # todo : finish filling default value and verify dicom correspondence
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
                ["InjectedRadioactivityUnits", (), "MBq"],
                ["InjectedMass", (), "n/a"],
                ["InjectedMassUnits", (), "n/a"],
                [
                    "SpecificRadioactivity",
                    ("RadiopharmaceuticalSpecificActivity",),
                    "n/a",
                ],
                ["SpecificRadioactivityUnits", (), "Bq/micromole"],
                ["ModeOfAdministration", (), "n/a"],  # todo : corresp ?
                # todo : required if bolus-infusion (below)
                ["InfusionRadioactivity", (), np.nan],  # todo : corresp ?
                ["InfusionStart", ("RadiopharmaceuticalStartTime",), np.nan],
                ["InfusionSpeed", (), np.nan],  # todo : corresp ?
                ["InfusionSpeedUnits", (), "n/a"],  # todo : corresp ?
                ["InjectedVolume", ("RadiopharmaceuticalVolume",), np.nan],
                # Reconstruction
                ["AcquisitionMode", (), "n/a"],  # todo : corresp ?
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

    def _translate_to_dict(self) -> dict:
        return {
            data.BIDSname: data.Value for _, data in self.meta_collection.iterrows()
        }

    def _build_dict(self, dcm_dir: Path, study: StudyName) -> dict:
        self._update_from_dcm(dcm_dir)
        # self._update_not_from_dcm
        # todo
        return self._translate_to_dict()

    def _get_admin_mode_from_study(self, study: StudyName) -> None:
        # todo : do we know for others ?
        if study == StudyName.ADNI:
            self.not_from_dcm["ModeOfAdministration"] = "bolus-infusion"

    def _get_injected_mass(self) -> None:
        # todo : "InjectedMass" : "InjectedRadioactivity" / "SpecificRadioactivity" -> need both else n/a
        # todo : Units : mole
        pass

    def _update_not_from_dcm(self, study: StudyName) -> None:
        self._get_admin_mode_from_study(study)
        self._get_injected_mass()

    def write_json(self, dcm_dir: Path, json_path: Path, study: StudyName):
        json_dict = self._build_dict(dcm_dir, study)
        with open(json_path, "w") as f:
            f.write(json.dumps(json_dict, skipkeys=True, indent=4))
