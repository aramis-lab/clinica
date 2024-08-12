import json
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import IO, List, MutableSequence, NewType, Optional, Union

from attrs import define, fields
from cattr.gen import make_dict_structure_fn, make_dict_unstructure_fn, override
from cattr.preconf.json import make_converter
from packaging.version import Version

from clinica.iotools.bids_dataset_description import get_bids_version
from clinica.utils.bids import BIDS_VERSION
from clinica.utils.exceptions import ClinicaBIDSError, ClinicaCAPSError
from clinica.utils.inputs import DatasetType
from clinica.utils.stream import cprint, log_and_raise, log_and_warn

__all__ = [
    "CAPS_VERSION",
    "CAPSDatasetDescription",
    "VersionComparisonPolicy",
    "are_versions_compatible",
    "write_caps_dataset_description",
    "build_caps_dataset_description",
]


CAPS_VERSION = Version("1.0.0")

IsoDate = NewType("IsoDate", datetime)


class VersionComparisonPolicy(str, Enum):
    """Defines the different ways we can compare version numbers in Clinica.

    STRICT: version numbers have to match exactly.
    MINOR : version numbers have to have the same major and minor numbers.
    MAJOR: version numbers only need to share the same major number.
    """

    STRICT = "strict"
    MINOR = "minor"
    MAJOR = "major"


def are_versions_compatible(
    version1: Union[str, Version],
    version2: Union[str, Version],
    policy: Optional[Union[str, VersionComparisonPolicy]] = None,
) -> bool:
    """Returns whether the two provided versions are compatible or not depending on the policy.

    Parameters
    ----------
    version1 : str or Version
        The first version number to compare.

    version2 : str or Version
        The second version number to compare.

    policy : str or VersionComparisonPolicy, optional
        The policy under which to compare version1 with version2.
        By default, a strict policy is used, meaning that version
        numbers have to match exactly.

    Returns
    -------
    bool :
        True if version1 is 'compatible' with version2, False otherwise.
    """
    if isinstance(version1, str):
        version1 = Version(version1)
    if isinstance(version2, str):
        version2 = Version(version2)
    if policy is None:
        policy = VersionComparisonPolicy.STRICT
    else:
        policy = VersionComparisonPolicy(policy)
    if policy == VersionComparisonPolicy.STRICT:
        return version1 == version2
    if policy == VersionComparisonPolicy.MINOR:
        return version1.major == version2.major and version1.minor == version2.minor
    if policy == VersionComparisonPolicy.MAJOR:
        return version1.major == version2.major


@define
class CAPSProcessingDescription:
    """This class models a CAPS processing pipeline metadata.

    Attributes
    ----------
    name : str
        The name of the processing pipeline.
        Example: 't1-linear'.

    date : IsoDate
        The date at which the processing pipeline has been run.
        More precisely, this is the date at which the dataset_description.json
        file is written to disk, which precedes the date at which the pipeline
        finishes processing.

    author : str
        This is the name of the user who ran this processing pipeline.

    machine : str
        This is the name of the machine of which the processing pipeline was run.

    input_path : str
        This is the path to the input dataset.
    """

    name: str
    date: IsoDate
    author: str
    machine: str
    input_path: str

    @classmethod
    def from_values(cls, name: str, input_path: str):
        return cls(
            name,
            _get_current_timestamp(),
            _get_username(),
            _get_machine_name(),
            input_path,
        )

    @classmethod
    def from_dict(cls, values: dict):
        return cls(
            values["Name"],
            values["Date"],
            values["Author"],
            values["Machine"],
            values["InputPath"],
        )

    def match(
        self,
        processing_name: Optional[str] = None,
        processing_input_path: Optional[str] = None,
    ) -> bool:
        return not (
            (processing_name is not None and self.name != processing_name)
            or (
                processing_input_path is not None
                and (self.input_path != processing_input_path)
            )
        )

    def write(self, to: IO[str]):
        json.dump(converter.unstructure(self), to, indent=4)

    def __str__(self):
        return json.dumps(converter.unstructure(self))

    @classmethod
    def from_file(cls, json_file: Path):
        with open(json_file, "r") as fp:
            content = json.load(fp)
        return converter.structure(content, CAPSProcessingDescription)


@define
class CAPSDatasetDescription:
    """Model representing a CAPS dataset description.

    Attributes
    ----------
    name : str
        The name of the CAPS dataset.

    bids_version : str
        The version number of the BIDS specifications used.

    caps_version : str
        The version number of the CAPS specifications used.

    dataset_type : DatasetType
        The dataset type.

    processing : List of CAPSProcessingDescription
        The list of processing pipelines that have been run.
    """

    name: str
    bids_version: Version = BIDS_VERSION
    caps_version: Version = CAPS_VERSION
    dataset_type: DatasetType = DatasetType.DERIVATIVE
    processing: MutableSequence[CAPSProcessingDescription] = []

    def write(self, to: IO[str]):
        json.dump(converter.unstructure(self), to, indent=4)

    def __str__(self):
        return json.dumps(converter.unstructure(self))

    def has_processing(
        self,
        processing_name: Optional[str] = None,
        processing_input_path: Optional[str] = None,
    ) -> bool:
        return any(
            [
                processing.match(
                    processing_name=processing_name,
                    processing_input_path=processing_input_path,
                )
                for processing in self.processing
            ]
        )

    def get_processing(
        self,
        processing_name: Optional[str] = None,
        processing_input_path: Optional[str] = None,
    ) -> List[CAPSProcessingDescription]:
        return [
            processing
            for processing in self.processing
            if processing.match(
                processing_name=processing_name,
                processing_input_path=processing_input_path,
            )
        ]

    def delete_processing(
        self,
        processing_name: Optional[str] = None,
        processing_input_path: Optional[str] = None,
    ):
        for processing in self.processing:
            if processing.match(
                processing_name=processing_name,
                processing_input_path=processing_input_path,
            ):
                self.processing.remove(processing)

    def add_processing(
        self,
        processing_name: str,
        processing_input_path: str,
    ):
        new_processing = CAPSProcessingDescription.from_values(
            processing_name, processing_input_path
        )
        existing_processings = self.get_processing(
            processing_name, processing_input_path
        )
        if existing_processings:
            existing_processing = existing_processings[0]
            log_and_warn(
                (
                    f"The CAPS dataset '{self.name}' already has a processing named {processing_name} "
                    f"with an input folder set to {processing_input_path}:\n"
                    f"{existing_processing}\nIt will be overwritten with the following:\n{new_processing}"
                ),
                UserWarning,
            )
            self.delete_processing(existing_processing.name)
        self.processing.append(new_processing)

    @classmethod
    def from_values(
        cls,
        name: Optional[str] = None,
        bids_version: Optional[Version] = None,
        caps_version: Optional[Version] = None,
        processing: Optional[List[CAPSProcessingDescription]] = None,
    ):
        return cls(
            name or _generate_random_name(),
            bids_version or BIDS_VERSION,
            caps_version or CAPS_VERSION,
            DatasetType.DERIVATIVE,
            processing or [],
        )

    @classmethod
    def from_file(cls, json_file: Path):
        with open(json_file, "r") as fp:
            content = json.load(fp)
        return converter.structure(content, CAPSDatasetDescription)

    @classmethod
    def from_dict(cls, values: dict):
        processing = []
        if "Processing" in values:
            processing = [
                CAPSProcessingDescription.from_dict(processing)
                for processing in values["Processing"]
            ]
        return cls(
            values["Name"],
            values["BIDSVersion"],
            values["CAPSVersion"],
            DatasetType(values["DatasetType"]),
            processing,
        )

    def is_compatible_with(
        self, other, policy: Optional[Union[str, VersionComparisonPolicy]] = None
    ) -> bool:
        return are_versions_compatible(
            self.bids_version, other.bids_version, policy=policy
        ) and are_versions_compatible(
            self.caps_version, other.caps_version, policy=policy
        )


def _get_username() -> str:
    import os
    import pwd

    return pwd.getpwuid(os.getuid()).pw_name


def _get_machine_name() -> str:
    import platform

    return platform.node()


def _get_current_timestamp() -> IsoDate:
    return IsoDate(datetime.now())


def _generate_random_name() -> str:
    import uuid

    return str(uuid.uuid4())


def _rename(name: str) -> str:
    """Rename attributes following the specification for the JSON file.

    Basically pascal case with known acronyms such as CAPS fully capitalized.
    """
    return "".join(
        word.upper() if word in ("bids", "caps") else word.capitalize()
        for word in name.split("_")
    )


# Register a JSON converter for the CAPS dataset description model.
converter = make_converter()

# Unstructuring hooks first
converter.register_unstructure_hook(Version, lambda dt: str(dt))
converter.register_unstructure_hook(IsoDate, lambda dt: dt.isoformat())
caps_processing_field_renaming = {
    a.name: override(rename=_rename(a.name)) for a in fields(CAPSProcessingDescription)
}
caps_processing_field_renaming_unstructure_hook = make_dict_unstructure_fn(
    CAPSProcessingDescription,
    converter,
    **caps_processing_field_renaming,
)
converter.register_unstructure_hook(
    CAPSProcessingDescription,
    caps_processing_field_renaming_unstructure_hook,
)
caps_dataset_description_field_renaming = {
    a.name: override(rename=_rename(a.name)) for a in fields(CAPSDatasetDescription)
}
caps_dataset_field_renaming_unstructure_hook = make_dict_unstructure_fn(
    CAPSDatasetDescription,
    converter,
    **caps_dataset_description_field_renaming,
)
converter.register_unstructure_hook(
    CAPSDatasetDescription,
    caps_dataset_field_renaming_unstructure_hook,
)

# And structuring hooks
converter.register_structure_hook(Version, lambda ts, _: Version(ts))
converter.register_structure_hook(IsoDate, lambda ts, _: datetime.fromisoformat(ts))
caps_processing_field_renaming_structure_hook = make_dict_structure_fn(
    CAPSProcessingDescription,
    converter,
    **caps_processing_field_renaming,
)
converter.register_structure_hook(
    CAPSProcessingDescription,
    caps_processing_field_renaming_structure_hook,
)
caps_dataset_field_renaming_structure_hook = make_dict_structure_fn(
    CAPSDatasetDescription,
    converter,
    **caps_dataset_description_field_renaming,
)
converter.register_structure_hook(
    CAPSDatasetDescription,
    caps_dataset_field_renaming_structure_hook,
)


def write_caps_dataset_description(
    input_dir: Path,
    output_dir: Path,
    processing_name: str,
    dataset_name: Optional[str] = None,
    bids_version: Optional[str] = None,
    caps_version: Optional[str] = None,
) -> None:
    """Write `dataset_description.json` at the root of the CAPS directory.

    Parameters
    ----------
    input_dir : Path
        The path to the folder of the input dataset.
        It can be a BIDS dataset or a CAPS dataset.

    output_dir : Path
        The path to the folder of the output dataset.
        This has to be a CAPS dataset, and this is where
        the requested dataset_description.json file will be written.

    processing_name : str
        The name of the processing performed. By default, pipelines of
        Clinica will set this as the name of the pipeline, but any name
        is possible.

    dataset_name : str, optional
        The name of the CAPS dataset. If not specified, a random identifier will
        be generated. If a dataset_description.json file already exists, the
        existing name will be kept.

    bids_version : str, optional
        The version of the BIDS specifications used.
        By default, this will be set as the BIDS version currently supported by Clinica.

    caps_version : str, optional
        The version of the CAPS specifications used.
        By default, this will be set as the CAPS version currently supported by Clinica.
    """
    description = build_caps_dataset_description(
        input_dir,
        output_dir,
        processing_name,
        dataset_name=dataset_name,
        bids_version=bids_version,
        caps_version=caps_version,
    )
    with open(output_dir / "dataset_description.json", "w") as f:
        description.write(to=f)


def build_caps_dataset_description(
    input_dir: Path,
    output_dir: Path,
    processing_name: str,
    dataset_name: Optional[str] = None,
    bids_version: Optional[str] = None,
    caps_version: Optional[str] = None,
) -> CAPSDatasetDescription:
    """Generate the CAPSDatasetDescription for a given CAPS dataset.

    Parameters
    ----------
    input_dir : Path
        The path to the folder of the input dataset.
        It can be a BIDS dataset or a CAPS dataset.

    output_dir : Path
        The path to the folder of the output dataset.
        This has to be a CAPS dataset, and this is where
        the requested dataset_description.json file will be written.

    processing_name : str
        The name of the processing performed. By default, pipelines of
        Clinica will set this as the name of the pipeline, but any name
        is possible.

    dataset_name : str, optional
        The name of the CAPS dataset. If not specified, a random identifier will
        be generated. If a dataset_description.json file already exists, the
        existing name will be kept.

    bids_version : str, optional
        The version of the BIDS specifications used.
        By default, this will be set as the BIDS version currently supported by Clinica.

    caps_version : str, optional
        The version of the CAPS specifications used.
        By default, this will be set as the CAPS version currently supported by Clinica.

    Returns
    -------
    CAPSDatasetDescription :
        The CAPSDatasetDescription generated.
    """
    bids_version_from_input_dir = None
    try:
        bids_version_from_input_dir = get_bids_version(input_dir)
    except ClinicaBIDSError:
        log_and_warn(
            (
                f"Unable to retrieve the BIDS version from the input folder {input_dir}."
                f"Please verify your input dataset. Clinica will assume a BIDS version of {BIDS_VERSION}."
            ),
            UserWarning,
        )
    if (
        bids_version is not None
        and bids_version_from_input_dir is not None
        and bids_version != bids_version_from_input_dir
    ):
        log_and_raise(
            f"The input dataset {input_dir} has BIDS specifications following "
            f"version {bids_version_from_input_dir}, while the BIDS specifications version "
            f"asked for the CAPS creation is {bids_version}. "
            "Please make sure the versions are the same before processing.",
            ClinicaBIDSError,
        )
    new_description = CAPSDatasetDescription.from_values(
        dataset_name, bids_version_from_input_dir, caps_version
    )
    if (output_dir / "dataset_description.json").exists():
        cprint(
            (
                f"The CAPS dataset '{dataset_name}', located at {output_dir}, already "
                "contains a 'dataset_description.json' file."
            ),
            lvl="info",
        )
        previous_description = CAPSDatasetDescription.from_file(
            output_dir / "dataset_description.json"
        )
        if not previous_description.is_compatible_with(
            new_description, VersionComparisonPolicy.STRICT
        ):
            msg = (
                f"Impossible to write the 'dataset_description.json' file in {output_dir} "
                "because it already exists and it contains incompatible metadata."
            )
            log_and_raise(msg, ClinicaCAPSError)
        if previous_description.name != new_description.name:
            log_and_warn(
                (
                    f"The existing CAPS dataset, located at {output_dir} has a "
                    f"name '{previous_description.name}' different from the new "
                    f"name '{new_description.name}'. The old name will be kept."
                ),
                UserWarning,
            )
            new_description.name = previous_description.name
        for processing in previous_description.processing:
            if not processing.match(
                processing_name=processing_name, processing_input_path=str(input_dir)
            ):
                new_description.processing.append(processing)
    new_description.add_processing(processing_name, str(input_dir))
    return new_description
