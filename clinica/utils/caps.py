import datetime
import json
from pathlib import Path
from typing import IO, List, MutableSequence, Optional

from attrs import define, fields
from cattr.gen import make_dict_structure_fn, make_dict_unstructure_fn, override
from cattr.preconf.json import make_converter

from clinica.utils.bids import BIDS_VERSION
from clinica.utils.exceptions import ClinicaBIDSError, ClinicaCAPSError
from clinica.utils.inputs import DatasetType
from clinica.utils.stream import log_and_raise, log_and_warn

__all__ = [
    "CAPS_VERSION",
    "CAPSDatasetDescription",
    "write_caps_dataset_description",
    "build_caps_dataset_description",
]


CAPS_VERSION = "1.0.0"


def _get_username() -> str:
    import os
    import pwd

    return pwd.getpwuid(os.getuid()).pw_name


def _get_machine_name() -> str:
    import platform

    return platform.node()


def _get_current_timestamp() -> datetime.datetime:
    return datetime.datetime.now()


def _generate_random_name() -> str:
    import uuid

    return str(uuid.uuid4())


def _get_bids_version(dataset_folder: Path):
    """Returns the BIDS version number of a BIDS or CAPS dataset."""
    try:
        with open(dataset_folder / "dataset_description.json", "r") as fp:
            bids_metadata = json.load(fp)
        return bids_metadata["BIDSVersion"]
    except FileNotFoundError:
        log_and_raise(
            (
                f"File {dataset_folder / 'dataset_description.json'} is missing "
                "while it is mandatory for a BIDS/CAPS dataset."
            ),
            ClinicaBIDSError,
        )
    except KeyError:
        log_and_raise(
            (
                f"File {dataset_folder / 'dataset_description.json'} is missing a "
                "'BIDSVersion' key while it is mandatory."
            ),
            ClinicaBIDSError,
        )
    except json.JSONDecodeError as e:
        log_and_raise(
            f"File {dataset_folder / 'dataset_description.json'} is not formatted correctly:\n{e}.",
            ClinicaBIDSError,
        )


@define
class CAPSProcessingDescription:
    """This class models a CAPS processing pipeline metadata.

    Attributes
    ----------
    name : str
        The name of the processing pipeline.
        Example: 't1-linear'.

    date : datetime
        The date at which the processing pipeline has been run.
        More precisely, this is the date at which the dataset_description.json
        file is written to disk, which precedes the date at which the pipeline
        finishes processing.

    author : str
        This is the name of the user who ran this processing pipeline.

    machine : str
        This is the name of the machine of which the processing pipeline was run.

    processing_path : str
        This is the path to the processing folder(s) relative to the root of the
        CAPS dataset.

    input_path : str
        This is the path to the input dataset.
    """

    name: str
    date: datetime.datetime
    author: str
    machine: str
    processing_path: str
    input_path: str

    @classmethod
    def from_values(cls, name: str, processing_path: str, input_path: str):
        return cls(
            name,
            _get_current_timestamp(),
            _get_username(),
            _get_machine_name(),
            processing_path,
            input_path,
        )

    @classmethod
    def from_dict(cls, values: dict):
        return cls(
            values["Name"],
            values["Date"],
            values["Author"],
            values["Machine"],
            values["ProcessingPath"],
            values["InputPath"],
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
    bids_version: str = BIDS_VERSION
    caps_version: str = CAPS_VERSION
    dataset_type: DatasetType = DatasetType.DERIVATIVE
    processing: MutableSequence[CAPSProcessingDescription] = []

    def write(self, to: IO[str]):
        json.dump(converter.unstructure(self), to, indent=4)

    def __str__(self):
        return json.dumps(converter.unstructure(self))

    def has_processing(self, processing_name: str) -> bool:
        return any(
            [processing.name == processing_name for processing in self.processing]
        )

    def get_processing(
        self, processing_name: str
    ) -> Optional[CAPSProcessingDescription]:
        for processing in self.processing:
            if processing.name == processing_name:
                return processing
        return None

    def delete_processing(self, processing_name: str):
        for processing in self.processing:
            if processing.name == processing_name:
                self.processing.remove(processing)

    def add_processing(
        self,
        processing_name: str,
        processing_output_path: str,
        processing_input_path: str,
    ):
        new_processing = CAPSProcessingDescription.from_values(
            processing_name, processing_output_path, processing_input_path
        )
        if (existing_processing := self.get_processing(processing_name)) is not None:
            log_and_warn(
                (
                    f"The CAPS dataset '{self.name}' already has a processing named {processing_name}:\n"
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
        bids_version: Optional[str] = None,
        caps_version: Optional[str] = None,
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

    def is_compatible_with(self, other) -> bool:
        if self.bids_version != other.bids_version:
            return False
        if self.caps_version != other.caps_version:
            return False
        return True


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
converter.register_unstructure_hook(datetime.datetime, lambda dt: dt.isoformat())
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
converter.register_structure_hook(
    datetime.datetime, lambda ts, _: datetime.datetime.fromisoformat(ts)
)
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
    processing_output_path: str,
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

    processing_output_path : str
        The path to the subfolder(s) in which the results of the processing
        will be stored, relative to the root of the CAPS dataset (defined as
        output_dir). If there are multiple folders, use a regexp.
        For example, for t1-linear: 'subjects/*/*/t1-linear'.

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
        processing_output_path,
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
    processing_output_path: str,
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

    processing_output_path : str
        The path to the subfolder(s) in which the results of the processing
        will be stored, relative to the root of the CAPS dataset (defined as
        output_dir). If there are multiple folders, use a regexp.
        For example, for t1-linear: 'subjects/*/*/t1-linear'.

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
    from clinica.utils.stream import cprint, log_and_raise

    bids_version_from_input_dir = None
    try:
        bids_version_from_input_dir = _get_bids_version(input_dir)
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
    new_desc = CAPSDatasetDescription.from_values(
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
        previous_desc = CAPSDatasetDescription.from_file(
            output_dir / "dataset_description.json"
        )
        if not previous_desc.is_compatible_with(new_desc):
            msg = (
                f"Impossible to write the 'dataset_description.json' file in {output_dir} "
                "because it already exists and it contains incompatible metadata."
            )
            log_and_raise(msg, ClinicaCAPSError)
        if previous_desc.name != new_desc.name:
            log_and_warn(
                (
                    f"The existing CAPS dataset, located at {output_dir} has a name '{previous_desc.name}' different "
                    f"from the new name '{new_desc.name}'. The old name will be kept."
                ),
                UserWarning,
            )
            new_desc.name = previous_desc.name
        for processing in previous_desc.processing:
            new_desc.processing.append(processing)
    new_desc.add_processing(processing_name, processing_output_path, str(input_dir))
    return new_desc
