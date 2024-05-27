import json
from enum import Enum
from pathlib import Path
from typing import IO, Optional

from attrs import define, fields
from cattr.gen import make_dict_unstructure_fn, override
from cattr.preconf.json import make_converter

from clinica.utils.bids import BIDS_VERSION
from clinica.utils.exceptions import ClinicaCAPSError
from clinica.utils.inputs import DatasetType

CAPS_VERSION = "1.0.0"


@define
class CAPSDatasetDescription:
    """Model representing a CAPS dataset description."""

    name: str
    bids_version: str = BIDS_VERSION
    caps_version: str = CAPS_VERSION
    dataset_type: DatasetType = DatasetType.DERIVATIVE

    def write(self, to: IO[str]):
        json.dump(converter.unstructure(self), to, indent=4)

    @classmethod
    def from_values(
        cls,
        name: str,
        bids_version: Optional[str] = None,
        caps_version: Optional[str] = None,
    ):
        return cls(
            name,
            bids_version or BIDS_VERSION,
            caps_version or CAPS_VERSION,
            DatasetType.DERIVATIVE,
        )

    @classmethod
    def from_file(cls, json_file: Path):
        parsed = json.loads(json_file.read_text())
        try:
            return cls(
                parsed["Name"],
                parsed["BidsVersion"],
                parsed["CAPSVersion"],
                DatasetType(parsed["DatasetType"]),
            )
        except KeyError:
            raise ClinicaCAPSError(
                f"CAPS dataset_description.json file {json_file} is not valid and "
                "cannot be parsed as a CAPSDatasetDescription. "
                "Please verify that the file is well formatted."
            )

    def is_compatible_with(self, other) -> bool:
        if self.bids_version != other.bids_version:
            return False
        if self.caps_version != other.caps_version:
            return False
        return True


def _rename(name: str) -> str:
    """Rename attributes following the specification for the JSON file.

    Basically pascal case with known acronyms such as BIDS fully capitalized.
    """
    return "".join(
        word.upper() if word == "caps" else word.capitalize()
        for word in name.split("_")
    )


# Register a JSON converter for the CAPS dataset description model.
converter = make_converter()

converter.register_unstructure_hook(
    CAPSDatasetDescription,
    make_dict_unstructure_fn(
        CAPSDatasetDescription,
        converter,
        **{
            a.name: override(rename=_rename(a.name))
            for a in fields(CAPSDatasetDescription)
        },
    ),
)


def write_caps_dataset_description(
    name: str,
    caps_dir: Path,
    bids_version: Optional[str] = None,
    caps_version: Optional[str] = None,
) -> None:
    """Write `dataset_description.json` at the root of the CAPS directory."""
    from clinica.utils.stream import cprint

    new_desc = CAPSDatasetDescription.from_values(name, bids_version, caps_version)
    if (caps_dir / "dataset_description.json").exists():
        cprint(
            f"The CAPS dataset {name} already contains a dataset_description.json file.",
            lvl="info",
        )
        previous_desc = CAPSDatasetDescription.from_file(
            caps_dir / "dataset_description.json"
        )
        if not previous_desc.is_compatible_with(new_desc):
            msg = (
                f"Impossible to write the dataset_description.json file in {caps_dir} "
                "because it already exists and it contains incompatible metadata."
            )
            cprint(msg, lvl="error")
            raise ClinicaCAPSError(msg)
        if previous_desc.name != new_desc.name:
            new_desc.name = f"{previous_desc.name} + {new_desc.name}"
    with open(caps_dir / "dataset_description.json", "w") as f:
        new_desc.write(to=f)
