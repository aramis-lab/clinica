import json
from pathlib import Path
from typing import IO

from attrs import define, fields
from cattr.gen import make_dict_unstructure_fn, override
from cattr.preconf.json import make_converter
from packaging.version import InvalidVersion, Version

from clinica.utils.bids import BIDS_VERSION
from clinica.utils.exceptions import ClinicaBIDSError
from clinica.utils.inputs import DatasetType
from clinica.utils.stream import log_and_raise

__all__ = [
    "BIDSDatasetDescription",
    "get_bids_version",
]


@define
class BIDSDatasetDescription:
    """Model representing a BIDS dataset description.

    See https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
    """

    name: str
    bids_version: Version = BIDS_VERSION
    dataset_type: DatasetType = DatasetType.RAW

    def write(self, to: IO[str]):
        import json

        json.dump(converter.unstructure(self), to, indent=4)


def _rename(name: str) -> str:
    """Rename attributes following the specification for the JSON file.

    Basically pascal case with known acronyms such as BIDS fully capitalized.
    """
    return "".join(
        word.upper() if word == "bids" else word.capitalize()
        for word in name.split("_")
    )


# Register a JSON converter for the BIDS dataset description model.
converter = make_converter()
converter.register_unstructure_hook(Version, lambda dt: str(dt))
converter.register_unstructure_hook(
    BIDSDatasetDescription,
    make_dict_unstructure_fn(
        BIDSDatasetDescription,
        converter,
        **{
            a.name: override(rename=_rename(a.name))
            for a in fields(BIDSDatasetDescription)
        },
    ),
)


def get_bids_version(dataset_folder: Path) -> Version:
    """Returns the BIDS version number of a BIDS or CAPS dataset."""
    try:
        with open(dataset_folder / "dataset_description.json", "r") as fp:
            bids_metadata = json.load(fp)
        return Version(bids_metadata["BIDSVersion"])
    except InvalidVersion as e:
        log_and_raise(
            (
                f"File {dataset_folder / 'dataset_description.json'} has a "
                f"BIDS version number not properly formatted:\n{e}"
            ),
            ClinicaBIDSError,
        )
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
