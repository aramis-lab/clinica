from enum import Enum
from typing import IO

from attrs import define, fields
from cattr.gen import make_dict_unstructure_fn, override
from cattr.preconf.json import make_converter

BIDS_VERSION = "1.7.0"


@define
class BIDSDatasetDescription:
    """Model representing a BIDS dataset description.

    See https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
    """

    class DatasetType(str, Enum):
        raw = "raw"
        derivative = "derivative"

    name: str
    bids_version: str = BIDS_VERSION
    dataset_type: DatasetType = DatasetType.raw

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
