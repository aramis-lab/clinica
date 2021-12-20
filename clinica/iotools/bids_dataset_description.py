from enum import Enum
from typing import IO

import attr
from cattr.preconf.json import make_converter
from cattr.gen import make_dict_unstructure_fn, override

BIDS_VERSION = "1.6.0"


@attr.define
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
    return "".join(
        word.upper() if word == "bids" else word.capitalize()
        for word in name.split("_")
    )


converter = make_converter()

converter.register_unstructure_hook(
    BIDSDatasetDescription,
    make_dict_unstructure_fn(
        BIDSDatasetDescription,
        converter,
        **{
            a.name: override(rename=_rename(a.name))
            for a in attr.fields(BIDSDatasetDescription)
        },
    ),
)
