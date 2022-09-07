from enum import Enum
from typing import IO

from attrs import define, fields
from cattr.gen import make_dict_unstructure_fn, override
from cattr.preconf.json import make_converter

BIDS_VERSION = "1.7.0"


@define
class BIDSReadme:
    """Model representing a BIDS ReadMe.

    See
    """

    class DatasetType(str, Enum):
        raw = "raw"
        derivative = "derivative"

    name: str
    bids_version: str = BIDS_VERSION
    dataset_type: DatasetType = DatasetType.raw

    def write(self, to: IO[str]):
        import clinica

        to.write(
            f"This BIDS directory was generated with Clinica v{clinica.__version__}.\n"
            f"More information on https://www.clinica.run\n"
        )


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
    BIDSReadme,
    make_dict_unstructure_fn(
        BIDSReadme,
        converter,
        **{a.name: override(rename=_rename(a.name)) for a in fields(BIDSReadme)},
    ),
)
