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

    name: str
    bids_version: str = BIDS_VERSION

    def write(self, to: IO[str], readme_dict):
        import clinica

        # datadict_link = {
        #     "HABS": "https://habs.mgh.harvard.edu",
        # }
        # datadict_description = {
        #     "HABS": "The overall goal of the Harvard Aging Brain Study (HABS) is to elucidate the earliest changes in molecular, functional and structural imaging markers that signal the transition from normal cognition to progressive cognitive decline along the trajectory of preclinical Alzheimer’s Disease.",
        # }

        to.write(
            f"This BIDS directory was generated with Clinica v{clinica.__version__}.\n"
            f"More information on https://www.clinica.run\n"
            f"\n"
            f"Study: {self.name}\n"
            f"\n"
            f"{readme_dict['desc']}\n\n"
            f"Find more about it and about the data user agreement: {readme_dict['link']}"
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
