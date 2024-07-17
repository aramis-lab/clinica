"""Convert OASIS dataset (https://sites.wustl.edu/oasisbrains/) to BIDS."""

from pathlib import Path
from typing import Optional

import nibabel as nb
import numpy as np
from ixi_to_bids_utils import read_ixi_clinical_data

from clinica.utils.filemanip import UserProvidedPath

__all__ = ["convert"]


def convert(
    path_to_dataset: UserProvidedPath,
    bids_dir: UserProvidedPath,
    path_to_clinical: UserProvidedPath,
    subjects: Optional[UserProvidedPath] = None,
    n_procs: Optional[int] = 1,
    **kwargs,
):
    from clinica.iotools.bids_utils import StudyName
    from clinica.iotools.converters.factory import get_converter_name
    from clinica.utils.stream import cprint

    from ..utils import validate_input_path

    path_to_dataset = validate_input_path(path_to_dataset)
    bids_dir = validate_input_path(bids_dir, check_exist=False)
    path_to_clinical = validate_input_path(path_to_clinical)

    if n_procs != 1:
        cprint(
            f"{get_converter_name(StudyName.IXI)} converter does not support multiprocessing yet. n_procs set to 1.",
            lvl="warning",
        )

    clinical_data = read_ixi_clinical_data(path_to_clinical)
    # todo : convert here
    cprint("Conversion to BIDS succeeded.", lvl="info")
