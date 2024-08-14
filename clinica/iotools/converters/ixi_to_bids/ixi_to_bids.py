"""Convert IXI dataset (https://brain-development.org/ixi-dataset/) to BIDS."""

from pathlib import Path
from typing import Optional

import nibabel as nb
import numpy as np

from clinica.iotools.bids_utils import write_modality_agnostic_files
from clinica.iotools.converters.ixi_to_bids.ixi_to_bids_utils import (
    define_participants,
    read_clinical_data,
    write_ixi_subject_data,
    write_scans,
    write_sessions,
)
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
    if subjects:
        subjects = validate_input_path(subjects)

    if n_procs != 1:
        cprint(
            f"{get_converter_name(StudyName.IXI)} converter does not support multiprocessing yet. n_procs set to 1.",
            lvl="warning",
        )

    clinical_data = read_clinical_data(path_to_clinical)
    participants = define_participants(path_to_dataset, clinical_data, subjects)

    for participant in participants:
        cprint(f"Converting IXI subject {participant} to BIDS")
        write_ixi_subject_data(
            bids_dir=bids_dir, participant=participant, path_to_dataset=path_to_dataset
        )
        write_sessions(
            bids_dir=bids_dir, participant=participant, clinical_data=clinical_data
        )
        write_scans(bids_dir=bids_dir, participant=participant)

    readme_data = {
        "link": "https://brain-development.org/ixi-dataset/",
        "desc": (
            "IXI is the nickname for the Information eXtraction from Images project, "
            "which issued a dataset of nearly 600 images from healthy subjects. The MR"
            "acquisition protocol includes T1,T2, PD weighted, MRA and diffusion-weighted"
            "images. Three hospitals in London were involved in data collection."
        ),
    }

    write_modality_agnostic_files(
        study_name=StudyName.IXI, readme_data=readme_data, bids_dir=bids_dir
    )

    cprint("Conversion to BIDS finished.", lvl="info")
