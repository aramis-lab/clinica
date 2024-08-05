"""Convert OASIS dataset (https://sites.wustl.edu/oasisbrains/) to BIDS."""

from pathlib import Path
from typing import Optional

import nibabel as nb
import numpy as np
from iotools.bids_utils import write_modality_agnostic_files
from ixi_to_bids_utils import *

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

    subjects = subjects if subjects else get_subjects_list_from_data(path_to_dataset)
    subjects = filter_subjects_list(subjects, clinical_data)

    # todo : treat DTI different (load and merge)

    image_data = get_img_data_df(path_to_dataset)

    for subject in subjects:
        write_subject(select_subject_data(image_data, subject), bids_dir)

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

    cprint("Conversion to BIDS succeeded.", lvl="info")
