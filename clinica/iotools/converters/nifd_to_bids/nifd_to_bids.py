"""Convert the NIFD dataset into BIDS."""

from os import PathLike
from typing import List


def convert_images(
    path_to_dataset: PathLike,
    bids_dir: PathLike,
    path_to_clinical: PathLike,
) -> List[PathLike]:
    """Convert the entire dataset in BIDS.

    Scans available files in the path_to_dataset,
    identifies the patients that have images described by the JSON file,
    converts the image with the highest quality for each category.
    """
    import clinica.iotools.bids_utils as bids

    from .nifd_utils import (
        dataset_to_bids,
        read_clinical_data,
        read_imaging_data,
        write_bids,
    )

    clinical_data = read_clinical_data(path_to_clinical)
    imaging_data = read_imaging_data(path_to_dataset)

    participants, sessions, scans = dataset_to_bids(
        imaging_data=imaging_data, clinical_data=clinical_data
    )

    written = write_bids(
        to=bids_dir,
        participants=participants,
        sessions=sessions,
        scans=scans,
    )
    bids.write_modality_agnostic_files(study_name="NIFD", bids_dir=bids_dir)
    return written
