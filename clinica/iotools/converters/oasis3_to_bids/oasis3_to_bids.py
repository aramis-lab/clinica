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
    from .oasis3_utils import (
        dataset_to_bids,
        intersect_data,
        read_clinical_data,
        read_imaging_data,
        write_bids,
    )

    # read the clinical data files
    df_pet, df_mri, df_subject, df_pup, df_adrc = read_clinical_data(path_to_clinical)

    # makes a df of the imaging data
    imaging_data = read_imaging_data(path_to_dataset)

    # intersect the data
    imaging_data, df_small = intersect_data(imaging_data, df_mri, df_subject, df_adrc)

    # build the tsv
    participants, sessions, scans = dataset_to_bids(imaging_data, df_small)

    written = write_bids(
        to=bids_dir, participants=participants, sessions=sessions, scans=scans
    )

    return written
