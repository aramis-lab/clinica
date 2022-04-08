"""Convert the UKB dataset into BIDS."""

from os import PathLike
from typing import List


def convert_images(
    path_to_dataset: PathLike,
    bids_dir: PathLike,
    path_to_clinical: PathLike,
) -> List[PathLike]:
    """Convert the entire dataset to BIDS.

    Scans available files in the path_to_dataset,
    identifies the patients that have images described by the JSON file,
    converts the image with the highest quality for each category.
    """
    from .ukb_utils import (
        complete_clinical,
        dataset_to_bids,
        find_clinical_data,
        intersect_data,
        read_imaging_data,
        write_bids,
    )

    # read the clinical data files
    df_clinical = find_clinical_data(path_to_clinical)

    # makes a df of the imaging data
    imaging_data = read_imaging_data(path_to_dataset)

    # intersect the data
    df_clinical = intersect_data(imaging_data, df_clinical)

    # complete clinical data
    df_clinical = complete_clinical(df_clinical)

    # build the tsv
    result = dataset_to_bids(df_clinical)

    written = write_bids(
        to=bids_dir,
        participants=result["participants"],
        sessions=result["sessions"],
        scans=result["scans"],
        dataset_directory=path_to_dataset,
    )

    return written
