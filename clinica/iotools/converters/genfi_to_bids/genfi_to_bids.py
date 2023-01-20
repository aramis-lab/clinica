"""Convert the GENFI dataset into BIDS."""

from os import PathLike
from typing import Optional


def convert_images(
    path_to_dataset: PathLike,
    bids_dir: PathLike,
    path_to_clinical: Optional[PathLike],
    gif: bool,
) -> None:
    """Convert the entire dataset to BIDS.

    Scans available files in the path_to_dataset,
    identifies the patients that have images described by the JSON file,
    converts the image with the highest quality for each category.

    Parameters
    ----------
    path_to_dataset: PathLike
        Path to the raw images

    bids_dir: PathLike
        Path to directory where the bids will be written

    path_to_clinical: PathLike, optional
        Path to the clinical data associated with the dataset.
        If None, the clinical data won't be converted.

    gif: bool
        If True, indicates the user wants to have the values of the gif parcellation
    """
    import os

    import clinica.iotools.bids_utils as bids

    from .genfi_to_bids_utils import (
        complete_clinical_data,
        complete_imaging_data,
        dataset_to_bids,
        find_clinical_data,
        intersect_data,
        read_imaging_data,
        write_bids,
    )

    # read the clinical data files
    if path_to_clinical:
        df_demographics, df_imaging, df_clinical = find_clinical_data(path_to_clinical)
    # makes a df of the imaging data
    imaging_data = read_imaging_data(path_to_dataset)

    # complete the data extracted
    imaging_data = complete_imaging_data(imaging_data)

    # complete clinical data
    if path_to_clinical:
        df_clinical_complete = complete_clinical_data(
            df_demographics, df_imaging, df_clinical
        )

    # intersect the data
    if path_to_clinical:
        df_complete = intersect_data(imaging_data, df_clinical_complete)
    else:
        df_complete = imaging_data
    # build the tsv
    results = dataset_to_bids(df_complete, gif)

    write_bids(
        to=bids_dir,
        participants=results["participants"],
        sessions=results["sessions"],
        scans=results["scans"],
    )
    # fetch data for readme from markdown
    readme_path = os.path.join(
        os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        ),
        "docs",
        "Converters",
        "GENFItoBIDS.md",
    )
    try:
        readme_data = {
            "link": bids.parse_url(readme_path)[0],
            "desc": bids.parse_description(readme_path, 4, 5),
        }
    except IndexError:
        raise ValueError("Could not parse information for dataset.")
    bids.write_modality_agnostic_files(
        study_name="GENFI", readme_data=readme_data, bids_dir=bids_dir
    )
