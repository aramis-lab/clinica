"""Convert the UKB dataset into BIDS."""

from os import PathLike


def convert_images(
    path_to_dataset: PathLike,
    bids_dir: PathLike,
    path_to_clinical: PathLike,
):
    """Convert the entire dataset to BIDS.

    Scans available files in the path_to_dataset,
    identifies the patients that have images described by the JSON file,
    converts the image with the highest quality for each category.
    """

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
    df_demographics, df_imaging, df_clinical = find_clinical_data(path_to_clinical)
    # makes a df of the imaging data
    imaging_data = read_imaging_data(path_to_dataset)

    # complete the data extracted
    imaging_data = complete_imaging_data(imaging_data)

    # complete clinical data
    df_clinical_complete = complete_clinical_data(
        df_demographics, df_imaging, df_clinical
    )

    # intersect the data
    df_complete = intersect_data(imaging_data, df_clinical_complete)

    # build the tsv
    results = dataset_to_bids(df_complete)

    write_bids(
        to=bids_dir,
        participants=results["participants"],
        sessions=results["sessions"],
        scans=results["scans"],
        dataset_directory=path_to_dataset,
    )
    readme_data = {
        "link": "https://www.genfi.org",
        "desc": (
            "UK Biobank is a large-scale biomedical database and research resource, containing in-depth genetic and "
            "health information from half a million UK participants. The database is regularly augmented with "
            "additional data and is globally accessible to approved researchers undertaking vital research into the "
            "most common and life-threatening diseases. It is a major contributor to the advancement of modern "
            "medicine it and has led to the discovery of several scientific advances and numerous treatments to "
            "improve human health."
        ),
    }
    bids.write_modality_agnostic_files(
        study_name="GENFI", readme_data=readme_data, bids_dir=bids_dir
    )
