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

    path_to_clinical: PathLike
        Path to the clinical data associated with the dataset

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
    # import markdown
    # path_to_markdown = os.path.join(
    #     os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))),
    #     "docs", "Converters",
    #     "GENFItoBIDS.md",)

    # f = open(path_to_markdown, 'r')
    # a = markdown.markdown( f.read())

    # print("markdown content:", markdown.markdown( f.read()) )
    readme_data = {
        "link": "https://www.genfi.org",
        # "desc":
        "desc": (
            "The Genetic Frontotemporal dementia Initiative (GENFI) is a group of research centres across Europe "
            "and Canada with expertise in familial FTD, and is co-ordinated by Professor Jonathan Rohrer at "
            "University College London. GENFI is the largest genetic FTD consortium to date and currently consists "
            "of sites across the UK, Netherlands, Belgium, France, Spain, Portugal, Italy, Germany, Sweden, "
            "Denmark, Finland and Canada. The aim of the study is to understand more about genetic FTD, "
            "particularly in those who have mutations in the progranulin (GRN), microtubule-associated protein tau "
            "(MAPT) and chromosome 9 open reading frame 72 (C9orf72) genes. GENFI investigates both people who have "
            "developed symptoms and also people who have a risk of developing symptoms in the future because they "
            "carry an abnormal genetic mutation. By studying these individuals who are destined to develop the "
            "disease later in life we can understand the development from the very earliest changes. The key "
            "objectives of GENFI are therefore to develop markers which help identify the disease at its earliest "
            "stage as well as markers that allow the progression of the disease to be tracked. We are now "
            "collaborating closely with other similar studies around the world through the FTD Prevention "
            "Initiative. Through this worldwide initiative we are working with pharmaceutical companies to help "
            "design clinical trials for genetic FTD."
        ),
    }
    bids.write_modality_agnostic_files(
        study_name="GENFI", readme_data=readme_data, bids_dir=bids_dir
    )
