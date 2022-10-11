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
    readme_data = {
        "link": "https://ida.loni.usc.edu/home/projectPage.jsp?project=NIFD&page=HOME&subPage=OVERVIEW_PR#",
        "desc": (
            "NIFD is the nickname for the frontotemporal lobar degeneration neuroimaging initiative "
            "(FTLDNI, AG032306), which was funded by the NIA and NINDS to characterize longitudinal clinical and "
            "imaging changes in FTLD.The imaging and clinical methods are the same for NIFD and for the 4-Repeat "
            "Tauopathy Neuroimaging Initiative (4RTNI), which is also available for download from LONI. Controls for "
            "NIFD are the same controls as those collected for 4RTNI."
        ),
    }
    bids.write_modality_agnostic_files(
        study_name="NIFD", readme_data=readme_data, bids_dir=bids_dir
    )
    return written
