"""Convert the NIFD dataset into BIDS."""

from pathlib import Path
from typing import Optional, Union

__all__ = ["convert"]


def convert(
    path_to_dataset: Path,
    bids_dir: Path,
    path_to_clinical: Path,
    subjects: Optional[Union[str, Path]] = None,
    n_procs: Optional[int] = 1,
    **kwargs,
):
    """Convert the entire dataset in BIDS.

    Scans available files in the path_to_dataset,
    identifies the patients that have images described by the JSON file,
    converts the image with the highest quality for each category.
    """
    from clinica.iotools.bids_utils import StudyName, write_modality_agnostic_files
    from clinica.iotools.converters.factory import get_converter_name
    from clinica.utils.check_dependency import ThirdPartySoftware, check_software
    from clinica.utils.stream import cprint

    from ..utils import validate_input_path
    from .nifd_utils import (
        dataset_to_bids,
        read_clinical_data,
        read_imaging_data,
        write_bids,
    )

    path_to_dataset = validate_input_path(path_to_dataset)
    bids_dir = validate_input_path(bids_dir, check_exist=False)
    path_to_clinical = validate_input_path(path_to_clinical)
    check_software(ThirdPartySoftware.DCM2NIIX)
    if subjects:
        cprint(
            (
                f"Subject filtering is not yet implemented in {get_converter_name(StudyName.NIFD)} converter. "
                "All subjects available will be converted."
            ),
            lvl="warning",
        )
    if n_procs != 1:
        cprint(
            f"{get_converter_name(StudyName.NIFD)} converter does not support multiprocessing yet. n_procs set to 1.",
            lvl="warning",
        )
    clinical_data = read_clinical_data(path_to_clinical)
    imaging_data = read_imaging_data(path_to_dataset)
    participants, sessions, scans = dataset_to_bids(
        imaging_data=imaging_data,
        clinical_data=clinical_data,
    )
    write_bids(
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
    write_modality_agnostic_files(
        study_name=StudyName.NIFD,
        readme_data=readme_data,
        bids_dir=bids_dir,
    )
    cprint("Conversion to BIDS succeeded.", lvl="info")
