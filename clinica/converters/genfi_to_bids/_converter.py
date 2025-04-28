"""Convert the GENFI dataset into BIDS."""

from pathlib import Path
from typing import Optional

from clinica.utils.filemanip import UserProvidedPath

__all__ = ["convert"]


def convert(
    path_to_dataset: UserProvidedPath,
    bids_dir: UserProvidedPath,
    path_to_clinical: Optional[UserProvidedPath] = None,
    gif: Optional[bool] = False,
    path_to_clinical_tsv: Optional[UserProvidedPath] = None,
    subjects: Optional[UserProvidedPath] = None,
    n_procs: Optional[int] = 1,
    **kwargs,
) -> None:
    """Convert the entire dataset to BIDS.

    Scans available files in the path_to_dataset,
    identifies the patients that have images described by the JSON file,
    converts the image with the highest quality for each category.

    Parameters
    ----------
    path_to_dataset: Path
        The path to the raw images.

    bids_dir: Path
        The path to directory where the bids will be written.

    path_to_clinical: Path, optional
        The path to the clinical data associated with the dataset.
        If None, the clinical data won't be converted.

    gif: bool, optional
        If True, indicates the user wants to have the values of the gif parcellation

    path_to_clinical_tsv: Path, optional
        The path to a TSV file containing the additional data the user wants to have in the BIDS output.
        If None, no additional data will be added.

    subjects : str or Path, optional
        The path to a file defining a subset of subjects to be converted.
        By default, all subjects available will be converted.

    n_procs : int, optional
        The requested number of processes.
        If specified, it should be between 1 and the number of available CPUs.
        Default=1.
    """
    from clinica.utils.check_dependency import ThirdPartySoftware, check_software
    from clinica.utils.stream import cprint

    from .._utils import validate_input_path, write_modality_agnostic_files
    from ..factory import get_converter_name
    from ..study_models import StudyName
    from ._utils import (
        merge_imaging_and_clinical_data,
        parse_clinical_data,
        parse_imaging_data,
        prepare_dataset_to_bids_format,
        write_bids,
    )

    path_to_dataset = validate_input_path(path_to_dataset)
    bids_dir = validate_input_path(bids_dir, check_exist=False)
    if path_to_clinical:
        path_to_clinical = validate_input_path(path_to_clinical)
    if path_to_clinical_tsv:
        path_to_clinical_tsv = validate_input_path(path_to_clinical_tsv)
    check_software(ThirdPartySoftware.DCM2NIIX)
    _check_clinical_path_inputs(path_to_clinical_tsv, path_to_clinical)
    if subjects:
        cprint(
            (
                f"Subject filtering is not yet implemented in {get_converter_name(StudyName.GENFI)} converter. "
                "All subjects available will be converted."
            ),
            lvl="warning",
        )
    if n_procs != 1:
        cprint(
            f"{get_converter_name(StudyName.GENFI)} converter does not support multiprocessing yet. n_procs set to 1.",
            lvl="warning",
        )
    imaging_data = parse_imaging_data(path_to_dataset)
    if path_to_clinical:
        clinical_data = parse_clinical_data(path_to_clinical)
        imaging_data = merge_imaging_and_clinical_data(imaging_data, clinical_data)
    results = prepare_dataset_to_bids_format(imaging_data, gif, path_to_clinical_tsv)
    write_bids(
        to=bids_dir,
        participants=results["participants"],
        sessions=results["sessions"],
        scans=results["scans"],
        source=path_to_dataset,
    )
    write_modality_agnostic_files(
        study_name=StudyName.GENFI,
        readme_data={
            "link": _get_link(),
            "desc": _get_description(),
        },
        bids_dir=bids_dir,
    )
    cprint("Conversion to BIDS succeeded.", lvl="info")


def _check_clinical_path_inputs(path_to_clinical_tsv: Path, path_to_clinical: Path):
    """Check that if a clinical tsv is given, a path to the clinical data is given as well."""
    from clinica.converters.factory import get_converter_name
    from clinica.converters.study_models import StudyName
    from clinica.utils.stream import cprint

    if path_to_clinical_tsv and not path_to_clinical:
        msg = (
            f"The {get_converter_name(StudyName.GENFI)} converter is unable to convert the clinical data because "
            "the path to these data was not provided while a TSV file with additional "
            f"data was given ({path_to_clinical_tsv}). You can either use the appropriate "
            "option from the clinica command line interface to provide the missing path, "
            "or chose to not convert clinical data at all."
        )
        cprint(msg, lvl="error")
        raise ValueError(msg)


def _get_link() -> str:
    return "https://www.genfi.org"


def _get_description() -> str:
    return (
        "The Genetic Frontotemporal dementia Initiative (GENFI) is a group of research "
        "centres across Europe and Canada with expertise in familial FTD, and is "
        "co-ordinated by Professor Jonathan Rohrer at University College London. "
        "GENFI is the largest genetic FTD consortium to date and currently consists "
        "of sites across the UK, Netherlands, Belgium, France, Spain, Portugal, Italy, "
        "Germany, Sweden, Denmark, Finland and Canada. The aim of the study is to "
        "understand more about genetic FTD, particularly in those who have mutations "
        "in the progranulin (GRN), microtubule-associated protein tau (MAPT) and "
        "chromosome 9 open reading frame 72 (C9orf72) genes. GENFI investigates both "
        "people who have developed symptoms and also people who have a risk of developing "
        "symptoms in the future because they carry an abnormal genetic mutation. "
        "By studying these individuals who are destined to develop the disease later in "
        "life we can understand the development from the very earliest changes. "
        "The key objectives of GENFI are therefore to develop markers which help identify "
        "the disease at its earliest stage as well as markers that allow the progression of "
        "the disease to be tracked. We are now collaborating closely with other similar "
        "studies around the world through the FTD Prevention Initiative. "
        "Through this worldwide initiative we are working with pharmaceutical companies to "
        "help design clinical trials for genetic FTD."
    )
