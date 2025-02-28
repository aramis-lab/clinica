"""
Convert the AIBL dataset (https://www.aibl.csiro.au/) into BIDS.
"""
from pathlib import Path
from typing import Optional

from clinica.utils.filemanip import UserProvidedPath

__all__ = ["convert"]


def convert(
    input_dataset: UserProvidedPath,
    output_dataset: UserProvidedPath,
    input_clinical_data: UserProvidedPath,
    overwrite: bool = False,
    clinical_data_only: bool = False,
    subjects: Optional[UserProvidedPath] = None,
    n_procs: Optional[int] = 1,
    **kwargs,
) -> None:
    """Convert the AIBL dataset (https://www.aibl.csiro.au/) in a BIDS dataset.

    Parameters
    ----------
    input_dataset : Path
        The path to the input AIBL dataset.

    output_dataset : Path
        The path to the BIDS directory in which to write the output.

    input_clinical_data : Path
        The path to the clinical CSV files associated with the input dataset.

    overwrite : bool, optional
        Overwrites previously written nifti and json files.
        Default=False.

    clinical_data_only : bool, optional
        If True, only convert the clinical data without the imaging data.
        If False, everything is converted.
        Default=False.

    subjects : str or Path, optional
        The path to a file defining a subset of subjects to be converted.
        By default, all subjects available will be converted.

    n_procs : int, optional
        The requested number of processes.
        If specified, it should be between 1 and the number of available CPUs.
        Default=1.
    """
    from clinica.converters.factory import get_converter_name
    from clinica.converters.study_models import StudyName
    from clinica.utils.check_dependency import ThirdPartySoftware, check_software
    from clinica.utils.stream import cprint

    from .._utils import validate_input_path

    input_dataset = validate_input_path(input_dataset)
    output_dataset = validate_input_path(output_dataset, check_exist=False)
    input_clinical_data = validate_input_path(input_clinical_data)
    check_software(ThirdPartySoftware.DCM2NIIX)
    if subjects:
        cprint(
            (
                f"Subject filtering is not yet implemented in {get_converter_name(StudyName.AIBL)} converter. "
                "All subjects available will be converted."
            ),
            lvl="warning",
        )
    output_dataset.mkdir(parents=True, exist_ok=True)
    if not clinical_data_only:
        _convert_images(
            input_dataset,
            input_clinical_data,
            output_dataset,
            overwrite,
            n_procs=n_procs,
        )
    _convert_clinical_data(input_clinical_data, output_dataset)


def _convert_images(
    input_dataset: Path,
    input_clinical_data: Path,
    output_dataset: Path,
    overwrite: bool = False,
    n_procs: Optional[int] = 1,
) -> None:
    """Conversion of the AIBL imaging data in BIDS.

    Parameters
    ----------
    input_dataset : Path
        The path to the input AIBL dataset.

    input_clinical_data : Path
        The path to the clinical CSV files associated with the input dataset.

    output_dataset : Path
        The path to the BIDS directory in which to write the output.

    overwrite : bool, optional
        Overwrites previously written nifti and json files.
        Default=False.

    n_procs : int, optional
        The requested number of processes.
        If specified, it should be between 1 and the number of available CPUs.
        Default=1.
    """
    from clinica.converters.aibl_to_bids.utils import Modality, paths_to_bids

    created_files = [
        paths_to_bids(
            input_dataset,
            input_clinical_data,
            output_dataset,
            modality,
            overwrite=overwrite,
            n_procs=n_procs,
        )
        for modality in Modality
    ]
    _warn_about_missing_files(created_files)


def _warn_about_missing_files(files: list[list[Optional[Path]]]):
    from clinica.utils.stream import cprint

    missing_files = []
    for files_for_given_modality in files:
        for file in files_for_given_modality:
            if file is not None and not file.exists():
                missing_files.append(file)
    if missing_files:
        msg = "The following file were not converted:\n" + "\n".join(
            (str(f) for f in missing_files)
        )
        cprint(msg=msg, lvl="warning")


def _convert_clinical_data(input_clinical_data: Path, output_dataset: Path) -> None:
    """Conversion of the AIBL clinical data in BIDS.

    Parameters
    ----------
    input_clinical_data : Path
        The path to the clinical CSV files associated with the input dataset.

    output_dataset : Path
        The path to the BIDS directory in which to write the output.
    """
    from clinica.converters._utils import write_modality_agnostic_files
    from clinica.converters.aibl_to_bids.utils import (
        create_participants_tsv_file,
        create_scans_tsv_file,
        create_sessions_tsv_file,
    )
    from clinica.converters.study_models import StudyName
    from clinica.utils.stream import cprint

    clinical_specifications_folder = Path(__file__).parents[1] / "specifications"
    if not clinical_specifications_folder.exists():
        msg = (
            f"{clinical_specifications_folder} folder cannot be found ! "
            "This is an internal folder of Clinica."
        )
        cprint(msg, lvl="error")
        raise FileNotFoundError(msg)

    cprint("Creating modality agnostic files...", lvl="info")
    readme_data = {
        "link": "http://adni.loni.usc.edu/study-design/collaborative-studies/aibl/",
        "desc": (
            "The Australian Imaging, Biomarker & Lifestyle Flagship Study of Ageing (AIBL) seeks to discover which "
            "biomarkers, cognitive characteristics, and health and lifestyle factors determine the development of AD. "
            "Although AIBL and ADNI have many of the same goals, there are differences between the two projects."
        ),
    }
    write_modality_agnostic_files(
        study_name=StudyName.AIBL,
        readme_data=readme_data,
        bids_dir=output_dataset,
    )
    cprint("Creating participants.tsv...", lvl="info")
    create_participants_tsv_file(
        output_dataset,
        clinical_specifications_folder,
        input_clinical_data,
        delete_non_bids_info=True,
    )
    cprint("Creating sessions files...", lvl="info")
    create_sessions_tsv_file(
        output_dataset, input_clinical_data, clinical_specifications_folder
    )
    cprint("Creating scans files...", lvl="info")
    create_scans_tsv_file(
        output_dataset, input_clinical_data, clinical_specifications_folder
    )
