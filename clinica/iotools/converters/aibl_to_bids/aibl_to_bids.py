"""
Convert the AIBL dataset (https://www.aibl.csiro.au/) into BIDS.
"""
from pathlib import Path
from typing import Optional


def convert(
    input_dataset: Path,
    input_clinical_data: Path,
    output_dataset: Path,
    overwrite: bool = False,
    clinical_data_only: bool = False,
    n_procs: Optional[int] = 1,
) -> None:
    """Convert the AIBL dataset (https://www.aibl.csiro.au/) in a BIDS dataset.

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

    clinical_data_only : bool, optional
        If True, only convert the clinical data without the imaging data.
        If False, everything is converted.
        Default=False.

    n_procs : int, optional
        The requested number of processes.
        If specified, it should be between 1 and the number of available CPUs.
        Default=1.
    """
    from clinica.utils.check_dependency import check_dcm2niix

    check_dcm2niix()

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
    from os.path import exists

    from clinica.iotools.converters.aibl_to_bids.utils import Modality, paths_to_bids
    from clinica.utils.stream import cprint

    list_of_created_files = [
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
    missing_files = []
    for modality_list in list_of_created_files:
        for file in modality_list:
            if not exists(str(file)):
                missing_files.append(file)
    if missing_files:
        msg = "The following file were not converted:\n" + "\n".join(missing_files)
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
    # clinical specifications in BIDS
    from os.path import join, realpath, split

    import clinica.iotools.bids_utils as bids
    from clinica.iotools.converters.aibl_to_bids.utils import (
        create_participants_tsv_file,
        create_scans_tsv_file,
        create_sessions_tsv_file,
    )
    from clinica.utils.stream import cprint

    clinical_spec_path = join(
        split(realpath(__file__))[0], "../../data/clinical_specifications"
    )
    # if not exists(clinical_spec_path):
    #    raise FileNotFoundError(
    #        f"{clinical_spec_path} file not found ! This is an internal file of Clinica."
    #    )

    cprint("Creating modality agnostic files...")
    readme_data = {
        "link": "http://adni.loni.usc.edu/study-design/collaborative-studies/aibl/",
        "desc": (
            "The Australian Imaging, Biomarker & Lifestyle Flagship Study of Ageing (AIBL) seeks to discover which "
            "biomarkers, cognitive characteristics, and health and lifestyle factors determine the development of AD. "
            "Although AIBL and ADNI have many of the same goals, there are differences between the two projects."
        ),
    }
    bids.write_modality_agnostic_files(
        study_name="AIBL", readme_data=readme_data, bids_dir=output_dataset
    )

    cprint("Creating participants.tsv...")
    create_participants_tsv_file(
        output_dataset,
        clinical_spec_path,
        input_clinical_data,
        delete_non_bids_info=True,
    )

    cprint("Creating sessions files...")
    create_sessions_tsv_file(output_dataset, input_clinical_data, clinical_spec_path)

    cprint("Creating scans files...")
    create_scans_tsv_file(output_dataset, input_clinical_data, clinical_spec_path)
