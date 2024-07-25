from os import PathLike
from typing import Iterable, Optional, Union

__all__ = ["center_nifti"]


def center_nifti(
    bids_directory: Union[str, PathLike],
    output_bids_directory: Union[str, PathLike],
    modalities: Optional[Iterable[str]] = None,
    center_all_files: bool = False,
    overwrite_existing_files: bool = False,
):
    """Center NIfTI files in a BIDS dataset.

    Parameters
    ----------
    bids_directory : PathLike
        The path to the BIDS directory.

    output_bids_directory : PathLike
        The path to

    modalities : Iterable of str, optional
        The modalities

    center_all_files : bool, optional
        Whether to center all file or not.
        Default=False.

    overwrite_existing_files : bool, optional
        If True and if the output BIDS directory already contain files,
        they might be overwritten. If False, the output BIDS has to be empty
        or non-existing otherwise a ClinicaExistingDatasetError will be raised.

    Notes
    -----
    This tool is mainly useful as a preprocessing step of SPM. In some cases, SPM is not able to segment T1 volumes
    because their respective center is not aligned with the origin of the world coordinate system. By default, only
    images detected as problematic are converted from INPUT_BIDS_DIRECTORY to OUTPUT_BIDS_DIRECTORY, whilst the others
    are copied verbatim.
    """
    import time
    from pathlib import Path

    from clinica.iotools.utils.data_handling import (
        center_all_nifti,
        write_list_of_files,
    )
    from clinica.utils.exceptions import ClinicaExistingDatasetError
    from clinica.utils.stream import cprint, log_and_raise

    bids_directory = Path(bids_directory)
    output_bids_directory = Path(output_bids_directory)
    if output_bids_directory.exists():
        files = [
            file.name
            for file in output_bids_directory.iterdir()
            if not file.name.startswith(".")
        ]
        if files and not overwrite_existing_files:
            raise ClinicaExistingDatasetError(output_bids_directory)
    cprint("Clinica is now centering the requested images.", lvl="info")

    centered_files = center_all_nifti(
        bids_directory,
        output_bids_directory,
        modalities,
        center_all_files,
    )
    # Write list of created files
    timestamp = time.strftime("%Y%m%d-%H%M%S", time.localtime(time.time()))
    log_file = output_bids_directory / f"centered_nifti_list_{timestamp}.txt"
    if not write_list_of_files(centered_files, log_file):
        log_and_raise(f"Could not create log file {log_file}.", IOError)

    cprint(
        f"{len(centered_files)} NIfTI files/images of BIDS folder:\n"
        f"\t{bids_directory}\n"
        f"for the modalities {modalities} have been centered in output folder:\n"
        f"\t{output_bids_directory}",
        lvl="info",
    )
    if log_file.is_file():
        cprint(
            f"The list of centered NIfTI files is available here: {log_file}.",
            lvl="info",
        )
    cprint(
        "Please note that the rest of the input BIDS folder has also been copied to the output folder.",
        lvl="info",
    )
