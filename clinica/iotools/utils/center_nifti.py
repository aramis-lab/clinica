from os import PathLike
from typing import Iterable, Optional, Union

__all__ = ["center_nifti"]


def center_nifti(
    bids_directory: Union[str, PathLike],
    output_bids_directory: Union[str, PathLike],
    modalities: Optional[Iterable[str]] = None,
    centering_threshold: int = 50,
    overwrite_existing_files: bool = False,
):
    """Center NIfTI files in a BIDS dataset.

    Parameters
    ----------
    bids_directory : PathLike
        The path to the BIDS directory.

    output_bids_directory : PathLike
        The path to the output BIDS directory that contains the centered images.

    modalities : Iterable of str, optional
        The modalities that will be processed.

    centering_threshold : int, optional
        Threshold above which images are centered. See the documentation for the reason behind the default.
        Default=50 (mm).

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
    from clinica.utils.stream import cprint, log_and_raise

    cprint(
        f"Clinica is now centering the requested images with a threshold of {centering_threshold} mm.",
        lvl="info",
    )

    centered_files = center_all_nifti(
        bids_directory,
        output_bids_directory,
        modalities,
        centering_threshold,
        overwrite_existing_files,
    )

    cprint(
        f"{len(centered_files)} NIfTI images of BIDS folder:\n"
        f"\t{bids_directory}\n"
        f"have been centered in the output folder:\n"
        f"\t{output_bids_directory}\n"
        f"for modalities:\n"
        f"\t{modalities}\n"
        f"Please note that the rest of the input BIDS folder has also been copied to the output folder.",
        lvl="info",
    )

    # Write list of created files
    timestamp = time.strftime("%Y%m%d-%H%M%S", time.localtime(time.time()))
    log_file = Path(output_bids_directory) / f"centered_nifti_list_{timestamp}.txt"
    if not write_list_of_files(
        [f.relative_to(output_bids_directory) for f in centered_files], log_file
    ):
        log_and_raise(f"Could not create log file {log_file}.", IOError)

    if log_file.is_file():
        cprint(
            f"The list of centered NIfTI files is available here: {log_file}.",
            lvl="info",
        )
