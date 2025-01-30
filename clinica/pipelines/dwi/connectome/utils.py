"""This module contains utilities used by the DWIConnectome pipeline."""

from pathlib import Path
from typing import List, Optional

__all__ = [
    "get_luts",
    "get_conversion_luts",
    "get_containers",
    "get_caps_filenames",
    "print_begin_pipeline",
    "print_end_pipeline",
    "convert_flirt_to_mrtrix_transformation",
]


def get_luts() -> List[str]:
    from clinica.utils.check_dependency import get_freesurfer_home

    return [
        str(get_freesurfer_home() / "FreeSurferColorLUT.txt"),
        str(get_freesurfer_home() / "FreeSurferColorLUT.txt"),
    ]


def get_conversion_luts() -> List[str]:
    from pathlib import Path

    path_to_mappings = (
        Path(__file__).resolve().parent.parent.parent.parent / "resources" / "mappings"
    )
    resulting_paths = []
    for filename in ("fs_default.txt", "fs_a2009s.txt"):
        file_path = path_to_mappings / filename
        if not file_path.is_file():
            file_path = _download_mrtrix3_file(filename, path_to_mappings)
        resulting_paths.append(str(file_path))
    return resulting_paths


def _download_mrtrix3_file(filename: str, path_to_mappings: Path) -> Path:
    from clinica.utils.inputs import RemoteFileStructure, fetch_file
    from clinica.utils.stream import cprint

    try:
        return fetch_file(
            RemoteFileStructure(
                filename=filename,
                url="https://raw.githubusercontent.com/MRtrix3/mrtrix3/master/share/mrtrix3/labelconvert/",
                checksum=_get_checksum_for_filename(filename),
            ),
            path_to_mappings,
        )
    except IOError as err:
        error_msg = f"Unable to download required MRTRIX mapping ({filename}) for processing: {err}"
        cprint(msg=error_msg, lvl="error")
        raise IOError(error_msg)


def _get_checksum_for_filename(filename: str) -> str:
    if filename == "fs_default.txt":
        return "bfebee26de22dc4cd03d5ee3f26524b046cce232679e1ba1bc26f18180d491f1"
    if filename == "fs_a2009s.txt":
        return "ae9660f2a9fb44b7d828dcf1f390ce81ed600471810af89042ba011c7a2a675f"
    raise ValueError(f"File name {filename} is not supported.")


def get_containers(subjects: List[str], sessions: List[str]) -> List[str]:
    from pathlib import Path

    return [
        str(Path("subjects") / subject / session / "dwi")
        for subject, session in zip(subjects, sessions)
    ]


def get_caps_filenames(dwi_file: str) -> tuple:
    import re

    error_msg = f"Input filename {dwi_file} is not in a CAPS compliant format."
    if (
        m := re.search(
            r"/(sub-[a-zA-Z0-9]+_ses-[a-zA-Z0-9]+.*_desc-preproc*)_dwi", dwi_file
        )
    ) is None:
        raise ValueError(error_msg)
    source_file_caps = m.group(1)
    if (
        m := re.search(
            r"/(sub-[a-zA-Z0-9]+_ses-[a-zA-Z0-9]+.*)_space-[a-zA-Z0-9]+_desc-preproc_dwi",
            dwi_file,
        )
    ) is None:
        raise ValueError(error_msg)
    source_file_bids = m.group(1)

    response = f"{source_file_caps}_model-CSD_responseFunction.txt"
    fod = f"{source_file_caps}_model-CSD_diffmodel.nii.gz"
    tracts = f"{source_file_caps}_model-CSD_tractography.tck"
    nodes = [
        f"{source_file_caps}_atlas-{atlas}_parcellation.nii.gz"
        for atlas in ("desikan", "destrieux")
    ]
    connectomes = [
        f"{source_file_bids}_model-CSD_atlas-{atlas}_connectivity.tsv"
        for atlas in ("desikan", "destrieux")
    ]

    return response, fod, tracts, nodes, connectomes


def print_begin_pipeline(in_bids_or_caps_file: str) -> None:
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    print_begin_image(get_subject_id(in_bids_or_caps_file))


def print_end_pipeline(in_bids_or_caps_file: str, final_file: str) -> None:
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(in_bids_or_caps_file))


def _check_file_presence(
    file: Path,
) -> None:
    if not file.is_file():
        raise FileNotFoundError(f"The file {file} was not found.")


def _get_transformation_cmd(
    source_image: Path,
    reference_image: Path,
    flirt_matrix: Path,
    mrtrix_matrix: Path,
) -> str:
    for file in (source_image, reference_image, flirt_matrix):
        _check_file_presence(file)

    return f"transformconvert {flirt_matrix} {source_image} {reference_image} flirt_import {mrtrix_matrix}"


def convert_flirt_to_mrtrix_transformation(
    source_image: Path,
    reference_image: Path,
    flirt_matrix: Path,
    name_output_matrix: Optional[str] = None,
) -> Path:
    """Convert flirt matrix to mrtrix matrix.

    This function converts a transformation matrix produced by FSL's flirt
    command into a format usable by MRtrix. The output of this function
    is usually for the mrtransform command.

    Parameters
    ----------
    source_image : Path
        File containing the source image used in FSL flirt with the -in flag.

    reference_image : Path
        File containing the reference image used in FSL flirt with the -ref flag.

    flirt_matrix : Path
        File containing the transformation matrix obtained by FSL flirt.

    name_output_matrix : Path, optional
        Name of the output matrix. Defaults to "mrtrix_matrix.mat".

    Returns
    -------
    mrtrix_matrix : Path
        Transformation matrix in MRtrix format.
    """
    import os

    from clinica.utils.check_dependency import ThirdPartySoftware, check_software

    check_software(ThirdPartySoftware.MRTRIX)

    name_output_matrix = name_output_matrix or "mrtrix_matrix.mat"
    mrtrix_matrix = Path(name_output_matrix)
    mrtrix_matrix = mrtrix_matrix.resolve()

    cmd = _get_transformation_cmd(
        source_image, reference_image, flirt_matrix, mrtrix_matrix
    )
    os.system(cmd)

    return mrtrix_matrix
