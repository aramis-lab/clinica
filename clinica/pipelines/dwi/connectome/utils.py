"""This module contains utilities used by the DWIConnectome pipeline."""

from pathlib import Path
from typing import List

__all__ = [
    "get_luts",
    "get_conversion_luts",
    "get_containers",
    "get_caps_filenames",
    "print_begin_pipeline",
    "print_end_pipeline",
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
        return "a8d561694887a1ca8d9df223aa5ef861b6c79d43ce9ed93835b9ce8aadc331b1"
    if filename == "fs_a2009s.txt":
        return "40b0d4d77bde7e1d265439347af5b30cc973748c1a88d203d7044cb35b3863e1"
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
