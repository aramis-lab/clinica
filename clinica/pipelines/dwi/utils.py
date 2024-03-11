"""This module contains utilities common to all DWI pipelines."""

from collections import namedtuple
from pathlib import Path
from typing import Dict, Tuple

__all__ = ["DWIDataset", "rename_files"]


DWIDataset = namedtuple("DWIDataset", "dwi b_values b_vectors")


def rename_files(in_caps_dwi: Path, mapping: Dict[Path, str]) -> Tuple[str, ...]:
    """Rename files provided.

    The new files are symbolic links to old files.
    For this reason, the old files still exists after renaming.

    Parameters
    ----------
    in_caps_dwi : Path
        A DWI file from the CAPS folder.
        This is used only to extract the BIDS identifier.

    mapping : dict<Path, str>
        Mapping between original file paths and suffixes for
        new file names.

    Returns
    -------
    tuple :
        New file names.
    """
    from nipype.interfaces.utility import Rename

    bids_id = _extract_bids_identifier_from_filename(str(in_caps_dwi))
    renamed_files = []
    for original_file, suffix in mapping.items():
        rename = Rename()
        rename.inputs.in_file = str(original_file)
        rename.inputs.format_string = str(original_file.parent / f"{bids_id}{suffix}")
        renamed_files.append(rename.run().outputs.out_file)

    return tuple(renamed_files)


def _extract_bids_identifier_from_filename(caps_dwi_filename: str) -> str:
    """Extract BIDS identifier from a DWI CAPS filename.

    Parameters
    ----------
    caps_dwi_filename : str
        DWI file name for which to extract the bids identifier.

    Returns
    -------
    str :
        The corresponding BIDS identifier.

    Examples
    --------
    >>> _extract_bids_identifier_from_filename("sub-01_ses-M000_dwi_space-b0_preproc.bval")
    'sub-01_ses-M000'
    >>> _extract_bids_identifier_from_filename("sub-01_ses-M000_dwi.bvec")
    'sub-01_ses-M000'
    >>> _extract_bids_identifier_from_filename("foo/bar/sub-01_ses-M000_dwi_baz.foo.bar")
    'sub-01_ses-M000'
    >>> _extract_bids_identifier_from_filename("foo/bar/sub-01_ses-M000_space-b0_desc-preproc_dwi.bval")
    'sub-01_ses-M000_space-b0_desc-preproc'
    """
    import re

    m = re.search(r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+).*_dwi", caps_dwi_filename)
    if not m:
        raise ValueError(
            f"Could not extract the BIDS identifier from the DWI input filename {caps_dwi_filename}."
        )
    identifier = m.group(0).rstrip("_dwi")

    return identifier
