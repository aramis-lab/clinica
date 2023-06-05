"""This module contains some functions to be used together with some Nipype functions.

In particular, this module currently contains functions used for Nipype DataSink.
"""


def container_from_filename(bids_or_caps_filename: str) -> str:
    """Extract container from BIDS or CAPS file.

    Parameters
    ----------
    bids_or_caps_filename : str
        Full path to BIDS or CAPS filename.

    Returns
    -------
    str :
        Container path of the form "subjects/<participant_id>/<session_id>".

    Examples
    --------
    >>> from clinica.utils.nipype import container_from_filename
    >>> container_from_filename('/path/to/bids/sub-CLNC01/ses-M000/anat/sub-CLNC01_ses-M000_T1w.nii.gz')
    'subjects/sub-CLNC01/ses-M000'
    >>> container_from_filename('caps/subjects/sub-CLNC01/ses-M000/dwi/preprocessing/sub-CLNC01_ses-M000_preproc.nii')
    'subjects/sub-CLNC01/ses-M000'
    """
    import os
    import re

    m = re.search(r"(sub-[a-zA-Z0-9]+)/(ses-[a-zA-Z0-9]+)", bids_or_caps_filename)
    if not m:
        raise ValueError(
            f"Input filename {bids_or_caps_filename} is not in a BIDS or CAPS compliant format."
            "It does not contain the participant and session ID."
        )
    subject = m.group(1)
    session = m.group(2)
    return os.path.join("subjects", subject, session)


def fix_join(path, *paths):
    """Fix joined path.

    This workaround function is used in pipelines like DWIPreprocessing* or PETVolume. In the workflow.connect part,
    you can use some function that are used as string, causing an import error
    """
    import os

    return os.path.join(path, *paths)
