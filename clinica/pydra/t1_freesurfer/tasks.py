import typing as ty

import pydra


@pydra.mark.task
@pydra.mark.annotate({"return": {"subject_id": "str"}})
def compute_freesurfer_subject_id(participant_id: str, session_id: ty.Optional[str]):
    """Compute FreeSurfer's subject ID from BIDS.

    Unlike BIDS, FreeSurfer's layout does not encode the concept of sessions.

    The latter is emulated using a compound subject ID, like `sub-P01_ses-M00`
    for BIDS participant `sub-P01` at imaging session `ses-M00`.

    If the BIDS dataset is single-session,
    then use the participant ID as FreeSurfer's subject ID.

    Parameters
    ----------
    participant_id : str
        BIDS subject identifier.
    session_id : str, optional
        BIDS session identifier.

    Returns
    -------
    str
        FreeSurfer's subject identifier.
    """
    return f"{participant_id}_{session_id}" if session_id else participant_id
