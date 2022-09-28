from os import PathLike

from pydra.mark import annotate, task


@task
@annotate({"return": {"out_status": None}})
def task_volume_location_in_world_coordinate_system(
    nifti_input: PathLike,
    bids_dir: PathLike,
    modality: str = "t1w",
    skip_question: bool = False,
) -> bool:
    """Check if images are centered around the origin of the world coordinate

    Parameters
    ----
    nifti_input: PathLike
        Path to nifti file

    bids_dir: PathLike
        path to bids directory associated with this check

    modality: str, optional
        the modality of the image. Default="t1w".

    skip_question: bool, optional
        if True, assume answer is yes. Default=False.

    Returns
    -------
    bool
        True if they are centered, False otherwise

    Warns
    ------
    If volume is not centered on origin of the world coordinate system
    """
    from clinica.iotools.utils.data_handling import (
        check_volume_location_in_world_coordinate_system,
    )

    nifti_list = [str(nifti_input)]

    return check_volume_location_in_world_coordinate_system(
        nifti_list, bids_dir, modality, skip_question
    )


@task
@annotate({"return": {"out_list": list}})
def wrap_list(in_list: list):
    """Wraps a list inside another list to comply with Pydra Nipype1Task requirements

    Parameters
    ----------
    in_list: list
        The input list to wrap

    Returns
    -------
    list
        The wrapped list
    """
    return [in_list]
