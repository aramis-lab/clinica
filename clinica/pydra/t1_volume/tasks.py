from os import PathLike

from pydra.mark import annotate, task


@task
@annotate({"return": {"out_status": None}})
def task_volume_location_in_world_coordinate_system(
    nifti_input: str,
    bids_dir: str,
    modality: str = "t1w",
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

    Returns
    -------
    bool
        True if they are centered, False otherwise

    Warns
    ------
    If volume is not centered on origin of the world coordinate system
    """
    from pathlib import Path

    from clinica.iotools.utils.data_handling import (
        are_images_centered_around_origin_of_world_coordinate_system,
    )

    nifti_list = [Path(nifti_input)]

    return are_images_centered_around_origin_of_world_coordinate_system(
        nifti_list, Path(bids_dir), modality
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


@task
@annotate({"return": {"image_files": list}})
def task_prepare_dartel_input_images(
    nifti_input: PathLike,
) -> bool:
    """Creates dartel images data structure for SPMCommand input

    Parameters
    ----
    nifti_input: PathLike
        Path to nifti file

    Returns
    -------
    list
        data structure with the right format for DARTELExistingTemplate(SPMCommand)
    """
    from clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_utils import (
        prepare_dartel_input_images_pydra,
    )

    return prepare_dartel_input_images_pydra(nifti_input)


@task
@annotate({"return": {"iteration_parameters": list}})
def task_create_iteration_parameters(
    template_input: PathLike,
) -> bool:
    """Creates iteration data structure for SPMCommand input

    Parameters
    ----
    template_input: PathLike
        Path to the template file

    Returns
    -------
    list
        data structure with the right format for DARTELExistingTemplate(SPMCommand) iteration parameter
    """
    from clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_utils import (
        create_iteration_parameters,
    )

    return create_iteration_parameters(template_input, None)
