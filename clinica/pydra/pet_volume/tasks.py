import typing as ty
from pathlib import PurePath

from pydra.mark import annotate, task


@task
@annotate({"return": {"suvr_pet_path": PurePath}})
def normalize_to_reference_task(pet_image: PurePath, region_mask: PurePath) -> PurePath:
    """Pydra task for normalizing the provided `pet_image`.

    .. note::
        Please refer to the documentation of function
        `clinica.pipelines.pet_volume.pet_volume_utils.normalize_to_reference`.
    """
    from clinica.pipelines.pet_volume.pet_volume_utils import normalize_to_reference

    return normalize_to_reference(pet_image, region_mask)


@task
@annotate({"return": {"out_mask": PurePath}})
def create_binary_mask_task(
    tissues: ty.List[PurePath],
    threshold: float = 0.3,
) -> PurePath:
    """Pydra task for creating a binary mask Nifti1Image from the list of tissues.

    .. note::
        Please refer to the documentation of function
        `clinica.pipelines.pet_volume.pet_volume_utils.create_binary_mask`.
    """
    from clinica.pipelines.pet_volume.pet_volume_utils import create_binary_mask

    return create_binary_mask(tissues, threshold)


@task
@annotate({"return": {"masked_image_path": PurePath}})
def apply_binary_mask_task(image: PurePath, binary_mask: PurePath) -> PurePath:
    """Pydra task for applying the provided `binary_mask` to the provided `image`.

    .. note::
        Please refer to the documentation of function
        `clinica.pipelines.pet_volume.pet_volume_utils.apply_binary_mask`.
    """
    from clinica.pipelines.pet_volume.pet_volume_utils import apply_binary_mask

    return apply_binary_mask(image, binary_mask)


@task
@annotate({"return": {"atlas_statistics": ty.List}})
def atlas_statistics_task(in_image: str, in_atlas_list: ty.List) -> ty.List:
    """Pydra task for generating regional measures from atlas_list in TSV files.

    .. note::
        Please refer to the documentation of function
        `clinica.pipelines.pet_volume.pet_volume_utils.atlas_statistics`.
    """
    from clinica.pipelines.pet_volume.pet_volume_utils import atlas_statistics

    return atlas_statistics(in_image, in_atlas_list)


@task
@annotate({"return": {"out_mask": PurePath}})
def create_pvc_mask_task(tissues: ty.List) -> PurePath:
    """Pydra task for creating a pvc mask from tissue list.

    .. note::
        Please refer to the documentation of function
        `clinica.pipelines.pet_volume.pet_volume_utils.create_pvc_mask`.
    """
    from clinica.pipelines.pet_volume.pet_volume_utils import create_pvc_mask

    return create_pvc_mask(tissues)


@task
@annotate({"return": {"pet_pvc_path": str}})
def pet_pvc_name_task(pet_image: PurePath, pvc_method: str) -> str:
    """Pydra task for building the name for the PET PVC interface.

    .. note::
        Please refer to the documentation of function
        `clinica.pipelines.pet_volume.pet_volume_utils.pet_pvc_name`.
    """
    from clinica.pipelines.pet_volume.pet_volume_utils import pet_pvc_name

    return pet_pvc_name(pet_image, pvc_method)


@task
@annotate({"return": {"psf_x": int, "psf_y": int, "psf_z": int}})
def get_psf_task(
    pvc_psf_tsv: PurePath,
    pet_filename: PurePath,
    pet_tracer: str,
) -> ty.Tuple[int]:
    """Returns the tuple (psf_x, psf_y, psf_z) for the subject at the
    session currently treated by the workflow.

    .. note::
        This function uses the name of the pet scan file to retrieve
        the subject and session IDs which are needed by the function
        `read_psf_information`.

    Parameters
    ----------
    pvc_psf_tsv : PurePath
        Path to the custom pvc_psf TSV file.
        See docstring of function `read_psf_information` to get more
        details on this file format.

    pet_filename : PurePath
        Path to the PET image for the subject and session treated by
        the workflow. See note above for details.

    pet_tracer : str
        Label of the PET tracer for which to retrieve psf.

    Returns
    -------
    (psf_x, psf_y, psf_z) : Tuple[int]
        The length 3 tuple for the subject session.
    """
    import re

    from clinica.utils.pet import read_psf_information

    m = re.search(r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)", str(pet_filename.name))
    if not m:
        raise ValueError(
            f"Could not extract subject and session IDs from filename {pet_filename}"
        )
    sub = m.group(1)
    ses = m.group(2)
    psf_x, psf_y, psf_z = read_psf_information(pvc_psf_tsv, [sub], [ses], pet_tracer)[0]
    return psf_x, psf_y, psf_z
