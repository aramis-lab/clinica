import os

from nibabel.nifti1 import Nifti1Image
from pydra.mark import annotate, task


@task
@annotate({"return": {"transforms_list": list}})
def concatenate_transforms_task(
    pet_to_t1w_tranform: str, t1w_to_mni_tranform: str
) -> list:
    """Pydra task for concatenating two input transformation files into a list.

    .. note::
        Please refer to the documentation of function
        `clinica.pipelines.pet_linear.pet_linear_utils.concatenate_transforms`.

    """
    from clinica.pipelines.pet.linear.utils import concatenate_transforms

    return concatenate_transforms(t1w_to_mni_tranform, pet_to_t1w_tranform)


@task
@annotate({"return": {"output_img": str}})
def suvr_normalization_task(
    input_img: os.PathLike,
    norm_img: os.PathLike,
    ref_mask: os.PathLike,
) -> str:
    """Pydra task for normalizing the input image according to the reference region.

    .. note::
        Please refer to the documentation of function
        `clinica.pipelines.pet_linear.pet_linear_utils.suvr_normalization`.

    """
    from pathlib import Path

    from clinica.pipelines.pet.linear.utils import perform_suvr_normalization

    return str(
        perform_suvr_normalization(Path(input_img), Path(norm_img), Path(ref_mask))
    )


@task
@annotate({"return": {"output_img": os.PathLike}})
def crop_nifti_task(input_img: os.PathLike, output_dir: os.PathLike) -> os.PathLike:
    """Pydra task for cropping the input image based on the reference.

    .. note::
        Please refer to the documentation of function
        `clinica.pipelines.pet_linear.pet_linear_utils.crop_nifti`.

    """
    from pathlib import Path

    from clinica.utils.image import crop_nifti

    return crop_nifti(Path(input_img), Path(output_dir))
