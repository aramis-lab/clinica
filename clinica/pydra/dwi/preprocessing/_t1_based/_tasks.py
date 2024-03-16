from os import PathLike

from pydra.mark import annotate, task


@task
@annotate({"return": {"output_image": PathLike}})
def merge_nifti_images_in_time_dimension_task(
    image1: PathLike, image2: PathLike
) -> Path:
    from clinica.utils.image import merge_nifti_images_in_time_dimension

    return merge_nifti_images_in_time_dimension((image1, image2))
