def perform_suvr_normalization_task(
    pet_image_path: str,
    normalizing_image_path: str,
    reference_mask_path: str,
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.linear.utils import perform_suvr_normalization

    return str(
        perform_suvr_normalization(
            Path(pet_image_path),
            Path(normalizing_image_path),
            Path(reference_mask_path),
        )
    )


def remove_mni_background_task(
    input_image: str,
    mni_mask_path: str,
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.linear.utils import remove_mni_background

    return str(
        remove_mni_background(
            Path(input_image),
            Path(mni_mask_path),
        )
    )


def clip_task(
    input_pet: str,
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.linear.utils import clip_img

    return str(clip_img(Path(input_pet)))


def rename_into_caps_task(
    pet_bids_image_filename: str,
    pet_preprocessed_image_filename: str,
    pet_to_mri_transformation_filename: str,
    suvr_reference_region: str,
    uncropped_image: bool,
    pet_filename_in_t1w_raw: str = None,
    output_dir: str = None,
) -> tuple:
    from pathlib import Path

    from clinica.pipelines.pet.linear.utils import rename_into_caps

    if pet_filename_in_t1w_raw:
        pet_filename_in_t1w_raw = Path(pet_filename_in_t1w_raw)
    if output_dir:
        output_dir = Path(output_dir)
    (
        pet_filename_caps,
        transformation_filename_caps,
        pet_filename_in_t1w_caps,
    ) = rename_into_caps(
        Path(pet_bids_image_filename),
        Path(pet_preprocessed_image_filename),
        Path(pet_to_mri_transformation_filename),
        suvr_reference_region,
        uncropped_image,
        pet_filename_in_t1w_raw,
        output_dir,
    )
    if pet_filename_in_t1w_caps:
        pet_filename_in_t1w_caps = str(pet_filename_in_t1w_caps)
    return (
        str(pet_filename_caps),
        str(transformation_filename_caps),
        pet_filename_in_t1w_caps,
    )
