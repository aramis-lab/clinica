"""This module contains Nipype tasks used in several pipelines."""


def crop_nifti_task(input_image: str, output_path: str) -> str:
    from pathlib import Path

    from clinica.utils.image import crop_nifti

    return str(crop_nifti(Path(input_image), Path(output_path)))


def get_filename_no_ext_task(filename: str) -> str:
    from pathlib import Path

    from clinica.utils.filemanip import get_filename_no_ext

    return get_filename_no_ext(Path(filename))


def run_n4biasfieldcorrection_task(
    input_image: str,
    bspline_fitting_distance: int,
    output_prefix=None,
    output_dir=None,
    save_bias=False,
    verbose=False,
) -> str:
    from pathlib import Path

    from clinica.pipelines.utils import run_n4biasfieldcorrection

    if output_dir:
        output_dir = Path(output_dir)

    return str(
        run_n4biasfieldcorrection(
            Path(input_image),
            bspline_fitting_distance,
            output_prefix,
            output_dir,
            save_bias,
            verbose,
        )
    )


def run_ants_registration_synquick_task(
    fixed_image: str,
    moving_image: str,
    random_seed: int,
    transform_type: str,
    output_prefix=None,
    output_dir=None,
    verbose: bool = False,
    return_inverse_transform: bool = False,
) -> tuple:
    from pathlib import Path

    from clinica.pipelines.utils import run_ants_registration_synquick

    if output_dir:
        output_dir = Path(output_dir)

    (
        warped_image_output_path,
        transformation_matrix_output_path,
        transformation_matrix_inverse_output_path,
    ) = run_ants_registration_synquick(
        Path(fixed_image),
        Path(moving_image),
        random_seed,
        transform_type,
        output_prefix,
        output_dir,
        verbose=verbose,
    )
    if return_inverse_transform:
        return (
            str(warped_image_output_path),
            str(transformation_matrix_output_path),
            str(transformation_matrix_inverse_output_path),
        )
    return str(warped_image_output_path), str(transformation_matrix_output_path)


def run_ants_registration_task(
    fixed_image: str,
    moving_image: str,
    random_seed: int,
    transform_type: str,
    output_prefix=None,
    output_dir=None,
    verbose: bool = False,
    shrink_factors=None,
    smoothing_sigmas=None,
    number_of_iterations=None,
    return_inverse_transform: bool = False,
) -> tuple:
    from pathlib import Path

    from clinica.pipelines.utils import run_ants_registration

    if output_dir:
        output_dir = Path(output_dir)

    (
        warped_image_output_path,
        transformation_matrix_output_path,
        transformation_matrix_inverse_output_path,
    ) = run_ants_registration(
        Path(fixed_image),
        Path(moving_image),
        random_seed,
        transform_type,
        output_prefix,
        output_dir,
        verbose=verbose,
        shrink_factors=shrink_factors,
        smoothing_sigmas=smoothing_sigmas,
        number_of_iterations=number_of_iterations,
    )

    if return_inverse_transform:
        return (
            str(warped_image_output_path),
            str(transformation_matrix_output_path),
            str(transformation_matrix_inverse_output_path),
        )
    return str(warped_image_output_path), str(transformation_matrix_output_path)


def run_ants_apply_transforms_task(
    reference_image: str,
    input_image: str,
    transforms: list,
    output_dir=None,
):
    from pathlib import Path

    from clinica.pipelines.utils import run_ants_apply_transforms

    if output_dir:
        output_dir = Path(output_dir)

    return str(
        run_ants_apply_transforms(
            Path(reference_image),
            Path(input_image),
            transforms,
            output_dir=output_dir,
        )
    )
