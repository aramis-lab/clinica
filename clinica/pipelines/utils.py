from enum import Enum
from pathlib import Path
from typing import List, Optional, Tuple, Union

__all__ = [
    "AntsRegistrationTransformType",
    "AntsRegistrationSynQuickTransformType",
    "run_n4biasfieldcorrection",
    "run_ants_registration",
    "run_ants_registration_synquick",
    "run_ants_apply_transforms",
]


class AntsRegistrationSynQuickTransformType(str, Enum):
    """The possible values for the transform type of AntsRegistrationSynQuick."""

    TRANSLATION = "antsRegistrationSyN[t]"
    RIGID = "antsRegistrationSyN[r]"
    SIMILARITY = "antsRegistrationSyN[s]"
    AFFINE = "antsRegistrationSyN[a]"


class AntsRegistrationTransformType(str, Enum):
    """The possible values for the transform type of AntsRegistration."""

    TRANSLATION = "Translation"
    RIGID = "Rigid"
    SIMILARITY = "Similarity"
    AFFINE = "Affine"
    SYN = "SyN"


def run_n4biasfieldcorrection(
    input_image: Path,
    bspline_fitting_distance: int,
    output_prefix: Optional[str] = None,
    output_dir: Optional[Path] = None,
    save_bias: bool = False,
    verbose: bool = False,
) -> Path:
    """Run n4biasfieldcorrection using antsPy.

    Parameters
    ----------
    input_image : Path
        The path to the input image.

    bspline_fitting_distance : int
        This is the 'spline_param' of n4biasfieldcorrection.

    output_prefix : str, optional
        The prefix to be put at the beginning of the output file names.
        Ex: 'sub-XXX_ses-MYYY'.

    output_dir : Path, optional
        The directory in which to write the output files.
        If not provided, these files will be written in the current directory.

    save_bias : bool, optional
        Whether to save the bias image or not.
        If set to True, the bias image is not returned but saved in the
        provided output_dir with a name of the form '{output_prefix}_bias_image.nii.gz'.
        Default=False.

    verbose : bool, optional
        Control the verbose mode of n4biasfieldcorrection. Set to True can be
        useful for debugging.
        Default=False.

    Returns
    -------
    bias_corrected_output_path : Path
        The path to the bias corrected image.
    """
    from clinica.utils.exceptions import ClinicaMissingDependencyError
    from clinica.utils.stream import cprint, log_and_raise

    try:
        import ants
    except ImportError:
        log_and_raise(
            "The package 'antsPy' is required to run antsRegistration in Python.",
            ClinicaMissingDependencyError,
        )

    output_prefix = output_prefix or ""
    bias_corrected_image = _call_n4_bias_field_correction(
        input_image, bspline_fitting_distance, save_bias=False, verbose=verbose
    )
    if save_bias:
        bias_image = _call_n4_bias_field_correction(
            input_image,
            bspline_fitting_distance,
            save_bias=True,
            verbose=verbose,
        )
        bias_output_path = (
            output_dir or Path.cwd()
        ) / f"{output_prefix}_bias_image.nii.gz"
        ants.image_write(bias_image, str(bias_output_path))
        cprint(f"Writing bias image to {bias_output_path}.", lvl="debug")
    bias_corrected_output_path = (
        output_dir or Path.cwd()
    ) / f"{output_prefix}_bias_corrected_image.nii.gz"
    cprint(
        f"Writing bias corrected image to {bias_corrected_output_path}.", lvl="debug"
    )
    ants.image_write(bias_corrected_image, str(bias_corrected_output_path))

    return bias_corrected_output_path


def _call_n4_bias_field_correction(
    input_image: Path,
    bspline_fitting_distance: int,
    save_bias: bool = False,
    verbose: bool = False,
) -> Path:
    import ants
    from ants.utils.bias_correction import n4_bias_field_correction

    return n4_bias_field_correction(
        ants.image_read(str(input_image)),
        spline_param=bspline_fitting_distance,
        return_bias_field=save_bias,
        verbose=verbose,
    )


def run_ants_registration_synquick(
    fixed_image: Path,
    moving_image: Path,
    random_seed: int,
    transform_type: Union[str, AntsRegistrationSynQuickTransformType],
    output_prefix: Optional[str] = None,
    output_dir: Optional[Path] = None,
    verbose: bool = False,
) -> Tuple[Path, Path, Path]:
    transform_type = AntsRegistrationSynQuickTransformType(transform_type)
    return _run_ants_registration(
        fixed_image,
        moving_image,
        random_seed,
        transform_type,
        output_prefix,
        output_dir,
        verbose,
    )


def run_ants_registration(
    fixed_image: Path,
    moving_image: Path,
    random_seed: int,
    transform_type: Union[str, AntsRegistrationTransformType],
    output_prefix: Optional[str] = None,
    output_dir: Optional[Path] = None,
    verbose: bool = False,
    shrink_factors: Optional[Tuple[int, ...]] = None,
    smoothing_sigmas: Optional[Tuple[int, ...]] = None,
    number_of_iterations: Optional[Tuple[int, ...]] = None,
) -> Tuple[Path, Path, Path]:
    transform_type = AntsRegistrationTransformType(transform_type)
    return _run_ants_registration(
        fixed_image,
        moving_image,
        random_seed,
        transform_type,
        output_prefix,
        output_dir,
        verbose,
        shrink_factors=shrink_factors,
        smoothing_sigmas=smoothing_sigmas,
        number_of_iterations=number_of_iterations,
    )


def _run_ants_registration(
    fixed_image: Path,
    moving_image: Path,
    random_seed: int,
    transform_type: Union[
        AntsRegistrationTransformType, AntsRegistrationSynQuickTransformType
    ],
    output_prefix: Optional[str] = None,
    output_dir: Optional[Path] = None,
    verbose: bool = False,
    shrink_factors: Optional[Tuple[int, ...]] = None,
    smoothing_sigmas: Optional[Tuple[int, ...]] = None,
    number_of_iterations: Optional[Tuple[int, ...]] = None,
) -> Tuple[Path, Path, Path]:
    """Run antsRegistration using antsPy.

    Parameters
    ----------
    fixed_image : Path
        The path to the fixed image.

    moving_image : Path
        The path to the moving image.

    random_seed : int
        The random seed to be used.

    transform_type : AntsRegistrationTransformType or AntsRegistrationSynQuickTransformType
        The type of transformation to be applied.

    output_prefix : str, optional
        The prefix to be put at the beginning of the output file names.
        Ex: 'sub-XXX_ses-MYYY'.

    output_dir : Path, optional
        The directory in which to write the output files.
        If not provided, these files will be written in the current directory.

    verbose : bool, optional
        Control the verbose mode of antsRegistration. Set to True can be
        useful for debugging.
        Default=False.

    Returns
    -------
    warped_image_output_path : Path
        The path to the warped nifti image generated by antsRegistration.

    transformation_matrix_output_path : Path
        The path to the transforms to move from moving to fixed image.
        This is a .mat file.

    Raises
    ------
    RuntimeError :
        If results cannot be extracted.
    """
    from clinica.utils.stream import log_and_raise

    registration_results = _call_ants_registration(
        fixed_image,
        moving_image,
        random_seed,
        transform_type,
        verbose=verbose,
        shrink_factors=shrink_factors,
        smoothing_sigmas=smoothing_sigmas,
        number_of_iterations=number_of_iterations,
    )
    try:
        warped_image = registration_results["warpedmovout"]
        transformation_matrix = registration_results["fwdtransforms"][-1]
        transformation_matrix_inverse = registration_results["invtransforms"][0]
    except (KeyError, IndexError):
        msg = (
            "Something went wrong when calling antsRegistration with the following parameters :\n"
            f"- fixed_image = {fixed_image}\n- moving_image = {moving_image}\n"
            f"- random_seed = {random_seed}\n- type_of_transformation='{transform_type.value}'\n"
        )
        log_and_raise(msg, RuntimeError)

    return _write_ants_registration_results(
        warped_image,
        transformation_matrix,
        transformation_matrix_inverse,
        output_prefix or "",
        output_dir,
    )


def _call_ants_registration(
    fixed_image: Path,
    moving_image: Path,
    random_seed: int,
    transform_type: Union[
        AntsRegistrationTransformType, AntsRegistrationSynQuickTransformType
    ],
    verbose: bool = False,
    shrink_factors: Optional[Tuple[int, ...]] = None,
    smoothing_sigmas: Optional[Tuple[int, ...]] = None,
    number_of_iterations: Optional[Tuple[int, ...]] = None,
) -> dict:
    from clinica.utils.exceptions import ClinicaMissingDependencyError
    from clinica.utils.stream import log_and_raise

    try:
        import ants
    except ImportError:
        log_and_raise(
            "The package 'antsPy' is required to run antsRegistration in Python.",
            ClinicaMissingDependencyError,
        )
    kwargs = {}
    if shrink_factors is not None:
        kwargs["aff_shrink_factors"] = shrink_factors
    if smoothing_sigmas is not None:
        kwargs["aff_smoothing_sigmas"] = smoothing_sigmas
    if number_of_iterations is not None:
        kwargs["aff_iterations"] = number_of_iterations

    return ants.registration(
        ants.image_read(str(fixed_image)),
        ants.image_read(str(moving_image)),
        type_of_transformation=transform_type.value,
        random_seed=random_seed,
        verbose=verbose,
        **kwargs,
    )


def _write_ants_registration_results(
    warped_image,
    transformation_matrix,
    transformation_matrix_inverse,
    output_prefix: str,
    output_dir: Optional[Path] = None,
) -> Tuple[Path, Path, Path]:
    import shutil

    import ants

    from clinica.utils.stream import cprint

    warped_image_output_path = (
        output_dir or Path.cwd()
    ) / f"{output_prefix}Warped.nii.gz"
    transformation_matrix_output_path = (
        output_dir or Path.cwd()
    ) / f"{output_prefix}0GenericAffine.mat"
    transformation_matrix_inverse_output_path = (
        output_dir or Path.cwd()
    ) / f"{output_prefix}inverse.mat"
    cprint(f"Writing warped image to {warped_image_output_path}.", lvl="debug")
    ants.image_write(warped_image, str(warped_image_output_path))
    shutil.copy(transformation_matrix, transformation_matrix_output_path)
    shutil.copy(
        transformation_matrix_inverse, transformation_matrix_inverse_output_path
    )

    return (
        warped_image_output_path,
        transformation_matrix_output_path,
        transformation_matrix_inverse_output_path,
    )


def run_ants_apply_transforms(
    fixed_image: Path,
    moving_image: Path,
    transformlist: List[str],
    output_dir: Optional[Path] = None,
) -> Path:
    import ants

    from clinica.utils.stream import cprint

    transformed_image = ants.apply_transforms(
        ants.image_read(str(fixed_image)),
        ants.image_read(str(moving_image)),
        transformlist=transformlist,
    )
    transformed_image_output_path = (output_dir or Path.cwd()) / "transformed.nii.gz"
    cprint(
        f"Writing transformed image to {transformed_image_output_path}.", lvl="debug"
    )
    ants.image_write(transformed_image, str(transformed_image_output_path))

    return transformed_image_output_path
