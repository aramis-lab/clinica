from pathlib import Path
from typing import Optional, Tuple


def get_substitutions_datasink_flair(bids_image_id: str) -> list:
    from clinica.pipelines.t1_linear.anat_linear_utils import (  # noqa
        _get_substitutions_datasink,
    )

    return _get_substitutions_datasink(bids_image_id, "FLAIR")


def get_substitutions_datasink_t1_linear(bids_image_id: str) -> list:
    from clinica.pipelines.t1_linear.anat_linear_utils import (  # noqa
        _get_substitutions_datasink,
    )

    return _get_substitutions_datasink(bids_image_id, "T1w")


def _get_substitutions_datasink(bids_image_id: str, suffix: str) -> list:
    """Return file name substitutions for renaming.

    Parameters
    ----------
    bids_image_id : str
        This is the original image BIDS file name without the extension.
        This will be used to get all the BIDS entities that shouldn't
        be modified (subject, session...).

    suffix : str
        The suffix to use for the new file.

    Returns
    -------
    substitutions : List of tuples of str
        List of length 3 containing the substitutions to perform.
    """
    if bids_image_id.endswith(f"_{suffix}"):
        bids_image_id_without_suffix = bids_image_id.removesuffix(f"_{suffix}")
    else:
        raise ValueError(
            f"bids image ID {bids_image_id} should end with provided {suffix}."
        )
    return [
        (
            f"{bids_image_id}Warped_cropped.nii.gz",
            f"{bids_image_id_without_suffix}_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_{suffix}.nii.gz",
        ),
        (
            f"{bids_image_id}0GenericAffine.mat",
            f"{bids_image_id_without_suffix}_space-MNI152NLin2009cSym_res-1x1x1_affine.mat",
        ),
        (
            f"{bids_image_id}Warped.nii.gz",
            f"{bids_image_id_without_suffix}_space-MNI152NLin2009cSym_res-1x1x1_{suffix}.nii.gz",
        ),
    ]


def print_end_pipeline(anat, final_file):
    """Display end message for <subject_id> when <final_file> is connected."""
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(anat))


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


def run_ants_registration(
    fixed_image: Path,
    moving_image: Path,
    random_seed: int,
    output_prefix: Optional[str] = None,
    output_dir: Optional[Path] = None,
    verbose: bool = False,
) -> Tuple[Path, Path]:
    """Run antsRegistration using antsPy.

    Parameters
    ----------
    fixed_image : Path
        The path to the fixed image.

    moving_image : Path
        The path to the moving image.

    random_seed : int
        The random seed to be used.

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
        fixed_image, moving_image, random_seed, verbose=verbose
    )
    try:
        warped_image = registration_results["warpedmovout"]
        transformation_matrix = registration_results["fwdtransforms"][-1]
    except (KeyError, IndexError):
        msg = (
            "Something went wrong when calling antsRegistration with the following parameters :\n"
            f"- fixed_image = {fixed_image}\n- moving_image = {moving_image}\n"
            f"- random_seed = {random_seed}\n- type_of_transformation='antsRegistrationSyN[a]'\n"
        )
        log_and_raise(msg, RuntimeError)

    return _write_ants_registration_results(
        warped_image, transformation_matrix, output_prefix or "", output_dir
    )


def _call_ants_registration(
    fixed_image: Path,
    moving_image: Path,
    random_seed: int,
    verbose: bool = False,
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
    return ants.registration(
        ants.image_read(str(fixed_image)),
        ants.image_read(str(moving_image)),
        type_of_transformation="antsRegistrationSyN[a]",
        random_seed=random_seed,
        verbose=verbose,
    )


def _write_ants_registration_results(
    warped_image,
    transformation_matrix,
    output_prefix: str,
    output_dir: Optional[Path] = None,
) -> Tuple[Path, Path]:
    import shutil

    import ants

    from clinica.utils.stream import cprint

    warped_image_output_path = (
        output_dir or Path.cwd()
    ) / f"{output_prefix}Warped.nii.gz"
    transformation_matrix_output_path = (
        output_dir or Path.cwd()
    ) / f"{output_prefix}0GenericAffine.mat"
    cprint(f"Writing warped image to {warped_image_output_path}.", lvl="debug")
    ants.image_write(warped_image, str(warped_image_output_path))
    shutil.copy(transformation_matrix, transformation_matrix_output_path)

    return warped_image_output_path, transformation_matrix_output_path
