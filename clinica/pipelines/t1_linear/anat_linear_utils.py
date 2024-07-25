from pathlib import Path
from typing import Optional


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
    import ants
    from ants.utils.bias_correction import n4_bias_field_correction

    from clinica.utils.stream import cprint

    output_prefix = output_prefix or ""
    input_image_ants = ants.image_read(str(input_image))
    bias_corrected_image = n4_bias_field_correction(
        input_image_ants, spline_param=bspline_fitting_distance, verbose=verbose
    )
    if save_bias:
        bias_image = n4_bias_field_correction(
            input_image_ants,
            spline_param=bspline_fitting_distance,
            return_bias_field=True,
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
