"""This module contains Nipype tasks used by the DWIPreprocessingUsingPhaseDiff pipeline.

Nipype tasks must be 'self-contained' such that only primitive type hints can be
used. The tasks are simple wrappers around a properly implemented Python function.
"""


def rename_into_caps_task(
    dwi_filename: str,
    dwi_preproc_filename: str,
    b_values_preproc_filename: str,
    b_vectors_preproc_filename: str,
    b0_brain_mask_filename: str,
    calibrated_magnitude_image_filename: str,
    calibrated_field_map_image_filename: str,
    calibrated_smoothed_field_map_image_filename: str,
) -> tuple:
    from pathlib import Path

    from clinica.pipelines.dwi.preprocessing.fmap.utils import rename_into_caps

    return rename_into_caps(
        Path(dwi_filename),
        Path(dwi_preproc_filename),
        Path(b_values_preproc_filename),
        Path(b_vectors_preproc_filename),
        Path(b0_brain_mask_filename),
        Path(calibrated_magnitude_image_filename),
        Path(calibrated_field_map_image_filename),
        Path(calibrated_smoothed_field_map_image_filename),
    )


def convert_phase_difference_to_hertz_task(
    phase_diff_filename: str,
    delta_echo_time: float,
    working_dir: str = None,
) -> str:
    """Wrapper for 'convert_phase_difference_to_hertz' to be used by Nipype."""
    from pathlib import Path

    from clinica.pipelines.dwi.preprocessing.fmap.utils import (
        convert_phase_difference_to_hertz,
    )

    if working_dir:
        working_dir = Path(working_dir)
    return str(
        convert_phase_difference_to_hertz(
            Path(phase_diff_filename), delta_echo_time, working_dir
        )
    )


def demean_image_task(
    input_image: str, mask: str = None, working_dir: str = None
) -> str:
    """Wrapper for 'demean_image' to be used by Nipype."""
    from pathlib import Path

    from clinica.pipelines.dwi.preprocessing.fmap.utils import demean_image

    if working_dir:
        working_dir = Path(working_dir)
    if mask:
        mask = Path(mask)
    return str(demean_image(Path(input_image), mask, working_dir))


def convert_phase_difference_to_rads_task(
    phase_diff_filename: str, working_dir: str = None
) -> str:
    """Wrapper for 'convert_phase_difference_to_rads' to be used by Nipype."""
    from pathlib import Path

    from clinica.pipelines.dwi.preprocessing.fmap.utils import (
        convert_phase_difference_to_rads,
    )

    if working_dir:
        working_dir = Path(working_dir)
    return str(convert_phase_difference_to_rads(Path(phase_diff_filename), working_dir))
