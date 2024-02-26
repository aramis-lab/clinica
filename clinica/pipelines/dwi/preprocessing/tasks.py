"""This module contains Nipype tasks used in multiple DWI preprocessing pipelines.

Nipype tasks must be 'self-contained' such that only primitive type hints can be
used. The tasks are simple wrappers around a properly implemented Python function.
"""


def generate_index_file_task(
    b_values_filename: str,
    image_id=None,
    output_dir=None,
) -> str:
    """Wrapper around 'generate_index_file' for Nipype."""
    from pathlib import Path

    from .utils import generate_index_file

    return str(
        generate_index_file(
            Path(b_values_filename),
            image_id,
            Path(output_dir) if output_dir else None,
        )
    )


def generate_acq_file_task(
    dwi_filename: str,
    fsl_phase_encoding_direction: str,
    total_readout_time: str,
    image_id=None,
    output_dir=None,
) -> str:
    """Wrapper around 'generate_acq_file' for Nipype."""
    from pathlib import Path

    from .utils import generate_acq_file

    return str(
        generate_acq_file(
            Path(dwi_filename),
            fsl_phase_encoding_direction,
            total_readout_time,
            image_id,
            Path(output_dir) if output_dir else None,
        )
    )


def compute_average_b0_task(
    dwi_filename: str,
    b_value_filename: str,
    b_value_threshold: float = 5.0,
    squeeze: bool = False,
    out_file: str = None,
) -> str:
    """Nipype task for 'compute_average_b0'."""
    from .utils import check_file, compute_average_b0

    return str(
        compute_average_b0(
            check_file(dwi_filename),
            check_file(b_value_filename),
            b_value_threshold,
            squeeze,
            out_file,
        )
    )
