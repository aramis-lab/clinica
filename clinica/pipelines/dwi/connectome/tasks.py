"""This module contains Nipype tasks used by the connectome pipeline.

Nipype tasks must be 'self-contained' such that only primitive type hints can be
used. The tasks are simple wrappers around a properly implemented Python function.
"""


def convert_flirt_to_mrtrix_transformation_task(
    source_image: str,
    reference_image: str,
    flirt_matrix: str,
    name_output_matrix=None,
) -> str:
    """Adapter for Nipype"""

    from pathlib import Path

    from clinica.pipelines.dwi.connectome.utils import (
        convert_flirt_to_mrtrix_transformation,
    )

    if name_output_matrix is not None:
        name_output_matrix = Path(name_output_matrix)

    return str(
        convert_flirt_to_mrtrix_transformation(
            Path(source_image),
            Path(reference_image),
            Path(flirt_matrix),
            name_output_matrix,
        )
    )
