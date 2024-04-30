"""This module contains Nipype tasks used by the DWIPreprocessingUsingT1 pipeline.

Nipype tasks must be 'self-contained' such that only primitive type hints can be
used. The tasks are simple wrappers around a properly implemented Python function.
"""


def prepare_reference_b0_task(
    dwi_filename: str,
    b_values_filename: str,
    b_vectors_filename: str,
    b_value_threshold: float = 5.0,
    working_directory=None,
) -> tuple:
    """Task called be Nipype to execute 'prepare_reference_b0'."""
    from pathlib import Path

    from clinica.pipelines.dwi.preprocessing.t1.utils import prepare_reference_b0
    from clinica.pipelines.dwi.utils import DWIDataset

    if working_directory:
        working_directory = Path(working_directory)

    b0_reference_filename, reference_dwi_dataset = prepare_reference_b0(
        DWIDataset(
            dwi=dwi_filename,
            b_values=b_values_filename,
            b_vectors=b_vectors_filename,
        ),
        b_value_threshold,
        working_directory,
    )

    return str(b0_reference_filename), *(str(_) for _ in reference_dwi_dataset)


def change_itk_transform_type_task(input_affine_file: str) -> str:
    """Task called be Nipype to execute 'change_itk_transform_type'."""
    from pathlib import Path

    from clinica.pipelines.dwi.preprocessing.t1.utils import change_itk_transform_type

    return str(change_itk_transform_type(Path(input_affine_file)))


def rotate_b_vectors_task(
    b_vectors_filename: str,
    matrix_filenames: list,
    output_dir: str = None,
) -> str:
    """Task called be Nipype to execute 'change_itk_transform_type'."""
    from pathlib import Path

    from clinica.pipelines.dwi.preprocessing.t1.utils import rotate_b_vectors

    return str(
        rotate_b_vectors(
            Path(b_vectors_filename),
            matrix_filenames,
            Path(output_dir) if output_dir else None,
        )
    )


def broadcast_matrix_filename_to_match_b_vector_length_task(
    matrix_filename: str, b_vectors_filename: str
) -> list:
    """Task called be Nipype to execute 'broadcast_matrix_filename_to_match_b_vector_length'."""
    from pathlib import Path

    from clinica.pipelines.dwi.preprocessing.t1.utils import (
        broadcast_matrix_filename_to_match_b_vector_length,
    )

    return [
        str(file_path)
        for file_path in broadcast_matrix_filename_to_match_b_vector_length(
            Path(matrix_filename), Path(b_vectors_filename)
        )
    ]


def rename_into_caps_task(
    dwi_filename: str,
    dwi_preproc_filename: str,
    b_values_preproc_filename: str,
    b_vectors_preproc_filename: str,
    b0_brain_mask_filename: str,
) -> tuple:
    """Task called be Nipype to execute 'rename_into_caps'."""
    from pathlib import Path

    from clinica.pipelines.dwi.preprocessing.t1.utils import rename_into_caps

    return rename_into_caps(
        Path(dwi_filename),
        Path(dwi_preproc_filename),
        Path(b_values_preproc_filename),
        Path(b_vectors_preproc_filename),
        Path(b0_brain_mask_filename),
    )


def merge_nifti_images_in_time_dimension_task(image1: str, image2: str) -> str:
    """Merges the two provided volumes in the time (4th) dimension."""
    from clinica.utils.image import merge_nifti_images_in_time_dimension

    return str(merge_nifti_images_in_time_dimension((image1, image2)))
