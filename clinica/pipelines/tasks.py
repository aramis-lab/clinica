"""This module contains Nipype tasks used in several pipelines."""


def crop_nifti_task(input_image: str, output_path: str) -> str:
    from pathlib import Path

    from clinica.utils.image import crop_nifti

    return str(crop_nifti(Path(input_image), Path(output_path)))


def get_filename_no_ext_task(filename: str) -> str:
    from pathlib import Path

    from clinica.utils.filemanip import get_filename_no_ext

    return get_filename_no_ext(Path(filename))
