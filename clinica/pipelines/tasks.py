"""This module contains Nipype tasks used in several pipelines."""


def crop_nifti_task(input_image: str, reference_image: str) -> str:
    from pathlib import Path

    from clinica.utils.image import crop_nifti

    return str(crop_nifti(Path(input_image), Path(reference_image)))
