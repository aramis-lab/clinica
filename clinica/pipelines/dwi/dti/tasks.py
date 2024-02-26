"""This module contains Nipype tasks used by the DWIDTI pipeline.

Nipype tasks must be 'self-contained' such that only primitive type hints can be
used. The tasks are simple wrappers around a properly implemented Python function.
"""


def compute_statistics_on_atlases_task(
    registered_map: str, name_map: str, dwi_preprocessed_file: str
) -> list:
    from pathlib import Path

    from clinica.pipelines.dwi.dti.utils import compute_statistics_on_atlases

    return compute_statistics_on_atlases(
        Path(registered_map),
        name_map,
        Path(dwi_preprocessed_file),
    )


def get_caps_filenames_task(caps_dwi_filename: str) -> tuple:
    from pathlib import Path

    from clinica.pipelines.dwi.dti.utils import get_caps_filenames

    names = get_caps_filenames(Path(caps_dwi_filename))

    return (
        names["bids_source"],
        names["diffmodel"],
        names["FA"],
        names["MD"],
        names["AD"],
        names["RD"],
        names["DECFA"],
    )


def rename_into_caps_task(
    in_caps_dwi: str,
    in_norm_fa: str,
    in_norm_md: str,
    in_norm_ad: str,
    in_norm_rd: str,
    in_b_spline_transform: str,
    in_affine_matrix: str,
) -> tuple:
    from pathlib import Path

    from clinica.pipelines.dwi.dti.utils import rename_into_caps

    return rename_into_caps(
        Path(in_caps_dwi),
        Path(in_norm_fa),
        Path(in_norm_md),
        Path(in_norm_ad),
        Path(in_norm_rd),
        Path(in_b_spline_transform),
        Path(in_affine_matrix),
    )
