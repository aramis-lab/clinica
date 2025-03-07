"""This module contains utilities used by the DWIDTI pipeline."""

from pathlib import Path
from typing import Dict, List, Tuple

from clinica.utils.dwi import DTIBasedMeasure

__all__ = [
    "compute_statistics_on_atlases",
    "get_caps_filenames",
    "rename_into_caps",
    "print_begin_pipeline",
    "print_end_pipeline",
    "get_ants_transforms",
]


def compute_statistics_on_atlases(
    registered_map: Path, name_map: str, dwi_preprocessed_file: Path
) -> List[Path]:
    """Computes a list of statistics files for each atlas.

    Parameters
    ----------
    registered_map : Path
        Map already registered on atlases.

    name_map : str
        Name of the registered map in CAPS format.

    dwi_preprocessed_file : Path
        The preprocessed DWI file name which contains the entities to be
        used for building the statistics file names.

    Returns
    -------
    list of str :
        List of paths leading to the statistics TSV files.
    """
    from pathlib import Path

    from clinica.bids import BIDSFileName
    from clinica.utils.atlas import atlas_factory
    from clinica.utils.statistics import statistics_on_atlas

    atlas_statistics_list = []
    for atlas_name in ("JHUDTI81", "JHUTracts0", "JHUTracts25"):
        atlas = atlas_factory(atlas_name)
        source = BIDSFileName.from_name(dwi_preprocessed_file)
        source.update_entity("space", atlas.name)
        source.update_entity("res", atlas.spatial_resolution)
        source.update_entity("map", name_map)
        source.suffix = "statistics"
        source.extension = ".tsv"
        out_atlas_statistics = (Path.cwd() / source.name).resolve()
        statistics_on_atlas(registered_map, atlas, out_atlas_statistics)
        atlas_statistics_list.append(out_atlas_statistics)

    return atlas_statistics_list


def get_caps_filenames(caps_dwi_filename: Path) -> Dict[str, str]:
    """Prepare some filenames with CAPS naming convention."""
    import re

    m = re.search(
        r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+).*_space-[a-zA-Z0-9]+_desc-preproc",
        str(caps_dwi_filename),
    )
    if not m:
        raise ValueError(
            f"Input filename {caps_dwi_filename} is not in a CAPS compliant format."
        )

    caps_prefix = m.group(0)
    results = {
        "bids_source": f"{m.group(1)}_{m.group(2)}",
        "diffmodel": f"{caps_prefix}_model-DTI_diffmodel.nii.gz",
    }
    results.update(
        {
            measure.value: f"{caps_prefix}_{measure.value}.nii.gz"
            for measure in DTIBasedMeasure
        }
    )
    results["DECFA"] = f"{caps_prefix}_DECFA.nii.gz"

    return results


def rename_into_caps(
    in_caps_dwi: Path,
    in_norm_fa: Path,
    in_norm_md: Path,
    in_norm_ad: Path,
    in_norm_rd: Path,
    in_b_spline_transform: Path,
    in_affine_matrix: Path,
) -> Tuple[str, ...]:
    """Rename different outputs of the pipelines into CAPS format.

    Parameters
    ----------
    in_caps_dwi : Path
    in_norm_fa : Path
    in_norm_md : Path
    in_norm_ad : Path
    in_norm_rd : Path
    in_b_spline_transform : Path
    in_affine_matrix : Path

    Returns
    -------
    tuple of str :
        The different outputs with CAPS naming convention
    """
    from clinica.pipelines.dwi.utils import rename_files

    return rename_files(
        in_caps_dwi,
        {
            in_norm_fa: "_space-MNI152Lin_res-1x1x1_FA.nii.gz",
            in_norm_md: "_space-MNI152Lin_res-1x1x1_MD.nii.gz",
            in_norm_ad: "_space-MNI152Lin_res-1x1x1_AD.nii.gz",
            in_norm_rd: "_space-MNI152Lin_res-1x1x1_RD.nii.gz",
            in_b_spline_transform: "_space-MNI152Lin_res-1x1x1_deformation.nii.gz",
            in_affine_matrix: "_space-MNI152Lin_res-1x1x1_affine.mat",
        },
    )


def print_begin_pipeline(in_bids_or_caps_file: str):
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    print_begin_image(get_subject_id(in_bids_or_caps_file))


def print_end_pipeline(in_bids_or_caps_file: str, final_file_1: str, final_file_2: str):
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(in_bids_or_caps_file))


def get_ants_transforms(
    in_affine_transformation: str, in_bspline_transformation: str
) -> list:
    """Combine transformations for antsApplyTransforms interface."""
    return [in_bspline_transformation, in_affine_transformation]
