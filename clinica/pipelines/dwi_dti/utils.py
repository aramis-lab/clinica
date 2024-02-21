def statistics_on_atlases(
    in_registered_map: str, name_map: str, dwi_preprocessed_file: str
) -> list:
    """Computes a list of statistics files for each atlas.

    Parameters
    ----------
    in_registered_map : str
        Map already registered on atlases.

    name_map : str
        Name of the registered map in CAPS format.

    dwi_preprocessed_file : str
        The preprocessed DWI file name which contains the entities to be
        used for building the statistics file names.

    Returns
    -------
    list of str :
        List of paths leading to the statistics TSV files.
    """
    from pathlib import Path

    from clinica.utils.atlas import (
        JHUDTI811mm,
        JHUTracts01mm,
        JHUTracts251mm,
    )
    from clinica.utils.bids import BIDSFileName
    from clinica.utils.statistics import statistics_on_atlas

    atlas_statistics_list = []
    for atlas in (JHUDTI811mm(), JHUTracts01mm(), JHUTracts251mm()):
        source = BIDSFileName.from_name(dwi_preprocessed_file)
        source.update_entity("space", atlas.get_name_atlas())
        source.update_entity("res", atlas.get_spatial_resolution())
        source.update_entity("map", name_map)
        source.suffix = "statistics"
        source.extension = ".tsv"
        out_atlas_statistics = str((Path.cwd() / source.name).resolve())
        statistics_on_atlas(in_registered_map, atlas, out_atlas_statistics)
        atlas_statistics_list.append(out_atlas_statistics)

    return atlas_statistics_list


def get_caps_filenames(caps_dwi_filename: str):
    """Prepare some filenames with CAPS naming convention."""
    import re

    m = re.search(
        r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+).*_space-[a-zA-Z0-9]+_desc-preproc",
        caps_dwi_filename,
    )
    if not m:
        raise ValueError(
            f"Input filename {caps_dwi_filename} is not in a CAPS compliant format."
        )

    caps_prefix = m.group(0)
    bids_source = f"{m.group(1)}_{m.group(2)}"

    out_dti = f"{caps_prefix}_model-DTI_diffmodel.nii.gz"
    out_fa = f"{caps_prefix}_FA.nii.gz"
    out_md = f"{caps_prefix}_MD.nii.gz"
    out_ad = f"{caps_prefix}_AD.nii.gz"
    out_rd = f"{caps_prefix}_RD.nii.gz"
    out_evec = f"{caps_prefix}_DECFA.nii.gz"

    return bids_source, out_dti, out_fa, out_md, out_ad, out_rd, out_evec


def rename_into_caps(
    in_caps_dwi: str,
    in_norm_fa: str,
    in_norm_md: str,
    in_norm_ad: str,
    in_norm_rd: str,
    in_b_spline_transform: str,
    in_affine_matrix: str,
) -> tuple:
    """Rename different outputs of the pipelines into CAPS format.

    Parameters
    ----------
    in_caps_dwi : str
    in_norm_fa : str
    in_norm_md : str
    in_norm_ad : str
    in_norm_rd : str
    in_b_spline_transform : str
    in_affine_matrix : str

    Returns
    -------
    tuple :
        The different outputs with CAPS naming convention
    """
    from clinica.utils.dwi import rename_files

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