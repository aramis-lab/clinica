def statistics_on_atlases(in_registered_map, name_map, prefix_file=None):
    """Computes a list of statistics files for each atlas.

    Args:
        in_registered_map (str): Map already registered on atlases.
        name_map (str): Name of the registered map in CAPS format.
        prefix_file (Opt[str]):
            <prefix_file>_space-<atlas_name>_map-<name_map>_statistics.tsv

    Returns:
        List of paths leading to the statistics TSV files.
    """
    import os

    from nipype.utils.filemanip import split_filename

    from clinica.utils.atlas import (
        AtlasAbstract,
        JHUDTI811mm,
        JHUTracts01mm,
        JHUTracts251mm,
    )
    from clinica.utils.statistics import statistics_on_atlas

    in_atlas_list = [JHUDTI811mm(), JHUTracts01mm(), JHUTracts251mm()]

    atlas_statistics_list = []
    for atlas in in_atlas_list:
        if not isinstance(atlas, AtlasAbstract):
            raise TypeError("Atlas element must be an AtlasAbstract type")

        if prefix_file:
            filename = (
                f"{prefix_file}_space-{atlas.get_name_atlas()}"
                f"_res-{atlas.get_spatial_resolution()}_map-{name_map}_statistics.tsv"
            )
        else:
            _, base, _ = split_filename(in_registered_map)
            filename = (
                f"{base}_space-{atlas.get_name_atlas()}"
                f"_res-{atlas.get_spatial_resolution()}_map-{name_map}_statistics.tsv"
            )

        out_atlas_statistics = os.path.abspath(os.path.join(os.getcwd(), filename))

        statistics_on_atlas(in_registered_map, atlas, out_atlas_statistics)
        atlas_statistics_list.append(out_atlas_statistics)

    return atlas_statistics_list


def get_caps_filenames(caps_dwi_filename: str):
    """Prepare some filenames with CAPS naming convention."""
    import re

    m = re.search(
        r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+).*_dwi_space-[a-zA-Z0-9]+",
        caps_dwi_filename,
    )
    if not m:
        raise ValueError(
            f"Input filename {caps_dwi_filename} is not in a CAPS compliant format."
        )

    caps_prefix = m.group(0)
    bids_source = f"{m.group(1)}_{m.group(2)}_dwi"

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


def print_begin_pipeline(in_bids_or_caps_file):
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    print_begin_image(get_subject_id(in_bids_or_caps_file))


def print_end_pipeline(in_bids_or_caps_file, final_file_1, final_file_2):
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(in_bids_or_caps_file))


def get_ants_transforms(in_affine_transformation, in_bspline_transformation):
    """Combine transformations for antsApplyTransforms interface."""
    return [in_bspline_transformation, in_affine_transformation]
