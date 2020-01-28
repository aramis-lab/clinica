# coding: utf8


def statistics_on_atlases(in_registered_map, name_map, prefix_file=None):
    """
    Computes a list of statistics files for each atlas.

    Args:
        in_registered_map (str): Map already registered on atlases.
        name_map (str): Name of the registered map in CAPS format.
        prefix_file (Opt[str]):
            <prefix_file>_space-<atlas_name>_map-<name_map>_statistics.tsv

    Returns:
        List of paths leading to the statistics TSV files.
    """
    from os import getcwd
    from os.path import abspath, join
    from nipype.utils.filemanip import split_filename
    from clinica.utils.atlas import (AtlasAbstract, JHUDTI811mm,
                                     JHUTracts01mm, JHUTracts251mm)
    from clinica.utils.statistics import statistics_on_atlas

    in_atlas_list = [JHUDTI811mm(),
                     JHUTracts01mm(), JHUTracts251mm()]

    atlas_statistics_list = []
    for atlas in in_atlas_list:
        if not isinstance(atlas, AtlasAbstract):
            raise TypeError("Atlas element must be an AtlasAbstract type")

        if prefix_file is None:
            _, base, _ = split_filename(in_registered_map)
            filename = '%s_space-%s_res-%s_map-%s_statistics.tsv' % \
                       (base, atlas.get_name_atlas(),
                        atlas.get_spatial_resolution(), name_map)
        else:
            filename = '%s_space-%s_res-%s_map-%s_statistics.tsv' % \
                       (prefix_file, atlas.get_name_atlas(),
                        atlas.get_spatial_resolution(), name_map)

        out_atlas_statistics = abspath(join(getcwd(), filename))

        statistics_on_atlas(in_registered_map, atlas, out_atlas_statistics)
        atlas_statistics_list.append(out_atlas_statistics)

    return atlas_statistics_list


def extract_bids_identifier_from_caps_filename(caps_dwi_filename):
    """Extract BIDS identifier from CAPS filename"""
    import re

    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+).*_dwi',
                  caps_dwi_filename)

    if m is None:
        raise ValueError('Input filename is not in a CAPS compliant format.')
    bids_identifier = m.group(0)

    return bids_identifier


def get_caps_filenames(caps_dwi_filename):
    """Prepare some filenames with CAPS naming convention."""
    import re

    m = re.search(
        r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+).*_dwi_space-[a-zA-Z0-9]+',
        caps_dwi_filename)
    if m is None:
        raise ValueError('Input filename is not in a CAPS compliant format.')

    caps_prefix = m.group(0)
    bids_source = m.group(1) + '_' + m.group(2) + '_dwi'

    out_dti = caps_prefix + '_model-DTI_diffmodel.nii.gz'
    out_fa = caps_prefix + '_FA.nii.gz'
    out_md = caps_prefix + '_MD.nii.gz'
    out_ad = caps_prefix + '_AD.nii.gz'
    out_rd = caps_prefix + '_RD.nii.gz'

    return bids_source, out_dti, out_fa, out_md, out_ad, out_rd


def rename_into_caps(in_caps_dwi,
                     in_norm_fa, in_norm_md, in_norm_ad, in_norm_rd,
                     in_b_spline_transform, in_affine_matrix):
    """Rename some outputs of the pipelines into CAPS format."""
    from nipype.interfaces.utility import Rename
    from clinica.pipelines.dwi_dti.dwi_dti_utils import extract_bids_identifier_from_caps_filename

    bids_identifier = extract_bids_identifier_from_caps_filename(in_caps_dwi)

    # CAPS normalized FA
    rename_fa = Rename()
    rename_fa.inputs.in_file = in_norm_fa
    rename_fa.inputs.format_string = bids_identifier + "_space-MNI152Lin_res-1x1x1_FA.nii.gz"
    out_caps_fa = rename_fa.run()
    # CAPS normalized MD
    rename_md = Rename()
    rename_md.inputs.in_file = in_norm_md
    rename_md.inputs.format_string = bids_identifier + "_space-MNI152Lin_res-1x1x1_MD.nii.gz"
    out_caps_md = rename_md.run()
    # CAPS normalized AD
    rename_ad = Rename()
    rename_ad.inputs.in_file = in_norm_ad
    rename_ad.inputs.format_string = bids_identifier + "_space-MNI152Lin_res-1x1x1_AD.nii.gz"
    out_caps_ad = rename_ad.run()
    # CAPS normalized RD
    rename_rd = Rename()
    rename_rd.inputs.in_file = in_norm_rd
    rename_rd.inputs.format_string = bids_identifier + "_space-MNI152Lin_res-1x1x1_RD.nii.gz"
    out_caps_rd = rename_rd.run()
    # CAPS B-spline transform
    rename_b_spline = Rename()
    rename_b_spline.inputs.in_file = in_b_spline_transform
    rename_b_spline.inputs.format_string = bids_identifier + "_space-MNI152Lin_res-1x1x1_deformation.nii.gz"
    out_caps_b_spline_transform = rename_b_spline.run()
    # CAPS Affine matrix
    rename_affine = Rename()
    rename_affine.inputs.in_file = in_affine_matrix
    rename_affine.inputs.format_string = bids_identifier + "_space-MNI152Lin_res-1x1x1_affine.mat"
    out_caps_affine_matrix = rename_affine.run()

    return out_caps_fa.outputs.out_file, out_caps_md.outputs.out_file,\
        out_caps_ad.outputs.out_file, out_caps_rd.outputs.out_file, \
        out_caps_b_spline_transform.outputs.out_file, \
        out_caps_affine_matrix.outputs.out_file


def print_begin_pipeline(in_bids_or_caps_file):
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    print_begin_image(get_subject_id(in_bids_or_caps_file))


def print_end_pipeline(in_bids_or_caps_file, final_file_1, final_file_2):
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(in_bids_or_caps_file))


def apply_ants_registration_syn_quick_transformation(
        in_image,
        in_reference_image,
        in_affine_transformation,
        in_bspline_transformation,
        name_output_image=None):
    """
    Apply a transformation obtained with antsRegistrationSyNQuick.sh.

    This function applies a rigid & deformable B-Spline syn transformation
    which has been estimated previously with antsRegistrationSyNQuick script.

    Args:
        in_image (str): File containing the input image to be transformed.
        in_reference_image (str): File defining the spacing, origin, size,
            and direction of the output warped image.
        in_affine_transformation (str): File containing the transformation
            matrix obtained by antsRegistrationSyNQuick (expected file:
            [Prefix]0GenericAffine.mat).
        in_bspline_transformation (str): File containing the transformation
            matrix obtained by antsRegistrationSyNQuick (expected file:
            [Prefix]1Warp.nii.gz).
        name_output_image (Optional[str]): Name of the output image
            (default=deformed_image.nii.gz).

    Returns:
        out_deformed_image (str): File containing the deformed image according
            to in_affine_transformation and in_bspline_transformation
            transformations.
    """
    import os
    from clinica.utils.check_dependency import check_ants
    check_ants()

    assert(os.path.isfile(in_image))
    assert(os.path.isfile(in_affine_transformation))
    assert(os.path.isfile(in_bspline_transformation))

    if name_output_image is None:
        out_deformed_image = os.path.abspath('deformed_image.nii.gz')
    else:
        out_deformed_image = os.path.abspath(name_output_image)

    cmd = 'antsApplyTransforms -d 3 -e 0 -i %s -o %s -t %s -t %s -r %s ' \
          '--interpolation Linear' \
          % (in_image, out_deformed_image, in_bspline_transformation,
             in_affine_transformation, in_reference_image)
    os.system(cmd)

    return out_deformed_image
