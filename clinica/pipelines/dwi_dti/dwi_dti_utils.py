# coding: utf8


def erode_mask(in_mask, npass=6, nthreads=2):
    """
    Erode the mask.

    This function is used for the erosion of the brain mask thus preventing
    voxels near the edge of the brain from potentially being erroneously
    selected as single-fibre voxels during the estimation of the response
    function.

    Args:
        in_mask (str): Binary mask.
        npass (Optional[int]): Number of times to repeatedly apply the erosion
            (default=6 which is used in MRtrix community).
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 1 disables multi-threading).

    Returns:
        An eroded mask.
    """
    import os.path as op
    import os

    assert(op.isfile(in_mask))

    out_eroded_mask = op.abspath('eroded_mask.nii.gz')
    cmd = 'maskfilter -npass %s -nthreads %s %s erode %s' \
          % (npass, nthreads, in_mask, out_eroded_mask)
    os.system(cmd)
    return out_eroded_mask


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


def dwi_container_from_filename(dwi_filename):
    import re
    from os.path import join
    m = re.search(r'(sub-[a-zA-Z0-9]+)/(ses-[a-zA-Z0-9]+)', dwi_filename)

    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format.'
            ' It does not contain the subject and session information.')

    subject = m.group(1)
    session = m.group(2)

    return join('subjects', subject, session)


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
    """
    Prepare some filenames with CAPS naming convention
    """
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
    """
    Rename different outputs of the pipelines into CAPS format.

    Returns:
        The different outputs with CAPS naming convention
    """
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
    """
    """
    from clinica.utils.stream import cprint
    import re
    import datetime
    from colorama import Fore

    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)',
                  in_bids_or_caps_file)
    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format.')
    now = datetime.datetime.now().strftime('%H:%M:%S')

    cprint('%s[%s]%s Running pipeline for %s...' % (
        Fore.BLUE, now, Fore.RESET, m.group(0)))


def print_end_pipeline(in_bids_or_caps_file, final_file_1, final_file_2):
    """
    """
    from clinica.utils.stream import cprint
    import re
    import datetime
    from colorama import Fore

    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)',
                  in_bids_or_caps_file)
    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format.')
    now = datetime.datetime.now().strftime('%H:%M:%S')

    cprint('%s[%s]%s ...%s has completed.' % (
        Fore.GREEN, now, Fore.RESET, m.group(0)))
