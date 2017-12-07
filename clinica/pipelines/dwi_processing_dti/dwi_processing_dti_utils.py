# coding: utf8


def convert_nifti_to_mrtrix_format(in_dwi_nii, in_bvals, in_bvecs, nthreads=2):
    """
    Convert to mrtrix image format.

    This function converts the diffusion images into a non-compressed format
    (not strictly necessary, but will make subsequent processing faster), embed
    the diffusion gradient encoding information within the image header,
    re-arrange the data strides to make volume data contiguous in memory
    for each voxel, and convert to floating-point representation (makes data
    access faster in subsequent commands).

    Args:
        in_dwi_nii (str): File containing DWI dataset in NIfTI format.
        in_bvals (str): File containing B-Value table in FSL format.
        in_bvecs (str): File containing Diffusion Gradient table in FSL format.
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 1 disables multi-threading).

    Returns:
        out_dwi_mif (str): DWI dataset in MRtrix format (the diffusion gradient
            encoding information is within the image header).
    """
    import os.path as op
    import os
    from nipype.utils.filemanip import split_filename

    assert(op.isfile(in_dwi_nii))
    assert(op.isfile(in_bvals))
    assert(op.isfile(in_bvecs))

    _, source_file_dwi, _ = split_filename(in_dwi_nii)

    out_dwi_mif = op.abspath(source_file_dwi + '.mif')
    cmd = 'mrconvert %s %s -fslgrad %s %s -datatype float32 -stride 0,0,0,1 ' \
          '-nthreads %s' \
          % (in_dwi_nii, out_dwi_mif, in_bvecs, in_bvals, nthreads)
    os.system(cmd)
    return out_dwi_mif


def compute_maximum_harmonic_order(in_bvecs):
    """
    Compute the maximum harmonic order according to the b-vectors.

    This function estimates the maximum harmonic order according to the number
    of distinct DW directions used in the acquisition is greater than the
    number of parameters that need to be estimated. A table of required
    parameters can be found here:
    https://jdtournier.github.io/mrtrix-0.2/tractography/preprocess.html

    .. warning:: If the number of distinct directions is greater than 91,
        the maximum harmonic order is automatically set to 12.

    Args:
        in_bvecs (str): Diffusion Gradient table in FSL format.

    Returns:
        out_lmax (int): The maximum harmonic order.
    """
    import os
    import warnings

    # TODO: Check if b=0 is present, nb_of_directions=nb_of_directions-1;
    number_of_directions = os.system("$(sort " + in_bvecs + " | uniq | wc -l)")
    out_lmax = 0
    if 6 < number_of_directions <= 15:
        out_lmax = 2
    elif 15 < number_of_directions <= 28:
        out_lmax = 4
    elif 28 < number_of_directions <= 45:
        out_lmax = 6
    elif 45 < number_of_directions <= 66:
        out_lmax = 8
    else:
        warnings.warn("The number of distinct directions is > to 66, lmax is "
                      "automatically set to 8", UserWarning)
        out_lmax = 8

    return out_lmax


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


def dwi_to_tensor(in_dwi_mif, in_b0_mask, nthreads=2, prefix_file=None):
    """
    Perform diffusion tensor estimation from DWI dataset.

    Diffusion (kurtosis) tensor estimation using iteratively reweighted linear
    least squares estimator.

    Args:
        in_dwi_mif (str): DWI dataset in MRtrix image format.
        in_b0_mask (str): Binary mask of the b0 image. Only perform computation
            within this specified binary brain mask image.
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 1 disables multi-threading).

    Returns:
        out_dti (str):
    """
    import os.path as op
    import os
    from nipype.utils.filemanip import split_filename

    assert(op.isfile(in_dwi_mif))
    assert(op.isfile(in_b0_mask))

    _, source_file_dwi, _ = split_filename(in_dwi_mif)

    if prefix_file is None:
        out_dti = op.abspath(source_file_dwi + '_model-dti_diffmodel.nii.gz')
    else:
        out_dti = op.abspath(prefix_file + '_model-dti_diffmodel.nii.gz')

    cmd = 'dwi2tensor -mask %s %s %s -nthreads %s' \
          % (in_b0_mask, in_dwi_mif, out_dti, nthreads)
    os.system(cmd)

    return out_dti


def tensor_to_metrics(in_dti, in_b0_mask, nthreads=2, prefix_file=None):
    """
    Generate maps of tensor-derived parameters.

    Args:
        in_dti (str): Tensor.
        in_b0_mask (str): Binary mask of the b0 image. Only perform computation
            within this specified binary brain mask image.
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 1 disables multi-threading).
        prefix_file (Optional[str]): prefix_name for outputs (Output:
            <prefix_file>_[fa|md|ad|rd|].nii.gz)

    Returns:
        out_fa (str): The tensor-derived parameter fractional anisotropy.
        out_md (str): The tensor-derived parameter mean diffusivity (also
            called mean apparent diffusion).
        out_ad (str): The tensor-derived parameter axial diffusivity.
        out_rd (str): The tensor-derived parameter radial diffusivity.
        out_ev (str): The directionally-encoded colour (DEC) fractional
            anisotropy map This corresponds to the
            first eigenvector modulated by the FA.
    """
    import os.path as op
    import os
    from nipype.utils.filemanip import split_filename

    assert(op.isfile(in_dti))
    assert(op.isfile(in_b0_mask))

    if prefix_file is None:
        out_fa = op.abspath('fa_map_from_dti.nii.gz')
    else:
        out_fa = op.abspath(prefix_file + '_fa.nii.gz')
    cmd = 'tensor2metric %s -nthreads %s -fa %s' % \
          (in_dti, nthreads, out_fa)
    os.system(cmd)

    if prefix_file is None:
        out_md = op.abspath('md_map_from_dti.nii.gz')
    else:
        out_md = op.abspath(prefix_file + '_md.nii.gz')
    cmd = 'tensor2metric %s -nthreads %s -adc %s' % \
          (in_dti, nthreads, out_md)
    os.system(cmd)

    if prefix_file is None:
        out_ad = op.abspath('ad_map_from_dti.nii.gz')
    else:
        out_ad = op.abspath(prefix_file + '_ad.nii.gz')
    cmd = 'tensor2metric %s -nthreads %s -ad %s' % \
          (in_dti, nthreads, out_ad)
    os.system(cmd)

    if prefix_file is None:
        out_rd = op.abspath('rd_map_from_dti.nii.gz')
    else:
        out_rd = op.abspath(prefix_file + '_rd.nii.gz')
    cmd = 'tensor2metric %s -nthreads %s -rd %s' % \
          (in_dti, nthreads, out_rd)
    os.system(cmd)

    if prefix_file is None:
        out_ev = op.abspath('dec_fa_map_from_dti.nii.gz')
    else:
        out_ev = op.abspath(prefix_file + '_decfa.nii.gz')
    cmd = 'tensor2metric %s -nthreads %s -vector %s' % \
          (in_dti, nthreads, out_ev)
    os.system(cmd)

    return out_fa, out_ad, out_md, out_rd, out_ev


def estimate_response(in_dwi_mif, in_b0_mask, lmax=None, algorithm='tax',
                      tmpdir=None, nthreads=2):
    """
    Estimate response function(s) for spherical deconvolution.

    This function generates an appropriate response function (chosen by
    `algorithm`) from the image data for spherical deconvolution.

    Args:
        in_dwi_mif (str): DWI dataset in MRtrix image format.
        in_b0_mask (str): Binary mask of the b0 image. Only perform computation
            within this specified binary brain mask image. It is advised to
            erode the mask (we assume that an eroded mask is given)
        algorithm (Optional[str]): Used algorithm to derive the response
            function. Currenly, choices are:
            - tax (default): Use the Tax et al. (2014) recursive calibration
                algorithm
            - tournier: Use the Tournier et al. (2013) iterative RF selection
                algorithm
        lmax (Optional[int]): The maximum harmonic degree(s) of response
            function estimation
        tmpdir (Optional[str]): Path where the temporary results are stored.
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 1 disables multi-threading).

    Returns:
        out_response_function (str): Text file containing response function
            coefficients.
    """
    import os.path as op
    import os
    import tempfile

    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()

    assert(op.isfile(in_dwi_mif))
    assert(op.isfile(in_b0_mask))

    if algorithm == 'tax':
        out_response_function = op.abspath('out_response_function_tax.txt')
    elif algorithm == 'tournier':
        out_response_function = op.abspath(
            'out_response_function_tournier.txt')
    else:
        raise ValueError(
            'Invalid choice of algorithm in estimate_response function')

    cmd = 'dwi2response ' + algorithm + ' -mask ' + in_b0_mask
    if lmax is not None:
        cmd = cmd + ' -lmax ' + str(lmax)
    cmd = '%s %s %s -nthreads %s -tempdir %s' % \
          (cmd, in_dwi_mif, out_response_function, nthreads, tmpdir)
    os.system(cmd)

    return out_response_function


def estimate_fod(in_dwi_mif, in_b0_mask, in_response_function_coefficients,
                 lmax=None, nthreads=2):
    """
    Estimate FOD.

    This function performs non-negativity constrained spherical deconvolution
    in order to estimate the fiber orientation distribution.

    Args:
        in_dwi_mif (str): DWI dataset in MRtrix image format.
        in_b0_mask (str): Binary mask of the b0 image. Only perform computation
            within this specified binary brain mask image.
        in_response_function_coefficients (str): Text file containing the
            diffusion-weighted signal response function coefficients for a
            single fibre population,
        lmax (int): Maximum harmonic order for the output series. By default,
            the program will use the highest possible `lmax` given the number
            of DWI, up to a maximum of 8.
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 1 disables multi-threading).

    Returns:
        out_sh_coefficients_image (str): Spherical harmonics coefficients
            image.

    References:
     Tournier, J.-D.; Calamante, F. & Connelly, A. Robust determination of the
     fibre orientation distribution in diffusion MRI: Non-negativity
     constrained super-resolved spherical deconvolution. NeuroImage, 2007, 35,
     1459-1472

     Tournier, J.-D.; Calamante, F., Gadian, D.G. & Connelly, A. Direct
     estimation of the fiber orientation density function from
     diffusion-weighted MRI data using spherical deconvolution.NeuroImage,
     2004, 23, 1176-1185
    """
    import os.path as op
    import os

    assert(op.isfile(in_dwi_mif))
    assert(op.isfile(in_b0_mask))
#    assert(op.isfile(in_response_function_coefficients))

    out_sh_coefficients_image = op.abspath('sh_coefficients_image.mif')

    cmd = 'dwi2fod csd -mask %s %s %s %s' % \
          (in_b0_mask, in_dwi_mif, in_response_function_coefficients,
           out_sh_coefficients_image)
    if lmax is not None:
        cmd = cmd + ' -lmax ' + str(lmax)
    os.system(cmd)

    return out_sh_coefficients_image


def streamlines_tractography(
        in_source,
        in_white_matter_binary_mask,
        algorithm='iFOD2',
        number_of_tracks='100K',
        fod_threshold=None,
        step_size=None,
        angle=None,
        nthreads=2):
    """
    Perform streamlines tractography.

    This function performs whole-brain tractograms, using the deterministic or
    probabilistic streamlines algorithm which can be coupled with the fibre
    orientations produced using CSD. The choice is given by the `algorithm`
    parameter.

    Args:
        in_source (str): Image containing the source data. According to
            `algorithm`, the type of data can be:
            - iFOD1, iFOD2 (default), Nulldist2: the spherical harmonics
            coefficients image resulting from CSD
            - Tensor_Det, Tensor_Prob: the DWI image in MRtrix format.
        in_white_matter_binary_mask (str): Binary mask of the white matter
            segmentation. Seed streamlines will be entirely generated at random
            within this mask.
        algorithm (Optional[str]): #TODO
        number_of_tracks (Optional[str]): Desired number of tracks. The program
            will continue to generate tracks until this number of tracks have
            been selected and written to the output file; set to 0 to ignore
            limit. ``number_of_tracks` expects values in string format e.g.
            100K for 100000.
        fod_threshold (Optional[float]): FA or FOD amplitude cutoff for
            terminating tracks (default is 0.1) # TODO: -cutoff
        step_size (Optional[float]): Step size of the algorithm in mm (default
            is 0.1 x voxelsize; for iFOD2: 0.5 x voxelsize). # TODO: -step
        angle (Optional[int]): Maximum angle between successive steps
            (default is 90deg x step_size / voxelsize). # TODO: -angle
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 1 disables multi-threading).

    Returns:
        out_tracks (str): File containing the tracks generated in *.tck format.
    """
    import os.path as op
    import os

#    assert(op.isfile(in_source))
    assert(op.isfile(in_white_matter_binary_mask))
    if algorithm not in \
            ('iFOD1', 'iFOD2', 'Nulldist2', 'Tensor_Det', 'Tensor_Prob'):
        raise ValueError(
            'Invalid choice of algorithm in streamlines_tractography function')

    out_tracks = op.abspath('out_tracks_' + number_of_tracks + '.tck')

    cmd = 'tckgen -algorithm %s -number %s -seed_image %s %s %s -nthreads %s' \
          % (algorithm, number_of_tracks, in_white_matter_binary_mask,
             in_source, out_tracks, nthreads)
    # + ' -rk4'
    # cmd = cmd + ' -step ' + str(step_size)

    os.system(cmd)

    return out_tracks


def tcksift(in_tracks, in_fod):
    """
    Perform filtering of tractograms.

    This function filters a whole-brain fibre-tracking data set such that the
    streamline densities match the FOD lobe integrals.

    Args:
        in_tracks (str): Input track file.
        in_fod (str): Input image containing the spherical harmonics of the
            fibre orientation distributions.

    Returns:
        out_tracks (str): Output filtered tracks file in *.tck format.
    """
    import os.path as op
    import os

    assert(op.isfile(in_tracks))
    assert(op.isfile(in_fod))
    out_tracks = op.abspath('out_tracks_sift1.tck')

    cmd = 'tcksift ' + in_tracks + ' ' + in_fod + ' ' + out_tracks
    os.system(cmd)

    return out_tracks


def statistics_on_atlases(in_registered_map, name_map, prefix_file=None):
    """
    Computes a list of statistics files for each atlas.

    Args:
        in_registered_map (str): Map already registered on atlases in Nifti format.
        name_map (str): Name of the registered map in CAPS format.
        prefix_file (Opt[str]): <prefix_file>_space-<atlas_name>_map-<name_map>_statistics.tsv

    Returns:
        List of paths leading to the statistics TSV files.
    """
    from os import getcwd
    from os.path import abspath, join
    from nipype.utils.filemanip import split_filename
    from clinica.utils.atlas import (AtlasAbstract, JHUDTI811mm,
                                     JHUTracts01mm, JHUTracts251mm,
                                     JHUTracts501mm)
    from clinica.utils.statistics import statistics_on_atlas

#    atlas_classes = AtlasAbstract.__subclasses__()

    in_atlas_list = [JHUDTI811mm(),
                     JHUTracts01mm(), JHUTracts251mm()]

    atlas_statistics_list = []
    for atlas in in_atlas_list:
        if not isinstance(atlas, AtlasAbstract):
            raise TypeError("Atlas element must be an AtlasAbstract type")

        if prefix_file is None:
            _, base, _ = split_filename(in_registered_map)
            out_atlas_statistics = abspath(
                join(getcwd(), '%s_space-%s_res-%s_map-%s_statistics.tsv' %
                     (base, atlas.get_name_atlas(),
                      atlas.get_spatial_resolution(), name_map)))
        else:
            out_atlas_statistics = abspath(
                join(getcwd(), '%s_space-%s_res-%s_map-%s_statistics.tsv' %
                     (prefix_file, atlas.get_name_atlas(),
                      atlas.get_spatial_resolution(), name_map)))

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

    from clinica.utils.stream import cprint
    cprint("BIDS identifier: " + bids_identifier)
    m = re.search(
        r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+).*_dwi_space-[a-zA-Z0-9]+',
        caps_dwi_filename)
    caps_identifier = m.group(0)
    cprint("CAPS identifier: " + caps_identifier)

    return bids_identifier


def extract_caps_identifier_from_caps_filename(caps_dwi_filename):
    """Extract CAPS identifier from CAPS filename"""
    import re

    m = re.search(
        r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+).*_dwi_space-[a-zA-Z0-9]+',
        caps_dwi_filename)

    if m is None:
        raise ValueError('Input filename is not in a CAPS compliant format.')

    caps_identifier = m.group(0)

    from clinica.utils.stream import cprint
    cprint("CAPS identifier: " + caps_identifier)

    return caps_identifier


def rename_into_caps(in_caps_dwi,
                     in_norm_fa, in_norm_md, in_norm_ad, in_norm_rd,
                     in_b_spline_transform, in_affine_matrix):
    """
    Rename the outputs of the pipelines into CAPS format namely:
    <source_file>_space-T1w_preproc[.nii.gz|bval|bvec]

    Args:

    Returns:
        The different outputs in CAPS format
    """
    from nipype.utils.filemanip import split_filename
    from nipype.interfaces.utility import Rename
    import os

    from clinica.pipelines.dwi_processing_dti.dwi_processing_dti_utils import extract_bids_identifier_from_caps_filename

    bids_identifier = extract_bids_identifier_from_caps_filename(in_caps_dwi)

    # Extract base path from fname:
    base_dir_norm_fa, _, _ = split_filename(in_norm_fa)
    base_dir_norm_md, _, _ = split_filename(in_norm_md)
    base_dir_norm_ad, _, _ = split_filename(in_norm_ad)
    base_dir_norm_rd, _, _ = split_filename(in_norm_rd)
    base_dir_b_spline_transform, _, _ = split_filename(in_b_spline_transform)
    base_dir_affine_matrix, _, _ = split_filename(in_affine_matrix)

    # Rename into CAPS FA:
    rename_fa = Rename()
    rename_fa.inputs.in_file = in_norm_fa
    rename_fa.inputs.format_string = os.path.join(
        base_dir_norm_fa, bids_identifier + "_space-MNI152Lin_res-1x1x1_fa.nii.gz")
    out_caps_fa = rename_fa.run()

    # Rename into CAPS MD:
    rename_md = Rename()
    rename_md.inputs.in_file = in_norm_md
    rename_md.inputs.format_string = os.path.join(
        base_dir_norm_md, bids_identifier + "_space-MNI152Lin_res-1x1x1_md.nii.gz")
    out_caps_md = rename_md.run()

    # Rename into CAPS AD:
    rename_ad = Rename()
    rename_ad.inputs.in_file = in_norm_ad
    rename_ad.inputs.format_string = os.path.join(
        base_dir_norm_ad, bids_identifier + "_space-MNI152Lin_res-1x1x1_ad.nii.gz")
    out_caps_ad = rename_ad.run()

    # Rename into CAPS RD:
    rename_rd = Rename()
    rename_rd.inputs.in_file = in_norm_rd
    rename_rd.inputs.format_string = os.path.join(
        base_dir_norm_rd, bids_identifier + "_space-MNI152Lin_res-1x1x1_rd.nii.gz")
    out_caps_rd = rename_rd.run()

    # Rename into CAPS B-spline transform:
    rename_b_spline = Rename()
    rename_b_spline.inputs.in_file = in_b_spline_transform
    rename_b_spline.inputs.format_string = os.path.join(
        base_dir_b_spline_transform, bids_identifier + "_space-MNI152Lin_res-1x1x1_deformation.nii.gz")
    out_caps_b_spline_transform = rename_b_spline.run()

    # Rename into CAPS Affine Matrix:
    rename_affine = Rename()
    rename_affine.inputs.in_file = in_affine_matrix
    rename_affine.inputs.format_string = os.path.join(
        base_dir_affine_matrix, bids_identifier + "_space-MNI152Lin_res-1x1x1_affine.mat")
    out_caps_affine_matrix = rename_affine.run()

    from clinica.utils.stream import cprint
    cprint("Renamed files:")
    cprint(out_caps_fa.outputs.out_file)
    cprint(out_caps_md.outputs.out_file)
    cprint(out_caps_ad.outputs.out_file)
    cprint(out_caps_rd.outputs.out_file)
    cprint(out_caps_b_spline_transform.outputs.out_file)
    cprint(out_caps_affine_matrix.outputs.out_file)

    return out_caps_fa.outputs.out_file, out_caps_md.outputs.out_file,\
        out_caps_ad.outputs.out_file, out_caps_rd.outputs.out_file, \
        out_caps_b_spline_transform.outputs.out_file, \
        out_caps_affine_matrix.outputs.out_file
