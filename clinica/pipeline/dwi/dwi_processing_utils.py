#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains functions used for the tractography pipeline."""


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

    assert(op.isfile(in_dwi_nii))
    assert(op.isfile(in_bvals))
    assert(op.isfile(in_bvecs))

    out_dwi_mif = op.abspath('dwi.mif')
    cmd = 'mrconvert ' + in_dwi_nii + ' ' + out_dwi_mif + ' -fslgrad ' + in_bvecs + ' ' + in_bvals + ' -datatype float32 -stride 0,0,0,1 -nthreads ' + str(nthreads)
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

    # TODO : Check if b=0 is present, number_of_directions=number_of_directions-1;
    number_of_directions = os.system("$(sort bvec | uniq | wc -l)")
    out_lmax = 0
    if number_of_directions > 6 and number_of_directions <= 15:
        out_lmax = 2
    elif number_of_directions > 15 and number_of_directions <= 28:
        out_lmax = 4
    elif number_of_directions > 28 and number_of_directions <= 45:
        out_lmax = 6
    elif number_of_directions > 45 and number_of_directions <= 66:
        out_lmax = 8
    else:
        warnings.warn("The number of distinct directions is > to 66, lmax is automatically set to 8", UserWarning)
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
        out_eroded_mask : FILE
          Output. Eroded mask.
    """
    import os.path as op
    import os

    assert(op.isfile(in_mask))

    out_eroded_mask = op.abspath('eroded_mask.nii.gz')
    cmd = 'maskfilter -npass ' + str(npass) + ' -nthreads ' + str(nthreads) + ' ' + in_mask + ' erode ' + out_eroded_mask
    os.system(cmd)
    return out_eroded_mask



def dwi_to_tensor(in_dwi_mif, in_b0_mask, nthreads=2):
    """
    Perform diffusion tensor estimation from DWI dataset.

    Diffusion (kurtosis) tensor estimation using iteratively reweighted linear least squares estimator.

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

    assert(op.isfile(in_dwi_mif))
    assert(op.isfile(in_b0_mask))

    out_dti = op.abspath('dti.mif')

    cmd = 'dwi2tensor -mask ' + in_b0_mask + ' ' + in_dwi_mif + ' ' + out_dti + ' -nthreads ' + str(nthreads)
    os.system(cmd)

    return out_dti



def tensor_to_metric(in_dti, in_b0_mask, metric='fa', nthreads=2):
    """
    Generate maps of tensor-derived parameters.

    ..warning:: Only FA and MD are available for the moment.

    Args:
        in_dti (str): Tensor.
        in_b0_mask (str): Binary mask of the b0 image. Only perform computation
            within this specified binary brain mask image.
        metric (str): Tensor-derived parameter. Currenly, choices are:
            - adc/md: mean apparent diffusion coefficient (also called
            mean diffusivity)
            - fa(default): fractional anisotropy
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 1 disables multi-threading).

    Returns:
        out_metric (str): The tensor-derived parameter `metric`.
    """
    import os.path as op
    import os

    assert(op.isfile(in_dti))
    assert(op.isfile(in_b0_mask))

    cmd = 'tensor2metric -mask ' + in_b0_mask + ' ' + in_dti + ' -nthreads ' + str(nthreads)
    if metric in ("ADC", "adc", "MD", "md"):
        out_metric = op.abspath('md.nii.gz')
        cmd = cmd + ' -adc ' + out_metric
    elif metric in ("FA", "fa"):
        out_metric = op.abspath('fa.nii.gz')
        cmd = cmd + ' -fa ' + out_metric
    else:
        raise ValueError('Invalid choice of metric in tensor2metric function')
    os.system(cmd)

    return out_metric



def tensor_to_metrics(in_dti, in_b0_mask, nthreads=2):
    """
    Generate maps of tensor-derived parameters.

    Args:
        in_dti (str): Tensor.
        in_b0_mask (str): Binary mask of the b0 image. Only perform computation
            within this specified binary brain mask image.
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 1 disables multi-threading).

    Returns:
        out_fa (str): The tensor-derived parameter fractional anisotropy.
        out_md (str): The tensor-derived parameter mean diffusivity (also called mean apparent diffusion).
        out_ad (str): The tensor-derived parameter axial diffusivity.
        out_rd (str): The tensor-derived parameter radial diffusivity.
        out_ev (str): The directionally-encoded colour (DEC) fractional anisotropy map This corresponds to the first eigenvector modulated by the FA.
    """
    import os.path as op
    import os

    assert(op.isfile(in_dti))
    assert(op.isfile(in_b0_mask))

    out_fa = op.abspath('fa_map_from_dti.nii.gz')
    cmd = 'tensor2metric -mask ' + in_b0_mask + ' ' + in_dti + ' -nthreads ' + str(nthreads) + ' -fa ' + out_fa
    os.system(cmd)

    out_md = op.abspath('fa_map_from_dti.nii.gz')
    cmd = 'tensor2metric -mask ' + in_b0_mask + ' ' + in_dti + ' -nthreads ' + str(nthreads) + ' -adc ' + out_md
    os.system(cmd)

    out_ad = op.abspath('ad_map_from_dti.nii.gz')
    cmd = 'tensor2metric -mask ' + in_b0_mask + ' ' + in_dti + ' -nthreads ' + str(nthreads) + ' -ad ' + out_ad
    os.system(cmd)

    out_rd = op.abspath('rd_map_from_dti.nii.gz')
    cmd = 'tensor2metric -mask ' + in_b0_mask + ' ' + in_dti + ' -nthreads ' + str(nthreads) + ' -rd ' + out_rd
    os.system(cmd)

    out_ev = op.abspath('dec_fa_map_from_dti.nii.gz')
    cmd = 'tensor2metric -mask ' + in_b0_mask + ' ' + in_dti + ' -nthreads ' + str(nthreads) + ' -vector ' + out_ev
    os.system(cmd)

    return out_fa, out_ad, out_md, out_rd, out_ev


def estimate_response(in_dwi_mif, in_b0_mask, lmax=None, algorithm='tax', tmpdir='/tmp/', nthreads=2):
    """
    Estimate response function(s) for spherical deconvolution.

    This function generates an appropriate response function (chosen by
    `algorithm`) from the image data for spherical deconvolution.

    Args:
        in_dwi_mif (str): DWI dataset in MRtrix image format.
        in_b0_mask (str): Binary mask of the b0 image. Only perform computation within this specified binary brain
            mask image. It is advised to erode the mask (we assume that an eroded mask is given)
        algorithm (Optional[str]): Used algorithm to derive the response function. Currenly, choices are:
            - tax (default): Use the Tax et al. (2014) recursive calibration algorithm
            - tournier: Use the Tournier et al. (2013) iterative RF selection algorithm
        lmax (Optional[int]): The maximum harmonic degree(s) of response function estimation
        tmpdir (Optional[str]): Path where the temporary results are stored.
        nthreads (Optional[int]): Number of threads used in this function (default=2, 1 disables multi-threading).

    Returns:
        out_response_function (str): Text file containing response function
            coefficients.
    """
    import os.path as op
    import os

    assert(op.isfile(in_dwi_mif))
    assert(op.isfile(in_b0_mask))

    if algorithm == 'tax':
        out_response_function = op.abspath('out_response_function_tax.txt')
    elif algorithm == 'tournier':
        out_response_function = op.abspath('out_response_function_tournier.txt')
    else:
        raise ValueError('Invalid choice of algorithm in estimate_response function')

    cmd = 'dwi2response ' + algorithm + ' -mask ' + in_b0_mask
    if lmax is not None:
        cmd = cmd + ' -lmax ' + str(lmax)
    cmd = cmd + ' ' + in_dwi_mif + ' ' +  out_response_function + ' -nthreads ' + str(nthreads)
    os.system(cmd)

    return out_response_function



def estimate_fod(in_dwi_mif, in_b0_mask, in_response_function_coefficients, lmax=None, nthreads=2):
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

    cmd = 'dwi2fod csd -mask ' + in_b0_mask + ' ' +  in_dwi_mif + ' ' + in_response_function_coefficients + ' ' + out_sh_coefficients_image
    if lmax is not None:
        cmd = cmd + ' -lmax ' + str(lmax)
    os.system(cmd)

    return out_sh_coefficients_image


def streamlines_tractography(
        in_source, in_white_matter_binary_mask, algorithm='iFOD2', number_of_tracks='100K',
        fod_threshold=None, step_size=None, angle=None, nthreads=2):
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
            (default is 90deg x stepsize / voxelsize). # TODO: -angle
        nthreads (Optional[int]): Number of threads used in this function
            (default=2, 1 disables multi-threading).

    Returns:
        out_tracks (str): File containing the tracks generated in *.tck format.
    """
    import os.path as op
    import os

#    assert(op.isfile(in_source))
    assert(op.isfile(in_white_matter_binary_mask))
    if algorithm not in ('iFOD1', 'iFOD2', 'Nulldist2', 'Tensor_Det', 'Tensor_Prob'):
        raise ValueError('Invalid choice of algorithm in streamlines_tractography function')

    out_tracks = op.abspath('out_tracks_' + number_of_tracks + '.tck')

    cmd = 'tckgen -algorithm ' + algorithm + ' -number ' + number_of_tracks + ' -seed_image ' + in_white_matter_binary_mask + \
        ' ' + in_source + ' ' + out_tracks +  ' -nthreads ' + str(nthreads) #+ ' -rk4'
    # cmd = cmd + ' -step ' + str(step_size)

    os.system(cmd)

    return out_tracks


def tcksift(in_tracks, in_fod):
    """
    Perform filtering of tractograms.

    This function filters a whole-brain fibre-tracking data set such that the streamline
    densities match the FOD lobe integrals.

    Args:
        in_source (str): Input track file.
        in_fod (str): Input image containing the spherical harmonics of the fibre orientation distributions.

    Returns:
        out_tracks (str): Output filtered tracks file in *.tck format.
    """
    import os.path as op
    import os

    assert(op.isfile(in_tracks))
    assert(op.isfile(in_fod))
    out_tracks = op.abspath('out_tracks_sift1.tck')

    cmd = 'tcksift ' + in_tracks+ ' ' + in_fod + ' ' + out_tracks
    os.system(cmd)

    return out_tracks
