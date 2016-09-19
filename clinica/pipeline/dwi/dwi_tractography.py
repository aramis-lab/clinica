#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains the tractography pipeline."""

def whole_brain_tractography_pipeline(
                datasink_directory, working_directory=None, max_harmonic_order=None,
                tractography_algorithm='iFOD2', tractography_nb_of_tracks="100K",
                tractography_fod_treshold=None, tractography_step_size=None, tractography_angle=None,
                nthreads=2, name="whole_brain_tractography_pipeline"):
    """
    Perform single-shell tractography.

    This pipeline performs a whole-brain single-shell tractography on a
    preprocessed DWI dataset. This Python implementation is using MRtrix3 and
    is based on the tutorial given by the MRtrix community during the ISMRM
    conference in 2015.

    .. warning :: This is not suited for multi-shell data at all.

    Args:
        datasink_directory (str): Directory where the results are stored.
        working_directory (Optional[str]): Directory where the temporary
            results are stored. If not specified, it is automatically
            generated (generally in /tmp/).
        max_harmonic_order (Optional[int]): Maximum harmonic order according to
            the b-vectors
        tractography_algorithm (Optional[str]): See streamlines_tractography
        tractography_nb_of_tracks (Optional[str]): See streamlines_tractography
        tractography_fod_treshold (Optional[float]): See streamlines_tractography
        step_size (Optional[int]): See streamlines_tractography
        angle (Optional[int]): See streamlines_tractography
        nthreads (Optional[int]): Number of threads used for the pipeline
            (default=2, 0 disables multi-threading).

    Inputnode:
        in_dwi_nii (str): File containing DWI dataset in NIfTI format.
        in_bvals (str): File containing B-Value table in FSL format.
        in_bvecs (str): File containing Diffusion Gradient table in FSL format.
        in_b0_mask (str): Binary mask of the b0 image. Only perform computation
            within this specified binary brain mask image.
        in_white_matter_binary_mask (str): Binary mask of the white matter
            segmentation. Seed streamlines will be entirely generated at random
            within this mask.

    Outputnode:
        out_dwi_mif (str): Preprocessed DWI in MRtrix format.
        out_dti (str): Tensor fitted to the DWI dataset.
        out_metrics (str): Maps of tensor-derived parameters namely fractional
            anisotropy, mean diffusivity (also called mean apparent diffusion),
            radial diffusivity and the first eigenvector modulated by the FA.
        out_eroded_mask (str): Eroded b0 mask (for debug puproses)
        out_response_function (str): Text file containing response function
            coefficients.
        out_sh_coefficients_image (str): File containing the spherical
            harmonics coefficients image
        out_tracks (str): File containing the generated tracks.

    Example:
        >>> from clinica.pipeline.dwi.dwi_tractography import whole_brain_tractography_pipeline
        >>> tractography = tractography_pipeline('/path/to/datasink/directory')
        >>> tractorgraphy.inputs.inputnode.in_dwi = 'subject_dwi.nii'
        >>> tractography.inputs.inputnode.in_bvecs = 'subject_dwi.bvecs'
        >>> tractography.inputs.inputnode.in_bvals = 'subject_dwi.bvals'
        >>> tractography.inputs.inputnode.in_b0_mask = 'subject_b0_mask.nii'
        >>> tractography.inputs.inputnode.white_matter_mask = 'subject_wm_mask.nii'
        >>> tractography.run()
    """
    import nipype.interfaces.io as nio
    import nipype.interfaces.utility as niu
    import nipype.pipeline.engine as pe
    from clinica.pipeline.dwi.dwi_tractography_utils import *
    from os.path import join
    import tempfile

    if working_directory is None:
        working_directory = tempfile.mkdtemp()

    inputnode = pe.Node(niu.IdentityInterface(
        fields=['in_dwi_nii', 'in_bvecs', 'in_bvals', 'in_b0_mask', 'in_white_matter_binary_mask']),
        name='inputnode')

    convert_nifti_to_mrtrix_format = pe.Node(interface=niu.Function(
        input_names=['in_dwi_nii', 'in_bvals', 'in_bvecs', 'nthreads'],
        output_names=['out_dwi_mif'],
        function=convert_nifti_to_mrtrix_format), name='convert_nifti_to_mrtrix_format')
    convert_nifti_to_mrtrix_format.inputs.nthreads = nthreads

    dwi_to_tensor = pe.Node(interface=niu.Function(
        input_names=['in_dwi_mif', 'in_b0_mask', 'nthreads'],
        output_names=['out_dti'],
        function=dwi_to_tensor), name='dwi_to_tensor')
    dwi_to_tensor.inputs.nthreads = nthreads

    tensor_to_metrics = pe.Node(interface=niu.Function(
        input_names=['in_dti', 'in_b0_mask', 'nthreads'],
        output_names=['out_fa', 'out_md', 'out_rd', 'out_ev'],
        function=tensor_to_metrics), name='tensor_to_metrics')

    erode_mask = pe.Node(interface=niu.Function(
        input_names=['in_mask', 'npass', 'nthreads'],
        output_names=['out_eroded_mask'], function=erode_mask), name='erode_mask')
    erode_mask.inputs.nthreads = nthreads

    estimate_response = pe.Node(interface=niu.Function(
        input_names=['in_dwi_mif', 'in_b0_mask', 'lmax', 'algorithm', 'tmpdir', 'nthreads'],
        output_names=['out_response_function'], function=estimate_response), name='estimate_response')
    estimate_response.inputs.lmax = max_harmonic_order
    estimate_response.inputs.tmpdir = working_directory
    estimate_response.inputs.nthreads = nthreads

    estimate_fod = pe.Node(interface=niu.Function(
        input_names=['in_dwi_mif', 'in_b0_mask', 'in_response_function_coefficients', 'lmax', 'nthreads'],
        output_names=['out_sh_coefficients_image'], function=estimate_fod), name='estimate_fod')
    estimate_fod.inputs.lmax = max_harmonic_order
    estimate_fod.inputs.nthreads = nthreads

    streamlines_tractography = pe.Node(interface=niu.Function(
        input_names=['in_source', 'in_white_matter_binary_mask', 'algorithm', 'number_of_tracks',
                     'fod_treshold', 'step_size', 'angle', 'nthreads'],
        output_names=['out_tracks'], function=streamlines_tractography), name='streamlines_tractography')
    streamlines_tractography.inputs.algorithm = tractography_algorithm
    streamlines_tractography.inputs.number_of_tracks = tractography_nb_of_tracks
    streamlines_tractography.inputs.fod_treshold = tractography_fod_treshold
    streamlines_tractography.inputs.step_size = tractography_step_size
    streamlines_tractography.inputs.angle = tractography_angle
    streamlines_tractography.inputs.nthreads = nthreads

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_dwi_mif', 'out_dti', 'out_metrics', 'out_fa', 'out_md', 'out_rd', 'out_ev',
                'out_response_function','out_sh_coefficients_image', 'out_tracks']),
        name='outputnode')

    datasink = pe.Node(nio.DataSink(), name='datasink')
    datasink.inputs.base_directory = join(datasink_directory, 'tractography/')

    wf = pe.Workflow(name=name)
    wf.connect([
        # Conversion to MRtrix format:
        (inputnode, convert_nifti_to_mrtrix_format, [('in_dwi_nii', 'in_dwi_nii'),
                                                     ('in_bvals', 'in_bvals'),
                                                     ('in_bvecs', 'in_bvecs')]),
        # Computation of the DTI:
        (inputnode,                      dwi_to_tensor, [('in_b0_mask', 'in_b0_mask')]),
        (convert_nifti_to_mrtrix_format, dwi_to_tensor, [('out_dwi_mif', 'in_dwi_mif')]),
        # Computation of the different metrics from the DTI:
        (inputnode,     tensor_to_metrics, [('in_b0_mask', 'in_b0_mask')]),
        (dwi_to_tensor, tensor_to_metrics, [('out_dti', 'in_dti')]),
        # Erosion of the b0 mask for the estimation of the response function:
        (inputnode, erode_mask, [('in_b0_mask','in_mask')]),
        # Estimation of the response function:
        (convert_nifti_to_mrtrix_format, estimate_response, [('out_dwi_mif', 'in_dwi_mif')]),
        (erode_mask,                     estimate_response, [('out_eroded_mask', 'in_b0_mask')]),
        # Estimation of the FOD using CSD:
        (convert_nifti_to_mrtrix_format, estimate_fod, [('out_dwi_mif', 'in_dwi_mif')]),
        (inputnode,                      estimate_fod, [('in_b0_mask', 'in_b0_mask')]),
        (estimate_response,              estimate_fod, [('out_response_function', 'in_response_function_coefficients')]),
        # Whole-brain tractography:
        (estimate_fod, streamlines_tractography, [('out_sh_coefficients_image', 'in_source')]),
        (inputnode,    streamlines_tractography, [('in_white_matter_binary_mask', 'in_white_matter_binary_mask')]),
        # Outputnode:
        (convert_nifti_to_mrtrix_format, outputnode, [('out_dwi_mif', 'out_dwi_mif')]),
        (dwi_to_tensor,                  outputnode, [('out_dti', 'out_dti')]),
        (tensor_to_metrics,              outputnode, [('out_fa', 'out_fa')]),
        (tensor_to_metrics,              outputnode, [('out_md', 'out_md')]),
        (tensor_to_metrics,              outputnode, [('out_rd', 'out_rd')]),
        (tensor_to_metrics,              outputnode, [('out_ev', 'out_ev')]),
        (erode_mask,                     outputnode, [('out_eroded_mask', 'out_eroded_mask')]),
        (estimate_response,              outputnode, [('out_response_function', 'out_response_function')]),
        (estimate_fod,                   outputnode, [('out_sh_coefficients_image', 'out_sh_coefficients_image')]),
        (streamlines_tractography,       outputnode, [('out_tracks', 'out_tracks')]),
        # Saving files with datasink:
        (convert_nifti_to_mrtrix_format, datasink, [('out_dwi_mif', 'out_dwi_mif')]),
        (dwi_to_tensor,                  datasink, [('out_dti', 'out_dti')]),
        (tensor_to_metrics,              datasink, [('out_fa', 'out_metrics.fa')]),
        (tensor_to_metrics,              datasink, [('out_md', 'out_metrics.md')]),
        (tensor_to_metrics,              datasink, [('out_rd', 'out_metrics.rd')]),
        (tensor_to_metrics,              datasink, [('out_ev', 'out_metrics.ev')]),
        (erode_mask,                     datasink, [('out_eroded_mask', 'out_eroded_mask')]),
        (estimate_response,              datasink, [('out_response_function', 'out_response_function')]),
        (estimate_fod,                   datasink, [('out_sh_coefficients_image', 'out_sh_coefficients_image')]),
        (streamlines_tractography,       datasink, [('out_tracks', 'out_tracks')])
    ])

    return wf
