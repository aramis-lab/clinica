# coding: utf8


def ants_combine_transform(in_file, transforms_list, reference):
    import os
    import os.path as op

    out_warp = op.abspath('out_warp.nii.gz')

    transforms = ""
    for trans in transforms_list:
        transforms += " " + trans
    cmd = 'antsApplyTransforms -o [out_warp.nii.gz,1] -i %s -r %s -t %s' \
          % (in_file, reference, transforms)
    os.system(cmd)

    return out_warp


def dwi_container_from_filename(bids_dwi_filename):
    """ Generate subjects/sub-<participant_id>/ses-<session_id> folder from BIDS filename."""
    import re
    from os.path import join
    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)_', bids_dwi_filename)

    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format. '
            + 'It does not contain the subject and session information.')

    subject = m.group(1)
    session = m.group(2)

    return join('subjects', subject, session)


def rename_into_caps(in_bids_dwi,
                     fname_dwi, fname_bval, fname_bvec, fname_brainmask):
    """
    Rename the outputs of the pipelines into CAPS format namely:
    <source_file>_space-T1w_preproc[.nii.gz|bval|bvec]

    Args:
        in_bids_dwi (str): Input BIDS DWI to extract the <source_file>
        fname_dwi (str): Preprocessed DWI.
        fname_bval (str): Preprocessed DWI.
        fname_bvec (str): Preprocessed DWI.
        fname_brainmask (str): B0 mask.

    Returns:
        The different outputs in CAPS format
    """
    from nipype.utils.filemanip import split_filename
    from nipype.interfaces.utility import Rename
    import os

    # Extract <source_file> in format sub-CLNC01_ses-M00[_acq-label]_dwi
    _, source_file_dwi, _ = split_filename(in_bids_dwi)

    # Extract base path from fname:
    base_dir_dwi, _, _ = split_filename(fname_dwi)
    base_dir_bval, _, _ = split_filename(fname_bval)
    base_dir_bvec, _, _ = split_filename(fname_bvec)
    base_dir_brainmask, _, _ = split_filename(fname_brainmask)

    # Rename into CAPS DWI:
    rename_dwi = Rename()
    rename_dwi.inputs.in_file = fname_dwi
    rename_dwi.inputs.format_string = os.path.join(
        base_dir_dwi, source_file_dwi + "_space-T1w_preproc.nii.gz")
    out_caps_dwi = rename_dwi.run()

    # Rename into CAPS bval:
    rename_bval = Rename()
    rename_bval.inputs.in_file = fname_bval
    rename_bval.inputs.format_string = os.path.join(
        base_dir_bval, source_file_dwi + "_space-T1w_preproc.bval")
    out_caps_bval = rename_bval.run()

    # Rename into CAPS bvec:
    rename_bvec = Rename()
    rename_bvec.inputs.in_file = fname_bvec
    rename_bvec.inputs.format_string = os.path.join(
        base_dir_bvec, source_file_dwi + "_space-T1w_preproc.bvec")
    out_caps_bvec = rename_bvec.run()

    # Rename into CAPS DWI:
    rename_brainmask = Rename()
    rename_brainmask.inputs.in_file = fname_brainmask
    rename_brainmask.inputs.format_string = os.path.join(
        base_dir_brainmask, source_file_dwi + "_space-T1w_brainmask.nii.gz")
    out_caps_brainmask = rename_brainmask.run()

    return out_caps_dwi.outputs.out_file, out_caps_bval.outputs.out_file, \
        out_caps_bvec.outputs.out_file, out_caps_brainmask.outputs.out_file


def change_itk_transform_type(input_affine_file):
    """
    This function takes in the affine.txt produced by the c3d_affine_tool
    (which converted an FSL FLIRT affine.mat into the affine.txt)
    it then modifies the 'Transform Type' of this affine.txt so that it is
    compatible with the antsApplyTransforms tool and produces a new affine
    file titled 'updated_affine.txt'
    """
    import os

    new_file_lines = []

    with open(input_affine_file) as f:
        for line in f:
            if 'Transform:' in line:
                if 'MatrixOffsetTransformBase_double_3_3' in line:
                    transform_line = 'Transform: AffineTransform_double_3_3\n'
                    new_file_lines.append(transform_line)
            else:
                new_file_lines.append(line)

    updated_affine_file = os.path.join(os.getcwd(), 'updated_affine.txt')

    with open(updated_affine_file, 'wt') as f:
        for line in new_file_lines:
            f.write(line)

    return updated_affine_file


def expend_matrix_list(in_matrix, in_bvec):
    import numpy as np

    bvecs = np.loadtxt(in_bvec).T
    out_matrix_list = [in_matrix]

    out_matrix_list = out_matrix_list * len(bvecs)

    return out_matrix_list


def ants_warp_image_multi_transform(fix_image, moving_image, ants_warp_affine):
    import os
    import os.path as op

    out_warp = op.abspath('warped_epi.nii.gz')

    cmd = ('WarpImageMultiTransform 3 %s %s -R %s %s %s %s' %
           (moving_image, out_warp, fix_image, ants_warp_affine[0], ants_warp_affine[1], ants_warp_affine[2]))
    os.system(cmd)

    return out_warp


def rotate_bvecs(in_bvec, in_matrix):
    """
    Rotates the input bvec file accordingly with a list of matrices.

    Notes:
        The input affine matrix transforms points in the destination image to their corresponding coordinates
        in the original image. Therefore, this matrix should be inverted first, as we want to know
        the target position of :math:`\\vec{r}`.
    """
    import os
    import numpy as np

    name, fext = os.path.splitext(os.path.basename(in_bvec))
    if fext == '.gz':
        name, _ = os.path.splitext(name)
    out_file = os.path.abspath('%s_rotated.bvec' % name)
    bvecs = np.loadtxt(in_bvec).T   # Warning, bvecs.txt are not in the good configuration, need to put '.T'
    new_bvecs = []

    if len(bvecs) != len(in_matrix):
        raise RuntimeError(('Number of b-vectors (%d) and rotation '
                            'matrices (%d) should match.') % (len(bvecs),
                                                              len(in_matrix)))

    for bvec, mat in zip(bvecs, in_matrix):
        if np.all(bvec == 0.0):
            new_bvecs.append(bvec)
        else:
            invrot = np.linalg.inv(np.loadtxt(mat))[:3, :3]
            newbvec = invrot.dot(bvec)
            new_bvecs.append((newbvec / np.linalg.norm(newbvec)))

    np.savetxt(out_file, np.array(new_bvecs).T, fmt='%0.15f')
    return out_file


def ants_combin_transform(fix_image, moving_image, ants_warp_affine):
    import os
    import os.path as op

    out_warp_field = op.abspath('out_warp_field.nii.gz')
    out_warped = op.abspath('out_warped.nii.gz')

    cmd = 'antsApplyTransforms -o [out_warp_field.nii.gz,1] -i ' + moving_image + ' -r ' + fix_image + ' -t ' + \
          ants_warp_affine[0] + ' ' + ants_warp_affine[1] + ' ' + ants_warp_affine[2]
    os.system(cmd)

    cmd1 = 'antsApplyTransforms -o out_warped.nii.gz -i ' + moving_image + ' -r ' + fix_image + ' -t ' + \
           ants_warp_affine[0] + ' ' + ants_warp_affine[1] + ' ' + ants_warp_affine[2]
    os.system(cmd1)

    return out_warp_field, out_warped


def create_jacobian_determinant_image(imageDimension, deformationField, outputImage):
    import os
    import os.path as op

    outputImage = op.abspath(outputImage)

    cmd = 'CreateJacobianDeterminantImage ' + str(imageDimension) + ' ' + deformationField + ' ' + outputImage
    os.system(cmd)

    return outputImage


def init_input_node(t1w, dwi, bvec, bval, dwi_json):
    """Extract "sub-<participant_id>_ses-<session_label>" from input node and print begin message."""
    from clinica.utils.filemanip import get_subject_id, extract_metadata_from_json
    from clinica.utils.dwi import check_dwi_volume, bids_dir_to_fsl_dir
    from clinica.utils.ux import print_begin_image

    # Check the image IDs for each file are the same
    image_id = get_subject_id(t1w)

    # Check that the number of DWI, bvec & bval are the same
    check_dwi_volume(dwi, bvec, bval)

    # Read metadata from DWI JSON file:
    [total_readout_time, phase_encoding_direction] = \
        extract_metadata_from_json(dwi_json, ['TotalReadoutTime', 'PhaseEncodingDirection'])
    phase_encoding_direction = bids_dir_to_fsl_dir(phase_encoding_direction)

    # Print begin message
    print_begin_image(image_id,
                      ['TotalReadoutTime', 'PhaseEncodingDirection'],
                      [total_readout_time, phase_encoding_direction])

    return (image_id, t1w, dwi, bvec, bval,
            total_readout_time, phase_encoding_direction)


def print_end_pipeline(image_id, final_file):
    """Display end message for `image_id` when `final_file` is connected."""
    from clinica.utils.ux import print_end_image
    print_end_image(image_id)


def prepare_reference_b0(in_dwi, in_bval, in_bvec, low_bval=5, working_directory=None):
    """
    This function prepare the data for further corrections. It co-registers
    the B0 images and then average it in order to obtain only
    one average B0 images.

    Args:
        in_dwi (str): Input DWI file.
        in_bvec (str): Vector file of the diffusion directions of
            the dwi dataset.
        in_bval (str): B-values file.
        low_bval (optional[int]):
        working_directory (str): temporary folder results where the results are stored

    Returns:
        out_reference_b0 (str): Average of the B0 images or the only B0 image.
        out_b0_mask (str): Binary mask obtained from the average of
            the B0 images.
        out_b0_dwi_merge (str): Average of B0 images merged to the DWIs.
        out_updated_bval (str): Updated gradient values table.
        out_updated_bvec (str): Updated gradient vectors table.
    """
    from clinica.utils.dwi import (insert_b0_into_dwi, b0_dwi_split,
                                   count_b0s, b0_average)
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows import b0_flirt_pipeline
    import hashlib

    import os.path as op

    import tempfile

    # Count the number of b0s
    nb_b0s = count_b0s(in_bval=in_bval, low_bval=low_bval)
    # cprint('Found %s b0 for %s' % (nb_b0s, in_dwi))

    # Split dataset into two datasets: the b0 and the b>low_bval datasets
    [extracted_b0, out_split_dwi, out_split_bval, out_split_bvec] = \
        b0_dwi_split(
            in_dwi=in_dwi, in_bval=in_bval, in_bvec=in_bvec, low_bval=low_bval)

    if nb_b0s == 1:
        # The reference b0 is the extracted b0
        # cprint('Only one b0 for %s' % in_dwi)
        out_reference_b0 = extracted_b0
    elif nb_b0s > 1:
        # Register the b0 onto the first b0
        b0_flirt = b0_flirt_pipeline(num_b0s=nb_b0s)
        b0_flirt.inputs.inputnode.in_file = extracted_b0
        if working_directory is None:
            working_directory = tempfile.mkdtemp()
        tmp_dir = op.join(working_directory, hashlib.md5(in_dwi.encode()).hexdigest())
        b0_flirt.base_dir = tmp_dir
        b0_flirt.run()
        # BUG: Nipype does allow to extract the output after running the
        # workflow: we need to 'guess' where the output will be generated
        # out_node = b0_flirt.get_node('outputnode')
        registered_b0s = op.abspath(op.join(
            tmp_dir, 'b0_coregistration', 'concat_ref_moving',
            'merged_files.nii.gz'))
        # cprint('B0 s will be averaged (file = ' + registered_b0s + ')')
        # Average the b0s to obtain the reference b0
        out_reference_b0 = b0_average(in_file=registered_b0s)
    else:
        raise ValueError(
            'The number of b0s should be strictly positive (b-val file =%s).'
            % in_bval)

    # Merge datasets such that bval(DWI) = (0 b1 ... bn)
    [out_b0_dwi_merge, out_updated_bval, out_updated_bvec] = \
        insert_b0_into_dwi(in_b0=out_reference_b0, in_dwi=out_split_dwi,
                           in_bval=out_split_bval, in_bvec=out_split_bvec)

    return out_reference_b0, out_b0_dwi_merge, out_updated_bval, out_updated_bvec