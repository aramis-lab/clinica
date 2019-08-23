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
    """ Generate subjects/sub-<participant_id>/ses-<session_id> folder
    from BIDS filename."""
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

    cmd = 'WarpImageMultiTransform 3 ' + moving_image + ' ' + out_warp + ' -R ' + fix_image + ' ' + ants_warp_affine[0] + ' ' + ants_warp_affine[1] + ' ' + ants_warp_affine[2]
    os.system(cmd)

    return out_warp


def rotate_bvecs(in_bvec, in_matrix):
    """
    Rotates the input bvec file accordingly with a list of matrices.
    .. note:: the input affine matrix transforms points in the destination
      image to their corresponding coordinates in the original image.
      Therefore, this matrix should be inverted first, as we want to know
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


def init_input_node(t1w, dwi, bvec, bval, total_readout_time, phase_encoding_direction):
    """Extract "sub-<participant_id>_ses-<session_label>" from input node and print begin message."""
    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.stream import cprint
    import datetime
    from colorama import Fore
    from clinica.utils.io import get_subject_id

    id_t1w = get_subject_id(t1w)
    id_dwi = get_subject_id(dwi)
    id_bvec = get_subject_id(bvec)
    id_bval = get_subject_id(bval)

    image_id = list(set([id_t1w, id_dwi, id_bvec, id_bval]))

    if not len(image_id) == 1:
        raise ClinicaException('<image_id> from input files mismatch (found: %s)' % image_id)

    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s Running pipeline for %s '
           '(TotalReadoutTime = %s, PhaseEncodingDirection = %s)' %
           (Fore.BLUE, now, Fore.RESET, image_id[0].replace('_', '|'),
            total_readout_time, phase_encoding_direction))
    return (image_id[0], t1w, dwi, bvec, bval,
            total_readout_time, phase_encoding_direction)


def print_end_pipeline(image_id, final_file):
    """Display end message for `image_id` when `final_file` is connected."""
    from clinica.utils.stream import cprint
    import datetime
    from colorama import Fore

    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s ...%s has completed.' % (
        Fore.GREEN, now, Fore.RESET, image_id.replace('_', '|')))


# def convert_eddy_2_hmc_ecc_flirt(eddy_parameters):
#
#     import os.path as op
#     import pandas as pd
#     import numpy as np
#
#     df_eddy = pd.read_csv(eddy_parameters, sep='  ', header=None)
#     df_hmc = df_eddy.iloc[:, 0:6]
#     df_ecc = df_eddy.iloc[:, 6:]
#
#     hmc_affine_list = []
#     ecc_affine_list = []
#
#     # convert the df_hmc into 4*4 affine matrix
#     num_vols = df_hmc.shape[0]
#     for i in range(num_vols):
#
#         # HMC
#         hmc_affine = np.ones([4, 4])
#         # for translation
#         hmc_affine[0, 3], hmc_affine[1, 3], hmc_affine[2, 3] = df_hmc.iloc[i, :][0], df_hmc.iloc[i, :][1], df_hmc.iloc[i, :][2]
#         hmc_affine[3, 0], hmc_affine[3, 1], hmc_affine[3, 2] = 0, 0, 0
#         # for rotation
#         Rx = np.array([[1, 0, 0],
#                        [0, np.cos(df_hmc.iloc[i, :][3]), -np.sin(df_hmc.iloc[i, :][3])],
#                        [0, np.sin(df_hmc.iloc[i, :][3]), np.cos(df_hmc.iloc[i, :][3])]])
#
#         Ry = np.array([[np.cos(df_hmc.iloc[i, :][4]), 0, np.sin(df_hmc.iloc[i, :][4])],
#                        [0, 1, 0],
#                        [-np.sin(df_hmc.iloc[i, :][4]), 0, np.cos(df_hmc.iloc[i, :][4])]])
#
#         Rz = np.array([[np.cos(df_hmc.iloc[i, :][5]), -np.sin(df_hmc.iloc[i, :][5]), 0],
#                        [np.sin(df_hmc.iloc[i, :][5]), np.cos(df_hmc.iloc[i, :][5]), 0],
#                        [0, 0, 1]])
#
#         R = np.multiply(np.multiply(Rz, Ry), Rx)
#
#         hmc_affine[0:3, 0:3] = R
#
#         np.savetxt(op.abspath("hmc_vol" + str(i) + '_affine.mat'), hmc_affine, delimiter='  ')
#         hmc_affine_list.append(op.abspath("hmc_vol" + str(i) + '_affine.mat'))
#
#         # ecc
#         # TODO, just set it into an identity matrix
#         ecc_affine = np.array([[1, 0, 0, 0],
#                                [0, 1, 0, 0],
#                                [0, 0, 1, 0],
#                                [0, 0, 0, 1]])
#
#         np.savetxt(op.abspath("ecc_vol" + str(i) + '_affine.mat'), ecc_affine, delimiter='  ')
#         ecc_affine_list.append(op.abspath("ecc_vol" + str(i) + '_affine.mat'))
#
#     return hmc_affine_list, ecc_affine_list
