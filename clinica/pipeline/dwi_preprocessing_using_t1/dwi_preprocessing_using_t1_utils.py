# coding: utf8


def ants_registration_syn_quick(fix_image, moving_image):
    import os
    import os.path as op

    image_warped = op.abspath('SyN_QuickWarped.nii.gz')
    affine_matrix = op.abspath('SyN_Quick0GenericAffine.mat')
    warp = op.abspath('SyN_Quick1Warp.nii.gz')
    inverse_warped = op.abspath('SyN_QuickInverseWarped.nii.gz')
    inverse_warp = op.abspath('SyN_Quick1InverseWarp.nii.gz')

    cmd = 'antsRegistrationSyNQuick.sh -t br -d 3 -f %s  -m %s -o SyN_Quick' \
          % (fix_image, moving_image)
    os.system(cmd)

    return image_warped, affine_matrix, warp, inverse_warped, inverse_warp


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


def dwi_container_from_filename(dwi_filename):
    import re
    from os.path import join
    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)_', dwi_filename)

    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format. '
            + 'It does not contain the subject and session information.')

    subject = m.group(1)
    session = m.group(2)

    return join('subjects', subject, session)
