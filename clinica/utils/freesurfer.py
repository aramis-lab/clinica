#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains FreeSurfer utilities."""




def freesurfer_volume_to_native_volume(freesurfer_volume, native_volume, name_output_volume=None):
    """
    Convert FreeSurfer volume in native space.

    This function converts any volume in FreeSurfer's conformed space
    (1x1x1mm voxel size, 256x256x256 dimension) into a volume in native space.

    For further details: https://surfer.nmr.mgh.harvard.edu/fswiki/FsAnat-to-NativeAnat

    Args:
        freesurfer_volume (str): Volume in FreeSurfer's conformed space (e.g. aparc+aseg.mgz containing
            the Desikan parcellation)
        native_volume (str): Volume in native space (You should choose ${SUBJECTS_DIR}/subject_id/mri/rawavg.mgz).
        name_output_volume (Optional[str]): Name of the output matrix (default=volume_in_native_space.nii.gz).

    Returns:
        out_volume (str): volume in native space (the file is saved here:
            ${SUBJECTS_DIR}/subject_id/native_space/label.nii)

    Example:
        >>> from clinica.utils.freesurfer import freesurfer_volume_to_native_volume
        >>> freesurfer_volume_to_native_volume(bert/mri/rawavg.mgz, bert/mri/aparc+aseg.mgz, 'aparc-in-native-space.nii')
    """
    import os
    import os.path as op

    assert(op.isfile(freesurfer_volume))
    assert(op.isfile(native_volume))

    if name_output_volume is None:
        out_volume = op.abspath('volume_in_native_space.nii.gz')
    else:
        out_volume = op.abspath(name_output_volume);

    cmd = 'mri_vol2vol --regheader --no-save-reg --mov ' + freesurfer_volume + ' --targ ' + native_volume + ' --o ' + out_volume
    os.system(cmd)

    return out_volume
