#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains FreeSurfer utilities."""


def label_to_volume_in_native_space(subject_id, label="aparc+aseg"):
    """
    Convert labels to a volume.

    This function converts a label or a set of labels (e.g. aseg.mgz,
    aparc+aseg.mgz, or aparc.a2009s+aseg.mgz where aparc (resp. aparc.a2009s)
    denotes the Desikan (resp. Destrieux) parcellation) into a volume in native
    space i.e. with the same header of subject_id originary image, not the
    FreeSurfer's conformed space.

    Args:
        subject_id (str): Folder located in ${SUBJECTS_DIR} containing
            a FreeSurfer segmentation.
        label (str): Name of the parcellation located in
            ${SUBJECTS_DIR}/subject_id/mri/label.mgz.

    Returns:
        out_volume (str): volume in native space (the file is saved here:
            ${SUBJECTS_DIR}/subject_id/native_space/label.nii)

    Example:
        >>> from clinica.utils.freesurfer import label_to_volume_in_native_space
        >>> label_to_volume_in_native_space(bert, aparc.a2009s+aseg)
    """
    #import numpy as np
    #import nibabel as nib
    import os

    subjects_dir = os.environ['SUBJECTS_DIR']
    assert(os.path.isdir(subjects_dir))

    path_to_subject_id = os.path.join(subjects_dir, subject_id)

    label_file = os.path.join(path_to_subject_id, 'mri', label, '.mgz')
    rawavg_file = os.path.join(path_to_subject_id, 'mri', 'rawavg.mgz')

    if not os.path.exists(os.path.join(path_to_subject_id, 'native_space')):
        os.makedirs(os.path.join(path_to_subject_id, 'native_space'))

    assert(os.path.isdir(path_to_subject_id))

    out_volume = os.path.join(path_to_subject_id, 'native_space', label, '.nii')
    cmd = 'mri_vol2vol --regheader --no-save-reg --mov ' + label_file + ' --targ ' + rawavg_file + ' --o ' + out_volume
    cmd = 'mri_label2vol --seg ' + label_file + '--temp rawavg.mgz --o ' + out_volume + ' --regheader ' + label_file
    os.system(cmd)

    return out_volume


def freesurfer_volume_to_native_volume(freesurfer_volume, native_volume, name_output_volume=None):
    """
    Convert FreeSurfer volume in native space.

    This function converts any volume in FreeSurfer's conformed space
    (1x1x1mm voxel size, 256x256x256 dimension) into a volume in native space.

    For further details: https://surfer.nmr.mgh.harvard.edu/fswiki/FsAnat-to-NativeAnat

    Args:
        freesurfer_volume (str): Volume in FreeSurfer's conformed space (e.g.
            aparc+aseg.mgz containing the Desikan parcellation)
        native_volume (str): Volume in native space (You should choose
            ${SUBJECTS_DIR}/subject_id/mri/rawavg.mgz).
        name_output_volume (Optional[str]): Name of the output matrix
            (default=volume_in_native_space.nii.gz).

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
