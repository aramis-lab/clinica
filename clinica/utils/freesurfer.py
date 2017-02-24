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


def fs_caps2reconall(caps_dir, dest_dir, subjects_visits_tsv):
    """
    This function is to transfer caps recon-all output structure to standard freesurfer recon-all output structure.
    :param caps_dir:
    :param dest_dir:
    :param subjects_visits_tsv:
    :param analysis_series_id:
    :return:
    """
    import os, csv
    from shutil import copytree

    subject_list = []
    session_list = []
    with open(subjects_visits_tsv, 'rb') as tsvin:
        tsv_reader = csv.reader(tsvin, delimiter='\t')

        for row in tsv_reader:
            if row[0] == 'participant_id':
                continue
            else:
                subject_list.append(row[0])
                session_list.append(row[1])

    output_path = os.path.expanduser(caps_dir)  # change the relative path to be absolute path
    caps_dir = os.path.join(output_path, 'subjects')

    for i in range(len(subject_list)):
        copytree(os.path.join(caps_dir, subject_list[i], session_list[i], 't1/freesurfer-cross-sectional', subject_list[i] + '_' + session_list[i]), os.path.join(dest_dir, subject_list[i] + '_' + session_list[i]))

def write_statistics_summary(subject_dir, subject_id, output_dir):
    """
    To write statistics summary for all the subjects after reconall
    :param subject_dir:
    :param subject_id:
    :param output_dir:
    :return:
    """
    ## TODO check if summary_tsv exits, if yes, add new info for new subjects, not to overwrite it.

    import os, errno

    # name all the 26 tsv output files.
    all_seg_volume = 'measure_all-seg.tsv'
    aseg_volume = 'measure_aseg-volume.tsv'

    aparc_desikan_lh_volume = 'hemisphere-lh_parcellation-desikan_measure-volume.tsv'
    aparc_desikan_rh_volume = 'hemisphere-rh_parcellation-desikan_measure-volume.tsv'
    aparc_desikan_lh_thickness = 'hemisphere-lh_parcellation-desikan_measure-thickness.tsv'
    aparc_desikan_rh_thickness = 'hemisphere-rh_parcellation-desikan_measure-thickness.tsv'
    aparc_desikan_lh_area = 'hemisphere-lh_parcellation-desikan_measure-area.tsv'
    aparc_desikan_rh_area = 'hemisphere-rh_parcellation-desikan_measure-area.tsv'
    aparc_desikan_lh_meancurv = 'hemisphere-lh_parcellation-desikan_measure-meancurv.tsv'
    aparc_desikan_rh_meancurv = 'hemisphere-rh_parcellation-desikan_measure-meancurv.tsv'

    aparc_destrieux_lh_volume = 'hemisphere-lh_parcellation-destrieux_measure-volume.tsv'
    aparc_destrieux_rh_volume = 'hemisphere-rh_parcellation-destrieux_measure-volume.tsv'
    aparc_destrieux_lh_thickness = 'hemisphere-lh_parcellation-destrieux_measure-thickness.tsv'
    aparc_destrieux_rh_thickness = 'hemisphere-rh_parcellation-destrieux_measure-thickness.tsv'
    aparc_destrieux_lh_area = 'hemisphere-lh_parcellation-destrieux_measure-area.tsv'
    aparc_destrieux_rh_area = 'hemisphere-rh_parcellation-destrieux_measure-area.tsv'
    aparc_destrieux_lh_meancurv = 'hemisphere-lh_parcellation-destrieux_measure-meancurv.tsv'
    aparc_destrieux_rh_meancurv = 'hemisphere-rh_parcellation-destrieux_measure-meancurv.tsv'

    aparc_BA_lh_volume = 'hemisphere-lh_parcellation-BA_measure-volume.tsv'
    aparc_BA_rh_volume = 'hemisphere-rh_parcellation-BA_measure-volume.tsv'
    aparc_BA_lh_thickness = 'hemisphere-lh_parcellation-BA_measure-thickness.tsv'
    aparc_BA_rh_thickness = 'hemisphere-rh_parcellation-BA_measure-thickness.tsv'
    aparc_BA_lh_area = 'hemisphere-lh_parcellation-BA_measure-area.tsv'
    aparc_BA_rh_area = 'hemisphere-rh_parcellation-BA_measure-area.tsv'
    aparc_BA_lh_meancurv = 'hemisphere-lh_parcellation-BA_measure-meancurv.tsv'
    aparc_BA_rh_meancurv = 'hemisphere-rh_parcellation-BA_measure-meancurv.tsv'

    output_path = os.path.expanduser(output_dir)

    cs_dir = os.path.join(output_path, 'subjects')
    if not os.path.isdir(cs_dir):
        print("ERROR: directory freesurfer-cross-sectional does not exist, it should be CAPS directory after running recon_all_pipeline!!!")
    else:
        pass
    dest_dir = cs_dir + '/regional_measures_summary'
    try:
        os.makedirs(dest_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST: # if dest_dir exists, go on, if its other error, raise
            raise

    # fetch the paths for all the 26 tsv files.
    all_seg_volume_tsv = os.path.join(dest_dir, all_seg_volume)
    aseg_volume_tsv = os.path.join(dest_dir, aseg_volume)
    # DESIKAN atlas(?h.aparc.stats)
    aparc_desikan_lh_volume_tsv = os.path.join(dest_dir, aparc_desikan_lh_volume)
    aparc_desikan_rh_volume_tsv = os.path.join(dest_dir, aparc_desikan_rh_volume)
    aparc_desikan_lh_thickness_tsv = os.path.join(dest_dir, aparc_desikan_lh_thickness)
    aparc_desikan_rh_thickness_tsv = os.path.join(dest_dir, aparc_desikan_rh_thickness)
    aparc_desikan_lh_area_tsv = os.path.join(dest_dir, aparc_desikan_lh_area)
    aparc_desikan_rh_area_tsv = os.path.join(dest_dir, aparc_desikan_rh_area)
    aparc_desikan_lh_meancurv_tsv = os.path.join(dest_dir, aparc_desikan_lh_meancurv)
    aparc_desikan_rh_meancurv_tsv = os.path.join(dest_dir, aparc_desikan_rh_meancurv)
    # DESTRIEUX atals
    aparc_destrieux_lh_volume_tsv = os.path.join(dest_dir, aparc_destrieux_lh_volume)
    aparc_destrieux_rh_volume_tsv = os.path.join(dest_dir, aparc_destrieux_rh_volume)
    aparc_destrieux_lh_thickness_tsv = os.path.join(dest_dir, aparc_destrieux_lh_thickness)
    aparc_destrieux_rh_thickness_tsv = os.path.join(dest_dir, aparc_destrieux_rh_thickness)
    aparc_destrieux_lh_area_tsv = os.path.join(dest_dir, aparc_destrieux_lh_area)
    aparc_destrieux_rh_area_tsv = os.path.join(dest_dir, aparc_destrieux_rh_area)
    aparc_destrieux_lh_meancurv_tsv = os.path.join(dest_dir, aparc_destrieux_lh_meancurv)
    aparc_destrieux_rh_meancurv_tsv = os.path.join(dest_dir, aparc_destrieux_rh_meancurv)
    # Brodmann Area atlas
    aparc_BA_lh_volume_tsv = os.path.join(dest_dir, aparc_BA_lh_volume)
    aparc_BA_rh_volume_tsv = os.path.join(dest_dir, aparc_BA_rh_volume)
    aparc_BA_lh_thickness_tsv = os.path.join(dest_dir, aparc_BA_lh_thickness)
    aparc_BA_rh_thickness_tsv = os.path.join(dest_dir, aparc_BA_rh_thickness)
    aparc_BA_lh_area_tsv = os.path.join(dest_dir, aparc_BA_lh_area)
    aparc_BA_rh_area_tsv = os.path.join(dest_dir, aparc_BA_rh_area)
    aparc_BA_lh_meancurv_tsv = os.path.join(dest_dir, aparc_BA_lh_meancurv)
    aparc_BA_rh_meancurv_tsv = os.path.join(dest_dir, aparc_BA_rh_meancurv)

    # get the cmd string for the command line wrappers
    subjects = ''
    for i in xrange(len(subject_dir)):
        subject_path = (os.path.join(subject_dir[i], subject_id[i]))
        subjects += subject_path + ' '

    cmd_all_seg = 'asegstats2table --subjects ' + subjects + '--meas volume --statsfile wmparc.stats --all-seg --tablefile ' + all_seg_volume_tsv
    os.system(cmd_all_seg)
    cmd_aseg = 'asegstats2table --subjects ' + subjects + '--meas volume --tablefile ' + aseg_volume_tsv
    os.system(cmd_aseg)

    cmd_aparc_desikan_lh_volume = 'aparcstats2table --subjects ' + subjects + '--hemi lh --meas volume --tablefile ' + aparc_desikan_lh_volume_tsv
    os.system(cmd_aparc_desikan_lh_volume)
    cmd_aparc_desikan_rh_volume = 'aparcstats2table --subjects ' + subjects + '--hemi rh --meas volume --tablefile ' + aparc_desikan_rh_volume_tsv
    os.system(cmd_aparc_desikan_rh_volume)
    cmd_parc_desikan_lh_thickness = 'aparcstats2table --subjects ' + subjects + '--hemi lh --meas thickness --tablefile ' + aparc_desikan_lh_thickness_tsv
    os.system(cmd_parc_desikan_lh_thickness)
    cmd_parc_desikan_rh_thickness = 'aparcstats2table --subjects ' + subjects + '--hemi rh --meas thickness --tablefile ' + aparc_desikan_rh_thickness_tsv
    os.system(cmd_parc_desikan_rh_thickness)
    cmd_aparc_desikan_lh_area = 'aparcstats2table --subjects ' + subjects + '--hemi lh --meas area --tablefile ' + aparc_desikan_lh_area_tsv
    os.system(cmd_aparc_desikan_lh_area)
    cmd_aparc_desikan_rh_area = 'aparcstats2table --subjects ' + subjects + '--hemi rh --meas area --tablefile ' + aparc_desikan_rh_area_tsv
    os.system(cmd_aparc_desikan_rh_area)
    cmd_aparc_desikan_lh_meancurv = 'aparcstats2table --subjects ' + subjects + '--hemi lh --meas meancurv --tablefile ' + aparc_desikan_lh_meancurv_tsv
    os.system(cmd_aparc_desikan_lh_meancurv)
    cmd_aparc_desikan_rh_meancurv = 'aparcstats2table --subjects ' + subjects + '--hemi rh --meas meancurv --tablefile ' + aparc_desikan_rh_meancurv_tsv
    os.system(cmd_aparc_desikan_rh_meancurv)

    cmd_aparc_destrieux_lh_volume = 'aparcstats2table --subjects ' + subjects + '--hemi lh --parc aparc.a2009s --meas volume --tablefile ' + aparc_destrieux_lh_volume_tsv
    os.system(cmd_aparc_destrieux_lh_volume)
    cmd_aparc_destrieux_rh_volume = 'aparcstats2table --subjects ' + subjects + '--hemi rh --parc aparc.a2009s --meas volume --tablefile ' + aparc_destrieux_rh_volume_tsv
    os.system(cmd_aparc_destrieux_rh_volume)
    cmd_parc_destrieux_lh_thickness = 'aparcstats2table --subjects ' + subjects + '--hemi lh --parc aparc.a2009s --meas thickness --tablefile ' + aparc_destrieux_lh_thickness_tsv
    os.system(cmd_parc_destrieux_lh_thickness)
    cmd_parc_destrieux_rh_thickness = 'aparcstats2table --subjects ' + subjects + '--hemi rh --parc aparc.a2009s --meas thickness --tablefile ' + aparc_destrieux_rh_thickness_tsv
    os.system(cmd_parc_destrieux_rh_thickness)
    cmd_aparc_destrieux_lh_area = 'aparcstats2table --subjects ' + subjects + '--hemi lh --parc aparc.a2009s --meas area --tablefile ' + aparc_destrieux_lh_area_tsv
    os.system(cmd_aparc_destrieux_lh_area)
    cmd_aparc_destrieux_rh_area = 'aparcstats2table --subjects ' + subjects + '--hemi rh --parc aparc.a2009s --meas area --tablefile ' + aparc_destrieux_rh_area_tsv
    os.system(cmd_aparc_destrieux_rh_area)
    cmd_aparc_destrieux_lh_meancurv = 'aparcstats2table --subjects ' + subjects + '--hemi lh --parc aparc.a2009s --meas meancurv --tablefile ' + aparc_destrieux_lh_meancurv_tsv
    os.system(cmd_aparc_destrieux_lh_meancurv)
    cmd_aparc_destrieux_rh_meancurv = 'aparcstats2table --subjects ' + subjects + '--hemi rh --parc aparc.a2009s --meas meancurv --tablefile ' + aparc_destrieux_rh_meancurv_tsv
    os.system(cmd_aparc_destrieux_rh_meancurv)

    cmd_aparc_BA_lh_volume = 'aparcstats2table --subjects ' + subjects + '--hemi lh --parc BA --meas volume --tablefile ' + aparc_BA_lh_volume_tsv
    os.system(cmd_aparc_BA_lh_volume)
    cmd_aparc_BA_rh_volume = 'aparcstats2table --subjects ' + subjects + '--hemi rh --parc BA --meas volume --tablefile ' + aparc_BA_rh_volume_tsv
    os.system(cmd_aparc_BA_rh_volume)
    cmd_parc_BA_lh_thickness = 'aparcstats2table --subjects ' + subjects + '--hemi lh --parc BA --meas thickness --tablefile ' + aparc_BA_lh_thickness_tsv
    os.system(cmd_parc_BA_lh_thickness)
    cmd_parc_BA_rh_thickness = 'aparcstats2table --subjects ' + subjects + '--hemi rh --parc BA --meas thickness --tablefile ' + aparc_BA_rh_thickness_tsv
    os.system(cmd_parc_BA_rh_thickness)
    cmd_aparc_BA_lh_area = 'aparcstats2table --subjects ' + subjects + '--hemi lh --parc BA --meas area --tablefile ' + aparc_BA_lh_area_tsv
    os.system(cmd_aparc_BA_lh_area)
    cmd_aparc_BA_rh_area = 'aparcstats2table --subjects ' + subjects + '--hemi rh --parc BA --meas area --tablefile ' + aparc_BA_rh_area_tsv
    os.system(cmd_aparc_BA_rh_area)
    cmd_aparc_BA_lh_meancurv = 'aparcstats2table --subjects ' + subjects + '--hemi lh --parc BA --meas meancurv --tablefile ' + aparc_BA_lh_meancurv_tsv
    os.system(cmd_aparc_BA_lh_meancurv)
    cmd_aparc_BA_rh_meancurv = 'aparcstats2table --subjects ' + subjects + '--hemi rh --parc BA --meas meancurv --tablefile ' + aparc_BA_rh_meancurv_tsv
    os.system(cmd_aparc_BA_rh_meancurv)