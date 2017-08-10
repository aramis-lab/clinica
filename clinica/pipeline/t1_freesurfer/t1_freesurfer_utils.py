#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains useful functions used for t1_freesurfer pipeline"""

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Michael Bacci", "Junhao Wen"]
__license__ = "See LICENSE.txt file"
__version__ = "1.0.0"
__maintainer__ = "Junhao Wen"
__email__ = "junhao.Wen@inria.fr"
__status__ = "Development"


def bids_datagrabber(input_dir, subject_list, session_list):
    """
        Fetch t1 images from a BIDS directory based on subject_list and a session_list
    Args:
        input_dir: BIDS directory
        subject_list: a list containing all the participant_id
        session_list: a list containing all the session_id

    Returns: a list containing all the t1 images

    """
    from bids.grabbids.bids_layout import BIDSLayout

    bidslayout = BIDSLayout(input_dir)
    anat_t1 = []
    if not bidslayout.get(target='run', return_type='id', type='T1w'):
        print "There is just one run for T1w image of this analysis"
        for i in range(len(subject_list)):
            t1 = bidslayout.get(return_type='file',
                                            type='T1w',
                                            extensions=['nii|nii.gz'],
                                            session=session_list[i].replace('ses-', ''),
                                        subject=subject_list[i].replace('sub-', ''))
            anat_t1.append(t1)
    else:
        print "There are more than one runs for T1w image for this analysis"
        for i in range(len(subject_list)):
            t1 = bidslayout.get(return_type='file',
                                            type='T1w',
                                            extensions=['nii|nii.gz'],
                                            session=session_list[i].replace('ses-', ''),
                                        subject=subject_list[i].replace('sub-', ''),
                                     run='1')
            anat_t1.append(t1)


    if len(anat_t1) == 0:
        raise ValueError("you have to grap at least one image, but the result is empty, please check it out!")
    if len(anat_t1) != len(subject_list) or len(anat_t1) != len(session_list):
        raise ValueError("Pybids found some missing files, you should remove them out from your analysis!!!")

    return anat_t1

def get_dirs_check_reconalled(output_dir, subject_list, session_list):
    """
        Get the info from subjects_visits_tsv, like subject_dir, subject_id, subject_list, session_list
        Also, this func is also to check out the rerun of the dataset, if the subject result folder has been created, you should
        check out the result and decide if you are going to rerun it or just ignore it.

    Args:
        output_dir: CAPS directory to contain the output
        subject_list:  a list containing all the participant_id
        session_list: a list containing all the session_id

    Returns: the related lists based on the tsv files

    """
    import os, errno
    from copy import deepcopy as cp

    subject_id = list(subject_list[i] + '_' + session_list[i] for i in range(len(subject_list)))
    subject_id_without_reconalled = cp(subject_id)
    subject_list_without_reconalled = cp(subject_list)
    session_list_without_reconalled = cp(session_list)

    output_path = os.path.expanduser(output_dir)  # change the relative path to be absolute path
    output_dir = os.path.join(output_path, 'subjects')

    try:
        os.makedirs(output_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:  # if the error is not exist error, raise, otherwise, pass
            raise

    subject_dir = []
    subject_dir_without_reconalled = []

    for i in range(len(subject_list)):
        subject = output_dir + '/' + subject_list[i] + '/' + session_list[i] + '/' + 't1' + '/' + 'freesurfer_cross_sectional'
        try:
            os.makedirs(subject)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        subject_dir.append(subject)
        subject_path = os.path.join(subject, subject_id[i])
        subject_path_abs = os.path.expanduser(subject_path)
        if os.path.exists(subject_path_abs):
            print "Note: subject %s seems to be already recon-alled or being runing, please check out the result folder and" \
                  "recon-all sumarry log, in case that the processing of recon-all has been killed accidentally, please" \
                  "delete the result foder and rerun it; In case that the subject has been run successfully with recon-all," \
                  "just ignore this message and  continue to run the new-added or non-recon-alled subjects!!! " % subject_id[i]
            subject_id_without_reconalled.remove(subject_id[i])
            subject_list_without_reconalled.remove(subject_list[i])
            session_list_without_reconalled.remove(session_list[i])
        else:
            subject_dir_without_reconalled.append(subject)

    # # Make sure the subject_id, subject_list and session list have the same length


    return subject_dir, subject_id, subject_dir_without_reconalled, subject_id_without_reconalled, subject_list_without_reconalled, session_list_without_reconalled

def checkfov(t1_list, recon_all_args):
    """
        Verifying size of inputs and FOV of each T1 image

    Args:
        t1_list: a list containing all the t1 images
        recon_all_args: default the -qache flag

    Returns: a list containing -qcache + -cw256 or -qcache based on the FOV

    """
    import sys
    import nibabel as nib
    from nipype.utils.filemanip import filename_to_list

    output_flags = []
    num_t1 = len(t1_list)
    t1_list = filename_to_list(t1_list)

    if num_t1 == 0:
        print("ERROR: No T1's Given")
        sys.exit(-1)

    f = nib.load(t1_list[0])
    voxel_size = f.header.get_zooms()
    t1_size = f.header.get_data_shape()
    if (voxel_size[0] * t1_size[0] > 256) or (voxel_size[1] * t1_size[1]> 256) or (voxel_size[2] * t1_size[2]> 256):
        print("Setting MRI Convert to crop images to 256 FOV")
        optional_flag = '-cw256'
    else:
        print("No need to add -cw256 flag")
        optional_flag = ''
    flag = "{0} ".format(recon_all_args) + optional_flag
    output_flags.append(flag)

    return output_flags


def create_flags_str(input_flags):
    """
        Create a commandline string from a list of input flags

    Args:
        input_flags: the added flag for recon-all command

    Returns: converted string flag

    """
    output_str = ""
    for flag in input_flags:
        output_str += "{0} ".format(flag)
    output_str.strip()  # stripped from the beginning and the end of the string (default whitespace characters).

    return output_str

def log_summary(subject_list, session_list, subject_id, output_dir):
    """
        create a log file to summarize the recon-all result for all the subjects, the first step quality check

    Args:
        subject_list: a list containing all the path to the subject
        session_list: a list containing all the session_id
        subject_id: a list containing all the participant_id
        output_dir: CAPS directory

    Returns:

    """
    import os
    from datetime import datetime

    output_path = os.path.expanduser(output_dir)
    dest_dir = os.path.join(output_path, 'subjects')
    if not os.path.isdir(dest_dir):
        print("ERROR: directory subjects does not exist, it should be CAPS directory after running recon_all_pipeline!!!")
    else:
        pass
    log_name = os.path.join(dest_dir, 'recon_all_summary.log')
    input_logs = []

    for i in xrange(len(subject_list)):
        input_log = os.path.join(dest_dir, subject_list[i], session_list[i], 't1', 'freesurfer_cross_sectional', subject_id[i], 'scripts', 'recon-all-status.log' )
        input_logs.append(input_log)

    bad_log = 0
    search_query = 'recon-all -s'
    with open(log_name, 'w') as f1:
        line1 = datetime.now().isoformat()
        line2 = 'Quality check: recon-all output summary'
        f1.write("%s\n%s\n\n" % (line1, line2))
        for log in input_logs:
            with open(log, 'r') as f2:
                lines = f2.readlines()
                for line in lines:
                    if (line.startswith(search_query)) and ('without error' not in line):
                        f1.write(line)
                        bad_log += 1
                        break
                    elif line.startswith(search_query):
                            f1.write(line)
                    else:
                        pass
        line3 = 'Number of subjects: %s \nNumber of bad recon-all is: %s ' % (len(subject_list), bad_log)
        f1.write(line3)

def write_statistics_per_subject(subject_id, output_dir):
    """
        A custom function to write the statistical measures into tsv files for each subject

    Args:
        subject_id: the subject's particiapnt_id
        output_dir: CAPS directory

    Returns:

    """

    import os, errno

    subject_list = subject_id.split('_')[0]
    session_list = subject_id.split('_')[1]
    # name all the 26 tsv output files.
    all_seg_volume = subject_id + '_parcellation-wm_volume.tsv'
    aseg_volume = subject_id + '_segmentationVolumes.tsv'

    aparc_desikan_lh_volume = subject_id + '_hemi-left_parcellation-desikan_volume.tsv'
    aparc_desikan_rh_volume = subject_id + '_hemi-right_parcellation-desikan_volume.tsv'
    aparc_desikan_lh_thickness = subject_id + '_hemi-left_parcellation-desikan_thickness.tsv'
    aparc_desikan_rh_thickness = subject_id + '_hemi-right_parcellation-desikan_thickness.tsv'
    aparc_desikan_lh_area = subject_id + '_hemi-left_parcellation-desikan_area.tsv'
    aparc_desikan_rh_area = subject_id + '_hemi-right_parcellation-desikan_area.tsv'
    aparc_desikan_lh_meancurv = subject_id + '_hemi-left_parcellation-desikan_meancurv.tsv'
    aparc_desikan_rh_meancurv = subject_id + '_hemi-right_parcellation-desikan_meancurv.tsv'

    aparc_destrieux_lh_volume = subject_id + '_hemi-left_parcellation-destrieux_volume.tsv'
    aparc_destrieux_rh_volume = subject_id + '_hemi-right_parcellation-destrieux_volume.tsv'
    aparc_destrieux_lh_thickness = subject_id + '_hemi-left_parcellation-destrieux_thickness.tsv'
    aparc_destrieux_rh_thickness = subject_id + '_hemi-right_parcellation-destrieux_thickness.tsv'
    aparc_destrieux_lh_area = subject_id + '_hemi-left_parcellation-destrieux_area.tsv'
    aparc_destrieux_rh_area = subject_id + '_hemi-right_parcellation-destrieux_area.tsv'
    aparc_destrieux_lh_meancurv = subject_id + '_hemi-left_parcellation-destrieux_meancurv.tsv'
    aparc_destrieux_rh_meancurv = subject_id + '_hemi-right_parcellation-destrieux_meancurv.tsv'

    aparc_BA_lh_volume = subject_id + '_hemi-left_parcellation-ba_volume.tsv'
    aparc_BA_rh_volume = subject_id + '_hemi-right_parcellation-ba_volume.tsv'
    aparc_BA_lh_thickness = subject_id + '_hemi-left_parcellation-ba_thickness.tsv'
    aparc_BA_rh_thickness = subject_id + '_hemi-right_parcellation-ba_thickness.tsv'
    aparc_BA_lh_area = subject_id + '_hemi-left_parcellation-ba_area.tsv'
    aparc_BA_rh_area = subject_id + '_hemi-right_parcellation-ba_area.tsv'
    aparc_BA_lh_meancurv = subject_id + '_hemi-left_parcellation-ba_meancurv.tsv'
    aparc_BA_rh_meancurv = subject_id + '_hemi-right_parcellation-ba_meancurv.tsv'

    # subject_name = subject_list + '_' + session_list
    output_path = os.path.expanduser(output_dir)
    cs_dir = os.path.join(output_path, 'subjects', subject_list, session_list, 't1', 'freesurfer_cross_sectional')
    if not os.path.isdir(cs_dir):
        print("ERROR: directory freesurfer_cross_sectional does not exist, it should be CAPS directory after running recon_all_pipeline!!!")
    else:
        pass
    dest_dir = cs_dir + '/regional_measures'
    try:
        os.makedirs(dest_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST: # if dest_dir exists, go on, if its other error, raise
            raise
    subject = os.path.join(cs_dir, subject_id)

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
    cmd_all_seg = 'asegstats2table --subjects ' + subject + ' --meas volume --statsfile wmparc.stats --all-seg --tablefile ' + all_seg_volume_tsv
    os.system(cmd_all_seg)
    cmd_aseg = 'asegstats2table --subjects ' + subject + ' --meas volume --tablefile ' + aseg_volume_tsv
    os.system(cmd_aseg)

    cmd_aparc_desikan_lh_volume = 'aparcstats2table --subjects ' + subject + ' --hemi lh --meas volume --tablefile ' + aparc_desikan_lh_volume_tsv
    os.system(cmd_aparc_desikan_lh_volume)
    cmd_aparc_desikan_rh_volume = 'aparcstats2table --subjects ' + subject + ' --hemi rh --meas volume --tablefile ' + aparc_desikan_rh_volume_tsv
    os.system(cmd_aparc_desikan_rh_volume)
    cmd_parc_desikan_lh_thickness = 'aparcstats2table --subjects ' + subject + ' --hemi lh --meas thickness --tablefile ' + aparc_desikan_lh_thickness_tsv
    os.system(cmd_parc_desikan_lh_thickness)
    cmd_parc_desikan_rh_thickness = 'aparcstats2table --subjects ' + subject + ' --hemi rh --meas thickness --tablefile ' + aparc_desikan_rh_thickness_tsv
    os.system(cmd_parc_desikan_rh_thickness)
    cmd_aparc_desikan_lh_area = 'aparcstats2table --subjects ' + subject + ' --hemi lh --meas area --tablefile ' + aparc_desikan_lh_area_tsv
    os.system(cmd_aparc_desikan_lh_area)
    cmd_aparc_desikan_rh_area = 'aparcstats2table --subjects ' + subject + ' --hemi rh --meas area --tablefile ' + aparc_desikan_rh_area_tsv
    os.system(cmd_aparc_desikan_rh_area)
    cmd_aparc_desikan_lh_meancurv = 'aparcstats2table --subjects ' + subject + ' --hemi lh --meas meancurv --tablefile ' + aparc_desikan_lh_meancurv_tsv
    os.system(cmd_aparc_desikan_lh_meancurv)
    cmd_aparc_desikan_rh_meancurv = 'aparcstats2table --subjects ' + subject + ' --hemi rh --meas meancurv --tablefile ' + aparc_desikan_rh_meancurv_tsv
    os.system(cmd_aparc_desikan_rh_meancurv)

    cmd_aparc_destrieux_lh_volume = 'aparcstats2table --subjects ' + subject + ' --hemi lh --parc aparc.a2009s --meas volume --tablefile ' + aparc_destrieux_lh_volume_tsv
    os.system(cmd_aparc_destrieux_lh_volume)
    cmd_aparc_destrieux_rh_volume = 'aparcstats2table --subjects ' + subject + ' --hemi rh --parc aparc.a2009s --meas volume --tablefile ' + aparc_destrieux_rh_volume_tsv
    os.system(cmd_aparc_destrieux_rh_volume)
    cmd_parc_destrieux_lh_thickness = 'aparcstats2table --subjects ' + subject + ' --hemi lh --parc aparc.a2009s --meas thickness --tablefile ' + aparc_destrieux_lh_thickness_tsv
    os.system(cmd_parc_destrieux_lh_thickness)
    cmd_parc_destrieux_rh_thickness = 'aparcstats2table --subjects ' + subject + ' --hemi rh --parc aparc.a2009s --meas thickness --tablefile ' + aparc_destrieux_rh_thickness_tsv
    os.system(cmd_parc_destrieux_rh_thickness)
    cmd_aparc_destrieux_lh_area = 'aparcstats2table --subjects ' + subject + ' --hemi lh --parc aparc.a2009s --meas area --tablefile ' + aparc_destrieux_lh_area_tsv
    os.system(cmd_aparc_destrieux_lh_area)
    cmd_aparc_destrieux_rh_area = 'aparcstats2table --subjects ' + subject + ' --hemi rh --parc aparc.a2009s --meas area --tablefile ' + aparc_destrieux_rh_area_tsv
    os.system(cmd_aparc_destrieux_rh_area)
    cmd_aparc_destrieux_lh_meancurv = 'aparcstats2table --subjects ' + subject + ' --hemi lh --parc aparc.a2009s --meas meancurv --tablefile ' + aparc_destrieux_lh_meancurv_tsv
    os.system(cmd_aparc_destrieux_lh_meancurv)
    cmd_aparc_destrieux_rh_meancurv = 'aparcstats2table --subjects ' + subject + ' --hemi rh --parc aparc.a2009s --meas meancurv --tablefile ' + aparc_destrieux_rh_meancurv_tsv
    os.system(cmd_aparc_destrieux_rh_meancurv)

    cmd_aparc_BA_lh_volume = 'aparcstats2table --subjects ' + subject + ' --hemi lh --parc BA --meas volume --tablefile ' + aparc_BA_lh_volume_tsv
    os.system(cmd_aparc_BA_lh_volume)
    cmd_aparc_BA_rh_volume = 'aparcstats2table --subjects ' + subject + ' --hemi rh --parc BA --meas volume --tablefile ' + aparc_BA_rh_volume_tsv
    os.system(cmd_aparc_BA_rh_volume)
    cmd_parc_BA_lh_thickness = 'aparcstats2table --subjects ' + subject + ' --hemi lh --parc BA --meas thickness --tablefile ' + aparc_BA_lh_thickness_tsv
    os.system(cmd_parc_BA_lh_thickness)
    cmd_parc_BA_rh_thickness = 'aparcstats2table --subjects ' + subject + ' --hemi rh --parc BA --meas thickness --tablefile ' + aparc_BA_rh_thickness_tsv
    os.system(cmd_parc_BA_rh_thickness)
    cmd_aparc_BA_lh_area = 'aparcstats2table --subjects ' + subject + ' --hemi lh --parc BA --meas area --tablefile ' + aparc_BA_lh_area_tsv
    os.system(cmd_aparc_BA_lh_area)
    cmd_aparc_BA_rh_area = 'aparcstats2table --subjects ' + subject + ' --hemi rh --parc BA --meas area --tablefile ' + aparc_BA_rh_area_tsv
    os.system(cmd_aparc_BA_rh_area)
    cmd_aparc_BA_lh_meancurv = 'aparcstats2table --subjects ' + subject + ' --hemi lh --parc BA --meas meancurv --tablefile ' + aparc_BA_lh_meancurv_tsv
    os.system(cmd_aparc_BA_lh_meancurv)
    cmd_aparc_BA_rh_meancurv = 'aparcstats2table --subjects ' + subject + ' --hemi rh --parc BA --meas meancurv --tablefile ' + aparc_BA_rh_meancurv_tsv
    os.system(cmd_aparc_BA_rh_meancurv)

    print "Writing statistical data to tsv file for %s finished!" % subject
