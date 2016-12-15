#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains functions used for the recon_all_pipeline() and recon_all_statistics_pipeline()"""
import os

# Function for recon_all_pipeline()

def absolute_path(arg):
    """
    Transfer any path to absolute path

    :param arg:
    :return:
    """

    if arg[:1] == '~':
        return os.path.expanduser(arg)
    elif arg[:1] == '.':
        return os.getcwd()
    else:
        return os.path.join(os.getcwd(), arg)


def get_dirs(output_dir, subjects_visits_tsv, analysis_series_id):
    """
    Define and create the CAPS output for recon_all_pipeline()

    :param output_dir:
    :param subjects_visits_tsv:
    :param analysis_series_id:
    :return:
    """
    import os, csv, errno
    subject_list = []
    session_list = []
    subject_id = []
    with open(subjects_visits_tsv, 'rb') as tsvin:
        tsv_reader = csv.reader(tsvin, delimiter='\t')

        for row in tsv_reader:
            if row[0] == 'participant_id':
                continue
            else:
                subject_list.append(row[0])
                session_list.append(row[1])
                subject_id.append(row[0] + '_' + row[1])

    output_path = os.path.expanduser(output_dir)  # change the relative path to be absolute path
    output_base = 'analysis-series-' + analysis_series_id + '/subjects'
    if output_path[-1] == '/':
        output_dir = output_path + output_base
    else:
        output_dir = output_path + '/' + output_base
    try:
        os.makedirs(output_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:  # if the error is not exist error, raise, otherwise, pass
            raise

    subject_dir = []
    num_subject = len(subject_list)
    for i in xrange(num_subject):
        subject = output_dir + '/' + subject_list[i] + '/' + session_list[i] + '/' + 't1' + '/' + 'freesurfer-cross-sectional'
        try:
            os.makedirs(subject)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        subject_dir.append(subject)

    return subject_dir, subject_id, subject_list, session_list

def checkfov(t1_list, recon_all_args):
    """
    Verifying size of inputs and FOV of each T1 image

    Note:node2mapnode, so every subject is running in parallel, we dont have to check out if they have the same SIZE, but if you
    node2node, it will be serialized, so that we can compare their size.

    :param t1_list:
    :param recon_all_args:
    :return:
    """
    import sys
    import nibabel as nib
    from nipype.utils.filemanip import filename_to_list

    # # this is the old version for node2node connection.
    # output_flags = []
    # num_t1 = len(t1_list)
    # t1_list = filename_to_list(t1_list)
    #
    # if num_t1 == 0:
    #     print("ERROR: No T1's Given")
    #     sys.exit(-1)
    #
    # # shape = nib.load(t1_list[0]).shape
    # for t1 in t1_list:
    #     f = nib.load(t1)
    #     voxel_size = f.header.get_zooms()
    #     t1_size = f.header.get_data_shape()
    #     # not sure if we should constrain all the T1 file should have the same size
    #     # if t1_size != shape:
    #     #     print("ERROR: T1s not the same size. Cannot process {0} and {1} "
    #     #           "together".format(t1_list[0], t1))
    #     #     sys.exit(-1)
    #     if voxel_size[0] * t1_size[0] > 256 or voxel_size[1] * t1_size[1] or voxel_size[2] * t1_size[2]:
    #         print("Setting MRI Convert to crop images to 256 FOV")
    #         optional_flag = '-cw256'
    #     else:
    #         print("No need to add -cw256 flag")
    #         optional_flag = ''
    #     flag = "{0} ".format(recon_all_args) + optional_flag
    #     output_flags.append(flag)


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

    :param input_flags:
    :return:
    """
    output_str = ""
    for flag in input_flags:
        output_str += "{0} ".format(flag)
    output_str.strip()  # stripped from the beginning and the end of the string (default whitespace characters).

    return output_str

def log_summary(subject_list, session_list, subject_id, output_dir, analysis_series_id):
    """
    create the txt file to summarize the reconall result for all the subjects

    :param subject_list:
    :param session_list:
    :param subject_id:
    :param output_dir:
    :param analysis_series_id:
    :return:
    """
    ## TODO check if log file exits, if yes, add new info for new subjects, not to overwrite it.
    import os
    from datetime import datetime
    # from nipype import config, logging

    output_path = os.path.expanduser(output_dir)
    dest_dir = output_path + '/analysis-series-' + analysis_series_id + '/subjects'
    if not os.path.isdir(dest_dir):
        print("ERROR: directory subjects does not exist, it should be CAPS directory after running recon_all_pipeline!!!")
    else:
        pass
    log_name = os.path.join(dest_dir, 'recon_all_summary.log')
    input_logs = []

    for i in xrange(len(subject_list)):
        input_log = os.path.join(dest_dir, subject_list[i], session_list[i], 't1', 'freesurfer-cross-sectional', subject_id[i], 'scripts', 'recon-all-status.log' )
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

    # logging.basicConfig(filename=log_name, format='%(asctime)s %(levelname)s:%(message)s',
    #                     datefmt='%m/%d/%Y %I:%M', level=logging.DEBUG)
    # config.update_config({'logging': {'log_directory': dest_dir,
    #                                   'log_to_file': True, 'workflow_level': 'DEBUG'}})
    # logging.update_logging(config)
    # bad_log = 0
    # line1 = 'Quality check: recon-all output summary'
    # logging.info(line1)
    # for log in input_logs:
    #     with open(log, 'r') as f2:
    #         line = f2.readlines()[-1]
    #         if 'without error' in line:
    #             logging.info(line)
    #         else:
    #             logging.warning(line)
    #             bad_log += 1
    #         f2.close()
    # line2 = 'Number of subjects: %s \nNumber of bad recon-all is: %s ' % (len(subject_list), bad_log)
    # logging.info(line2)

def write_statistics(subject_dir, subject_id, analysis_series_id, output_dir):
    """

    :param subject_dir:
    :param subject_id:
    :param analysis_series_id:
    :param output_dir:
    :return:
    """
    ## TODO check if summary_tsv exits, if yes, add new info for new subjects, not to overwrite it.

    import os, errno

    # name all the 26 tsv output files.
    all_seg_volume = 'all_seg.tsv'
    aseg_volume = 'aseg_volume.tsv'

    aparc_desikan_lh_volume = 'lh_aparc_desikan_volume.tsv'
    aparc_desikan_rh_volume = 'rh_aparc_desikan_volume.tsv'
    aparc_desikan_lh_thickness = 'lh_aparc_desikan_thickness.tsv'
    aparc_desikan_rh_thickness = 'rh_aparc_desikan_thickness.tsv'
    aparc_desikan_lh_area = 'lh_aparc_desikan_area.tsv'
    aparc_desikan_rh_area = 'rh_aparc_desikan_area.tsv'
    aparc_desikan_lh_meancurv = 'lh_aparc_desikan_meancurv.tsv'
    aparc_desikan_rh_meancurv = 'rh_aparc_desikan_meancurv.tsv'

    aparc_destrieux_lh_volume = 'lh_aparc_destrieux_volume.tsv'
    aparc_destrieux_rh_volume = 'rh_aparc_destrieux_volume.tsv'
    aparc_destrieux_lh_thickness = 'lh_aparc_destrieux_thickness.tsv'
    aparc_destrieux_rh_thickness = 'rh_aparc_destrieux_thickness.tsv'
    aparc_destrieux_lh_area = 'lh_aparc_destrieux_area.tsv'
    aparc_destrieux_rh_area = 'rh_aparc_destrieux_area.tsv'
    aparc_destrieux_lh_meancurv = 'lh_aparc_destrieux_meancurv.tsv'
    aparc_destrieux_rh_meancurv = 'rh_aparc_destrieux_meancurv.tsv'

    aparc_BA_lh_volume = 'lh_aparc_BA_volume.tsv'
    aparc_BA_rh_volume = 'rh_aparc_BA_volume.tsv'
    aparc_BA_lh_thickness = 'lh_aparc_BA_thickness.tsv'
    aparc_BA_rh_thickness = 'rh_aparc_BA_thickness.tsv'
    aparc_BA_lh_area = 'lh_aparc_BA_area.tsv'
    aparc_BA_rh_area = 'rh_aparc_BA_area.tsv'
    aparc_BA_lh_meancurv = 'lh_aparc_BA_meancurv.tsv'
    aparc_BA_rh_meancurv = 'rh_aparc_BA_meancurv.tsv'

    output_path = os.path.expanduser(output_dir)

    cs_dir = output_path + '/analysis-series-' + analysis_series_id + '/subjects'
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

def write_statistics_per_subject(subject_id, analysis_series_id, output_dir):

    import os, errno

    subject_list = subject_id.split('_')[0]
    session_list = subject_id.split('_')[1]
    # name all the 26 tsv output files.
    all_seg_volume = subject_id + '_all_seg.tsv'
    aseg_volume = subject_id + '_aseg_volume.tsv'

    aparc_desikan_lh_volume = subject_id + '_lh_aparc_desikan_volume.tsv'
    aparc_desikan_rh_volume = subject_id + '_rh_aparc_desikan_volume.tsv'
    aparc_desikan_lh_thickness = subject_id + '_lh_aparc_desikan_thickness.tsv'
    aparc_desikan_rh_thickness = subject_id + '_rh_aparc_desikan_thickness.tsv'
    aparc_desikan_lh_area = subject_id + '_lh_aparc_desikan_area.tsv'
    aparc_desikan_rh_area = subject_id + '_rh_aparc_desikan_area.tsv'
    aparc_desikan_lh_meancurv = subject_id + '_lh_aparc_desikan_meancurv.tsv'
    aparc_desikan_rh_meancurv = subject_id + '_rh_aparc_desikan_meancurv.tsv'

    aparc_destrieux_lh_volume = subject_id + '_lh_aparc_destrieux_volume.tsv'
    aparc_destrieux_rh_volume = subject_id + '_rh_aparc_destrieux_volume.tsv'
    aparc_destrieux_lh_thickness = subject_id + '_lh_aparc_destrieux_thickness.tsv'
    aparc_destrieux_rh_thickness = subject_id + '_rh_aparc_destrieux_thickness.tsv'
    aparc_destrieux_lh_area = subject_id + '_lh_aparc_destrieux_area.tsv'
    aparc_destrieux_rh_area = subject_id + '_rh_aparc_destrieux_area.tsv'
    aparc_destrieux_lh_meancurv = subject_id + '_lh_aparc_destrieux_meancurv.tsv'
    aparc_destrieux_rh_meancurv = subject_id + '_rh_aparc_destrieux_meancurv.tsv'

    aparc_BA_lh_volume = subject_id + '_lh_aparc_BA_volume.tsv'
    aparc_BA_rh_volume = subject_id + '_rh_aparc_BA_volume.tsv'
    aparc_BA_lh_thickness = subject_id + '_lh_aparc_BA_thickness.tsv'
    aparc_BA_rh_thickness = subject_id + '_rh_aparc_BA_thickness.tsv'
    aparc_BA_lh_area = subject_id + '_lh_aparc_BA_area.tsv'
    aparc_BA_rh_area = subject_id + '_rh_aparc_BA_area.tsv'
    aparc_BA_lh_meancurv = subject_id + '_lh_aparc_BA_meancurv.tsv'
    aparc_BA_rh_meancurv = subject_id + '_rh_aparc_BA_meancurv.tsv'

    # subject_name = subject_list + '_' + session_list
    output_path = os.path.expanduser(output_dir)

    cs_dir = output_path + '/analysis-series-' + analysis_series_id + '/subjects/' + subject_list + '/' + session_list + '/t1/freesurfer-cross-sectional'
    if not os.path.isdir(cs_dir):
        print("ERROR: directory freesurfer-cross-sectional does not exist, it should be CAPS directory after running recon_all_pipeline!!!")
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

