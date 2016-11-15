#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains functions used for the recon_all_pipeline() and recon_all_statistics_pipeline()"""
import os, errno, csv

# Function for recon_all_pipeline()

def absolute_path(arg):
    """Transfer any path to absolute path"""

    if arg[:1] == '~':
        return os.path.expanduser(arg)
    elif arg[:1] == '.':
        return os.getcwd()
    else:
        return os.path.join(os.getcwd(), arg)


def CAPS_output(output_dir, subjects_visits_tsv, analysis_series_id):
    """Define and create the CAPS output for recon_all_pipeline()"""
    subject_list = []
    session_list = []
    subject_id = []
    with open(subjects_visits_tsv, 'rb') as tsvin:
        tsv_reader = csv.reader(tsvin, delimiter='\t')

        for row in tsv_reader:
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
    for i in range(num_subject):
        subject = output_dir + '/' + subject_list[i] + '/' + session_list[i] + '/' + 't1' + '/' + 'freesurfer-cross-sectional'
        try:
            os.makedirs(subject)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        subject_dir.append(subject)

    return subject_dir, subject_id, subject_list, session_list

def checkfov(t1_list, recon_all_args):
    """Verifying size of inputs and FOV of each T1 image

    Note:node2mapnode, so every subject is running in parallel, we dont have to check out if they have the same SIZE, but if you
    node2node, it will be serialized, so that we can compare their size.
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
    if voxel_size[0] * t1_size[0] > 256 or voxel_size[1] * t1_size[1] or voxel_size[2] * t1_size[2]:
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
    """
    output_str = ""
    for flag in input_flags:
        output_str += "{0} ".format(flag)
    output_str.strip()  # stripped from the beginning and the end of the string (default whitespace characters).

    return output_str

# Function for recon_all_statistics_pipeline()

def get_vars(subjects_visits_tsv):
    """fetch some necessary vars for this pipeline"""
    subject_list = []
    session_list = []
    with open(subjects_visits_tsv, 'rb') as tsvin:
        tsv_reader = csv.reader(tsvin, delimiter='\t')

        for row in tsv_reader:
            subject_list.append(row[0])
            session_list.append(row[1])

    return subject_list, session_list



def freesurferstatas_to_tsv(subject,
                            all_seg_volume_tsv,
                            aseg_volume_tsv,
                            aparc_desikan_lh_volume_tsv,
                            aparc_desikan_rh_volume_tsv,
                            aparc_desikan_lh_thickness_tsv,
                            aparc_desikan_rh_thickness_tsv,
                            aparc_desikan_lh_area_tsv,
                            aparc_desikan_rh_area_tsv,
                            aparc_desikan_lh_meancurv_tsv,
                            aparc_desikan_rh_meancurv_tsv,
                            aparc_destrieux_lh_volume_tsv,
                            aparc_destrieux_rh_volume_tsv,
                            aparc_destrieux_lh_thickness_tsv,
                            aparc_destrieux_rh_thickness_tsv,
                            aparc_destrieux_lh_area_tsv,
                            aparc_destrieux_rh_area_tsv,
                            aparc_destrieux_lh_meancurv_tsv,
                            aparc_destrieux_rh_meancurv_tsv,
                            aparc_BA_lh_volume_tsv,
                            aparc_BA_rh_volume_tsv,
                            aparc_BA_lh_thickness_tsv,
                            aparc_BA_rh_thickness_tsv,
                            aparc_BA_lh_area_tsv,
                            aparc_BA_rh_area_tsv,
                            aparc_BA_lh_meancurv_tsv,
                            aparc_BA_rh_meancurv_tsv):
    """
    wrapper FS command line freesurfer2table.

    This script will generate text/ascii tables of freesurfer aseg stats data, aseg.stats.
    This can then be easily imported into a spreadsheet and/or stats program.
    """

    cmd_all_seg = 'asegstats2table -s ' + subject + ' --meas volume --statsfile wmparc.stats --all-seg --tablefile ' + all_seg_volume_tsv
    os.system(cmd_all_seg)
    cmd_aseg = 'asegstats2table -s ' + subject + ' --meas volume --tablefile ' + aseg_volume_tsv
    os.system(cmd_aseg)

    cmd_aparc_desikan_lh_volume = 'aparcstats2table -s ' + subject + ' --hemi lh --meas volume --tablefile ' + aparc_desikan_lh_volume_tsv
    os.system(cmd_aparc_desikan_lh_volume)
    cmd_aparc_desikan_rh_volume = 'aparcstats2table -s ' + subject + ' --hemi rh --meas volume --tablefile ' + aparc_desikan_rh_volume_tsv
    os.system(cmd_aparc_desikan_rh_volume)
    cmd_parc_desikan_lh_thickness = 'aparcstats2table -s ' + subject + ' --hemi lh --meas thickness --tablefile ' + aparc_desikan_lh_thickness_tsv
    os.system(cmd_parc_desikan_lh_thickness)
    cmd_parc_desikan_rh_thickness = 'aparcstats2table -s ' + subject + ' --hemi rh --meas thickness --tablefile ' + aparc_desikan_rh_thickness_tsv
    os.system(cmd_parc_desikan_rh_thickness)
    cmd_aparc_desikan_lh_area = 'aparcstats2table -s ' + subject + ' --hemi lh --meas area --tablefile ' + aparc_desikan_lh_area_tsv
    os.system(cmd_aparc_desikan_lh_area)
    cmd_aparc_desikan_rh_area = 'aparcstats2table -s ' + subject + ' --hemi rh --meas area --tablefile ' + aparc_desikan_rh_area_tsv
    os.system(cmd_aparc_desikan_rh_area)
    cmd_aparc_desikan_lh_meancurv = 'aparcstats2table -s ' + subject + ' --hemi lh --meas meancurv --tablefile ' + aparc_desikan_lh_meancurv_tsv
    os.system(cmd_aparc_desikan_lh_meancurv)
    cmd_aparc_desikan_rh_meancurv = 'aparcstats2table -s ' + subject + ' --hemi rh --meas meancurv --tablefile ' + aparc_desikan_rh_meancurv_tsv
    os.system(cmd_aparc_desikan_rh_meancurv)

    cmd_aparc_destrieux_lh_volume = 'aparcstats2table -s ' + subject + ' --hemi lh --parc aparc.a2009s --meas volume --tablefile ' + aparc_destrieux_lh_volume_tsv
    os.system(cmd_aparc_destrieux_lh_volume)
    cmd_aparc_destrieux_rh_volume = 'aparcstats2table -s ' + subject + ' --hemi rh --parc aparc.a2009s --meas volume --tablefile ' + aparc_destrieux_rh_volume_tsv
    os.system(cmd_aparc_destrieux_rh_volume)
    cmd_parc_destrieux_lh_thickness = 'aparcstats2table -s ' + subject + ' --hemi lh --parc aparc.a2009s --meas thickness --tablefile ' + aparc_destrieux_lh_thickness_tsv
    os.system(cmd_parc_destrieux_lh_thickness)
    cmd_parc_destrieux_rh_thickness = 'aparcstats2table -s ' + subject + ' --hemi rh --parc aparc.a2009s --meas thickness --tablefile ' + aparc_destrieux_rh_thickness_tsv
    os.system(cmd_parc_destrieux_rh_thickness)
    cmd_aparc_destrieux_lh_area = 'aparcstats2table -s ' + subject + ' --hemi lh --parc aparc.a2009s --meas area --tablefile ' + aparc_destrieux_lh_area_tsv
    os.system(cmd_aparc_destrieux_lh_area)
    cmd_aparc_destrieux_rh_area = 'aparcstats2table -s ' + subject + ' --hemi rh --parc aparc.a2009s --meas area --tablefile ' + aparc_destrieux_rh_area_tsv
    os.system(cmd_aparc_destrieux_rh_area)
    cmd_aparc_destrieux_lh_meancurv = 'aparcstats2table -s ' + subject + ' --hemi lh --parc aparc.a2009s --meas meancurv --tablefile ' + aparc_destrieux_lh_meancurv_tsv
    os.system(cmd_aparc_destrieux_lh_meancurv)
    cmd_aparc_destrieux_rh_meancurv = 'aparcstats2table -s ' + subject + ' --hemi rh --parc aparc.a2009s --meas meancurv --tablefile ' + aparc_destrieux_rh_meancurv_tsv
    os.system(cmd_aparc_destrieux_rh_meancurv)

    cmd_aparc_BA_lh_volume = 'aparcstats2table -s ' + subject + ' --hemi lh --parc BA --meas volume --tablefile ' + aparc_BA_lh_volume_tsv
    os.system(cmd_aparc_BA_lh_volume)
    cmd_aparc_BA_rh_volume = 'aparcstats2table -s ' + subject + ' --hemi rh --parc BA --meas volume --tablefile ' + aparc_BA_rh_volume_tsv
    os.system(cmd_aparc_BA_rh_volume)
    cmd_parc_BA_lh_thickness = 'aparcstats2table -s ' + subject + ' --hemi lh --parc BA --meas thickness --tablefile ' + aparc_BA_lh_thickness_tsv
    os.system(cmd_parc_BA_lh_thickness)
    cmd_parc_BA_rh_thickness = 'aparcstats2table -s ' + subject + ' --hemi rh --parc BA --meas thickness --tablefile ' + aparc_BA_rh_thickness_tsv
    os.system(cmd_parc_BA_rh_thickness)
    cmd_aparc_BA_lh_area = 'aparcstats2table -s ' + subject + ' --hemi lh --parc BA --meas area --tablefile ' + aparc_BA_lh_area_tsv
    os.system(cmd_aparc_BA_lh_area)
    cmd_aparc_BA_rh_area = 'aparcstats2table -s ' + subject + ' --hemi rh --parc BA --meas area --tablefile ' + aparc_BA_rh_area_tsv
    os.system(cmd_aparc_BA_rh_area)
    cmd_aparc_BA_lh_meancurv = 'aparcstats2table -s ' + subject + ' --hemi lh --parc BA --meas meancurv --tablefile ' + aparc_BA_lh_meancurv_tsv
    os.system(cmd_aparc_BA_lh_meancurv)
    cmd_aparc_BA_rh_meancurv = 'aparcstats2table -s ' + subject + ' --hemi rh --parc BA --meas meancurv --tablefile ' + aparc_BA_rh_meancurv_tsv
    os.system(cmd_aparc_BA_rh_meancurv)

    print "Writing statistical data to tsv file for %s finished!" % subject


def write_statistics(subject_list, session_list, analysis_series_id, output_dir):
    """wrapper for FS command line freesurferstats2table
    Note:node2mapnode, so every subject is running in parallel, this will be faster than node2node which is serialized.
    """

    import os, errno
    from clinica.pipeline.t1.t1_freesurfer_utils import freesurferstatas_to_tsv

    all_seg_volume = subject_list + '_' + session_list + '_all_seg.tsv'
    aseg_volume = subject_list + '_' + session_list + '_aseg_volume.tsv'

    aparc_desikan_lh_volume = subject_list + '_' + session_list + '_lh_aparc_desikan_volume.tsv'
    aparc_desikan_rh_volume = subject_list + '_' + session_list + '_rh_aparc_desikan_volume.tsv'
    aparc_desikan_lh_thickness = subject_list + '_' + session_list + '_lh_aparc_desikan_thickness.tsv'
    aparc_desikan_rh_thickness = subject_list + '_' + session_list + '_rh_aparc_desikan_thickness.tsv'
    aparc_desikan_lh_area = subject_list + '_' + session_list + '_lh_aparc_desikan_area.tsv'
    aparc_desikan_rh_area = subject_list + '_' + session_list + '_rh_aparc_desikan_area.tsv'
    aparc_desikan_lh_meancurv = subject_list + '_' + session_list + '_lh_aparc_desikan_meancurv.tsv'
    aparc_desikan_rh_meancurv = subject_list + '_' + session_list + '_rh_aparc_desikan_meancurv.tsv'

    aparc_destrieux_lh_volume = subject_list + '_' + session_list + '_lh_aparc_destrieux_volume.tsv'
    aparc_destrieux_rh_volume = subject_list + '_' + session_list + '_rh_aparc_destrieux_volume.tsv'
    aparc_destrieux_lh_thickness = subject_list + '_' + session_list + '_lh_aparc_destrieux_thickness.tsv'
    aparc_destrieux_rh_thickness = subject_list + '_' + session_list + '_rh_aparc_destrieux_thickness.tsv'
    aparc_destrieux_lh_area = subject_list + '_' + session_list + '_lh_aparc_destrieux_area.tsv'
    aparc_destrieux_rh_area = subject_list + '_' + session_list + '_rh_aparc_destrieux_area.tsv'
    aparc_destrieux_lh_meancurv = subject_list + '_' + session_list + '_lh_aparc_destrieux_meancurv.tsv'
    aparc_destrieux_rh_meancurv = subject_list + '_' + session_list + '_rh_aparc_destrieux_meancurv.tsv'

    aparc_BA_lh_volume = subject_list + '_' + session_list + '_lh_aparc_BA_volume.tsv'
    aparc_BA_rh_volume = subject_list + '_' + session_list + '_rh_aparc_BA_volume.tsv'
    aparc_BA_lh_thickness = subject_list + '_' + session_list + '_lh_aparc_BA_thickness.tsv'
    aparc_BA_rh_thickness = subject_list + '_' + session_list + '_rh_aparc_BA_thickness.tsv'
    aparc_BA_lh_area = subject_list + '_' + session_list + '_lh_aparc_BA_area.tsv'
    aparc_BA_rh_area = subject_list + '_' + session_list + '_rh_aparc_BA_area.tsv'
    aparc_BA_lh_meancurv = subject_list + '_' + session_list + '_lh_aparc_BA_meancurv.tsv'
    aparc_BA_rh_meancurv = subject_list + '_' + session_list + '_rh_aparc_BA_meancurv.tsv'

    subject_name = subject_list + '_' + session_list
    output_path = os.path.expanduser(output_dir)
    cs_dir = output_path + '/' + 'analysis-series-' + analysis_series_id + '/subjects' + '/' + subject_list + '/' + session_list + '/t1/freesurfer-cross-sectional'
    if not os.path.isdir(cs_dir):
        print("ERROR: directory freesurfer-cross-sectional does not exist, it should be CAPS directory after running recon_all_pipeline!!!")
    else:
        pass
    dest_dir = cs_dir + '/statistics'
    try:
        os.makedirs(dest_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST: # if dest_dir exists, go on, if its other error, raise
            raise

    subject = os.path.join(cs_dir, subject_name)
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

    freesurferstatas_to_tsv(subject,
                            all_seg_volume_tsv,
                            aseg_volume_tsv,
                            aparc_desikan_lh_volume_tsv,
                            aparc_desikan_rh_volume_tsv,
                            aparc_desikan_lh_thickness_tsv,
                            aparc_desikan_rh_thickness_tsv,
                            aparc_desikan_lh_area_tsv,
                            aparc_desikan_rh_area_tsv,
                            aparc_desikan_lh_meancurv_tsv,
                            aparc_desikan_rh_meancurv_tsv,
                            aparc_destrieux_lh_volume_tsv,
                            aparc_destrieux_rh_volume_tsv,
                            aparc_destrieux_lh_thickness_tsv,
                            aparc_destrieux_rh_thickness_tsv,
                            aparc_destrieux_lh_area_tsv,
                            aparc_destrieux_rh_area_tsv,
                            aparc_destrieux_lh_meancurv_tsv,
                            aparc_destrieux_rh_meancurv_tsv,
                            aparc_BA_lh_volume_tsv,
                            aparc_BA_rh_volume_tsv,
                            aparc_BA_lh_thickness_tsv,
                            aparc_BA_rh_thickness_tsv,
                            aparc_BA_lh_area_tsv,
                            aparc_BA_rh_area_tsv,
                            aparc_BA_lh_meancurv_tsv,
                            aparc_BA_rh_meancurv_tsv)

