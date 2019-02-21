# coding: utf8

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Junhao Wen"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Junhao Wen"
__email__ = "Junhao.Wen@inria.fr"
__status__ = "Development"


def bids_datagrabber(input_dir, subject_list, session_list):
    """
    Fetch t1 images from a BIDS directory based on subject_list and a
    session_list

    Args:
        input_dir: BIDS directory
        subject_list: a list containing all the participant_id
        session_list: a list containing all the session_id

    Returns: a list containing all the t1 images
    """
    from bids.grabbids.bids_layout import BIDSLayout
    from clinica.utils.stream import cprint

    bidslayout = BIDSLayout(input_dir)
    anat_t1 = []
    missing_subject_session = []
    if not bidslayout.get(target='run', return_type='id', type='T1w'):
        for i in range(len(subject_list)):
            t1 = bidslayout.get(return_type='file',
                                type='T1w',
                                extensions=['nii|nii.gz'],
                                session=session_list[i].replace('ses-', ''),
                                subject=subject_list[i].replace('sub-', ''))
            if len(t1) == 0:
                missing_subject_session.append([subject_list[i], session_list[i]])
            else:
                anat_t1.append(t1)
    else:
        cprint("There are more than one run for T1w image for this analysis")
        for i in range(len(subject_list)):
            t1 = bidslayout.get(return_type='file',
                                type='T1w',
                                extensions=['nii|nii.gz'],
                                session=session_list[i].replace('ses-', ''),
                                subject=subject_list[i].replace('sub-', ''),
                                run='1')
            if len(t1) == 0:
                missing_subject_session.append([subject_list[i], session_list[i]])
            else:
                anat_t1.append(t1)

    # Check if pybids works well to find all the T1 images
    if len(missing_subject_session) > 0:
        error_string = 'Please verify there is no error in your TSV file. Clinica could not find T1w for those ' + str(len(missing_subject_session)) + ' subjects - session :'
        for e in missing_subject_session:
            error_string += '\n' + e[0] + ' with session ' + e[1]
        raise IOError(error_string)
    if len(anat_t1) != len(subject_list) or len(anat_t1) != len(session_list):
        raise ValueError('Found ' + str(len(anat_t1)) + '  T1w but there are ' + str(len(subject_list)) + ' subjects.')

    return anat_t1


def get_dirs_check_reconalled(output_dir, subject_list, session_list):
    """
    Get the info from subjects_visits_tsv, like subject_dir, subject_id,
    subject_list, session_list. Also, this function checks out the rerun of
    the dataset, if the subject result folder has been created, you should
    check out the result and decide if you are going to rerun it or just
    ignore it.

    Args:
        output_dir: CAPS directory to contain the output
        subject_list:  a list containing all the participant_id
        session_list: a list containing all the session_id

    Returns: the related lists based on the tsv files
    """
    import os
    import errno
    from copy import deepcopy as cp
    import subprocess
    from clinica.utils.stream import cprint

    # subject_id, subject_list and session_list
    subject_id = list(subject_list[i] + '_' + session_list[i] for i in range(len(subject_list)))
    subject_id_without_reconalled = cp(subject_id)
    subject_list_without_reconalled = cp(subject_list)
    session_list_without_reconalled = cp(session_list)

    # output_path is the path to CAPS
    output_path = os.path.expanduser(output_dir)  # change the relative path to be absolute path
    output_dir = os.path.join(output_path, 'subjects')

    try:
        os.makedirs(output_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:  # if the error is not exist error, raise, otherwise, pass
            raise

    # subject_dir is the real path to FreeSurfer output path
    subject_dir = []
    subject_dir_without_reconalled = []

    for i in range(len(subject_list)):
        subject = os.path.join(output_dir,
                               subject_list[i],
                               session_list[i],
                               't1',
                               'freesurfer_cross_sectional')
        try:
            os.makedirs(subject)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        # add the FS path into a list with all subjects and check the recon-all.log file to see if this subject have been run recon-all successfully.
        subject_dir.append(subject)
        subject_path = os.path.join(subject, subject_id[i])
        subject_path_abs = os.path.expanduser(subject_path)
        # check the recon-all.log
        log_file = os.path.join(subject_path_abs, 'scripts', 'recon-all.log')
        if os.path.isfile(log_file):
            last_line = subprocess.check_output(['tail', '-1', log_file])
            if b'finished without error' in last_line:
                cprint("Skipping %s (FreeSurfer segmentation without error was found)" % subject_id[i])
                subject_id_without_reconalled.remove(subject_id[i])
                subject_list_without_reconalled.remove(subject_list[i])
                session_list_without_reconalled.remove(session_list[i])
            else:
                subject_dir_without_reconalled.append(subject)
        else:
            subject_dir_without_reconalled.append(subject)

    return subject_dir, subject_id, subject_dir_without_reconalled, subject_id_without_reconalled, subject_list_without_reconalled, session_list_without_reconalled


def check_fov(t1_list, recon_all_args):
    """
    Check size of inputs and field of view (FOV) of each T1 image.

    Args:
        t1_list: a list containing all the t1 images
        recon_all_args: default the -qcache flag

    Returns:
        A list containing "<recon_all_args> -cw256" or "<recon_all_args>"
        based on the FOV.
    """
    import sys
    import nibabel as nib
    from nipype.utils.filemanip import filename_to_list
    from clinica.utils.stream import cprint

    output_flags = []
    num_t1 = len(t1_list)
    t1_list = filename_to_list(t1_list)

    if num_t1 == 0:
        cprint("ERROR: No T1's Given")
        sys.exit(-1)

    f = nib.load(t1_list[0])
    voxel_size = f.header.get_zooms()
    t1_size = f.header.get_data_shape()
    if (voxel_size[0] * t1_size[0] > 256) or \
            (voxel_size[1] * t1_size[1] > 256) or \
            (voxel_size[2] * t1_size[2] > 256):
        # cprint("Setting MRI Convert to crop images to 256 FOV")
        optional_flag = '-cw256'
    else:
        # cprint("No need to add -cw256 flag")
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


def write_tsv_files(subject_id, output_dir):
    """
    Generate statistics TSV files for a given subject.
    """
    import os
    from clinica.utils.freesurfer import generate_regional_measures

    participant_id = subject_id.split('_')[0]
    session_id = subject_id.split('_')[1]

    path_segmentation = os.path.join(
        os.path.expanduser(output_dir),
        'subjects',
        participant_id,
        session_id,
        't1',
        'freesurfer_cross_sectional'
    )

    generate_regional_measures(path_segmentation, subject_id)


def print_begin_pipeline(subject_id):
    """
    Display begin message for a given subject.

    Args:
        subject_id: Subject ID (e.g. sub-CLNC01_ses-M0)
    """
    from clinica.utils.stream import cprint
    import datetime
    from colorama import Fore

    now = datetime.datetime.now().strftime('%H:%M:%S')

    cprint('%s[%s]%s Running pipeline for %s...' % (
        Fore.BLUE, now, Fore.RESET, subject_id))


def print_end_pipeline(subject_id, final_file):
    """
    Display end message for a given subject.

    Args:
        subject_id: Subject ID (e.g. sub-CLNC01_ses-M0)
        final_file: One of the last file generated by the pipeline.
            Connect a final file to this argument in order to display the end
            message.
    """
    from clinica.utils.stream import cprint
    import datetime
    from colorama import Fore

    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s ...%s has completed.' % (
        Fore.GREEN, now, Fore.RESET, subject_id))
