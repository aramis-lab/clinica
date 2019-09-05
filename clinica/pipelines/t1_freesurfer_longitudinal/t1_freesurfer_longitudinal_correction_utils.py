# coding: utf8

__author__ = "Alexis Guyot"
__copyright__ = "Copyright 2016-2019, The Aramis Lab Team"
__credits__ = ["Alexis Guyot"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Alexis Guyot"
__email__ = "alexis.guyot@icm-institute.org"
__status__ = "Development"


def check_template_reconalled(
        caps_dir, unique_subject_list, persubject_session_list2):
    """Check if the template creation has been run for each subject.

    Args:
        caps_dir (string): CAPS subdirectory containing the subjects
            (absolute path)
        unique_subject_list (list of string): list of unique
            participant_id
        persubject_session_list2 (list of string): list of per-subject
            list of sessions

    Returns:
        N/A
    """
    import os
    import subprocess
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as template_utils

    unique_subject_number = len(unique_subject_list)

    # check recon_all template creation for each subject
    for subject_index in range(unique_subject_number):
        subject = unique_subject_list[subject_index]
        session_list = persubject_session_list2[subject_index]
        # get template location in the CAPS directory
        caps_template_target = template_utils.get_capstemplate_path(
            caps_dir, subject, session_list)
        # check if template subfolder exists
        if not os.path.isdir(caps_template_target):
            error_msg = 'Error: {0} does not exist.'.format(
                caps_template_target)
            raise IOError(error_msg)
        # Check that the recon-all template creation run was successful
        log_file = os.path.join(caps_template_target, 'scripts', 'recon-all.log')
        if os.path.isfile(log_file):
            last_line = subprocess.check_output(['tail', '-1', log_file])
            if b'finished without error' not in last_line:
                error_msg = 'No template found for subject {0}.'.format(
                    subject)
                error_msg = '{0} Re-run t1-freesurfer-template.'.format(
                    error_msg)
                raise IOError(error_msg)


def check_reconalled(caps_dir, subject_list, session_list):
    """Check cross-sectional and -base recon-all run previously

    Checks if t1-freesurfer-cross-sectional ('recon-all -all') and if
    t1-freesurfer-template ('recon-all -base') have been correctly run
    beforehand.
    Will crash if any of those have failed (i.e., data are missing in
    the CAPS folder)

    Args:
        caps_dir (string): CAPS directory to contain the output
        subject_list (list): a list containing all the participant_id
        session_list (list): a list containing all the session_id

    Returns:
        N/A
    """
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as template_utils
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_utils as utils

    # Check if t1-freesurfer-cross-sectional has been correctly run
    # beforehand. Will crash if it has not
    template_utils.check_xsectional_reconalled(
        caps_dir, subject_list, session_list)

    # get unique list of subjects and corresponding per-subject sessions
    [unique_subject_list,
     persubject_session_list2] = template_utils.get_unique_subjects(
         subject_list, session_list)

    # Check if t1-freesurfer-template has been correctly run
    # beforehand. Will crash if it has not
    utils.check_template_reconalled(
        caps_dir, unique_subject_list, persubject_session_list2)


def check_processed(
        subject_list,
        session_list,
        caps_dir,
        template_unpcssd_sublist):
    """Check which subjects have (/not) been processed already

    Check the .tsv file (to force overwrite already-processed subjects)
    generated while running template pipeline to see which subjects have
    been processed already and which have not yet.
    Additionally, this checks two extra things:
    1. that the user has not provided a new .tsv file after executing
    the template pipeline for running the longitudinal-correction
    2. that the CAPS folder does not find subjects for which a
    longitudinal-correction was condcuted with no corresponding template
    in CAPS dir (logical inconsistency)

    Args:
        subject_list (list of string): list of participant_id
        session_list (list of string): list of session_id corresponding
            to a specific participant_id
        caps_dir (string): location of CAPS folder
        template_unpcssd_sublist (list of string): list of all the
            subjects for which no processing was found in the CAPS
            folder

    Returns:
        unpcssd_subject_list (list of strings): list of all the subjects
            for which no template/longitudinal-correction processing was
            found in the CAPS folder
        unpcssd_session_list (list of list of strings): list of all the
            sessions for which no template/longitudinal-correction
            processing was found in the CAPS folder
        unpcssd_capstarget_list (string): list, for all subjects for
            which no template/longitudinal-correction processing was
            found in the CAPS folder, of associated CAPS folder where
            the template will get stored
        all_sublist (list of strings): list of all the subjects
        all_seslist (list of strings): list of all sessions
        all_capstargetlist (string): list, for all subjects and
            corresponding sessions of associated CAPS folder where the
            template will get stored
        rundir (str): path to the directory where clinica is being
            launched. Used to create a .tsv file in running directory
            if the user wishes to re-run the pipeline on subjects that
            had already been processed.
    """
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as template_utils

    # Check processed/unprocessed subjects for longitudinal-correction
    [
        longcor_pcssd_sublist,
        dummy,
        dummy,
        longcor_unpcssd_sublist,
        longcor_unpcssd_seslist,
        longcor_unpcssd_capstargetlist,
        all_sublist,
        all_seslist,
        all_capstargetlist] = template_utils.check_caps_processing(
            'longitudinal-correction',
            subject_list,
            session_list,
            caps_dir)

    # Raise error if some subjects longitudinal-correction processed but
    # not template processed
    for longcor_pccsd_sub in longcor_pcssd_sublist:
        if longcor_pccsd_sub in template_unpcssd_sublist:
            error_msg = 'Error: longitudinal processing found for'
            error_msg = '{0} subject {1} '.format(error_msg, longcor_pccsd_sub)
            error_msg = '{0} but no corresponding template found.'.format(error_msg)
            raise ValueError(error_msg)

    # deduce list of unprocessed subjects/sessions (i.e., all
    # longitudinal-correction unprocessed subjects/sessions)
    unpcssd_sublist = longcor_unpcssd_sublist
    unpcssd_seslist = longcor_unpcssd_seslist
    unpcssd_capstargetlist = longcor_unpcssd_capstargetlist

    return [
        unpcssd_sublist,
        unpcssd_seslist,
        unpcssd_capstargetlist,
        all_sublist,
        all_seslist,
        all_capstargetlist]


def get_warning_msg(
        pcssd_capstargetlist,
        caps_dir,
        tsv_path,
        working_directory,
        n_procs):
    """Generate a warning message if processing in CAPS folder already

    Will generate a warning message to indicate to the user that there
    already is processing in the CAPS folder. Will also show the command
    to run if they wish to overwrite the existing CAPS data

    Args:
        processing_type (string): Either 'template' or 'longitudinal
            correction'. Used to build a warning message specific to the
            type of data that is being looked for in the CAPS folder.
        pcssd_capstargetlist (list of strings): list of CAPS subfolders
            which contain template / longitudinal correction data that
            have been processed already
        caps_dir (string): location of CAPS folder
        tsv_path (string): path to .tsv file provided by user. None if
            no subjects already processed.
            Used to display the command line to force re-running the
            pipeline if a previous run was detected in CAPS folder.
        working_directory (string): path to working directory
            provided by user. None if they did not provide any.
            Used to display the command line to force re-running the
            pipeline if a previous run was detected in CAPS folder.
        n_procs (int): number of processors provided by user
            Used to display the command line to force re-running the
            pipeline if a previous run was detected in CAPS folder.

    Returns:
        warning_msg (string): warning message
    """
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as template_utils

    # define command line to force overwriting data in CAPS folder
    cl_overwrite_caps = template_utils.get_cl_overwritecaps(
        caps_dir, tsv_path, working_directory, n_procs)

    # sanity check: warning message only if at least one subject
    # processed already
    if pcssd_capstargetlist:
        # define warning message
        if len(pcssd_capstargetlist) == 1:
            warning_msg = 'There is already a template directory in:'
        if len(pcssd_capstargetlist) > 1:
            warning_msg = 'There are already template directories in:'
        warning_msg = '{0}\n'.format(warning_msg)
        for pcssd_capstarget in pcssd_capstargetlist:
            warning_msg = '{0}- {1}\n'.format(warning_msg, pcssd_capstarget)
        if len(pcssd_capstargetlist) == 1:
            warning_msg = '{0} To overwrite this directory, run:\n'.format(
                warning_msg)
        if len(pcssd_capstargetlist) > 1:
            warning_msg = '{0} To overwrite these directories, run:\n'.format(
                warning_msg)
        # add command line to force overwriting data in CAPS folder to
        # warning message
        warning_msg = '{0} \'{1}\''.format(warning_msg, cl_overwrite_caps)
    else:
        warning_msg = None

    return warning_msg


def to_process(
        pcssd_capstargetlist,
        unpcssd_sublist,
        unpcssd_seslist,
        unpcssd_capstargetlist,
        all_sublist,
        all_seslist,
        all_capstargetlist,
        in_caps_dir,
        in_overwrite_caps,
        in_working_directory,
        in_n_procs,
        overwrite_tsv_path):
    """ Get list of subjects to be processed

    The list of subjects to be processed will depend upon whether the
    user provided the option to --force-overwrite already existing
    processing or not. If the option --force-overwrite was not passed
    but existing processing was found, a warning message will be
    generated to show the command to run in order to force-overwrite the
    existing processing.

    Args:
        pcssd_capstargetlist (string): list, for all subjects for
            which template/longitudinal-correction processing was
            found in the CAPS folder, of associated CAPS folder where
            the template will get stored
        unpcssd_sublist (list of strings): list of all the subjects
            for which no processing was found in the CAPS folder
        unpcssd_seslist (list of objects): list, for all subjects
            for which no processing was found in the CAPS folder, of the
            list of associated session
        unpcssd_capstargetlist (string): list, for all subjects for
            which no processing was found in the CAPS folder, of
            associated CAPS folder where the template or
            longitudinal-correction processing will get stored
        all_sublist (list of strings): list of all the subjects
        all_seslist (list of objects): list, for all subjects of
            the list of associated sessions
        all_capstargetlist (string): list, for all subjects and
            corresponding sessions of associated CAPS folder where the
            template will get stored
        in_caps_dir (string): CAPS directory to contain the output
        in_overwrite_caps (string): Option provided by the user to state
            whether the data in the CAPS folder (i.e.,
            previously-computed template) should be overwritten or not
            by a new run. By default, do not overwrite. Will overwrite
            for the following values: 'true', 'True' and 'TRUE'
        in_working_directory (string): path to working directory
            provided by user. None if they did not provide any.
            Used to display the command line to force re-running the
            pipeline if a previous run was detected in CAPS folder.
        in_n_procs (int): number of processors provided by user
            Used to display the command line to force re-running the
            pipeline if a previous run was detected in CAPS folder.
        overwrite_tsv_path (str): path to the .tsv file to use if the user
            wishes to re-run the pipeline on subjects that had already
            been processed. None if no subjects already processed.

    Returns:
        topcss_sublist (list of strings): list of the subjects that will
            be processed. All if --force-overwrite, else only subjects
            that have not been processed already
        topcss_seslist2 (list of list): list of the list of sessions
            that correspond to each subject to be process
        topcss_capstargetlist (list of strings): list, for all subjects
            to be processed, of the associated CAPS folder where the
            template will get stored
        overwrite_warning (string): warning message to display in case
            template) processing has been found in the CAPS folder for
            any of the input subject/session
    """
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as template_utils
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_utils as utils

    # Convert overwrite-caps flag to boolean
    force_overwrite = template_utils.boolean_overwrite_caps(in_overwrite_caps)

    # Generate warning message if some subjects have already been
    # processed and --force-overwrite option not provided
    overwrite_warning = None
    some_processed = len(unpcssd_sublist) != len(all_sublist)
    if (some_processed) and (not force_overwrite):
        overwrite_warning = utils.get_warning_msg(
            pcssd_capstargetlist,
            in_caps_dir,
            overwrite_tsv_path,
            in_working_directory,
            in_n_procs)

    # check what subjects / sessions / CAPS target dirs to process
    # depending on --force-overwrite flag
    if force_overwrite:
        topcss_sublist = all_sublist
        topcss_seslist2 = all_seslist
        topcss_capstargetlist = all_capstargetlist
    else:
        topcss_sublist = unpcssd_sublist
        topcss_seslist2 = unpcssd_seslist
        topcss_capstargetlist = unpcssd_capstargetlist

    return [
        topcss_sublist,
        topcss_seslist2,
        topcss_capstargetlist,
        overwrite_warning]


def process_input_node(
        in_caps_dir,
        in_subject_list,
        in_session_list,
        in_unpcssd_sublist,
        in_pcssd_capstargetlist,
        in_overwrite_caps,
        in_working_directory,
        in_n_procs,
        in_overwrite_tsv):
    """ Carry out all processing for the input node

    Performs all the processing for the input node:
    1. Check if t1-freesurfer-cross-sectional has been correctly run
        beforehand
    2. Check if there are already template subfolders in the CAPS dir.

    Args:
        in_caps_dir (string): CAPS directory to contain the output
        in_subject_list (list): a list containing all the participant_id
        in_session_list (list): a list containing all the session_id
        in_unpcssd_sublist (list of string): list of participant_id, where
            each participant appears only once, wich have not already
            been (template) processed
        in_pcssd_capstargetlist (list of string): list of path to the
            template directory in CAPS folder for each participant, for
            the participants that have already been processed
        in_overwrite_caps (string): Option provided by the user to state
            whether the data in the CAPS folder (i.e.,
            previously-computed template) should be overwritten or not
            by a new run. By default, do not overwrite. Will overwrite
            for the following values: 'true', 'True' and 'TRUE'
        in_working_directory (string): path to working directory
            provided by user. None if they did not provide any.
            Used to display the command line to force re-running the
            pipeline if a previous run was detected in CAPS folder.
        in_n_procs (int): number of processors provided by user
            Used to display the command line to force re-running the
            pipeline if a previous run was detected in CAPS folder.
        in_overwrite_tsv (str): path to the .tsv file to use if the user
            wishes to re-run the pipeline on subjects that had already
            been processed. None if no subjects already processed.

    Returns:
        out_subject_list (list of string): list of participant_id
        out_session_list (list of list): list of session_id associated
            to any participant
        out_capstarget_list (list of string): list of path to the
            template directory in CAPS folder for each participant
        out_overwrite_warning (string): warning message to display in
            case (template or longitudinal-correction) processing has
            been found in the CAPS folder for any of the input
            subject/session
        out_overwrite_sublist (list of strings): list of subjects for
            which either template or longitudinal processing should be
            re-done if the user re-runs the pipeline using the flag
            '--force_overwrite True'
        out_overwrite_seslist (list of strings): list of sessions for
            which either template or longitudinal processing should be
            re-done if the user re-runs the pipeline using the flag
            '--force_overwrite True'
    """
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_utils as utils

    # check if we have the necessary recon-all data (cross-sectional
    # runs and created templates) prior to running longitudinal
    # processing. Will crash if we do not have this.
    utils.check_reconalled(in_caps_dir, in_subject_list, in_session_list)

    # Check target CAPS target dir for template /
    # longitudinal-correction processing before conducting any
    # longitudinal processing
    [
        unpcssd_sublist,
        unpcssd_seslist,
        unpcssd_capstargetlist,
        all_sublist,
        all_seslist,
        all_capstargetlist] = utils.check_processed(
            in_subject_list,
            in_session_list,
            in_caps_dir,
            in_unpcssd_sublist)

    # Deduce subjects to be processed next (depending on whether
    # '--force-overwrite' flag was provided or not)
    [
        out_subject_list,
        out_session_list,
        out_capstarget_list,
        out_overwrite_warning] = utils.to_process(
            in_pcssd_capstargetlist,
            unpcssd_sublist,
            unpcssd_seslist,
            unpcssd_capstargetlist,
            all_sublist,
            all_seslist,
            all_capstargetlist,
            in_caps_dir,
            in_overwrite_caps,
            in_working_directory,
            in_n_procs,
            in_overwrite_tsv)

    return [
        out_subject_list,
        out_session_list,
        out_capstarget_list,
        out_overwrite_warning]


def create_symlinks(
        in_caps_dir,
        in_subject_list,
        in_session_list,
        in_longsubdirnames_dict):
    """Create symlinks for recon-all execution

    t1-freesurfer-cross-sectional has created the following folders:
    [caps_dir]/subjects/sub-[sub_i]/ses-[ses_ij]/t1/...
    ...freesurfer_cross_sectional/sub-[sub_i]_ses-[ses_ij]
    and t1-freesurfer-template has created:
    [caps_dir]/subjects/sub-[sub_i]/long-[long_i]/...
    ...freesurfer_unbiased_template/sub-[sub_i]
    where i in range [1, subject_number], [long_i] corresponds to all
    seessions associated with [sub_i] and ij represents any of the
    potential sessions for subject i
    This function creates for all above folders the following symbolic
    link ./subjects/sub-[sub_i]_ses-[ses_ij]
    and
    ./subjects/sub-[sub_i] (template) in the working directory.
    All the symlinks are located in the same folder, disrespective of
    the subject (i.e. use a Node and not a MapNode). The path of that
    folder is returned so recon-all -base can later be run on all the
    sessions corresponding to subject i.

    Args:
        in_caps_dir (string): CAPS directory to contain the output
        in_subject_list (list of string): list of participant_id
        in_session_list (list of string): list of corresponding
            session_id
        in_longsubdirnames_dict (dictionary of strings): for each
            unique subject, name of longitudinal subofolder associated
            to the current set of visits

    Returns:
        out_subject_symlink_path (string): path to the folder
            containing all session symlinks and all template symlinks
            for all subject_id
    """
    import os
    import errno

    # change the relative path to be absolute path
    caps_path = os.path.expanduser(in_caps_dir)
    caps_dir = os.path.join(caps_path, 'subjects')

    # create subject-specific path in clinica's working directory
    # (later used as $SUBJECT_DIR)
    out_subject_symlink_path = os.path.abspath('./subjects')
    try:
        os.mkdir(out_subject_symlink_path)
    except OSError as oserror:
        if oserror.errno != errno.EEXIST:
            raise

    subject_number = len(in_subject_list)
    for subject_index in range(subject_number):
        # retrieve current {subject,session}
        subject = in_subject_list[subject_index]
        session = in_session_list[subject_index]

        # find all the sessions associated to the current subject
        long_subdir_name = in_longsubdirnames_dict[subject]

        # generate symlink for cross-sectional reconall run
        caps_subjectid_path = "{0}/{1}/{2}/t1/freesurfer_cross_sectional/{1}_{2}".format(
            caps_dir, subject, session)
        wd_subjectid_symlink_path = "{0}/{1}_{2}".format(
            out_subject_symlink_path, subject, session)
        if not os.path.exists(wd_subjectid_symlink_path):
            os.symlink(caps_subjectid_path, wd_subjectid_symlink_path)

        # generate symlink for the template associated with subject_id
        caps_template_path = "{0}/{1}/{2}/freesurfer_unbiased_template/{1}".format(
            caps_dir, subject, long_subdir_name)
        wd_template_symlink_path = "{0}/{1}".format(
            out_subject_symlink_path, subject)
        if not os.path.exists(wd_template_symlink_path):
            os.symlink(caps_template_path, wd_template_symlink_path)

    return out_subject_symlink_path


def get_reconalllong_flags(in_subject, in_session):
    """Get recon-all -long flags for each subject

    Creates the flags that are required to run the recon-all -long for
    each subject (i.e. run freesurfer longitudinal for each
    {subject;session}).

    Args:
        in_subject (string): participant_id (/!\ from self.subjects /!\)
        in_session (string): session_id (/!\ from self.sessions /!\)

    Returns:
        out_reconalllong_flags (string): all the flags to run
            recon-all -long for any particular subject and their
            (unique) associated session. There are as many
            'recon-all flags' as there are {subject, session}
    """
    template_id = in_subject
    out_reconalllong_flags = '-long {0}_{1} {2} -all'.format(
        in_subject, in_session, template_id)

    return out_reconalllong_flags


def check_reconall_longitudinal_single(subjects_dir,
                                       in_subject,
                                       in_session):
    """Check if longitudinal run successfully for any single subject

    Args:
        subjects_dir (string): CAPS subdirectory containing the subject
            (absolute path)
        in_subject (string): participant_id
        in_session (string): session_id

    Returns:
        out_longitudinal_result (Boolean): True if success, False
            otherwise
    """
    import os
    import subprocess

    # check if longitudinal subfolder exists
    subjectid_longitudinal = "{0}/{1}_{2}.long.{1}".format(
        subjects_dir, in_subject, in_session)
    if not os.path.isdir(subjectid_longitudinal):
        raise OSError(
            '{0} is not a valid folder.'.format(subjectid_longitudinal))

    # check the recon-all.log
    log_file = os.path.join(subjectid_longitudinal, 'scripts', 'recon-all.log')
    if os.path.isfile(log_file):
        last_line = subprocess.check_output(['tail', '-1', log_file])
        if b'finished without error' not in last_line:
            out_longitudinal_result = False
            ioerror_msg = 'Subject {0} has not been processed yet'.format(
                in_subject)
            ioerror_msg = '{0} Re-run t1-freesurfer-longitudinal.'.format(
                ioerror_msg)
            raise IOError(ioerror_msg)
        else:
            out_longitudinal_result = True
    else:
        raise IOError(
            'no recon-all.log script for subject {0}'.format(in_subject))

    return out_longitudinal_result


def run_reconalllong(in_subject,
                     in_session,
                     in_reconalllong_flags,
                     in_symlink_path):
    """Run recon-all -long for each subject

    Run freesurfer longitudinal for each {subject;session}). The
    recon-all command is built with previously computed arguments
    (flags).

    Args:
        in_subject (string): participant_id
        in_session (string): session_id
        in_reconalllong_flags (string): all the flags to run
            recon-all -long for any particular subject and their
            (unique) associated session
        in_symlink_path (string): paths to subject-specific folder
            containing all session symlinks for subject_id associated to
            in_reconalllong_flags.

    Returns:
        out_longitudinal_result (Boolean): True if success, False
            otherwise
    """
    import subprocess
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_utils as utils

    # run recon-all -long
    reconalllong_command = 'recon-all {0} -sd {1}'.format(
        in_reconalllong_flags, in_symlink_path)
    subprocess_run_reconalllong = subprocess.run(
        reconalllong_command,
        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if subprocess_run_reconalllong.returncode != 0:
        raise ValueError('recon-all -long failed, returned non-zero code')

    # check whether recon-all -long failed or not
    out_longitudinal_result = utils.check_reconall_longitudinal_single(
        in_symlink_path, in_subject, in_session)

    return out_longitudinal_result


def write_stats_files(
        in_subject, in_session, in_symlink_path, in_longitudinal_result):
    """Generate statistics TSV files for a given {subject,session}.

    Statistics are derived from recon-all's own .stats files and are
    meant to provide a summary of those.

    Args:
        in_subject (string): participant_id
        in_session (string): session_id
        in_symlink_path (string): paths to subject-specific folder
            containing all session symlinks for subject_id associated to
            in_reconalllong_flags. Used for 1) retrieving the path to
            recon-all generated .tsv files 2) write stats files in a
            subfolder of the folder containing the symlink to
            {subject,session}
        in_longitudinal_result (Boolean): Used to force longitudinal
            recon-all to be run before generating the .tsv files from
            the longitudinal processing

    Returns:
        out_stats_path (string): path to the folder containing
            all the .tsv stats files for {subject,session}
    """
    import os
    from clinica.utils.freesurfer import generate_regional_measures

    # get path to longitudinal folder in working dir (where recon-all
    # results are stored, including .stats files)
    seg_path = os.path.expanduser(in_symlink_path)

    # define output dir for the .tsv files
    subses_long_id = '{0}_{1}.long.{0}'.format(in_subject, in_session)
    subses_long_stats_subdir = "{0}_stats".format(subses_long_id)
    out_stats_path = os.path.join(seg_path, subses_long_stats_subdir)

    # generate the .tsv files from .stats files
    subses_id = '{0}_{1}'.format(in_subject, in_session)
    generate_regional_measures(
        seg_path, subses_id, output_dir=out_stats_path, longitudinal=True)

    return out_stats_path


def copy_to_caps(
        in_subject,
        in_session,
        in_subject_dir,
        in_caps_target,
        in_caps_dir,
        in_overwrite_warning,
        in_overwrite_caps,
        in_stats_path):
    """ Copy template folder to CAPS dir

    Move the longitudinal-correction folder created by recon-all
    longitudinal to its correct location in the CAPS tree structure.

    Args:
        in_subject (string): participant_id
        in_session (string): session_id
        in_subject_dir (string): location of the subjects (in the
            working directory)
        in_caps_target (string): path to CAPS target subfolder for the
            longitudinal correction connected to the current
            {subject,session}
        in_caps_dir (string): CAPS directory to contain the output
            Used to make sure the folder we are overwriting (if
            overwrite-caps option provided) are subolders of the input
            CAPS dir
        in_overwrite_warning (string): warning message if some subjects
            have already been processed. Shows the subjects/sessions
            that have been processed and the command line to run in case
            the user wishes to re-run clinica on the subjects that have
            been processed already
        in_overwrite_caps (string): Option provided by the user to state
            whether the data in the CAPS folder (i.e.,
            previously-computed template) should be overwritten or not
            by a new run. By default, do not overwrite. Will overwrite
            for the following values: 'true', 'True' and 'TRUE'
        in_stats_path (string): path to the folder (in working
            directory) that contains

    Returns:
        out_copy2_caps (boolean): flag to indicate whether or not the
            copy was successful
    """
    import os
    import shutil
    import warnings
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as template_utils

    # Convert overwrite-caps flag to boolean
    force_overwrite = template_utils.boolean_overwrite_caps(in_overwrite_caps)

    # define longitudinal-correction location in the working directory
    wd_longitudinal_location = os.path.join(
        in_subject_dir, '{0}_{1}.long.{0}'.format(in_subject, in_session))

    # define location of .tsv stats files in the CAPS directory,
    # following a similar pattern to that of in_caps_target:
    # the CAPS target is as follows:
    # [caps_dir]/subjects/[subject]/[session]/t1/[longsubdir]/...
    # ...freesurfer_longitudinal/...
    # ...sub-[subject]_ses-[session].long.sub-[subject]
    # The CAPS stats dir differs as it is:
    # [caps_dir]/subjects/[subject]/[session]/t1/[longsubdir]/...
    # ...freesurfer_longitudinal/regional_measures
    caps_target_stats = os.path.join(
        os.path.dirname(in_caps_target),
        'regional_measures')

    # Check if overwrite-caps option provided. If provided, delete the
    # target subdirectory in CAPS folder
    if force_overwrite:
        # check the target subdirectory exists in CAPS
        if os.path.exists(in_caps_target):
            # remove
            template_utils.safer_rmtree(in_caps_target, in_caps_dir)
        # check the stats subdirectory exists in CAPS
        if os.path.exists(caps_target_stats):
            # remove
            template_utils.safer_rmtree(caps_target_stats, in_caps_dir)

    # copy the longitudinal processing from working directory to CAPS
    # directory
    shutil.copytree(wd_longitudinal_location, in_caps_target)

    # copy the .tsv stats files from working directory to CAPS directory
    shutil.copytree(in_stats_path, caps_target_stats)
    out_copy2_caps = True

    # show warning message
    if in_overwrite_warning is not None:
        warnings.warn(in_overwrite_warning)

    return out_copy2_caps
