# coding: utf8

__author__ = "Alexis Guyot"
__copyright__ = "Copyright 2016-2019, The Aramis Lab Team"
__credits__ = ["Alexis Guyot"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Alexis Guyot"
__email__ = "alexis.guyot@icm-institute.org"
__status__ = "Development"


def get_capsxsectional_path(caps_dir, subject, session):
    """Get path to cross-sectional run in CAPS dir

    Args:
        caps_dir (string): location of CAPS folder
        subject (string): participant ID
        session_list (strings): session

    Returns:
        caps_xsectional_path (string): path to the cross-sectional
            subfolder in the CAPS directory
    """
    import os

    caps_xsectional_path = os.path.join(
        caps_dir,
        'subjects',
        subject,
        session,
        't1',
        'freesurfer_cross_sectional'
        )

    return caps_xsectional_path


def check_xsectional_reconalled(caps_dir, subject_list, session_list):
    """Check cross-sectional recon-all was run previously

    Checks if t1-freesurfer-cross-sectional has been
    correctly run beforehand.

    Args:
        caps_dir (string): CAPS directory to contain the output
        subject_list (list): a list containing all the participant_id
        session_list (list): a list containing all the session_id

    Returns:
        N/A
    """
    import os
    import subprocess
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils

    # subses_list, subject_list and session_list
    subject_number = len(subject_list)
    subses_list = list(
        '{0}_{1}'.format(
            subject_list[subject_index],
            session_list[subject_index]
            ) for subject_index in range(subject_number))

    # change the relative path to be absolute path
    abs_caps_path = os.path.abspath(caps_dir)
    subjects_path = os.path.join(abs_caps_path, 'subjects')

    if not os.path.isdir(subjects_path):
        error_msg = 'Error: {0} does not exist.'.format(subjects_path)
        raise IOError(error_msg)

    for subject_index in range(subject_number):
        subject = subject_list[subject_index]
        session = session_list[subject_index]
        subses = subses_list[subject_index]

        xsectional_path = utils.get_capsxsectional_path(
            caps_dir, subject, session)
        if not os.path.isdir(xsectional_path):
            error_msg = 'Error: {0} does not exist.'.format(xsectional_path)
            raise IOError(error_msg)

        # check the recon-all.log
        log_file = os.path.join(
            xsectional_path, subses, 'scripts', 'recon-all.log')
        if os.path.isfile(log_file):
            last_line = subprocess.check_output(['tail', '-1', log_file])
            if b'finished without error' not in last_line:
                error_msg = 'Subject {0} has not been processed yet.'.format(
                    subses)
                error_msg = '{0} Re-run t1-freesurfer-cross-sectional.'.format(
                    error_msg)
                raise IOError(error_msg)


def get_unique_subjects(in_subject_list, in_session_list):
    """Get unique participant IDs

    The function to read the .tsv file returns the following
    participant_id and session_id lists:
    participant1, participant1, ..., participant2, participant2, ...
    session1    , session2    , ..., session1    , session2    , ...
    This function returns a list where all participants are only selected
    once:
    participant1, participant2, ..., participant_n
    and for each participant, the list of corresponding session id
    eg.:
    participant1 -> [session1, session2]
    participant2 -> [session1]
    ...
    participant_n -> [session1, session2, session3]

    Args:
        in_subject_list (list of strings): list of participant_id
        in_session_list (list of strings): list of session_id

    Returns:
        out_unique_subject_list (list of strings): list of
            participant_id, where each participant appears only once
        out_persubject_session_list2 (list of list): list of list
            (list2) of session_id associated to any single participant
    """

    import numpy as np

    subject_array = np.array(in_subject_list)
    session_array = np.array(in_session_list)

    # The second returned element indicates for each participant_id the
    # element they correspond to in the 'unique' list. We will use this
    # to link each session_id in the repeated list of session_id to
    # their corresponding unique participant_id

    unique_subject_array, out_inverse_positions = np.unique(
        subject_array, return_inverse=True)
    out_unique_subject_list = unique_subject_array.tolist()

    subject_number = len(out_unique_subject_list)
    out_persubject_session_list2 = [
        session_array[
            out_inverse_positions == subject_index
            ].tolist() for subject_index in range(subject_number)]

    return out_unique_subject_list, out_persubject_session_list2


def sessionid_to_sessionlabel(sessionid):
    """Extract session label from session ID

    Example:
    if session_id='ses_M00', the session_label='M00'

    Args:
        sessionid (string): session ID

    Returns:
        sessionlabel (string): session ID
    """
    sessionlabel = sessionid[4:]

    return sessionlabel


def get_longsubdir_name(session_list):
    """Returns a subfolder name associated to a set of visits

    This will create a unique identifier for a subject and its
    corresponding sessions.
    Used to be able to run the longitudinal processing at different
    points in time, where new visits are successively added.
    eg.
    - 1. t=0 year, subject: sub01, sessions: ses-M00, ses-M01
        -> identifier: long-M00M01
    - 2. t=1 year, subject: sub01, sessions: ses-M00, ses-M01, ses-M02
        -> identifier: long-M00M01M02
    - 2. t=2 year, subject: sub01, sessions: ses-M00, ses-M01, ses-M02,
            ses-M03, ses-M04
        -> identifier: long-M00M01M02M03M04

    Args:
        session_list (list of strings):

    Returns:
        long_subdirname (string): name of the subfolder
            corresponding to all the sessions of a unique subject. This
            will be used as part of both:
            - the path to the unbiased template for a subject and its
                sessions
            - the corresponding path for the longitudinal correction
    """
    # get session label for each session ID
    sessionlabel_list = [
        sessionid_to_sessionlabel(session) for session in session_list]

    # concatenate all the sessions labels
    long_subdirname = ''.join(sessionlabel_list)
    long_subdirname = 'long-{0}'.format(long_subdirname)

    return long_subdirname


def get_capstemplate_path(caps_dir, subject, session_list):
    """Get path to template in CAPS dir

    Args:
        caps_dir (string): location of CAPS folder
        subject (string): participant ID
        session_list (list of strings): list of all the sessions
            corresponding to the participant

    Returns:
        caps_template_path (string): path to the template subfolder in the
            CAPS directory
    """
    import os
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils

    # get unique identifier corresponding to the subject and all
    # associated sessions. Will be used as part of a path in the CAPS
    # dir
    long_subdirname = utils.get_longsubdir_name(session_list)
    # define template location in the CAPS directory
    caps_template_path = os.path.join(
        caps_dir,
        'subjects',
        subject,
        long_subdirname,
        'freesurfer_unbiased_template',
        subject)

    return caps_template_path


def get_longsubdir_dict(
        in_subject_list, in_session_list):
    """Get subdir name for each subject

    For each unique subject in the provided subject list, return the
    name of the longitudinal subfolder that is associated to the
    current set of visits

    Args:
        in_subject_list (list of string): list of participant_id
        in_session_list (list of string): list of sessions

    Returns:
        longsubdir_dict (dictionary of strings): for each
            unique subject, name of longitudinal subofolder associated
            to the current set of visits
    """
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils

    # get list of unique subjects and for each individual subjects get
    # list of associated sessions
    [
        unique_subject_list,
        persubject_session_list2] = utils.get_unique_subjects(
            in_subject_list, in_session_list)

    # get longitudinal identifier for ech subject and store into
    # dictionary
    longsubdir_dict = dict()
    subject_number = len(unique_subject_list)
    for subject_index in range(subject_number):
        subject = unique_subject_list[subject_index]
        subject_session_list = persubject_session_list2[subject_index]
        longsubdir_dict[subject] = utils.get_longsubdir_name(
            subject_session_list)

    return longsubdir_dict


def get_capslongcor_path(caps_dir, subject, session, longsubdir_dict):
    """Get path to longitudinal correction subfolder in CAPS dir

    Args:
        caps_dir (string): location of CAPS folder
        subject (string): participant ID
        session (string): session ID
        longsubdir_dict (dictionary of strings): for each
            unique subject, name of longitudinal subofolder associated
            to the current set of visits

    Returns:
        caps_longcor_path (string): path to the longitudinal correction
            subfolder in the CAPS directory
    """
    import os

    # get longitudinal subfolder identifier for the current subject ID
    longsubdir = longsubdir_dict[subject]

    # get target CAPS subfolder
    capslongcor_path = os.path.join(
        caps_dir,
        'subjects',
        subject,
        session,
        't1',
        longsubdir,
        'freesurfer_longitudinal',
        '{0}_{1}.long.{0}'.format(subject, session))

    return capslongcor_path


def boolean_overwrite_caps(str_overwrite_caps):
    """Convert overwrite-caps flag to Boolean

    Args:
        str_overwrite_caps (string): any option in 'true', 'True',
            'TRUE', 'false', 'False', 'FALSE'

    Returns:
        overwrite_caps (Boolean): corresponding boolean value
    """
    if str_overwrite_caps in ['true', 'True', 'TRUE']:
        overwrite_caps = True
    elif str_overwrite_caps in ['false', 'False', 'FALSE']:
        overwrite_caps = False
    else:
        error_msg = 'Error: --force-overwrite should be \'true\' or \'false\''
        raise ValueError(error_msg)

    return overwrite_caps


def get_cl_overwritecaps(
        caps_dir,
        tsv_path,
        working_directory,
        n_procs):
    """Define command-line to force overwriting CAPS folder

    Output the command line a user will need to run if they want to
    force overwriting the CAPS folder, using the same parameters as
    were passed to the terminal for the current execution of the
    pipeline.

    Args:
        caps_dir (string): location of CAPS folder
        tsv_path (string): path to .tsv file provided by user
            Used to display the command line to force re-running the
            pipeline if a previous run was detected in CAPS folder.
            None if no subjects already processed.
        working_directory (string): path to working directory
            provided by user. None if they did not provide any.
            Used to display the command line to force re-running the
            pipeline if a previous run was detected in CAPS folder.
        n_procs (int): number of processors provided by user
            Used to display the command line to force re-running the
            pipeline if a previous run was detected in CAPS folder.

    Returns:
        cl_overwrite_caps (string): command line to force overwriting
            the CAPS folder. Parameters (e.g. CAPS dir, tsv file,
            working directory and n_procs) are the same as those
            provided by the user for the current execution of the
            pipeline.
    """
    # initialise
    cl_overwrite_caps = 'clinica run t1-freesurfer-longitudinal'
    # CAPS dir: mandatory argument
    cl_overwrite_caps = '{0} {1}'.format(cl_overwrite_caps, caps_dir)
    # TSV file: optional
    if tsv_path is not None:
        cl_overwrite_caps = '{0} -tsv {1}'.format(cl_overwrite_caps, tsv_path)
    # working directory: optional
    if working_directory is not None:
        cl_overwrite_caps = '{0} -wd {1}'.format(
            cl_overwrite_caps, working_directory)
    # number of processors: optional
    if n_procs is not None:
        cl_overwrite_caps = '{0} -np {1}'.format(cl_overwrite_caps, n_procs)
    # option to force overwrite: mandatory
    cl_overwrite_caps = '{0} --overwrite_caps True'.format(cl_overwrite_caps)

    return cl_overwrite_caps


def check_caps_processing(
        processing_type,
        subject_list,
        sesobj_list,
        caps_dir):
    """Check there is no processing done in CAPS already

    Check, for any of the 'template' or 'longitudinal-correction'
    processing type which subjects/sessions have been / have not been processed already.
    correspond folder for the current subject and its corresponding
    sessions.
    Three cases are possible:
    1. CAPS folder contains no processing
    2. All subjects/sessions have already been processed
    3. Some subjects/sessions (but not all) have been processed already.

    Args:
        processing_type (string): Either 'template' or 'longitudinal
            correction'. Used to build a warning message specific to the
            type of data that is being looked for in the CAPS folder.
        subject_list (list of string): list of participant_id. Will
            either be a list of unique subjects (template) or a list of
            all subjects for all sub-sess (longitudinal-correction)
        sesobj_list (list of objects): either a list of list of sessions
            for all unique participant_id (template) or a list of
            sessions (longitudinal_correction)
        caps_dir (string): location of CAPS folder

    Returns:
        pcssd_sublist (list of strings): list of all the subjects for
            which processing was found in the CAPS folder
        pcssd_sesobjlist (list of objects): list, for all subjects
            for which processing was found in the CAPS folder, of either
            the list of associated session (template) or of single
            sessions (longitudinal correction)
        unpcssd_capstargetlist (string): list, for all subjects for
            which processing was found in the CAPS folder, of
            associated CAPS folder where the template or
            longitudinal-correction processing will get stored
        unpcssd_sublist (list of strings): list of all the subjects
            for which no processing was found in the CAPS folder
        unpcssd_sesobjlist (list of objects): list, for all subjects
            for which no processing was found in the CAPS folder, of
            either the list of associated session (template) or of
            single sessions (longitudinal correction)
        unpcssd_capstargetlist (string): list, for all subjects for
            which no processing was found in the CAPS folder, of
            associated CAPS folder where the template or
            longitudinal-correction processing will get stored
        all_sublist (list of strings): list of all the subjects
        all_sesobjlist (list of objects): list, for all subjects of
            either the list of associated session (template) or of
            single sessions (longitudinal correction)
        all_capstargetlist (string): list, for all subjects and
            corresponding sessions of associated CAPS folder where the
            template will get stored
    """
    import os
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils

    # initialise lists of CAPS targets and list of processed /
    # unprocessed subjects
    pcssd_sublist = list()
    pcssd_sesobjlist = list()
    pcssd_capstargetlist = list()
    unpcssd_sublist = list()
    unpcssd_sesobjlist = list()
    unpcssd_capstargetlist = list()
    all_sublist = list()
    all_sesobjlist = list()
    all_capstargetlist = list()

    # if dealing with longitudinal-correction, get longitudinal
    # subfolder identifier for each unique subject prior to looping over
    # each sub/ses
    if processing_type == 'longitudinal-correction':
        longsubdir_dict = get_longsubdir_dict(subject_list, sesobj_list)

    # loop through all subjects and their associated sesobj
    subject_number = len(subject_list)
    for subject_index in range(subject_number):
        subject = subject_list[subject_index]
        sesobj = sesobj_list[subject_index]
        # get location in the CAPS directory
        if processing_type == 'template':
            caps_target = utils.get_capstemplate_path(
                caps_dir, subject, sesobj)
        elif processing_type == 'longitudinal-correction':
            caps_target = utils.get_capslongcor_path(
                caps_dir, subject, sesobj, longsubdir_dict)
        all_sublist.append(subject)
        all_sesobjlist.append(sesobj)
        all_capstargetlist.append(caps_target)
        # check if there is already a template in the CAPS directory
        if os.path.isdir(caps_target):
            # there is a template already. Update list of processed
            # subjects
            pcssd_sublist.append(subject)
            pcssd_sesobjlist.append(sesobj)
            pcssd_capstargetlist.append(caps_target)
        else:
            # there is no template folder in CAPS corresponding to the
            # subject. Update list of unprocessed subjects and
            # unprocessed sessions
            unpcssd_sublist.append(subject)
            unpcssd_sesobjlist.append(sesobj)
            unpcssd_capstargetlist.append(caps_target)

    return [
        pcssd_sublist,
        pcssd_sesobjlist,
        pcssd_capstargetlist,
        unpcssd_sublist,
        unpcssd_sesobjlist,
        unpcssd_capstargetlist,
        all_sublist,
        all_sesobjlist,
        all_capstargetlist]


def check_caps_template(
        subject_list,
        session_list,
        caps_dir):
    """Check there is not a template in CAPS already

    Check the CAPS folder does not already contain a template folder
    for the current subject and its corresponding sessions.
    Three cases are possible:
    1. CAPS folder contains no processed templates
    -> return target path for all subjects/sessions
    2. All templates have already been processed and sent to CAPS folder
    -> if force_overwrite option provided, return target path for all
        subjects/sessions, else exit
    3. Some templates (but not all) have been processed already.
    -> if force_overwrite option provided, return target path for all
        subjects/sessions, else return target path for the
        subjects/sessions that have not been processed yet

    Args:
        subject_list (list of string): list of participant_id
        session_list (list of string): list of session_id corresponding
            to a specific participant_id
        caps_dir (string): location of CAPS folder

    Returns:
        pcssd_sublist (list of strings): list of all the subjects for
            which processing was found in the CAPS folder
        pcssd_seslist2 (list of objects): list, for all subjects
            for which processing was found in the CAPS folder, of the
            list of associated sessions
        unpcssd_sublist (list of strings): list of all the subjects
            for which no processing was found in the CAPS folder
        unpcssd_seslist2 (list of objects): list, for all subjects
            for which no processing was found in the CAPS folder, of the
            list of associated session
        unpcssd_capstargetlist (string): list, for all subjects for
            which no processing was found in the CAPS folder, of
            associated CAPS folder where the template or
            longitudinal-correction processing will get stored
        all_sublist (list of strings): list of all the subjects
        all_seslist2 (list of objects): list, for all subjects of
            the list of associated sessions
        all_capstargetlist (string): list, for all subjects and
            corresponding sessions of associated CAPS folder where the
            template will get stored
    """
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils

    # get a unique list of subjects and associated per-subject list of
    # sessions
    unique_subject_list, persubject_session_list2 = utils.get_unique_subjects(
        subject_list, session_list)

    # Check processed/unprocessed subjects
    [
        pcssd_sublist,
        pcssd_seslist2,
        pcssd_capstargetlist,
        unpcssd_sublist,
        unpcssd_seslist2,
        unpcssd_capstargetlist,
        all_sublist,
        all_seslist2,
        all_capstargetlist] = utils.check_caps_processing(
            'template',
            unique_subject_list,
            persubject_session_list2,
            caps_dir)

    return [
        pcssd_sublist,
        pcssd_seslist2,
        pcssd_capstargetlist,
        unpcssd_sublist,
        unpcssd_seslist2,
        unpcssd_capstargetlist,
        all_sublist,
        all_seslist2,
        all_capstargetlist]


def to_process(
        pcssd_sublist,
        pcssd_seslist2,
        unpcssd_sublist,
        unpcssd_seslist2,
        unpcssd_capstargetlist,
        all_sublist,
        all_seslist2,
        all_capstargetlist,
        in_caps_dir,
        in_overwrite_caps,
        in_working_directory,
        in_n_procs,
        in_rundir):
    """ Get list of subjects to be processed

    The list of subjects to be processed will depend upon whether the
    user provided the option to --force-overwrite already existing
    processing or not. If the option --force-overwrite was not passed
    but existing processing was found, a .tsv file will be generated
    that contains the list of subjects/sessions already processed.

    Args:
        pcssd_sublist (list of strings): list of all the subjects for
            which processing was found in the CAPS folder
        pcssd_seslist2 (list of objects): list, for all subjects
            for which processing was found in the CAPS folder, of the
            list of associated sessions
        unpcssd_sublist (list of strings): list of all the subjects
            for which no processing was found in the CAPS folder
        unpcssd_seslist2 (list of objects): list, for all subjects
            for which no processing was found in the CAPS folder, of the
            list of associated session
        unpcssd_capstargetlist (string): list, for all subjects for
            which no processing was found in the CAPS folder, of
            associated CAPS folder where the template or
            longitudinal-correction processing will get stored
        all_sublist (list of strings): list of all the subjects
        all_seslist2 (list of objects): list, for all subjects of
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
        in_rundir (string): path to the directory where clinica is being
            launched. Used to create a .tsv file in running directory
            if the user wishes to re-run the pipeline on subjects that
            had already been processed.

    Returns:
        topcss_sublist (list of strings): list of the subjects that will
            be processed. All if --force-overwrite, else only subjects
            that have not been processed already
        topcss_seslist2 (list of list): list of the list of sessions
            that correspond to each subject to be processed
        topcss_capstargetlist (list of strings): list, for all subjects
            to be processed, of the associated CAPS folder where the
            template will get stored
        overwrite_tsv_path (string): path to the .tsv to use if the user
            wishes to overwrite the subject that have already been
            processed in the CAPS folder. None if no subjects already
            processed.
    """
    import os
    import csv
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils

    # Convert overwrite-caps flag to boolean
    force_overwrite = utils.boolean_overwrite_caps(in_overwrite_caps)

    # Generate .tsv file. Will be used to re-run with --force-overwrite
    # enabled.
    pcssd_sub_number = len(pcssd_sublist)
    if pcssd_sub_number > 0:
        # only write the .tsv file if at least one subject already
        # processed
        overwrite_tsv_path = os.path.join(in_rundir, 'overwrite.tsv')
        with open(overwrite_tsv_path, 'w') as overwrite_tsv_file:
            # create tsv file
            tsv_writer = csv.writer(
                overwrite_tsv_file, delimiter='\t')
            tsv_writer.writerow(['participant_id', 'session_id'])
            # loop through all processed subject/session
            for pcssd_sub_index in range(pcssd_sub_number):
                pcssd_sub = pcssd_sublist[pcssd_sub_index]
                pcssd_seslist = pcssd_seslist2[pcssd_sub_index]
                pcssd_seslist_number = len(pcssd_seslist)
                for pcssd_seslist_index in range(pcssd_seslist_number):
                    pcssd_ses = pcssd_seslist[pcssd_seslist_index]
                    tsv_writer.writerow([pcssd_sub, pcssd_ses])
    else:
        overwrite_tsv_path = None

    # if all subjects already processed and --force-overwrite option not
    # provided, exit
    if (not unpcssd_sublist) and (not force_overwrite):
        cl_overwrite_caps = utils.get_cl_overwritecaps(
            in_caps_dir,
            overwrite_tsv_path,
            in_working_directory,
            in_n_procs)
        error_msg = 'Error: all subjects are already processed. To overwrite'
        error_msg = '{0} the data in CAPS folder, please run:\n'.format(
            error_msg)
        error_msg = '{0} {1}'.format(error_msg, cl_overwrite_caps)
        raise ValueError(error_msg)

    # check what subjects / sessions / CAPS target dirs to process
    # depending on --force-overwrite flag
    if force_overwrite:
        topcss_sublist = all_sublist
        topcss_seslist2 = all_seslist2
        topcss_capstargetlist = all_capstargetlist
    else:
        topcss_sublist = unpcssd_sublist
        topcss_seslist2 = unpcssd_seslist2
        topcss_capstargetlist = unpcssd_capstargetlist

    return [
        topcss_sublist,
        topcss_seslist2,
        topcss_capstargetlist,
        overwrite_tsv_path]


def check_single_timepoint(subject_list, session_list2):
    """Check and warn in case subject with 1 timepoint detected

    Args:
        sublist (list of strings): list of the subjects that will be
            processed. All if --force-overwrite, else only subjects that
            have not been processed already
        seslist2 (list of list): list of the list of sessions that
            correspond to each subject to be processed

    Returns:
        N/A
    """
    import sys
    import os
    import warnings
    import colorama

    # Loop through all subjects and associated session lists
    single_timepoint_detected = False
    single_timepoint_subject_list = list()
    for subject, subses_list in zip(subject_list, session_list2):
        if len(subses_list) == 1:
            single_timepoint_detected = True
            single_timepoint_subject_list.append(subject)

    # Warn the user if a subject was found with a single time point
    if single_timepoint_detected:
        # Retrieve the system temporary folder location
        # Check operating system
        platform = sys.platform
        if platform.startswith('linux'):
            # if Linux: location is /tmp
            temp_folder_location = ' (/tmp)'
        elif platform == 'darwin':
            # if OS X: location is given by $TMPDIR
            temp_folder_location = ' ({0})'.format(os.environ['TMPDIR'])
        else:
            # other OS (e.g., other Unix, cygwin). Not supported by Clinica
            temp_folder_location = ''
        # Build the warning message
        warning_msg = '{0}'.format(colorama.Fore.RED)
        if len(single_timepoint_subject_list) == 1:
            warning_msg = '{0}Only one session was detected for the following subject:'.format(
                warning_msg)
        else:
            warning_msg = '{0}Only one session was detected for the following subjects:'.format(
                warning_msg)
        for single_timepoint_subject in single_timepoint_subject_list:
            warning_msg = '{0}\n    - {1}'.format(warning_msg, single_timepoint_subject)
        if len(single_timepoint_subject_list) == 1:
            warning_msg = '{0}\n Processing for this subject'.format(warning_msg)
        else:
            warning_msg = '{0}\n Processing for these subjects'.format(warning_msg)
        warning_msg = '{0} will be done inside the system temporary folder'.format(warning_msg)
        warning_msg = '{0}{1}.'.format(warning_msg, temp_folder_location)
        warning_msg = '{0}{1}'.format(warning_msg, colorama.Fore.RESET)
        # Show the warning
        warnings.warn(warning_msg)


def process_input_node(
        in_caps_dir,
        in_subject_list,
        in_session_list,
        in_overwrite_caps,
        in_working_directory,
        in_n_procs,
        in_rundir):
    """Carry out all processing for the input node

    Performs all the processing for the input node:
    1. Check if t1-freesurfer-cross-sectional has been correctly run
        beforehand
    2. Check if there are already template subfolders in the CAPS dir.

    Args:
        in_caps_dir (string): CAPS directory to contain the output
        in_subject_list (list): a list containing all the participant_id
        in_session_list (list): a list containing all the session_id
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
        in_rundir (str): path to the directory where clinica is being
            launched. Used to create a .tsv file in running directory
            if the user wishes to re-run the pipeline on subjects that
            had already been processed.

    Returns:
        out_subject_list (list of string): list of participant_id, where
            each participant appears only once
        out_session_list2 (list of list): list of list (list2) of
            session_id associated to any single participant
        out_capstarget_list (list of string): list of path to the
            template directory in CAPS folder for each participant
        out_unpcssd_sublist (list of string): list of participant_id,
            where each participant appears only once, wich have not
            already been (template) processed
        out_pcssd_capstargetlist (list of string): list of path to the
            template directory in CAPS folder for each participant, for
            the participants that have already been processed
        out_overwrite_tsv (str): same as in_rundir. Path to the .tsv to
            use if the user wishes to overwrite the subjects that have
            already been processed in the CAPS folder. None if no
            subjects already processed.
    """
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils

    # sanity check: check at least one subject
    if not in_subject_list:
        error_msg = 'Error: there should be at least one subject'
        error_msg = '{0} to conduct the longitudinal study'.format(error_msg)
        raise ValueError(error_msg)

    # check if we have the necessary (cross-sectional) recon-all data
    # prior to running the template creation. Will crash if we do not
    # have this.
    utils.check_xsectional_reconalled(
        in_caps_dir, in_subject_list, in_session_list)

    # Check target CAPS target dir for template creation before
    # conducting any processing
    [
        pcssd_sublist,
        pcssd_seslist2,
        out_pcssd_capstargetlist,
        out_unpcssd_sublist,
        unpcssd_seslist2,
        unpcssd_capstargetlist,
        all_sublist,
        all_seslist2,
        all_capstargetlist] = utils.check_caps_template(
            in_subject_list, in_session_list, in_caps_dir)

    # Deduce subjects to be processed next (depending on whether
    # '--force-overwrite' flag was provided or not)
    [
        out_subject_list,
        out_session_list2,
        out_capstarget_list,
        out_overwrite_tsv] = utils.to_process(
            pcssd_sublist,
            pcssd_seslist2,
            out_unpcssd_sublist,
            unpcssd_seslist2,
            unpcssd_capstargetlist,
            all_sublist,
            all_seslist2,
            all_capstargetlist,
            in_caps_dir,
            in_overwrite_caps,
            in_working_directory,
            in_n_procs,
            in_rundir)

    # Check if any subject can be found for which only one time point
    # has been found, in which case warn users of potential issues
    utils.check_single_timepoint(out_subject_list, out_session_list2)

    return [
        out_subject_list,
        out_session_list2,
        out_capstarget_list,
        out_unpcssd_sublist,
        out_pcssd_capstargetlist,
        out_overwrite_tsv]


def store_reconallbase_results(in_subject_list):
    """Create a folder to store all recon-all -base outputs

    This creates a 'store 'folder that will store, for each subject, the
    recon-all -base output. The function returns the path to the store
    folder.
    The folder is initially empty. At a later stage, a MapNode is run
    on each subject to carry out the recon-all -base in a temp folder
    and when completed, send the result to this store folder. This
    way, the map node can be re-run (in case it crashed) and won't have
    to re-compute the template if already done for any particular
    subject.

    Args:
        in_subject_list (list): list of participant_id. This is only
            used to connect the input_node to the store_reconall_node
            and serves no other purpose within the current function.

    Returns:
        out_workdirstore_path (string): path to the 'store' folder that
            will contain the output for recon-all -base on all subjects
            when the command is run in a subsequent node
    """
    import os
    import errno

    # create subject-specific path in Clinica's working directory
    # for the FreeSurfer 'recon-all' commands
    out_workdirstore_path = os.path.abspath("./subjects")
    try:
        os.mkdir(out_workdirstore_path)
    except OSError as oserror:
        if oserror.errno != errno.EEXIST:
            raise

    return out_workdirstore_path


def get_reconallbase_flags(in_subject, in_session_list):
    """Create reconall -base flags

    Creates the flags that are required to run the recon-all -base for
    each subject (i.e. to create the unbiased template for all time
    points for the subject).

    Args:
        in_subject (string): (unique) participant_id
        in_session_list (list of string): list of session_id
            corresponding to a specific participant_id containing all
            session symlinks for any particular subject_id

    Returns:
        out_reconallbase_flags (string): all the flags to run
            recon-all -base for any particular subject and their
            associated list of sessions (time points)
    """
    # template ID: subject name
    out_reconallbase_flags = '-base {0}'.format(in_subject)

    # timepoints
    for session in in_session_list:
        sub_ses = '{0}_{1}'.format(in_subject, session)
        out_reconallbase_flags = '{0} -tp {1}'.format(
            out_reconallbase_flags, sub_ses)

    # -all directive
    out_reconallbase_flags = '{0} -all'.format(out_reconallbase_flags)

    return out_reconallbase_flags


def check_reconall_base_single(subjects_dir, in_subject):
    """Check if recon-all base run successfully for any single subject

    Args:
        subjects_dir (string): CAPS subdirectory containing the subject
            (absolute path)
        in_subject (string): participant_id

    Returns:
        out_template_created (Boolean): True if success, False otherwise
    """
    import os
    import subprocess

    # check if longitudinal base subfolder exists
    subjectid_longitudinal = "{0}/{1}".format(subjects_dir, in_subject)
    if not os.path.isdir(subjectid_longitudinal):
        error_msg = 'Error: {0} does not exist'.format(subjectid_longitudinal)
        raise IOError(error_msg)

    # check the recon-all.log
    log_file = os.path.join(subjectid_longitudinal, 'scripts', 'recon-all.log')
    if os.path.isfile(log_file):
        last_line = subprocess.check_output(['tail', '-1', log_file])
        if b'finished without error' not in last_line:
            out_template_created = False
            error_msg = 'Subject {0} has not been processed yet'.format(
                in_subject)
            error_msg = '{0}. Re-run t1-freesurfer-template.'.format(error_msg)
            raise IOError(error_msg)
        else:
            out_template_created = True
    else:
        error_msg = 'Error: {0} does not exist'.format(log_file)
        raise IOError(error_msg)

    return out_template_created


def create_fssubdir_path(subject, session_list):
    """Create a subdir to FS recon-all

    This creates a subdir that will contain the outpute of FreeSurfer
    'recon-all -base' command. Distinguishes between two cases:
    1) only one session is linked to the subject
    - there currently (as of 22 Feb 2019) is a bug in Freesurfer
    recon-all -base, which in some cases (e.g., only one time point),
    will crash as it's trying to write lines too long for the shell to
    handle. This is caused by the path to FreeSurfer SUBJECT_DIR being
    too long itself. The current function works around this issue by
    checking if there only is one session associated to a subject, and
    in that case, putting the SUBJECT_DIR inside the system temporary
    folder (Linux: /tmp, OS X: $TMPDIR) so that its path is as short as
    possible.
    - This causes an issue as sytem temporary folders are automatically
    cleaned up (e.g., every time the computer is rebooted or at
    periodical instances). In that case, a crash within the node
    corresponding to the current function means we can only re-run the
    pipeline safely if we re-do the recon-all -base computation in a new
    temporary location
    2) More than one session are linked to the subject: creates the
    temporary folder inside the working directory.
    The subdir is eventually moved to another store location (in
    function run_reconallbase), so that folders for both case 1) and
    case 2) are ultimately stored in the same place

    Args:
        subject (string): subject ID
        session_list (list of string): list of sessions ID associated to
            the current subject

    Returns:
        fssubdir_path (string): path to the subdir containing the output
            of 'recon-all -base'.
    """
    import os
    import errno
    import tempfile

    if len(session_list) == 1:
        # create FS subdir in temporary folder
        fssubdir_path = tempfile.mkdtemp()
    else:
        # create FS subdir in working dir, current node
        fssubdir_path = os.path.join(os.path.abspath('./'), subject)
        # create folder if not exists (else recon-all computations will
        # be resumed in the existing folder)
        try:
            os.mkdir(fssubdir_path)
        except OSError as oserror:
            if oserror.errno != errno.EEXIST:
                raise

    return fssubdir_path


def run_reconallbase(
        in_caps_dir,
        in_subject,
        in_session_list,
        in_reconallbase_flags,
        in_workdirstore_path):
    """Run recon-all -base function

    Run freesurfer template creation for each unique subject and their
    associated list of sessions. The recon-all command is built with
    previously computed arguments (flags).

    Args:
        in_caps_dir (string): CAPS directory to contain the output
        in_subject (string): participant_id
        in_session_list (list of string): list of all the sessions
            corresponding to the subject
        in_reconallbase_flags (string): all the flags to run
            recon-all -base for any particular subject and their
            associated list of sessions (time points)
        in_workdirstore_path (string): path to the 'store' folder that
            contain the output for recon-all -base on all subjects
            after the command has been run

    Returns:
        out_template_created (Boolean): True if successful, False
            otherwise
    """
    import os
    import errno
    import subprocess
    import shutil
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils

    # check if current subject has been processed already
    # (i.e., if corresponding folder already stored to working dir)
    workdir_subject_path = "{0}/{1}".format(in_workdirstore_path, in_subject)
    if not os.path.isdir(workdir_subject_path):
        # Hasn't been processed, or crashed before being transferred to
        # the working directory. We re-run recon-all -base

        # create temporary folder (store path in .txt file for
        # traceability) and insert a FreeSurfer SUBJECT_DIR there
        fssubdir_path = utils.create_fssubdir_path(in_subject, in_session_list)
        storepath_filename = os.path.abspath(
            "{0}/{1}_store_path.txt".format(in_workdirstore_path, in_subject))
        # Note: we might want to check first if file already exists and
        #       send a warning accordingly
        with open(storepath_filename, 'w') as storepath_file:
            storepath_file.write(fssubdir_path)

        # add symlinks to all subject_section cross-sectional run in the
        # FreeSurfer SUBJECT_DIR (temporary folder)
        caps_path = os.path.expanduser(in_caps_dir)
        caps_dir = os.path.join(caps_path, 'subjects')
        for session in in_session_list:
            # retrieve path to {subject,session} cross-sectional reconall
            # runs
            caps_subses_path = "{0}/{1}/{2}/t1/freesurfer_cross_sectional/{1}_{2}".format(
                caps_dir, in_subject, session)
            # define symlink path
            symlink_subses_path = "{0}/{1}_{2}".format(
                fssubdir_path, in_subject, session)
            if not os.path.exists(symlink_subses_path):
                os.symlink(caps_subses_path, symlink_subses_path)

        # run recon-all -base in temp FreeSurfer SUBJECT_DIR
        reconallbase_command = 'recon-all {0} -sd {1}'.format(
            in_reconallbase_flags, fssubdir_path)
        subprocess_run_reconallbase = subprocess.run(
            reconallbase_command,
            shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if subprocess_run_reconallbase.returncode != 0:
            raise ValueError('recon-all -base failed, returned non-zero code')

        # check if completed successfully
        out_template_created = utils.check_reconall_base_single(
            fssubdir_path, in_subject)
        if out_template_created:
            # create subject-specific workdir store location for
            # recon-all -base results
            workdirstore_subject_path = os.path.join(
                in_workdirstore_path, in_subject)
            try:
                os.mkdir(workdirstore_subject_path)
            except OSError as oserror:
                if oserror.errno != errno.EEXIST:
                    raise
            # move the recon-all -base output to working directory
            list_filefolder = os.listdir(fssubdir_path)
            for entry in list_filefolder:
                shutil.move(
                    '{0}/{1}'.format(fssubdir_path, entry),
                    workdirstore_subject_path)
        else:
            # return an error
            error_msg = 'Error: recon-all -base did not complete.'
            error_msg = '{0} Template not created.'.format(error_msg)
            raise IOError(error_msg)

        # remove temporary folder
        # for now do nothing, the associated system temporary folder
        # should be empty any way...

    return out_template_created


def safer_rmtree(caps_target, caps_dir):
    """Only remove CAPS target if it's as subfolder of CAPS dir

    Will check if the CAPS target folder is a subfolder of CAPS dir
    prior to deleting it.
    This is to reduce the possibility of deleting files outside the CAPS
    dir in case of an accidental erroneous definition of the CAPS target
    in code.

    Args:
        caps_target (string): path to CAPS target subfolder for the
            template connected to the current subject
        caps_dir (string): CAPS directory to contain the output
            Used to make sure the folder we are overwriting (if
            overwrite-caps option provided) are subolders of the input
            CAPS dir

    Returns:
        N/A
    """
    import os
    import shutil

    # get absolute path + real (i.e., actual path for symlinks)
    caps_target_real = os.path.realpath(caps_target)
    caps_dir_real = os.path.realpath(caps_dir)
    caps_dir_real_sep = "{0}{1}".format(caps_dir_real, os.sep)

    # check if the path to image encompasses the path to folder
    if not caps_target_real.startswith(caps_dir_real_sep):
        error_msg = 'Error: CAPS target {0}'.format(caps_target)
        error_msg = '{0} is not a subfolder of caps_dir {1}.'.format(
            error_msg, caps_dir)
        error_msg = '{0} Cannot safely remove.'.format(error_msg)
        raise IOError(error_msg)

    # remove the caps target
    if os.path.exists(caps_target_real):
        shutil.rmtree(caps_target_real)
    else:
        error_msg = 'Error: folder {0} does not exist'.format(caps_target)
        raise IOError(error_msg)


def copy_to_caps(
        in_subject,
        in_subject_dir,
        in_caps_target,
        in_caps_dir,
        in_overwrite_caps,
        in_template_created):
    """Copy template folder to CAPS dir

    Move the template folder created by recon-all base to its correct
    location in the CAPS tree structure.

    Args:
        in_subject (string): (unique) participant_id
        in_subject_dir (string): location of the template folder (in the
            working directory)
        in_caps_target (string): path to CAPS target subfolder for the
            template connected to the current subject
        in_caps_dir (string): CAPS directory to contain the output
            Used to make sure the folder we are overwriting (if
            overwrite-caps option provided) are subolders of the input
            CAPS dir
        in_overwrite_caps (string): Option provided by the user to state
            whether the data in the CAPS folder (i.e.,
            previously-computed template) should be overwritten or not
            by a new run. By default, do not overwrite. Will overwrite
            for the following values: 'true', 'True' and 'TRUE'
        in_template_created (Boolean): Used to force template
            recon-all -base to be run before moving anything to the CAPS
            directory

    Returns:
        out_copy2_caps (boolean): flag to indicate whether or not the
            copy was successful
    """
    import os
    import shutil
    import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils

    # Convert overwrite-caps flag to boolean
    force_overwrite = utils.boolean_overwrite_caps(in_overwrite_caps)

    # define template location in the working directory
    # (repeating twice in_subject because in_subject_dir/in_subject
    # contains both in_subject template and cross-sectional
    # symlinks)
    wd_template_location = os.path.join(in_subject_dir, in_subject, in_subject)

    # Check if overwrite-caps option provided. If provided, delete the
    # target subdirectory in CAPS folder
    if force_overwrite:
        # check the target subdirectory exists in CAPS
        if os.path.exists(in_caps_target):
            # remove
            utils.safer_rmtree(in_caps_target, in_caps_dir)

    # copy the template from working directory to CAPS directory
    shutil.copytree(wd_template_location, in_caps_target)
    out_copy2_caps = True

    return out_copy2_caps


def sendto_longcorr(
        in_copy2_caps,
        in_unpcssd_sublist,
        in_pcssd_capstargetlist,
        in_overwrite_tsv):
    """Identity function

    Used to store the data from the template pipeline to be sent to the
    longitudinal pipeline.

    Args:
        in_copy2_caps (Boolean): Force the copy2caps_node to be run
            before sending data to longitudinal-correction (i.e., before
            starting the longitudinal-correction pipeline)
        in_unpcssd_sublist (list of string): list of participant_id,
            where each participant appears only once, wich have not
            already been (template) processed
        in_pcssd_capstargetlist (list of string): list of path to the
            template directory in CAPS folder for each participant, for
            the participants that have already been processed
        in_overwrite_tsv (str): path to the .tsv file to use if the user
            wishes to re-run the pipeline on subjects that had already
            been processed. None if no subjects already processed.

    Returns:
        out_unpcssd_sublist (list of string): same as in_unpcssd_sublist
        out_pcssd_capstargetlist (list of string): same as
            in_pcssd_capstargetlist
        out_overwrite_tsv (str): same as in_overwrite_tsv
    """
    out_unpcssd_sublist = in_unpcssd_sublist
    out_pcssd_capstargetlist = in_pcssd_capstargetlist
    out_overwrite_tsv = in_overwrite_tsv

    return out_unpcssd_sublist, out_pcssd_capstargetlist, out_overwrite_tsv
