# coding: utf8


def zip_nii(in_file, same_dir=False):
    from os import getcwd
    from os.path import abspath, join
    import gzip
    import shutil
    from nipype.utils.filemanip import split_filename
    from traits.trait_base import _Undefined

    if (in_file is None) or isinstance(in_file, _Undefined):
        return None

    if not isinstance(in_file, str):  # type(in_file) is list:
        return [zip_nii(f, same_dir) for f in in_file]

    orig_dir, base, ext = split_filename(str(in_file))

    # Already compressed
    if ext[-3:].lower() == ".gz":
        return in_file
    # Not compressed

    if same_dir:
        out_file = abspath(join(orig_dir, base + ext + '.gz'))
    else:
        out_file = abspath(join(getcwd(), base + ext + '.gz'))

    with open(in_file, 'rb') as f_in, gzip.open(out_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    return out_file


def unzip_nii(in_file):
    from nipype.utils.filemanip import split_filename
    from nipype.algorithms.misc import Gunzip
    from traits.trait_base import _Undefined

    if (in_file is None) or isinstance(in_file, _Undefined):
        return None

    if not isinstance(in_file, str):  # type(in_file) is list:
        return [unzip_nii(f) for f in in_file]

    _, base, ext = split_filename(in_file)

    # Not compressed
    if ext[-3:].lower() != ".gz":
        return in_file
    # Compressed
    gunzip = Gunzip(in_file=in_file)
    gunzip.run()
    return gunzip.aggregate_outputs().out_file


def fix_join(path, *paths):
    # This workaround is used in pipelines like DWIPreprocessingUsingT1
    # In the workflow.connect part, you can use some function that are used as string, causing an import error
    import os
    return os.path.join(path, *paths)


def save_participants_sessions(participant_ids, session_ids, out_folder, out_file=None):
    """
    Save <participant_ids> <session_ids> in <out_folder>/<out_file> TSV file.
    """
    import os
    import errno
    import pandas
    from clinica.utils.stream import cprint

    assert(len(participant_ids) == len(session_ids))

    try:
        os.makedirs(out_folder)
    except OSError as e:
        if e.errno != errno.EEXIST:  # EEXIST: folder already exists
            raise e

    if out_file:
        tsv_file = os.path.join(out_folder, out_file)
    else:
        tsv_file = os.path.join(out_folder, 'participants.tsv')

    try:
        data = pandas.DataFrame({
            'participant_id': participant_ids,
            'session_id': session_ids,
        })
        data.to_csv(tsv_file, sep='\t', index=False, encoding='utf-8')
    except Exception as e:
        cprint("Impossible to save %s with pandas" % out_file)
        raise e


def get_subject_id(bids_or_caps_file):
    """
    Extracts "sub-<participant_id>_ses-<session_label>" from BIDS or CAPS file
    """
    import re

    m = re.search(r'(sub-[a-zA-Z0-9]+)/(ses-[a-zA-Z0-9]+)', bids_or_caps_file)

    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format.'
            ' It does not contain the subject and session information.')

    subject_id = m.group(1) + '_' + m.group(2)

    return subject_id


def extract_image_ids(bids_or_caps_files):
    """Extract image IDs (e.g. ['sub-CLNC01_ses-M00', 'sub-CLNC01_ses-M18']  from `bids_or_caps_files`."""
    import re
    id_bids_or_caps_files = [re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)', file).group()
                             for file in bids_or_caps_files]
    return id_bids_or_caps_files


def check_bids_folder(bids_directory):
    """
    check_bids_folder function checks the following items:
        - bids_directory is a string
        - the provided path exists and is a directory
        - provided path is not a CAPS folder (BIDS and CAPS could be swapped by user). We simply check that there is
          not a folder called 'subjects' in the provided path (that exists in CAPS hierarchy)
        - provided folder is not empty
        - provided folder must contains at least one directory whose name starts with 'sub-'
    """
    from os.path import isdir, join
    from os import listdir
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaBIDSError

    assert isinstance(bids_directory, str), 'Argument you provided to check_bids_folder() is not a string.'

    if not isdir(bids_directory):
        raise ClinicaBIDSError(Fore.RED + '\n[Error] The BIDS directory you gave is not a folder.\n' + Fore.RESET
                               + Fore.YELLOW + '\nError explanations:\n' + Fore.RESET
                               + ' - Clinica expected the following path to be a folder:' + Fore.BLUE + bids_directory
                               + Fore.RESET + '\n'
                               + ' - If you gave relative path, did you run Clinica on the good folder?')

    if isdir(join(bids_directory, 'subjects')):
        raise ClinicaBIDSError(Fore.RED + '\n[Error] The BIDS directory (' + bids_directory + ') you provided seems to '
                               + 'be a CAPS directory due to the presence of a \'subjects\' folder.' + Fore.RESET)

    if len(listdir(bids_directory)) == 0:
        raise ClinicaBIDSError(Fore.RED + '\n[Error] The BIDS directory you provided  is empty. (' + bids_directory
                               + ').' + Fore.RESET)

    if len([item for item in listdir(bids_directory) if item.startswith('sub-')]) == 0:
        raise ClinicaBIDSError(Fore.RED + '\n[Error] Your BIDS directory does not contains a single folder whose name '
                               + 'starts with \'sub-\'. Check that your folder follow BIDS standard' + Fore.RESET)


def check_caps_folder(caps_directory):
    """
    check_caps_folder function checks the following items:
        - caps_directory is a string
        - the provided path exists and is a directory
        - provided path is not a BIDS folder (BIDS and CAPS could be swapped by user). We simply check that there is
          not a folder whose name starts with 'sub-' in the provided path (that exists in BIDS hierarchy)
    Keep in mind that CAPS folder can be empty
    """
    from os import listdir
    import os
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaCAPSError

    assert isinstance(caps_directory, str), 'Argument you provided to check_caps_folder() is not a string.'

    if not os.path.isdir(caps_directory):
        raise ClinicaCAPSError(Fore.RED + '\n[Error] The CAPS directory you gave is not a folder.\n' + Fore.RESET
                               + Fore.YELLOW + '\nError explanations:\n' + Fore.RESET
                               + ' - Clinica expected the following path to be a folder:' + Fore.BLUE + caps_directory
                               + Fore.RESET + '\n'
                               + ' - If you gave relative path, did you run Clinica on the good folder?')

    sub_folders = [item for item in listdir(caps_directory) if item.startswith('sub-')]
    if len(sub_folders) > 0:
        error_string = '\n[Error] Your CAPS directory contains at least one folder whose name ' \
                       + 'starts with \'sub-\'. Check that you did not swap BIDS and CAPS folders.\n' \
                       + ' Folder(s) found that match(es) BIDS architecture:\n'
        for dir in sub_folders:
            error_string += '\t' + dir + '\n'
        error_string += 'A CAPS directory has a folder \'subjects\' at its root, in which are stored the output ' \
                        + 'of the pipeline for each subject.'
        raise ClinicaCAPSError(error_string)


def extract_crash_files_from_log_file(filename):
    """Extract crash files (*.pklz) from `filename`.
    """
    import os
    import re

    assert(os.path.isfile(filename)),\
        'extract_crash_files_from_log_file: filename parameter is not a file (%s)' % filename

    log_file = open(filename, "r")
    crash_files = []
    for line in log_file:
        if re.match("(.*)crashfile:(.*)", line):
            crash_files.append(line.replace('\t crashfile:', '').replace('\n', ''))

    return crash_files


def read_participant_tsv(tsv_file):
    """Extract participant IDs and session IDs from TSV file.

    Raise:
        ClinicaException if tsv_file is not a file
        ClinicaException if participant_id or session_id column is missing from TSV file
    """
    import os
    import pandas as pd
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaException

    if not os.path.isfile(tsv_file):
        raise ClinicaException(
            "\n%s[Error] The TSV file you gave is not a file.%s\n"
            "\n%sError explanations:%s\n"
            " - Clinica expected the following path to be a file: %s%s%s\n"
            " - If you gave relative path, did you run Clinica on the good folder?" %
            (Fore.RED, Fore.RESET,
             Fore.YELLOW, Fore.RESET,
             Fore.BLUE, tsv_file, Fore.RESET)
        )
    ss_df = pd.io.parsers.read_csv(tsv_file, sep='\t')
    if 'participant_id' not in list(ss_df.columns.values):
        raise ClinicaException(
            "\n%s[Error] The TSV file does not contain participant_id column (path: %s)%s" %
            (Fore.RED, tsv_file, Fore.RESET)
        )
    if 'session_id' not in list(ss_df.columns.values):
        raise ClinicaException(
            "\n%s[Error] The TSV file does not contain session_id column (path: %s)%s" %
            (Fore.RED, tsv_file, Fore.RESET)
        )
    participants = list(ss_df.participant_id)
    sessions = list(ss_df.session_id)

    # Remove potential whitespace in participant_id or session_id
    return [sub.strip(' ') for sub in participants], [ses.strip(' ') for ses in sessions]


def get_subject_session_list(input_dir, ss_file=None, is_bids_dir=True, use_session_tsv=False):
    """Parses a BIDS or CAPS directory to get the subjects and sessions.

    This function lists all the subjects and sessions based on the content of
    the BIDS or CAPS directory or (if specified) on the provided
    subject-sessions TSV file.

    Args:
        input_dir: A BIDS or CAPS directory path.
        ss_file: A subjects-sessions file (.tsv format).
        is_bids_dir: Indicates if input_dir is a BIDS or CAPS directory
        use_session_tsv (boolean): Specify if the list uses the sessions listed in the sessions.tsv files

    Returns:
        subjects: A subjects list.
        sessions: A sessions list.
    """
    import os
    import tempfile
    from time import time, strftime, localtime
    import clinica.iotools.utils.data_handling as cdh

    if not ss_file:
        output_dir = tempfile.mkdtemp()
        timestamp = strftime('%Y%m%d_%H%M%S', localtime(time()))
        tsv_file = 'subjects_sessions_list_%s.tsv' % timestamp
        ss_file = os.path.join(output_dir, tsv_file)

        cdh.create_subs_sess_list(
            input_dir=input_dir,
            output_dir=output_dir,
            file_name=tsv_file,
            is_bids_dir=is_bids_dir,
            use_session_tsv=use_session_tsv)

    participant_ids, session_ids = read_participant_tsv(ss_file)
    return session_ids, participant_ids
