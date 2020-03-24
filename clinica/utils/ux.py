# coding: utf8


LINES_TO_DISPLAY = 25


def print_images_to_process(list_participant_id, list_session_id):
    """Print which images will be processed by the pipeline."""
    from .stream import cprint
    from clinica.utils.participant import get_unique_subjects
    unique_participants, sessions_per_participant = get_unique_subjects(list_participant_id, list_session_id)

    cprint('The pipeline will be run on the following %s image(s):' % len(list_participant_id))
    for i in range(0, min(len(unique_participants), LINES_TO_DISPLAY)):
        sessions_i_th_participant = ', '.join(s_id for s_id in sessions_per_participant[i])
        cprint("\t%s | %s," % (unique_participants[i], sessions_i_th_participant))

    if len(unique_participants) > LINES_TO_DISPLAY:
        cprint("\t...")
        sessions_last_participant = ', '.join(s_id for s_id in sessions_per_participant[-1])
        cprint("\t%s | %s" % (unique_participants[-1], sessions_last_participant))


def print_begin_image(image_id, list_keys=None, list_values=None):
    """Print begin run pipeline message for a given image `image_id`."""
    import datetime
    from colorama import Fore
    from .stream import cprint

    if list_keys is not None:
        assert(len(list_keys) == len(list_values))

    begin_message = 'Running pipeline for %s' % (image_id.replace('_', ' | '))
    if list_keys and list_values:
        begin_message += ' ('
        begin_message += ', '.join(key + ' = ' + value for key, value in zip(list_keys, list_values))
        begin_message += ')'
    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s %s' % (Fore.BLUE, now, Fore.RESET, begin_message))


def print_end_image(image_id):
    """Print end pipeline message for a given image `image_id`."""
    import datetime
    from colorama import Fore
    from .stream import cprint

    end_message = '%s has completed' % image_id.replace('_', ' | ')
    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s %s' % (Fore.GREEN, now, Fore.RESET, end_message))


def print_end_pipeline(cli_name, working_directory, working_directory_was_specified):
    """Print end pipeline message after its execution."""
    import os
    import datetime
    from colorama import Fore
    from .stream import cprint
    from os.path import abspath

    now = datetime.datetime.now().strftime('%H:%M:%S')
    if working_directory_was_specified:
        cprint('%s[%s]%s The %s pipeline has completed. You can now delete the working directory (%s).' %
               (Fore.GREEN, now, Fore.RESET, cli_name, os.path.join(working_directory, cli_name)))
    else:
        cprint('%s[%s]%s The %s pipeline has completed. Working directory was automatically deleted.' %
               (Fore.GREEN, now, Fore.RESET, cli_name))


def print_failed_images(cli_name, image_ids):
    """Print missing images in CAPS folder after a RuntimeError from Nipype."""
    import datetime
    from colorama import Fore
    from .filemanip import extract_subjects_sessions_from_filename
    from .stream import cprint
    from clinica.utils.participant import get_unique_subjects
    list_participant_id, list_session_id = extract_subjects_sessions_from_filename(image_ids)
    unique_participants, sessions_per_participant = get_unique_subjects(list_participant_id, list_session_id)

    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('\n%s[%s] The %s pipeline finished with errors.%s\n' %
           (Fore.RED, now, cli_name, Fore.RESET))
    cprint('%sCAPS outputs were not found for %s image(s):%s' %
           (Fore.RED, len(image_ids), Fore.RESET))
    for i in range(0, min(len(unique_participants), LINES_TO_DISPLAY)):
        sessions_i_th_participant = ', '.join(s_id for s_id in sessions_per_participant[i])
        cprint("\t%s%s | %s%s" % (Fore.RED, unique_participants[i], sessions_i_th_participant, Fore.RESET))

    if len(unique_participants) > LINES_TO_DISPLAY:
        cprint("\t...")
        sessions_last_participant = ', '.join(s_id for s_id in sessions_per_participant[-1])
        cprint("\t%s%s | %s%s" % (Fore.RED, unique_participants[-1], sessions_last_participant, Fore.RESET))


def print_crash_files_and_exit(log_file, working_directory):
    """Print command(s) to type in order to extract details after a Nipype RuntimeError and exit with an exception."""
    from colorama import Fore
    from .filemanip import extract_crash_files_from_log_file
    from .exceptions import ClinicaException
    from .stream import cprint

    cprint('%s\nError details can be found by opening the crash file(s) with the following command(s):%s' %
           (Fore.YELLOW, Fore.RESET))

    crash_files = extract_crash_files_from_log_file(log_file)
    for file in crash_files:
        cprint('%s- nipypecli crash %s%s' %
               (Fore.YELLOW, file,  Fore.RESET))

    cprint('%s\n'
           'If your pipeline crashed due to lack of space of network issues, '
           're-run the pipeline with the working directory (-wd %s).\n'
           'Known issues are displayed here: http://www.clinica.run/doc/InteractingWithClinica/#known-issues\n'
           'Otherwise, you can delete it.%s' %
           (Fore.YELLOW, working_directory, Fore.RESET))
    # Force the display of "Documentation can be found..."
    raise ClinicaException('')


def print_groups_in_caps_directory(caps_directory):
    """Print group IDs based on `caps_directory`/groups folder."""
    from .group import extract_group_ids
    from .stream import cprint

    group_ids = extract_group_ids(caps_directory)
    if group_ids == ['']:
        cprint('No group was found in CAPS directory')
    else:
        cprint("Groups that exist in your CAPS directory are %s." %
               ', '.join(g_id.replace('group-', '') for g_id in group_ids))
