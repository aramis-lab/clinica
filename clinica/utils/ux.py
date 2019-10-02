# coding: utf8


def print_images_to_process(list_participant_id, list_session_id):
    """Print which images will be processed by the pipeline."""
    from .stream import cprint
    images_to_process = ', '.join(participant_id + '|' + session_id
                                  for participant_id, session_id in zip(list_participant_id, list_session_id)
                                  )
    if len(list_participant_id) < 25:
        cprint('The pipeline will be run on the following %s image(s): %s' %
               (len(list_participant_id), images_to_process))
    else:
        cprint('The pipeline will be run on %s images (too many IDs to display)' %
               (len(list_participant_id)))


def print_begin_image(image_id, list_keys=None, list_values=None):
    """Print begin run pipeline message for a given image `image_id`."""
    import datetime
    from colorama import Fore
    from .stream import cprint

    assert(len(list_keys) == len(list_values))

    begin_message = 'Running pipeline for %s' % (image_id.replace('_', '|'))
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

    end_message = '%s has completed' % image_id.replace('_', '|')
    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s %s' % (Fore.GREEN, now, Fore.RESET, end_message))


def print_no_image_to_process():
    """Print end pipeline message when no image was found."""
    from colorama import Fore
    from .stream import cprint

    cprint('%sEither all the images were already run by the pipeline or no image was found to run the pipeline. '
           'The program will now exit.%s' % (Fore.BLUE, Fore.RESET))


def print_end_pipeline(cli_name, working_directory):
    """Print end pipeline message after its execution."""
    import datetime
    from colorama import Fore
    from .stream import cprint

    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s The %s pipeline has completed. You can now delete the working directory (%s).' %
           (Fore.GREEN, now, Fore.RESET, cli_name, working_directory))


def print_failed_images(cli_name, image_ids, log_file):
    """Print missing images in CAPS folder after a RuntimeError from Nipype."""
    import datetime
    from colorama import Fore
    from .stream import cprint
    from .io import extract_crash_files_from_log_file

    missing_caps = ', '.join(
        image_ids[i].split('_')[0][4:] + '|' + image_ids[i].split('_')[1][4:]
        for i in range(len(image_ids)))

    # Display summary of the pipeline
    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('\n%s[%s] The %s pipeline finished with errors.%s\n' %
           (Fore.RED, now, cli_name, Fore.RESET))
    cprint('%sCAPS outputs were not found for %s subject(s): %s%s\n' %
           (Fore.RED, len(image_ids), missing_caps, Fore.RESET))
    cprint('%sError details can be found either by opening the log file (%s) or '
           'by opening the crash file(s) with the following command(s):%s' %
           (Fore.YELLOW, log_file, Fore.RESET))

    crash_files = extract_crash_files_from_log_file(log_file)
    for file in crash_files:
        cprint('%s- nipypecli crash %s%s' %
               (Fore.YELLOW, file,  Fore.RESET))
