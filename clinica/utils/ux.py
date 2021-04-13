# coding: utf8

"""This module gathers formatted messages that are displayed when running Clinica.

These functions are mainly called by the pipelines.
"""

LINES_TO_DISPLAY = 25


def print_images_to_process(list_participant_id, list_session_id):
    """Print which images will be processed by the pipeline."""
    from clinica.utils.participant import get_unique_subjects

    from .stream import cprint

    unique_participants, sessions_per_participant = get_unique_subjects(
        list_participant_id, list_session_id
    )

    cprint(
        f"The pipeline will be run on the following {len(list_participant_id)} image(s):"
    )
    for i in range(0, min(len(unique_participants), LINES_TO_DISPLAY)):
        sessions_i_th_participant = ", ".join(
            s_id for s_id in sessions_per_participant[i]
        )
        cprint(f"\t{unique_participants[i]} | {sessions_i_th_participant},")

    if len(unique_participants) > LINES_TO_DISPLAY:
        cprint("\t...")
        sessions_last_participant = ", ".join(
            s_id for s_id in sessions_per_participant[-1]
        )
        cprint(f"\t{unique_participants[-1]} | {sessions_last_participant},")


def print_begin_image(image_id, list_keys=None, list_values=None):
    """Print begin run pipeline message for a given image `image_id`."""
    import datetime

    from colorama import Fore

    from .stream import cprint

    if list_keys is not None:
        assert len(list_keys) == len(list_values)

    begin_message = f"Running pipeline for {image_id.replace('_', ' | ')}"
    if list_keys and list_values:
        begin_message += " ("
        begin_message += ", ".join(
            f"{key} = {key_value}" for key, key_value in zip(list_keys, list_values)
        )
        begin_message += ")"
    now = datetime.datetime.now().strftime("%H:%M:%S")
    cprint(f"{Fore.BLUE}[{now}]{Fore.RESET} {begin_message}")


def print_end_image(image_id):
    """Print end pipeline message for a given image `image_id`."""
    import datetime

    from colorama import Fore

    from .stream import cprint

    end_message = f"{image_id.replace(' ', ' | ')} has completed"
    now = datetime.datetime.now().strftime("%H:%M:%S")
    cprint(f"{Fore.GREEN}[{now}]{Fore.RESET} {end_message}")


def print_end_pipeline(cli_name, working_directory, working_directory_was_specified):
    """Print end pipeline message after its execution."""
    import datetime
    import os

    from colorama import Fore

    from .stream import cprint

    now = datetime.datetime.now().strftime("%H:%M:%S")
    if working_directory_was_specified:
        cprint(
            f"{Fore.GREEN}[{now}]{Fore.RESET} The {cli_name} pipeline has completed. "
            f"You can now delete the working directory ({os.path.join(working_directory, cli_name)})."
        )
    else:
        cprint(
            f"{Fore.GREEN}[{now}]{Fore.RESET} The {cli_name} pipeline has completed. "
            f"Working directory was automatically deleted."
        )


def print_failed_images(cli_name, image_ids):
    """Print missing images in CAPS folder after a RuntimeError from Nipype."""
    import datetime

    from colorama import Fore

    from clinica.utils.participant import get_unique_subjects

    from .filemanip import extract_subjects_sessions_from_filename
    from .stream import cprint

    list_participant_id, list_session_id = extract_subjects_sessions_from_filename(
        image_ids
    )
    unique_participants, sessions_per_participant = get_unique_subjects(
        list_participant_id, list_session_id
    )

    now = datetime.datetime.now().strftime("%H:%M:%S")
    cprint(
        f"\n{Fore.RED}[{now}] The {cli_name} pipeline finished with errors.{Fore.RESET}\n"
    )
    cprint(
        f"{Fore.RED}CAPS outputs were not found for {len(image_ids)} image(s):{Fore.RESET}"
    )
    for i in range(0, min(len(unique_participants), LINES_TO_DISPLAY)):
        sessions_i_th_participant = ", ".join(
            s_id for s_id in sessions_per_participant[i]
        )
        cprint(
            f"\t{Fore.RED}{unique_participants[i]} | {sessions_i_th_participant}{Fore.RESET}"
        )

    if len(unique_participants) > LINES_TO_DISPLAY:
        cprint("\t...")
        sessions_last_participant = ", ".join(
            s_id for s_id in sessions_per_participant[-1]
        )
        cprint(
            f"\t{Fore.RED}{unique_participants[-1]} | {sessions_last_participant}{Fore.RESET}"
        )


def print_crash_files_and_exit(log_file, working_directory):
    """Print command(s) to type in order to extract details after a Nipype RuntimeError and exit with an exception."""
    from colorama import Fore

    from .exceptions import ClinicaException
    from .filemanip import extract_crash_files_from_log_file
    from .stream import cprint

    cprint(
        f"{Fore.YELLOW}\nError details can be found by opening the crash file(s) "
        f"with the following command(s):{Fore.RESET}"
    )

    crash_files = extract_crash_files_from_log_file(log_file)
    for file in crash_files:
        cprint(f"{Fore.YELLOW}- nipypecli crash {file}{Fore.RESET}")

    cprint(
        f"{Fore.YELLOW}\n"
        f"If your pipeline crashed due to lack of space of network issues, "
        f"re-run the pipeline with the working directory (-wd {working_directory}).\n"
        f"Known issues are displayed here: https://aramislab.paris.inria.fr/clinica/docs/public/latest/InteractingWithClinica/#known-issues\n"
        f"Otherwise, you can delete it.{Fore.RESET}"
    )
    # Force the display of "Documentation can be found..."
    raise ClinicaException("")


def print_groups_in_caps_directory(caps_directory):
    """Print group IDs based on `caps_directory`/groups folder."""
    from .group import extract_group_ids
    from .stream import cprint

    group_ids = extract_group_ids(caps_directory)
    if group_ids == [""]:
        cprint("No group was found in CAPS directory")
    else:
        found_groups = ", ".join(g_id.replace("group-", "") for g_id in group_ids)
        cprint(f"Groups that exist in your CAPS directory are {found_groups}.")
