# coding: utf8


def init_input_node(t1w):
    """
    Extract "sub-<participant_id>_ses-<session_label>" from input node
    and print begin message.
    """
    from clinica.utils.stream import cprint
    import datetime
    from colorama import Fore
    from clinica.utils.io import get_subject_id

    subject_id = get_subject_id(t1w)

    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s Running pipeline for %s...' %
           (Fore.BLUE, now, Fore.RESET, subject_id.replace('_', '|')))

    return subject_id, t1w


def check_flags(in_t1w, recon_all_args):
    """
    Check <recon_all_args> given by the user and
    add -cw256 if the FOV of <in_t1w> is greater than 256.
    """
    import nibabel as nib
    # from clinica.utils.stream import cprint

    f = nib.load(in_t1w)
    voxel_size = f.header.get_zooms()
    t1_size = f.header.get_data_shape()
    if (voxel_size[0] * t1_size[0] > 256) or \
            (voxel_size[1] * t1_size[1] > 256) or \
            (voxel_size[2] * t1_size[2] > 256):
        # cprint("Setting MRI Convert to crop images to 256 FOV for %s file." % in_t1w)
        optional_flag = '-cw256'
    else:
        # cprint("No need to add -cw256 flag for %s file." % in_t1w)
        optional_flag = ''
    flags = "{0} ".format(recon_all_args) + optional_flag

    return flags


def create_subjects_dir(output_dir, subject_id):
    """
    Create and return <subjects_dir>
    where <subjects_dir> = <output_dir>/<subject_id>.
    """
    import os
    subjects_dir = os.path.join(output_dir, subject_id)
    try:
        os.makedirs(subjects_dir)
    except OSError as e:
        if e.errno != 17:  # Errno 17: folder already exists
            raise e
    return subjects_dir


def write_tsv_files(subjects_dir, subject_id):
    """
    Generate statistics TSV files in <subjects_dir>/regional_measures folder
    for <subject_id>.
    """
    from clinica.utils.freesurfer import generate_regional_measures
    generate_regional_measures(subjects_dir, subject_id)
    return subject_id


def print_end_pipeline(subject_id):
    """
    Display end message for <subject_id> (e.g. sub-CLNC01_ses-M0).
    """
    from clinica.utils.stream import cprint
    import datetime
    from colorama import Fore

    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s %s has completed.' %
           (Fore.GREEN, now, Fore.RESET, subject_id.replace('_', '|')))


def save_to_caps(source_dir, subject_id, caps_dir, overwrite_caps=False):
    raise NotImplementedError('t1-freesurfer::save_to_caps - Not Implemented Yet')
    return subject_id
