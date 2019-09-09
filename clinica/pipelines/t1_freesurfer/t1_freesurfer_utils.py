# coding: utf8


def init_input_node(t1w):
    """Extract "sub-<participant_id>_ses-<session_label>" from input node and print begin message."""
    from clinica.utils.stream import cprint
    import datetime
    from colorama import Fore
    from clinica.utils.io import get_subject_id

    image_id = get_subject_id(t1w)

    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s Running pipeline for %s' %
           (Fore.BLUE, now, Fore.RESET, image_id.replace('_', '|')))

    return image_id, t1w


def check_flags(in_t1w, recon_all_args):
    """Check `recon_all_args` given by the user and add -cw256 if the FOV of `in_t1w` is greater than 256."""
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


def create_subjects_dir(output_dir, image_id):
    """Create and return `subjects_dir` where `subjects_dir` = `output_dir`/`image_id`."""
    import os
    import errno
    subjects_dir = os.path.join(output_dir, image_id)
    try:
        os.makedirs(subjects_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:  # EEXIST: folder already exists
            raise e
    return subjects_dir


def write_tsv_files(subjects_dir, image_id):
    """Generate statistics TSV files in `subjects_dir`/regional_measures folder for `image_id`."""
    from clinica.utils.freesurfer import generate_regional_measures
    generate_regional_measures(subjects_dir, image_id)
    return image_id


def save_to_caps(source_dir, image_id, caps_dir, overwrite_caps=False):
    """Save `source_dir`/`image_id`/ to CAPS folder.

    This function copies outputs of `source_dir`/`image_id`/ to
    `caps_dir`/subjects/<participant_id>/<session_id>/t1_freesurfer_cross_sectional/
    where `image_id` = <participant_id>_<session_id>.

    Raise:
        FileNotFoundError: If symbolic links in `source_dir`/`image_id` folder are not removed
        IOError: If the `source_dir`/`image_id` folder does not contain FreeSurfer segmentation.
    """
    import os
    import errno
    import subprocess
    from colorama import Fore
    from clinica.utils.stream import cprint
    import shutil
    import datetime

    participant_id = image_id.split('_')[0]
    session_id = image_id.split('_')[1]

    destination_dir = os.path.join(
        os.path.expanduser(caps_dir),
        'subjects',
        participant_id,
        session_id,
        't1',
        'freesurfer_cross_sectional'
    )

    log_file = os.path.join(os.path.expanduser(source_dir), image_id, image_id, 'scripts', 'recon-all.log')
    if os.path.isfile(log_file):
        try:
            os.unlink(os.path.join(os.path.expanduser(source_dir), image_id, 'fsaverage'))
            os.unlink(os.path.join(os.path.expanduser(source_dir), image_id, 'lh.EC_average'))
            os.unlink(os.path.join(os.path.expanduser(source_dir), image_id, 'rh.EC_average'))
        except FileNotFoundError as e:
            if e.errno != errno.ENOENT:
                raise e
        last_line = subprocess.check_output(['tail', '-1', log_file])
        if b'finished without error' in last_line:
            if os.path.isfile(os.path.join(destination_dir, image_id, 'scripts', 'recon-all.log')):
                if overwrite_caps:
                    raise NotImplementedError('Overwritten of CAPS folder in t1-freesurfer pipeline not implemented')
                else:
                    now = datetime.datetime.now().strftime('%H:%M:%S')
                    cprint('%s[%s] Previous run of FreeSurfer was found in CAPS folder for %s. The copy will be skipped.%s' %
                           (Fore.YELLOW, now, image_id.replace('_', '|'), Fore.RESET))
            else:
                shutil.copytree(os.path.join(source_dir, image_id), destination_dir, symlinks=True)
                now = datetime.datetime.now().strftime('%H:%M:%S')
                cprint('%s[%s]%s %s has completed.' %
                       (Fore.GREEN, now, Fore.RESET, image_id.replace('_', '|')))

        else:
            cprint("%sFreeSurfer segmentation for %s will not be saved: errors were found.%s" %
                   (Fore.RED, image_id, Fore.RESET))
    else:
        raise IOError('Unexpected error: the log file should be present (%s)' % log_file)
    return image_id
