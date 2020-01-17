# coding: utf8


def init_input_node(caps_dir, participant_id, list_session_ids, output_dir):
    """Initialize the pipeline.

    This function will create folders and symbolic links in SUBJECTS_DIR env variable for upcoming run of recon-all.

    Note(AG):
        There currently (as of 22 Feb 2019) is a bug in FreeSurfer recon-all -base, which in some cases (e.g., only one
        time point), will crash as it's trying to write lines too long for the shell to handle. This is caused by
        the path to FreeSurfer SUBJECT_DIR being too long itself.

    The current function works around this issue by checking if there only is one session associated to a subject, and
    in that case, putting the SUBJECT_DIR inside the system temporary folder so that its path is as short as possible.
    """
    import os
    import errno
    from clinica.utils.longitudinal import get_long_id
    from clinica.utils.ux import print_begin_image

    # Extract <image_id>
    long_id = get_long_id(list_session_ids)
    image_id = participant_id + '_' + long_id

    # Create SUBJECTS_DIR for recon-all (otherwise, the command won't run)
    subjects_dir = os.path.join(output_dir, image_id)
    try:
        os.makedirs(subjects_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:  # EEXIST: folder already exists
            raise e

    # Create symbolic links containing cross-sectional segmentation(s) in SUBJECTS_DIR so that recon-all can run
    for session_id in list_session_ids:
        cross_sectional_path = os.path.join(
            caps_dir,
            'subjects',
            participant_id,
            session_id,
            't1',
            'freesurfer_cross_sectional',
            participant_id + '_' + session_id
        )
        try:
            os.symlink(cross_sectional_path, os.path.join(subjects_dir, participant_id + '_' + session_id))
        except FileExistsError as e:
            if e.errno != errno.EEXIST:  # EEXIST: folder already exists
                raise e

    # Prepare arguments for recon-all.
    flags = ""
    for session_id in list_session_ids:
        flags += " -tp " + participant_id + "_" + session_id

    print_begin_image(image_id)

    return image_id, subjects_dir, flags


def run_recon_all_base(subjects_dir,
                       subject_id,
                       flags,
                       directive):
    """Run recon-all to create an unbiased template from all time points.

    Note:
    Unbiased template creation with FreeSurfer expects arguments to follow this syntax:
    recon-all -base <template_id> -tp <tp1_id> -tp <tp2_id> ... -all

    Currently, Nipype does not provide interface for "recon-all -base" case. As a result, ReconAll interface should be
    modified to handle this situation. In the meantime, the arguments of this function follows ReconAll.inputs names
    namely:
        - "-base <template_id>" is likely to be fed to ReconAll.inputs.subject_id or new input (e.g. template_id?)
        - "-tp <tp1_id> -tp <tp2_id> ..." is likely to  be fed to ReconAll.inputs.args
        - "-all" will be fed to ReconAll.inputs.directive

    See official documentation (https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalProcessing) for details.
    """
    import subprocess

    recon_all_base_command = 'recon-all -base {0} {1} -sd {2} {3}'.format(
        subject_id, flags, subjects_dir, directive)
    subprocess_run_recon_all_base = subprocess.run(
        recon_all_base_command,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL)
    if subprocess_run_recon_all_base.returncode != 0:
        raise ValueError('recon-all -base failed, returned non-zero code')

    return subject_id


def save_to_caps(source_dir, image_id, list_session_ids, caps_dir, overwrite_caps):
    """Save `subjects_dir` to CAPS folder.

    This function copies outputs of `subjects_dir`/ to
    `caps_dir`/subjects/<participant_id>/<long_id>/freesurfer_unbiased_template/
    The `subjects_dir`/ folder should contain the following elements:
        - fsaverage, lh.EC_average and rh.EC_average symbolic links
        - symbolic links to cross-sectional segmentations
        - `participant_id`_`long_id`/ folder containing the FreeSurfer segmentation
    Notes:
        We do not need to check the line "finished without error" in scripts/recon-all.log.
        If an error occurs, it will be detected by Nipype and the next nodes (i.e. save_to_caps will not be called).
    Raise:
        FileNotFoundError: If symbolic links in `source_dir`/`image_id` folder are not removed
    """
    import os
    import datetime
    import errno
    import shutil
    from colorama import Fore
    from clinica.utils.stream import cprint
    from clinica.utils.longitudinal import save_long_id
    from clinica.utils.ux import print_end_image

    participant_id = image_id.split('_')[0]
    long_id = image_id.split('_')[1]

    destination_dir = os.path.join(
        os.path.expanduser(caps_dir),
        'subjects',
        participant_id,
        long_id,
        'freesurfer_unbiased_template'
    )

    # Save <long_id>_sessions.tsv to retrieve sessions used to generate <long_id>
    sessions_tsv_path = os.path.join(
        os.path.expanduser(caps_dir),
        'subjects',
        participant_id,
        long_id
    )
    if not os.path.isfile(os.path.join(sessions_tsv_path, long_id + '_sessions.tsv')):
        save_long_id(list_session_ids, sessions_tsv_path, long_id + '_sessions.tsv')

    # Save FreeSurfer segmentation
    representative_file = os.path.join(image_id, 'mri', 'aparc+aseg.mgz')
    representative_source_file = os.path.join(os.path.expanduser(source_dir), image_id, representative_file)
    representative_destination_file = os.path.join(destination_dir, representative_file)
    if os.path.isfile(representative_source_file):
        # Remove symbolic links before the copy
        try:
            os.unlink(os.path.join(os.path.expanduser(source_dir), image_id, 'fsaverage'))
            os.unlink(os.path.join(os.path.expanduser(source_dir), image_id, 'lh.EC_average'))
            os.unlink(os.path.join(os.path.expanduser(source_dir), image_id, 'rh.EC_average'))
            for session_id in list_session_ids:
                os.unlink(os.path.join(os.path.expanduser(source_dir), image_id, participant_id + '_' + session_id))
        except FileNotFoundError as e:
            if e.errno != errno.ENOENT:
                raise e

        if os.path.isfile(representative_destination_file) and overwrite_caps:
            shutil.rmtree(destination_dir)
            shutil.copytree(os.path.join(source_dir, image_id), destination_dir, symlinks=True)
            print_end_image(image_id)

    else:
        now = datetime.datetime.now().strftime('%H:%M:%S')
        cprint('%s[%s] %s does not contain mri/aseg+aparc.mgz file. '
               'Copy will be skipped.%s' %
               (Fore.YELLOW, now, image_id.replace('_', ' | '), Fore.RESET))

    return image_id
