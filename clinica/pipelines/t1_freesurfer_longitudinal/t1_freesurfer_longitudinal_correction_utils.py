# coding: utf8


def init_input_node(caps_dir, participant_id, session_id, long_id, output_dir):
    """Initialize the pipeline."""
    import os
    import errno
    import datetime
    import platform
    from tempfile import mkdtemp
    from colorama import Fore
    from clinica.utils.longitudinal import read_sessions
    from clinica.utils.ux import print_begin_image
    from clinica.utils.stream import cprint

    # Extract <image_id>
    image_id = "{0}_{1}_{2}".format(participant_id, session_id, long_id)

    # Create SUBJECTS_DIR for recon-all (otherwise, the command won't run)
    if platform.system().lower().startswith("darwin"):
        # Special case: On macOS, 'recon-all -long' can failed if the $SUBJECTS_DIR is too long
        # To circumvent this issue, we create a sym link in $(TMP) so that $SUBJECTS_DIR is a short path
        subjects_dir = mkdtemp()
        now = datetime.datetime.now().strftime("%H:%M:%S")
        cprint(
            "%s[%s] Needs to create a $SUBJECTS_DIR folder in %s for %s (macOS case). %s"
            % (Fore.YELLOW, now, subjects_dir, image_id.replace("_", " | "), Fore.RESET)
        )
    else:
        subjects_dir = os.path.join(output_dir, image_id)

    try:
        os.makedirs(subjects_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:  # EEXIST: folder already exists
            raise e

    # Create symbolic link containing cross-sectional segmentation(s) in SUBJECTS_DIR so that recon-all can run
    for s_id in read_sessions(caps_dir, participant_id, long_id):
        cross_sectional_path = os.path.join(
            caps_dir,
            "subjects",
            participant_id,
            s_id,
            "t1",
            "freesurfer_cross_sectional",
            participant_id + "_" + s_id,
        )
        os.symlink(
            cross_sectional_path,
            os.path.join(subjects_dir, participant_id + "_" + s_id),
        )

    # Create symbolic links containing unbiased template in SUBJECTS_DIR so that recon-all can run
    template_path = os.path.join(
        caps_dir,
        "subjects",
        participant_id,
        long_id,
        "freesurfer_unbiased_template",
        participant_id + "_" + long_id,
    )
    os.symlink(
        template_path, os.path.join(subjects_dir, participant_id + "_" + long_id)
    )

    print_begin_image(image_id)

    return subjects_dir


def run_recon_all_long(subjects_dir, participant_id, session_id, long_id, directive):
    """Run recon-all to create a longitudinal correction of a time point.

    Note:
    Longitudinal correction with FreeSurfer expects arguments to follow this syntax:
    recon-all -long <tpN_id> <template_id> -all; e.g.: recon-all -long sub-CLNC01_ses-M00 sub-CLNC01_long-M00M18 -all

    Currently, Nipype does not provide interface for "recon-all -long" case. As a result, ReconAll interface should be
    modified to handle this situation. In the meantime, the arguments of this function follows ReconAll.inputs names
    namely:
        - "-long <tpN_id> <template_id>" is likely to be fed to ReconAll.inputs.args
        - "-all" will be fed to ReconAll.inputs.directive

    Folder containing the longitudinal correction has the following convention:
    <tpN_id>.long.<template_id>; e.g.: sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18
    which is automatically generated by FreeSurfer.

    This folder name is likely to be retrieved in ReconAll.outputs.subject_id.

    See official documentation (https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalProcessing) for details.
    """
    import subprocess

    # Prepare arguments for recon-all.
    flags = " -long {0}_{1} {0}_{2} ".format(participant_id, session_id, long_id)

    recon_all_long_command = "recon-all {0} -sd {1} {2}".format(
        flags, subjects_dir, directive
    )
    subprocess_run_recon_all_long = subprocess.run(
        recon_all_long_command,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if subprocess_run_recon_all_long.returncode != 0:
        raise ValueError("recon-all -long failed, returned non-zero code")

    subject_id = "{0}_{1}.long.{0}_{2}".format(participant_id, session_id, long_id)

    return subject_id


def write_tsv_files(subjects_dir, subject_id):
    """
    Generate statistics TSV files in `subjects_dir`/regional_measures folder for `subject_id`.

    Notes:
        We do not need to check the line "finished without error" in scripts/recon-all.log.
        If an error occurs, it will be detected by Nipype and the next nodes (i.e.
        write_tsv_files will not be called).
    """
    import os
    import datetime
    from colorama import Fore
    from clinica.utils.stream import cprint
    from clinica.utils.freesurfer import (
        generate_regional_measures,
        extract_image_id_from_longitudinal_segmentation,
    )

    image_id = extract_image_id_from_longitudinal_segmentation(subject_id)
    str_image_id = (
        image_id.participant_id + "_" + image_id.session_id + "_" + image_id.long_id
    )
    if os.path.isfile(os.path.join(subjects_dir, subject_id, "mri", "aparc+aseg.mgz")):
        generate_regional_measures(subjects_dir, subject_id)
    else:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        cprint(
            "%s[%s] %s does not contain mri/aseg+aparc.mgz file. "
            "Creation of regional_measures/ folder will be skipped.%s"
            % (Fore.YELLOW, now, str_image_id.replace("_", " | "), Fore.RESET)
        )
    return subject_id


def move_subjects_dir_to_source_dir(subjects_dir, source_dir, subject_id):
    """
    Move content of `subjects_dir`/`subject_id` to `source_dir`.

    This function will move content of `subject_id` if recon-all has run
    in $(TMP). This happens when FreeSurfer is run on macOS. Content of
    $(TMP)/`subject_id` is copied to `source_dir` before the deletion of $(TMP).

    Args:
        subjects_dir: $(TMP), if segmentation was performed on macOS,
            <base_dir>/<Pipeline.Name>/ReconAll/`subject_id` otherwise
        source_dir: <base_dir>/<Pipeline.Name>/ReconAll folder
        subject_id: Subject ID (e.g. "sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18")

    Returns:
        subject_id for node connection with Nipype
    """
    import os
    import shutil
    import datetime
    from colorama import Fore
    from clinica.utils.freesurfer import extract_image_id_from_longitudinal_segmentation
    from clinica.utils.stream import cprint

    image_id = extract_image_id_from_longitudinal_segmentation(subject_id)
    participant_id = image_id.participant_id
    session_id = image_id.session_id
    long_id = image_id.long_id
    str_image_id = (
        image_id.participant_id + "_" + image_id.session_id + "_" + image_id.long_id
    )

    if source_dir not in subjects_dir:
        shutil.copytree(
            src=os.path.join(subjects_dir, subject_id),
            dst=os.path.join(source_dir, str_image_id, subject_id),
            symlinks=True,
        )
        shutil.copytree(
            src=os.path.join(subjects_dir, "regional_measures"),
            dst=os.path.join(source_dir, str_image_id, "regional_measures"),
            symlinks=True,
        )
        shutil.rmtree(subjects_dir)
        now = datetime.datetime.now().strftime("%H:%M:%S")
        cprint(
            f"{Fore.YELLOW}[{now}] Segmentation of {subject_id.replace('_', ' | ')} "
            f"has moved to working directory and $SUBJECTS_DIR folder ({subjects_dir}) was deleted{Fore.RESET}"
        )

    return subject_id


def save_to_caps(source_dir, subject_id, caps_dir, overwrite_caps=False):
    """Save `source_dir`/`subject_id`/ to CAPS folder.

    This function copies FreeSurfer segmentation and regional_measures folder of `source_dir`/`image_id`/ to
    `caps_dir`/subjects/<participant_id>/<session_id>/<long_id>/t1_freesurfer_longitudinal/
    where `image_id` = <participant_id>_<session_id>_<long_id>.

    The `source_dir`/`image_id`/ folder should contain the following elements:
        - fsaverage, lh.EC_average and rh.EC_average symbolic links automatically generated by recon-all
        - symbolic links to cross-sectional segmentation(s) and unbiased template needed for recon-all
        - <participant_id>_<session_id>.long.<participant_id>_<long_id>/ folder containing the FreeSurfer segmentation
        - regional_measures/ folder containing TSV files

    Notes:
        We do not need to check the line "finished without error" in scripts/recon-all.log.
        If an error occurs, it will be detected by Nipype and the next nodes (i.e. save_to_caps will not be called).
    """
    import os
    import datetime
    import shutil
    from colorama import Fore
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_end_image
    from clinica.utils.freesurfer import extract_image_id_from_longitudinal_segmentation

    image_id = extract_image_id_from_longitudinal_segmentation(subject_id)
    participant_id = image_id.participant_id
    session_id = image_id.session_id
    long_id = image_id.long_id
    str_image_id = (
        image_id.participant_id + "_" + image_id.session_id + "_" + image_id.long_id
    )

    destination_dir = os.path.join(
        os.path.expanduser(caps_dir),
        "subjects",
        participant_id,
        session_id,
        "t1",
        long_id,
        "freesurfer_longitudinal",
    )

    # Save FreeSurfer segmentation
    representative_file = os.path.join(subject_id, "mri", "aparc+aseg.mgz")
    representative_source_file = os.path.join(
        os.path.expanduser(source_dir), str_image_id, representative_file
    )
    representative_destination_file = os.path.join(destination_dir, representative_file)
    if os.path.isfile(representative_source_file):
        if os.path.isfile(representative_destination_file):
            if overwrite_caps:
                shutil.rmtree(destination_dir)
                shutil.copytree(
                    src=os.path.join(source_dir, subject_id, subject_id),
                    dst=os.path.join(destination_dir, subject_id),
                    symlinks=True,
                )
                shutil.copytree(
                    src=os.path.join(source_dir, subject_id, "regional_measures"),
                    dst=os.path.join(destination_dir, "regional_measures"),
                    symlinks=True,
                )
        else:
            shutil.copytree(
                src=os.path.join(source_dir, str_image_id, subject_id),
                dst=os.path.join(destination_dir, subject_id),
                symlinks=True,
            )
            shutil.copytree(
                src=os.path.join(source_dir, str_image_id, "regional_measures"),
                dst=os.path.join(destination_dir, "regional_measures"),
                symlinks=True,
            )
        print_end_image(str_image_id)
    else:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        cprint(
            f"{Fore.YELLOW}[{now}] {subject_id.replace('_', ' | ')}  does not contain "
            f"mri/aseg+aparc.mgz file. Copy will be skipped.{Fore.RESET}"
        )

    return str_image_id


def get_processed_images(caps_directory, part_ids, sess_ids, long_ids):
    """
    Extract image IDs (e.g. ["sub-CLNC01_ses-M00_long-M00M18", "sub-CLNC01_ses-M18_long-M00M18"]) of outputs
    already processed by T1FreeSurferLongitudinalCorrection pipeline.
    """
    import os

    image_ids = []
    if os.path.isdir(caps_directory):
        for (participant_id, session_id, long_id) in zip(part_ids, sess_ids, long_ids):
            output_file = os.path.join(
                os.path.expanduser(caps_directory),
                "subjects",
                participant_id,
                session_id,
                "t1",
                long_id,
                "freesurfer_longitudinal",
                participant_id
                + "_"
                + session_id
                + ".long."
                + participant_id
                + "_"
                + long_id,
                "mri",
                "aparc+aseg.mgz",
            )
            if os.path.isfile(output_file):
                image_ids.append(participant_id + "_" + session_id + "_" + long_id)
    return image_ids
