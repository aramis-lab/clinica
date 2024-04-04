from pathlib import Path

__all__ = [
    "extract_participant_long_ids_from_filename",
    "init_input_node",
    "run_recon_all_base",
    "save_to_caps",
]


def extract_participant_long_ids_from_filename(
    caps_files: list[str],
) -> tuple[list[str], list[str]]:
    """TODO: Find a way to merge with utils/filemanip.py::extract_subjects_sessions_from_filename into one util."""
    import re

    caps_files = [
        re.search(r"(sub-[a-zA-Z0-9]+)_(long-[a-zA-Z0-9]+)", file).group()
        for file in caps_files
    ]
    split = [image_id.split("_") for image_id in caps_files]
    part_ids = [p_id[0] for p_id in split]
    long_ids = [l_id[1] for l_id in split]
    return part_ids, long_ids


def init_input_node(
    caps_dir: Path,
    participant_id: str,
    list_session_ids: list[str],
    output_dir: Path,
) -> tuple[str, Path, str]:
    """Initialize the pipeline.

    This function will create folders and symbolic links in SUBJECTS_DIR env variable for upcoming run of recon-all.

    Note (@alexis-g-icm):
        There currently (as of 22 Feb 2019) is a bug in FreeSurfer recon-all -base, which in some cases (e.g., only one
        time point), will crash as it's trying to write lines too long for the shell to handle. This is caused by
        the path to FreeSurfer SUBJECT_DIR being too long itself.

    The current function works around this issue by checking if there only is one session associated to a subject, and
    in that case, putting the SUBJECT_DIR inside the system temporary folder so that its path is as short as possible.
    """
    from tempfile import mkdtemp

    from clinica.compat import errno
    from clinica.utils.longitudinal import get_long_id
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_begin_image

    long_id = get_long_id(list_session_ids)
    image_id = f"{participant_id}_{long_id}"
    if len(list_session_ids) == 1:
        # Special case: When only one time point is used, 'recon-all -base' can failed
        # if the $SUBJECTS_DIR is too long ('Word too long.' error).
        # To circumvent this issue, we create a sym link in $(TMP) so that $SUBJECTS_DIR is a short path
        subjects_dir = Path(mkdtemp())
        cprint(
            msg=(
                f"{image_id.replace('_', ' | ')} has only one time point. "
                f"Needs to create a $SUBJECTS_DIR folder in {subjects_dir}"
            ),
            lvl="warning",
        )
    else:
        subjects_dir = output_dir / image_id
    subjects_dir.mkdir(exist_ok=True)
    # Create symbolic links containing cross-sectional segmentation(s) in SUBJECTS_DIR so that recon-all can run
    for session_id in list_session_ids:
        cross_sectional_path = (
            caps_dir
            / "subjects"
            / participant_id
            / session_id
            / "t1"
            / "freesurfer_cross_sectional"
            / f"{participant_id}_{session_id}"
        )
        try:
            (subjects_dir / f"{participant_id}_{session_id}").symlink_to(
                cross_sectional_path
            )
        except FileExistsError as e:
            if e.errno != errno.EEXIST:  # EEXIST: folder already exists
                raise e
    flags = ""
    for session_id in list_session_ids:
        flags += f" -tp {participant_id}_{session_id}"
    print_begin_image(image_id)

    return image_id, subjects_dir, flags


def run_recon_all_base(subjects_dir: str, subject_id: str, flags, directive):
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

    recon_all_base_command = "recon-all -base {0} {1} -sd {2} {3}".format(
        subject_id, flags, subjects_dir, directive
    )
    subprocess_run_recon_all_base = subprocess.run(
        recon_all_base_command,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if subprocess_run_recon_all_base.returncode != 0:
        raise ValueError("recon-all -base failed, returned non-zero code")

    return subject_id


def save_to_caps(
    source_dir: Path,
    image_id: str,
    list_session_ids: list[str],
    caps_dir: Path,
    overwrite_caps: bool = False,
) -> str:
    """Save `source_dir`/`image_id`/ to CAPS folder.

    This function copies FreeSurfer segmentation and regional_measures folder of `source_dir`/`image_id`/ to
    `caps_dir`/subjects/<participant_id>/<long_id>/freesurfer_unbiased_template/
    where `image_id` = <participant_id>_<long_id>.

    The `source_dir`/`image_id`/ folder should contain the following elements:
        - fsaverage, lh.EC_average and rh.EC_average symbolic links automatically generated by recon-all
        - symbolic links to cross-sectional segmentation(s) and unbiased template needed for recon-all
        - <participant_id>_<long_id>/ folder containing the FreeSurfer segmentation

    Notes:
        We do not need to check the line "finished without error" in scripts/recon-all.log.
        If an error occurs, it will be detected by Nipype and the next nodes (i.e. save_to_caps will not be called).
    """
    import shutil

    from clinica.utils.longitudinal import save_long_id
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_end_image

    participant_id = image_id.split("_")[0]
    long_id = image_id.split("_")[1]
    sessions_tsv_path = caps_dir.expanduser() / "subjects" / participant_id / long_id
    destination_dir = sessions_tsv_path / "freesurfer_unbiased_template"
    if not (sessions_tsv_path / f"{long_id}_sessions.tsv").is_file():
        save_long_id(list_session_ids, sessions_tsv_path, f"{long_id}_sessions.tsv")
    source_file = (
        source_dir.expanduser() / image_id / image_id / "mri" / "aparc+aseg.mgz"
    )
    destination_file = destination_dir / image_id / "mri" / "aparc+aseg.mgz"
    if source_file.is_file():
        if destination_file.is_file():
            if overwrite_caps:
                shutil.rmtree(destination_dir)
                shutil.copytree(
                    src=source_dir / image_id / image_id,
                    dst=destination_dir / image_id,
                    symlinks=True,
                )
        else:
            shutil.copytree(
                src=source_dir / image_id / image_id,
                dst=destination_dir / image_id,
                symlinks=True,
            )
        print_end_image(image_id)
    else:
        cprint(
            msg=(
                f"{image_id.replace('_', ' | ')}  does not contain "
                "mri/aseg+aparc.mgz file. Copy will be skipped."
            ),
            lvl="warning",
        )
    return image_id
