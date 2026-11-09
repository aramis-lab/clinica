from pathlib import Path
from typing import Optional

__all__ = [
    "extract_subject_session_longitudinal_ids_from_filename",
    "get_processed_images",
    "init_input_node",
    "read_part_sess_long_ids_from_tsv",
    "run_recon_all_long",
    "save_to_caps",
    "write_tsv_files",
]


def extract_subject_session_longitudinal_ids_from_filename(
    bids_or_caps_files: list[str],
) -> tuple[list[str], list[str], list[str]]:
    """Extract participant/session/longitudinal IDs from filename.

    Example:
        ['sub-CLNC01', 'sub-CLNC01'], ['ses-M000', 'ses-M018'], ['long-M000M018', 'long-M000M018'])

    TODO: Find a way to merge with utils/filemanip.py::extract_subjects_sessions_from_filename into one util
    """
    import re

    id_bids_or_caps_files = [
        re.search(
            r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)_(long-[a-zA-Z0-9]+)", file
        ).group()
        for file in bids_or_caps_files
    ]
    split = [image_id.split("_") for image_id in id_bids_or_caps_files]
    part_ids = [p_id[0] for p_id in split]
    sess_ids = [s_id[1] for s_id in split]
    long_ids = [l_id[2] for l_id in split]
    return part_ids, sess_ids, long_ids


def read_part_sess_long_ids_from_tsv(
    tsv_file: Path,
) -> tuple[list[str], list[str], list[str]]:
    """Extract participant, session and longitudinal from TSV file.

    TODO: Find a way to merge with utils/filemanip.py::read_participant_tsv into one util
    """
    import pandas as pd

    from clinica.utils.exceptions import ClinicaException

    if not tsv_file.is_file():
        raise ClinicaException(
            "The TSV file you gave is not a file.\n"
            "Error explanations:\n"
            f" - Clinica expected the following path to be a file: {tsv_file}\n"
            " - If you gave relative path, did you run Clinica on the good folder?"
        )
    df = pd.read_csv(tsv_file, sep="\t")
    if not {"participant_id", "session_id", "long_id"}.issubset(df.columns):
        raise ClinicaException(
            f"The TSV file ({tsv_file}) does not contain all required columns."
        )

    return (
        [sub.strip(" ") for sub in list(df.participant_id)],
        [ses.strip(" ") for ses in list(df.session_id)],
        [lng.strip(" ") for lng in list(df.long_id)],
    )


def init_input_node(
    caps_dir: Path, participant_id: str, session_id: str, long_id: str, output_dir: Path
) -> Path:
    """Initialize the pipeline."""
    import platform
    from tempfile import mkdtemp

    from clinica.utils.longitudinal import read_sessions
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_begin_image

    image_id = "_".join((participant_id, session_id, long_id))

    # Create SUBJECTS_DIR for recon-all (otherwise, the command won't run)
    if platform.system().lower().startswith("darwin"):
        # Special case: On macOS, 'recon-all -long' can fail if the $SUBJECTS_DIR is too long
        # To circumvent this issue, we create a sym link in $(TMP) so that $SUBJECTS_DIR is a short path
        subjects_dir = Path(mkdtemp())
        cprint(
            msg=(
                f"Needs to create a $SUBJECTS_DIR folder "
                f"in {subjects_dir} for {image_id.replace('_', ' | ')} (macOS case)."
            ),
            lvl="warning",
        )
    else:
        subjects_dir = output_dir / image_id
    subjects_dir.mkdir(parents=True, exist_ok=True)
    # Create symbolic link containing cross-sectional segmentation(s) in SUBJECTS_DIR so that recon-all can run
    for s_id in read_sessions(caps_dir, participant_id, long_id):
        cross_sectional_path = (
            caps_dir
            / "subjects"
            / participant_id
            / s_id
            / "t1"
            / "freesurfer_cross_sectional"
            / f"{participant_id}_{s_id}"
        )
        (subjects_dir / f"{participant_id}_{s_id}").symlink_to(cross_sectional_path)
    # Create symbolic links containing unbiased template in SUBJECTS_DIR so that recon-all can run
    template_path = (
        caps_dir
        / "subjects"
        / participant_id
        / long_id
        / "freesurfer_unbiased_template"
        / f"{participant_id}_{long_id}"
    )
    (subjects_dir / f"{participant_id}_{long_id}").symlink_to(template_path)

    print_begin_image(image_id)

    return subjects_dir


def run_recon_all_long(
    subjects_dir: str, participant_id: str, session_id: str, long_id: str, directive
) -> str:
    """Run recon-all to create a longitudinal correction of a time point.

    Notes
    -----
    Longitudinal correction with FreeSurfer expects arguments to follow this syntax:
    recon-all -long <tpN_id> <template_id> -all; e.g.: recon-all -long sub-CLNC01_ses-M000 sub-CLNC01_long-M000M018 -all

    Currently, Nipype does not provide interface for "recon-all -long" case. As a result, ReconAll interface should be
    modified to handle this situation. In the meantime, the arguments of this function follows ReconAll.inputs names
    namely:
        - "-long <tpN_id> <template_id>" is likely to be fed to ReconAll.inputs.args
        - "-all" will be fed to ReconAll.inputs.directive

    Folder containing the longitudinal correction has the following convention:
    <tpN_id>.long.<template_id>; e.g.: sub-CLNC01_ses-M000.long.sub-CLNC01_long-M000M018
    which is automatically generated by FreeSurfer.

    This folder name is likely to be retrieved in ReconAll.outputs.subject_id.

    See official documentation (https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalProcessing) for details.
    """
    import subprocess

    flags = f" -long {participant_id}_{session_id} {participant_id}_{long_id} "
    recon_all_long_command = f"recon-all {flags} -sd {subjects_dir} {directive}"
    result = subprocess.run(
        recon_all_long_command,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if result.returncode != 0:
        raise ValueError("recon-all -long failed, returned non-zero code")

    return f"{participant_id}_{session_id}.long.{participant_id}_{long_id}"


def write_tsv_files(subjects_dir: Path, subject_id: str) -> str:
    """Generate statistics TSV files in `subjects_dir`/regional_measures folder for `subject_id`.

    Notes
    -----
    We do not need to check the line "finished without error" in scripts/recon-all.log.
    If an error occurs, it will be detected by Nipype and the next nodes (i.e.
    write_tsv_files will not be called).
    """
    from clinica.pipelines.anatomical.freesurfer.utils import (
        extract_image_id_from_freesurfer_id,
        generate_regional_measures,
    )
    from clinica.utils.stream import cprint

    if (subjects_dir / subject_id / "mri" / "aparc+aseg.mgz").is_file():
        generate_regional_measures(subjects_dir, subject_id)
    else:
        image_id = " | ".join(extract_image_id_from_freesurfer_id(subject_id))
        cprint(
            msg=(
                f"{image_id} does not contain mri/aseg+aparc.mgz file. "
                "Creation of regional_measures/ folder will be skipped."
            ),
            lvl="warning",
        )
    return subject_id


def save_to_caps(
    source_dir: Path, subject_id: str, caps_dir: Path, overwrite_caps: bool = False
) -> str:
    """Save `source_dir`/`subject_id`/ to CAPS folder.

    This function copies FreeSurfer segmentation and regional_measures folder of `source_dir`/`image_id`/ to
    `caps_dir`/subjects/<participant_id>/<session_id>/<long_id>/t1_freesurfer_longitudinal/
    where `image_id` = <participant_id>_<session_id>_<long_id>.

    The `source_dir`/`image_id`/ folder should contain the following elements:
        - fsaverage, lh.EC_average and rh.EC_average symbolic links automatically generated by recon-all
        - symbolic links to cross-sectional segmentation(s) and unbiased template needed for recon-all
        - <participant_id>_<session_id>.long.<participant_id>_<long_id>/ folder containing the FreeSurfer segmentation
        - regional_measures/ folder containing TSV files

    Notes
    -----
    We do not need to check the line "finished without error" in scripts/recon-all.log.
    If an error occurs, it will be detected by Nipype and the next nodes (i.e. save_to_caps will not be called).
    """
    import shutil

    from clinica.pipelines.anatomical.freesurfer.utils import (
        extract_image_id_from_freesurfer_id,
    )
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_end_image

    image_id = extract_image_id_from_freesurfer_id(subject_id)
    str_image_id = "_".join(image_id)
    destination_dir = (
        caps_dir.expanduser()
        / "subjects"
        / image_id.participant_id
        / image_id.session_id
        / "t1"
        / image_id.long_id
        / "freesurfer_longitudinal"
    )
    source_file = (
        source_dir.expanduser() / str_image_id / subject_id / "mri" / "aparc+aseg.mgz"
    )
    destination_file = destination_dir / subject_id / "mri" / "aparc+aseg.mgz"
    if source_file.is_file():
        if destination_file.is_file():
            if overwrite_caps:
                shutil.rmtree(destination_dir)
                shutil.copytree(
                    src=source_dir / subject_id / subject_id,
                    dst=destination_dir / subject_id,
                    symlinks=True,
                )
                shutil.copytree(
                    src=source_dir / subject_id / "regional_measures",
                    dst=destination_dir / "regional_measures",
                    symlinks=True,
                )
        else:
            shutil.copytree(
                src=source_dir / str_image_id / subject_id,
                dst=destination_dir / subject_id,
                symlinks=True,
            )
            shutil.copytree(
                src=source_dir / str_image_id / "regional_measures",
                dst=destination_dir / "regional_measures",
                symlinks=True,
            )
        print_end_image(str_image_id)
    else:
        cprint(
            msg=(
                f"{subject_id.replace('_', ' | ')}  does not "
                "contain mri/aseg+aparc.mgz file. Copy will be skipped."
            ),
            lvl="warning",
        )

    return str_image_id


def get_processed_images(
    caps_directory: Path,
    part_ids: list,
    sess_ids: list,
    long_ids: list,
    atlas: Optional[str] = None,
) -> list[str]:
    """
    Extract image IDs (e.g. ["sub-CLNC01_ses-M000_long-M000+M018", "sub-CLNC01_ses-M018_long-M000+M018"]) of outputs
    already processed by T1FreeSurferLongitudinalCorrection pipeline.
    """
    image_ids = []
    if caps_directory.is_dir():
        for participant_id, session_id, long_id in zip(part_ids, sess_ids, long_ids):
            output_file_pre = (
                caps_directory.expanduser()
                / "subjects"
                / participant_id
                / session_id
                / "t1"
                / long_id
                / "freesurfer_longitudinal"
                / f"{participant_id}_{session_id}.long.{participant_id}_{long_id}"
            )
            if atlas:
                output_file = output_file_pre / "stat" / "rh." / atlas / ".stats"
            else:
                output_file = output_file_pre / "mri" / "aparc+aseg.mgz"
            if output_file.is_file():
                image_ids.append(f"{participant_id}_{session_id}_{long_id}")
    return image_ids
