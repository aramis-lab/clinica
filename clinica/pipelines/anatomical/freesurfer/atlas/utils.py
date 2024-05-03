from pathlib import Path

__all__ = ["compute_atlas", "write_tsv_files"]


def compute_atlas(
    caps_directory: Path,
    to_process_with_atlases: tuple,
    path_to_atlas: Path,
) -> tuple[Path, str, str]:
    import subprocess
    from pathlib import Path
    from tempfile import mkdtemp

    from clinica.compat import errno
    from clinica.utils.stream import cprint

    subject_dir = Path()
    image_id = ""
    atlas = ""

    for path in path_to_atlas.rglob(f"*{to_process_with_atlases[0]}_6p0.gcs"):
        hemisphere = Path(path.stem).stem
        atlas_name = Path(path.stem).suffix[1:].split("_")[0]
        if "long" not in to_process_with_atlases[1]:
            sub, ses = to_process_with_atlases[1].split("_")
            subject_dir = (
                caps_directory
                / "subjects"
                / sub
                / ses
                / "t1"
                / "freesurfer_cross_sectional"
            )
            substitute_dir = subject_dir
            freesurfer_cross = subject_dir / to_process_with_atlases[1]
            output_path_annot = (
                freesurfer_cross / "label" / f"{hemisphere}.{atlas_name}.annot"
            )
            sphere_reg = freesurfer_cross / "surf" / f"{hemisphere}.sphere.reg"
            output_path_stats = (
                freesurfer_cross / "stats" / f"{hemisphere}.{atlas_name}.stats"
            )
            subject_id = to_process_with_atlases[1]
        else:
            sub, ses, long = to_process_with_atlases[1].split("_")
            substitute_dir = Path(mkdtemp())
            substitute_dir.mkdir(exist_ok=True)
            subject_id = f"{sub}_{ses}.long.{sub}_{long}"
            subject_dir = (
                caps_directory
                / "subjects"
                / sub
                / ses
                / "t1"
                / long
                / "freesurfer_longitudinal"
            )
            freesurfer_cross = subject_dir / subject_id
            try:
                (substitute_dir / subjid).symlink_to(freesurfer_cross)
            except FileExistsError as e:
                if e.errno != errno.EEXIST:  # EEXIST: folder already exists
                    raise e
            substitute_dir_id = substitute_dir / subject_id
            output_path_annot = (
                substitute_dir_id / "label" / f"{hemisphere}.{atlas_name}.annot"
            )
            sphere_reg = substitute_dir_id / "surf" / f"{hemisphere}.sphere.reg"
            output_path_stats = (
                substitute_dir_id / "stats" / f"{hemisphere}.{atlas_name}.stats"
            )
        if not sphere_reg.is_file():
            cprint(
                f"The {hemisphere}.sphere.reg file appears to be missing. "
                f"The data for {to_process_with_atlases[1]} will not be "
                f"processed with {to_process_with_atlases[0]}",
                lvl="warning",
            )
        mris_ca_label_command = (
            f"mris_ca_label -sdir {substitute_dir} {subject_id} {hemisphere} "
            f"{freesurfer_cross}/surf/{hemisphere}.sphere.reg {path} {output_path_annot}"
        )
        subprocess.run(mris_ca_label_command, shell=True, capture_output=True)
        mris_anatomical_stats_command = (
            f"export SUBJECTS_DIR={substitute_dir}; "
            f"mris_anatomical_stats -a {output_path_annot} -f {output_path_stats} -b {subject_id} {hemisphere}"
        )
        subprocess.run(mris_anatomical_stats_command, shell=True, capture_output=True)
        image_id, atlas = (
            to_process_with_atlases[1],
            to_process_with_atlases[0],
        )
    return subject_dir, image_id, atlas


def write_tsv_files(subject_dir: Path, image_id: str, atlas: str) -> str:
    """
    Generate statistics TSV files in `subjects_dir`/regional_measures folder for `image_id`.

    Notes
    -----
    We do not need to check the line "finished without error" in scripts/recon-all.log.
    If an error occurs, it will be detected by Nipype and the next nodes (including
    write_tsv_files will not be called).
    """
    from clinica.pipelines.anatomical.freesurfer.utils import generate_regional_measures
    from clinica.utils.stream import cprint

    if "long" in image_id:
        sub, ses, long = image_id.split("_")
        folder = f"{sub}_{ses}.long.{sub}_{long}"
    else:
        folder = image_id
    if (subject_dir / folder / "mri" / "aparc+aseg.mgz").is_file():
        generate_regional_measures(subject_dir, folder, [atlas])
    else:
        cprint(
            msg=(
                f"{image_id.replace('_', ' | ')} does not contain "
                f"mri/aseg+aparc.mgz file. Creation of regional_measures/ folder will be skipped."
            ),
            lvl="warning",
        )
    return image_id
