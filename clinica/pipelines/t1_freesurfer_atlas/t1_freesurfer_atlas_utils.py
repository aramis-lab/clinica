# coding: utf8


from pathlib import Path

from clinica.utils.stream import cprint


def compute_atlas(
    caps_directory: str, to_process_with_atlases: tuple, path_to_atlas: str
):
    import os
    import subprocess
    from pathlib import Path
    from tempfile import mkdtemp

    from clinica.compat import errno
    from clinica.utils.stream import cprint

    subject_dir = ""
    image_id = ""
    atlas = ""

    # we iterate over left and right hemisphere
    for path in Path(path_to_atlas).rglob(
        "*" + to_process_with_atlases[0] + "_6p0.gcs"
    ):
        hemisphere = Path(path.stem).stem
        atlas_name = Path(path.stem).suffix[1:].split("_")[0]
        if "long" not in to_process_with_atlases[1]:
            sub, ses = to_process_with_atlases[1].split("_")
            substitute_dir = (
                Path(caps_directory)
                / "subjects"
                / sub
                / ses
                / ("t1/freesurfer_cross_sectional")
            )
            subject_dir = substitute_dir
            path_to_freesurfer_cross = Path(substitute_dir) / Path(
                to_process_with_atlases[1]
            )
            output_path_annot = (
                path_to_freesurfer_cross
                / "label"
                / (hemisphere + "." + atlas_name + ".annot")
            )
            sphere_reg = (
                path_to_freesurfer_cross / "surf" / (hemisphere + ".sphere.reg")
            )
            output_path_stats = (
                path_to_freesurfer_cross
                / "stats"
                / (hemisphere + "." + atlas_name + ".stats")
            )
            subjid = to_process_with_atlases[1]
        else:
            sub, ses, long = to_process_with_atlases[1].split("_")
            substitute_dir = mkdtemp()
            os.makedirs(substitute_dir, exist_ok=True)
            subjid = sub + "_" + ses + ".long." + sub + "_" + long
            subject_dir = (
                Path(caps_directory)
                / "subjects"
                / sub
                / ses
                / "t1"
                / long
                / "freesurfer_longitudinal"
            )
            path_to_freesurfer_cross = subject_dir / (
                sub + "_" + ses + ".long." + sub + "_" + long
            )
            try:
                os.symlink(
                    path_to_freesurfer_cross, os.path.join(substitute_dir, subjid)
                )
            except FileExistsError as e:
                if e.errno != errno.EEXIST:  # EEXIST: folder already exists
                    raise e
            substitute_dir_id = os.path.join(substitute_dir, subjid)

            output_path_annot = (
                Path(substitute_dir_id)
                / "label"
                / (hemisphere + "." + atlas_name + ".annot")
            )
            sphere_reg = Path(substitute_dir_id) / "surf" / (hemisphere + ".sphere.reg")
            output_path_stats = (
                Path(substitute_dir_id)
                / "stats"
                / (hemisphere + "." + atlas_name + ".stats")
            )
        if not os.path.isfile(sphere_reg):
            cprint(
                f"The {hemisphere}.sphere.reg file appears to be missing. The data for {to_process_with_atlases[1]} will not be processed with {to_process_with_atlases[0]}",
                lvl="warning",
            )
        mris_ca_label_command = f"mris_ca_label -sdir {substitute_dir} {subjid} {hemisphere} {path_to_freesurfer_cross}/surf/{hemisphere}.sphere.reg {path} {output_path_annot}"
        a = subprocess.run(mris_ca_label_command, shell=True, capture_output=True)
        mris_anatomical_stats_command = f"export SUBJECTS_DIR={substitute_dir}; mris_anatomical_stats -a {output_path_annot} -f {output_path_stats} -b {subjid} {hemisphere}"
        c = subprocess.run(
            mris_anatomical_stats_command, shell=True, capture_output=True
        )
        image_id, atlas = (
            to_process_with_atlases[1],
            to_process_with_atlases[0],
        )
    return subject_dir, image_id, atlas


def write_tsv_files(subject_dir: str, image_id: str, atlas: str) -> str:
    """
    Generate statistics TSV files in `subjects_dir`/regional_measures folder for `image_id`.

    Notes:
        We do not need to check the line "finished without error" in scripts/recon-all.log.
        If an error occurs, it will be detected by Nipype and the next nodes (including
        write_tsv_files will not be called).
    """
    import os

    from clinica.utils.freesurfer import generate_regional_measures
    from clinica.utils.stream import cprint

    if "long" not in image_id:
        folder = image_id
    else:
        sub, ses, long = image_id.split("_")
        folder = sub + "_" + ses + ".long." + sub + "_" + long
    if os.path.isfile(os.path.join(subject_dir, folder, "mri", "aparc+aseg.mgz")):
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
