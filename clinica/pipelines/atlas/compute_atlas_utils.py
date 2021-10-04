# coding: utf8


from clinica.utils.stream import cprint


def compute_atlases(
    caps_directory: str, to_process_with_atlases: tuple, path_to_atlas: str
) -> tuple(str, str, str):
    import os
    import subprocess
    from pathlib import Path

    from clinica.utils.stream import cprint

    subject_dir = ""
    image_id = ""
    atlas = ""
    if to_process_with_atlases != []:
        for path in Path(path_to_atlas).rglob(
            "*" + to_process_with_atlases[0][0] + "_6p0.gcs"
        ):
            hemisphere = os.path.split(path)[1].rsplit(".")[0]
            atlas_name = os.path.split(path)[1].rsplit(".")[1].split("_")[0]
            sub, ses = to_process_with_atlases[0][1].split("_")

            path_to_freesurfer_cross = (
                caps_directory
                + "/subjects/"
                + sub
                + "/"
                + ses
                + "/t1/freesurfer_cross_sectional/"
                + to_process_with_atlases[0][1]
            )
            output_path_annot = (
                path_to_freesurfer_cross
                + "/label/"
                + hemisphere
                + "."
                + atlas_name
                + ".annot"
            )
            sphere_reg = (
                path_to_freesurfer_cross + "/surf/" + hemisphere + ".sphere.reg"
            )
            if not os.path.isfile(sphere_reg):
                cprint(
                    f"The {hemisphere}.sphere.reg file appears to be missing. The data for {to_process_with_atlases[0][1]} will not be processed with {to_process_with_atlases[0][0]}",
                    lvl="warning",
                )
            command = f"mris_ca_label {to_process_with_atlases[0][1]} {hemisphere} {path_to_freesurfer_cross}/surf/{hemisphere}.sphere.reg {path} {output_path_annot}"
            a = subprocess.run(command, shell=True, capture_output=True)

            output_path_stats = (
                path_to_freesurfer_cross
                + "/stats/"
                + hemisphere
                + "."
                + atlas_name
                + ".stats"
            )
            command2 = f"mris_anatomical_stats -a {output_path_annot} -f {output_path_stats} -b {to_process_with_atlases[0][1]} {hemisphere}"
            c = subprocess.run(command2, shell=True, capture_output=True)

            image_id, atlas = (
                to_process_with_atlases[0][1],
                to_process_with_atlases[0][0],
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

    from clinica.utils.freesurfer import generate_regional_measures_alt
    from clinica.utils.stream import cprint

    if os.path.isfile(os.path.join(subject_dir, image_id, "mri", "aparc+aseg.mgz")):
        generate_regional_measures_alt(subject_dir, image_id, atlas)
    else:
        cprint(
            msg=(
                f"{image_id.replace('_', ' | ')} does not contain "
                f"mri/aseg+aparc.mgz file. Creation of regional_measures/ folder will be skipped."
            ),
            lvl="warning",
        )
    return image_id
