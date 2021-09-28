# coding: utf8


def init_input_node(t1w, recon_all_args, output_dir):
    """Initialize the pipeline.

    This function will:
        - Extract <image_id> (e.g. sub-CLNC01_ses-M00) T1w filename;
        - Check FOV of T1w;
        - Create SUBJECTS_DIR for recon-all (otherwise, the command won't run);
        - Print begin execution message.
    """
    import os

    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.freesurfer import check_flags
    from clinica.utils.ux import print_begin_image

    # Extract <image_id>
    image_id = get_subject_id(t1w)

    # Create SUBJECTS_DIR for recon-all (otherwise, the command won't run)
    subjects_dir = os.path.join(output_dir, image_id)
    os.makedirs(subjects_dir, exist_ok=True)

    return image_id, t1w, subjects_dir


def compute_atlases(caps_directory, to_process_with_atlases, path_to_atlas):
    import os
    import subprocess
    from pathlib import Path

    from clinica.utils.freesurfer import generate_regional_measures_alt

    print("\n\nqqch")
    print("ze caps: ", caps_directory)
    print("ze atlas path: ", path_to_atlas)
    print("ze current subject and session: ", "something smart")
    print("ze NOT dictionnary: ", to_process_with_atlases, "\n\n")

    # for sub_n_ses in to_process_with_atlases[0]:
    #     print("my favorite atlas: ", to_process_with_atlases[0],
    #     "\nandis sub_n_ses: ", sub_n_ses)
    subject_dir = ""
    image_id = ""
    atlas = ""
    print("to_process_with_atlases: ", to_process_with_atlases)
    if to_process_with_atlases != []:
        print(
            "ze atlas: ",
            to_process_with_atlases[0][0],
            "\nand ze subject: ",
            to_process_with_atlases[0][1],
        )

        for path in Path(path_to_atlas).rglob(
            "*" + to_process_with_atlases[0][0] + "_6p0.gcs"
        ):
            hemisphere = os.path.split(path)[1].rsplit(".")[0]
            atlas_name = os.path.split(path)[1].rsplit(".")[1].split("_")[0]
            sub, ses = to_process_with_atlases[0][1].split("_")
            print("sub: ", sub)
            print("ses: ", ses)

            atlas_file = hemisphere + "." + atlas_name + "_6p0.gcs"
            subject_dir = (
                caps_directory
                + "/subjects/"
                + sub
                + "/"
                + ses
                + "/t1/freesurfer_cross_sectional/"
            )
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

            command = f"mris_ca_label {to_process_with_atlases[0][1]} {hemisphere} {path_to_freesurfer_cross}/surf/{hemisphere}.sphere.reg {path} {output_path_annot}"
            # print("\n\nze command: ", command, "\n\n")
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
            print(
                "calling with: \n",
                "path_to_freesurfer_cross: ",
                path_to_freesurfer_cross,
                "\nsubnses: ",
                to_process_with_atlases[0][1],
                "\natlas: ",
                to_process_with_atlases[0][0],
            )

            image_id, atlas = (
                to_process_with_atlases[0][1],
                to_process_with_atlases[0][0],
            )
    return subject_dir, image_id, atlas


def write_tsv_files(subject_dir, image_id, atlas):
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
