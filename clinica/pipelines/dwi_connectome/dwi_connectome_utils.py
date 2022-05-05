def get_luts():
    import os

    from clinica.utils.exceptions import ClinicaException

    try:
        # For aparc+aseg.mgz file:
        default = os.path.join(os.environ["FREESURFER_HOME"], "FreeSurferColorLUT.txt")
        # For aparc.a2009s+aseg.mgz file:
        a2009s = os.path.join(os.environ["FREESURFER_HOME"], "FreeSurferColorLUT.txt")

        # TODO: Add custom Lausanne2008 LUTs here.
    except KeyError:
        raise ClinicaException("Could not find FREESURFER_HOME environment variable.")
    return [default, a2009s]


def get_conversion_luts_offline():

    from pathlib import Path

    from clinica.utils.stream import cprint

    root = Path.cwd()

    path_to_mappings = root / "clinica" / "resources" / "mappings"
    ref_fs_default = path_to_mappings / "fs_default.txt"
    ref_fs_a2009 = path_to_mappings / "fs_a2009.txt"

    if not (ref_fs_default.is_file()):
        cprint(
            msg=f"Unable to locate (fs_default.txt) for processing",
            lvl="error",
        )
    if not (ref_fs_a2009.is_file()):
        cprint(
            msg=f"Unable to locate (fs_a2009.txt) for processing",
            lvl="error",
        )

    return [ref_fs_default, ref_fs_a2009]


def get_conversion_luts():
    from pathlib import Path
    from os import pardir
    from os.path import abspath, join, dirname

    from clinica.utils.inputs import RemoteFileStructure, fetch_file
    from clinica.utils.stream import cprint

    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))

    path_to_mappings = root / "resources" / "mappings"

    url_mrtrix = "https://raw.githubusercontent.com/MRtrix3/mrtrix3/master/share/mrtrix3/labelconvert/"

    FILE1 = RemoteFileStructure(
        filename="fs_default.txt",
        url=url_mrtrix,
        checksum="6ee07088915fdbcf52b05147ddae86e5fcaf3efc63db5b0ba8f361637dfa11ef",
    )

    FILE2 = RemoteFileStructure(
        filename="fs_a2009s.txt",
        url=url_mrtrix,
        checksum="b472f09cfe92ac0b6694fb6b00a87baf15dd269566e4a92b8a151ff1080bf170",
    )

    ref_fs_default = path_to_mappings / Path(FILE1.filename)
    ref_fs_a2009 = path_to_mappings / Path(FILE2.filename)

    if not (ref_fs_default.is_file()):
        try:
            ref_fs_default = fetch_file(FILE1, path_to_mappings)
        except IOError as err:
            cprint(
                msg=f"Unable to download required MRTRIX mapping (fs_default.txt) for processing: {err}",
                lvl="error",
            )
    if not (ref_fs_a2009.is_file()):
        try:
            ref_fs_a2009 = fetch_file(FILE2, path_to_mappings)
        except IOError as err:
            cprint(
                msg=f"Unable to download required MRTRIX mapping (fs_a2009.txt) for processing: {err}",
                lvl="error",
            )

    return [ref_fs_default, ref_fs_a2009]


def get_containers(subjects, sessions):
    import os

    return [
        os.path.join("subjects", subjects[i], sessions[i], "dwi")
        for i in range(len(subjects))
    ]


def get_caps_filenames(dwi_file: str):
    import re

    m = re.search(r"/(sub-[a-zA-Z0-9]+_ses-[a-zA-Z0-9]+.*)_preproc", dwi_file)
    if not m:
        raise ValueError(
            f"Input filename {dwi_file} is not in a CAPS compliant format."
        )
    source_file_caps = m.group(1)

    m = re.search(
        r"/(sub-[a-zA-Z0-9]+_ses-[a-zA-Z0-9]+.*)_space-[a-zA-Z0-9]+_preproc", dwi_file
    )
    if not m:
        raise ValueError(
            f"Input filename {dwi_file} is not in a CAPS compliant format."
        )
    source_file_bids = m.group(1)

    response = f"{source_file_caps}_model-CSD_responseFunction.txt"
    fod = f"{source_file_caps}_model-CSD_diffmodel.nii.gz"
    tracts = f"{source_file_caps}_model-CSD_tractography.tck"
    nodes = [
        f"{source_file_caps}_atlas-desikan_parcellation.nii.gz",
        f"{source_file_caps}_atlas-destrieux_parcellation.nii.gz",
    ]
    # TODO: Add custom Lausanne2008 node files here.
    connectomes = [
        f"{source_file_bids}_model-CSD_atlas-desikan_connectivity.tsv",
        f"{source_file_bids}_model-CSD_atlas-destrieux_connectivity.tsv",
    ]
    # TODO: Add custom Lausanne2008 connectome files here.

    return response, fod, tracts, nodes, connectomes


def print_begin_pipeline(in_bids_or_caps_file: str) -> None:
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    print_begin_image(get_subject_id(in_bids_or_caps_file))


def print_end_pipeline(in_bids_or_caps_file: str, final_file: str) -> None:
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(in_bids_or_caps_file))
