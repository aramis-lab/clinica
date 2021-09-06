# coding: utf8


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


def get_conversion_luts():
    import os

    from clinica.utils.exceptions import ClinicaException

    try:
        # For aparc+aseg.mgz file:
        default = os.path.join(
            os.environ["MRTRIX_HOME"],
            "share",
            "mrtrix3",
            "labelconvert",
            "fs_default.txt",
        )
        # For aparc.a2009s+aseg.mgz file:
        a2009s = os.path.join(
            os.environ["MRTRIX_HOME"],
            "share",
            "mrtrix3",
            "labelconvert",
            "fs_a2009s.txt",
        )

        # TODO: Add custom Lausanne2008 conversion LUTs here.
    except KeyError:
        raise ClinicaException("Could not find MRTRIX_HOME environment variable.")
    return [default, a2009s]


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
