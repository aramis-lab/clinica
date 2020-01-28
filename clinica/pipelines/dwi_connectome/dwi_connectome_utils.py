# coding: utf8


def get_luts():
    import os
    from clinica.utils.exceptions import ClinicaException

    try:
        # For aparc+aseg.mgz file:
        default = os.path.join(os.environ['FREESURFER_HOME'],
                               'FreeSurferColorLUT.txt')
        # For aparc.a2009s+aseg.mgz file:
        a2009s = os.path.join(os.environ['FREESURFER_HOME'],
                              'FreeSurferColorLUT.txt')

        # TODO: Add custom Lausanne2008 LUTs here.
    except KeyError:
        raise ClinicaException('Could not find FREESURFER_HOME environment variable.')

    return [default, a2009s]


def get_conversion_luts():
    import os
    from clinica.utils.exceptions import ClinicaException

    try:
        # For aparc+aseg.mgz file:
        default = os.path.join(os.environ['MRTRIX_HOME'],
                               'share/mrtrix3/labelconvert/fs_default.txt')
        # For aparc.a2009s+aseg.mgz file:
        a2009s = os.path.join(os.environ['MRTRIX_HOME'],
                              'share/mrtrix3/labelconvert/fs_a2009s.txt')

        # TODO: Add custom Lausanne2008 conversion LUTs here.
    except KeyError:
        raise ClinicaException('Could not find MRTRIX_HOME environment variable.')

    return [default, a2009s]


def get_containers(subjects, sessions):

    return [
        'subjects/' + subjects[i] + '/' + sessions[i] + '/dwi'
        for i in range(len(subjects))
    ]


def get_caps_filenames(dwi_file):

    import re

    m = re.search(r'\/(sub-[a-zA-Z0-9]+_ses-[a-zA-Z0-9]+.*)_preproc', dwi_file)
    if m is None:
        raise ValueError('Input filename is not in a CAPS compliant format.')
    source_file_caps = m.group(1)

    m = re.search(r'\/(sub-[a-zA-Z0-9]+_ses-[a-zA-Z0-9]+.*)_space-[a-zA-Z0-9]+_preproc', dwi_file)
    if m is None:
        raise ValueError('Input filename is not in a CAPS compliant format.')
    source_file_bids = m.group(1)

    response = source_file_caps + '_model-CSD_responseFunction.txt'
    fod = source_file_caps + '_model-CSD_diffmodel.nii.gz'
    tracts = source_file_caps + '_model-CSD_tractography.tck'
    nodes = [source_file_caps + '_atlas-desikan_parcellation.nii.gz',
             source_file_caps + '_atlas-destrieux_parcellation.nii.gz']
    # TODO: Add custom Lausanne2008 node files here.
    connectomes = [source_file_bids + '_model-CSD_atlas-desikan_connectivity.tsv',
                   source_file_bids + '_model-CSD_atlas-destrieux_connectivity.tsv']
    # TODO: Add custom Lausanne2008 connectome files here.

    return response, fod, tracts, nodes, connectomes


def print_begin_pipeline(in_bids_or_caps_file):
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    print_begin_image(get_subject_id(in_bids_or_caps_file))


def print_end_pipeline(in_bids_or_caps_file, final_file):
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(in_bids_or_caps_file))


def convert_flirt_transformation_to_mrtrix_transformation(
        in_source_image,
        in_reference_image,
        in_flirt_matrix,
        name_output_matrix=None):
    """
    Convert flirt matrix to mrtrix matrix.

    This function converts a transformation matrix produced by FSL's flirt
    command into a format usable by MRtrix. The output of this function
    is usually for the mrtransform command.

    Args:
        in_source_image (str): File containing the source image used in
            FSL flirt with the -in flag.
        in_reference_image (str): File containing the reference image used in
            FSL flirt with the -ref flag.
        in_flirt_matrix (str): File containing the transformation matrix
            obtained by FSL flirt.
        name_output_matrix (Optional[str]): Name of the output matrix
            (default=deformed_image.nii.gz).

    Returns:
        out_mrtrix_matrix (str): Transformation matrix in MRtrix format.
    """
    import os
    from clinica.utils.check_dependency import check_mrtrix
    check_mrtrix()

    assert(os.path.isfile(in_source_image))
    assert(os.path.isfile(in_reference_image))
    assert(os.path.isfile(in_flirt_matrix))

    if name_output_matrix is None:
        out_mrtrix_matrix = os.path.abspath('mrtrix_matrix.mat')
    else:
        out_mrtrix_matrix = os.path.abspath(name_output_matrix)

    cmd = 'transformconvert %s %s %s flirt_import %s' \
          % (in_flirt_matrix, in_source_image, in_reference_image,
             out_mrtrix_matrix)
    os.system(cmd)

    return out_mrtrix_matrix
