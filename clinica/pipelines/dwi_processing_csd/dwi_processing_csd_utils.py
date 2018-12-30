# coding: utf8

"""Tractography - Clinica Utilities.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""


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
        default = os.path.join(os.environ['MRTRIX3_HOME'],
                               'share/mrtrix3/labelconvert/fs_default.txt')
        # For aparc.a2009s+aseg.mgz file:
        a2009s = os.path.join(os.environ['MRTRIX3_HOME'],
                              'share/mrtrix3/labelconvert/fs_a2009s.txt')

        # TODO: Add custom Lausanne2008 conversion LUTs here.
    except KeyError:
        raise ClinicaException('Could not find MRTRIX3_HOME environment variable.')

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

    source_file = m.group(1)

    response = source_file + '_model-CSD_diffmodel.txt'
    fod = source_file + '_model-CSD_FOD.mif'
    tracts = source_file + '_model-CSD_tractography.tck'
    nodes = [source_file + '_parcellation-desikan_node.nii.gz',
             source_file + '_parcellation-destrieux_node.nii.gz']
    # TODO: Add custom Lausanne2008 node files here.
    connectomes = [source_file + '_model-CSD_parcellation-desikan_connectivity.tsv',
                   source_file + '_model-CSD_parcellation-destrieux_connectivity.tsv']
    # TODO: Add custom Lausanne2008 connectome files here.

    return response, fod, tracts, nodes, connectomes
