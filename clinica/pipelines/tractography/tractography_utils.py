# coding: utf8

"""Tractography - Clinica Utilities.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""

def get_containers(subjects, sessions):

    return [
        'subjects/' + subjects[i] + '/' + sessions[i] + '/dwi  '
        for i in range(len(subjects))
    ]

def get_substitutions(subjects, sessions):
    return [
        [
            ('response/wm.txt', subjects[i] + '_response.txt'),
            ('fod/wm.mif', subjects[i] + '_fod.mif'),
            ('tracts/tracked.tck', subjects[i] + '_tracts.tck'),
            ('trait_added', ''),
        ]
        for i in range(len(subjects))
    ]


def get_caps_filenames(dwi_file):

    import re
    from clinica.utils.stream import cprint

    m = re.search(r'\/(sub-[a-zA-Z0-9]+_ses-[a-zA-Z0-9]+.*)_preproc', dwi_file)

    if m is None:
        raise ValueError('Input filename is not in a CAPS compliant format.')

    source_file = m.group(1)

    cprint(source_file)

    response = source_file + '_responsefunction.txt'
    fod =  source_file + '_fod.mif'
    tracts = source_file + '_tract.tck'

    return response, fod, tracts