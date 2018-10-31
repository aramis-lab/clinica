# coding: utf8

"""
Convert the AIBL dataset into the BIDS specification.
"""

__author__ = "Simona Bottani"
__copyright__ = "Copyright 2016-2018 The Aramis Lab Team"
__credits__ = ["Simona Bottani"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Development"


def convert_images(path_to_dataset, path_to_csv, bids_dir):

    # Conversion of the entire dataset in BIDS

    from clinica.utils.stream import cprint
    from clinica.iotools.converters.aibl_to_bids.aibl_utils import av45_paths_to_bids, flute_paths_to_bids, \
        pib_paths_to_bids
    from clinica.iotools.converters.aibl_to_bids.aibl_utils import t1_paths_to_bids

    for converter in [av45_paths_to_bids, flute_paths_to_bids, pib_paths_to_bids, t1_paths_to_bids] :
        try:
            converter(path_to_dataset, path_to_csv, bids_dir)
        except:
            cprint("Convertion failed ")


def convert_clinical_data(bids_dir, path_to_csv):
    # clinical specifications in BIDS
    import os
    from clinica.utils.stream import cprint
    from os.path import join, split, realpath
    from clinica.iotools.converters.aibl_to_bids.aibl_utils import create_participants_df_AIBL, \
        create_sessions_dict_AIBL

    clinical_spec_path = join(split(realpath(__file__))[0], '../../iotools/data/clinical_specifications.xlsx')
    try:
        create_participants_df_AIBL(bids_dir, clinical_spec_path, path_to_csv, delete_non_bids_info=True)
        create_sessions_dict_AIBL(bids_dir, path_to_csv, clinical_spec_path)
    except:
        cprint("Convertion clinical data failed")

