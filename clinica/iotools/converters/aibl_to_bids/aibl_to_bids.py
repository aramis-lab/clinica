# coding=utf-8
"""
Covert the AIBL dataset into the BIDS specification.


@author: Simona Bottani

"""

__author__ = "Simona Bottani"
__copyright__ = "Copyright 2011, The Aramis Lab Team"
__credits__ = ["Simona Bottani"]
__license__ = ""
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Development"


def convert_images(path_to_dataset, path_to_csv, bids_dir):

    # Conversion of the entire dataset in BIDS

    from clinica.iotools.converters.aibl_to_bids.aibl_utils import av45_paths_to_bids, flute_paths_to_bids, \
        pib_paths_to_bids
    from clinica.iotools.converters.aibl_to_bids.aibl_utils import t1_paths_to_bids

    av45_paths_to_bids(path_to_dataset, path_to_csv, bids_dir)
    flute_paths_to_bids(path_to_dataset, path_to_csv, bids_dir)
    pib_paths_to_bids(path_to_dataset, path_to_csv, bids_dir)
    t1_paths_to_bids(path_to_dataset, path_to_csv, bids_dir)


def convert_clinical_data(bids_dir, path_to_csv):
    # clinical specifications in BIDS
    import os
    from clinica.iotools.converters.aibl_to_bids.aibl_utils import create_participants_df_AIBL, \
        create_sessions_dict_AIBL

    CLINICA_HOME = os.environ.get('CLINICA_HOME', '')
    if not CLINICA_HOME:
        raise Exception('CLINICA_HOME variable from Clinica software is not set')

    clinical_spec_path = os.path.join(CLINICA_HOME, 'clinica', 'iotools', 'data', 'clinical_specifications.xlsx')
    create_participants_df_AIBL(bids_dir, clinical_spec_path, path_to_csv, delete_non_bids_info=True)
    create_sessions_dict_AIBL(bids_dir, path_to_csv, clinical_spec_path)
