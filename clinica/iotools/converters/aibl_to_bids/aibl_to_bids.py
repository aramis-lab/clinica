# coding: utf8

"""
Convert the AIBL dataset (http://www.aibl.csiro.au/) into BIDS.
"""

__author__ = "Simona Bottani"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Simona Bottani"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Development"


def convert_images(path_to_dataset, path_to_csv, bids_dir):

    # Conversion of the entire dataset in BIDS
    from clinica.utils.stream import cprint
    from os.path import exists
    from colorama import Fore

    from clinica.iotools.converters.aibl_to_bids.aibl_utils import paths_to_bids
    list_of_created_files = []

    for modality in ['t1', 'av45', 'flute', 'pib']:
        list_of_created_files.append(paths_to_bids(path_to_dataset,
                                                   path_to_csv,
                                                   bids_dir,
                                                   modality))

    error_string = ''
    for modality_list in list_of_created_files:
        for file in modality_list:
            if not exists(str(file)):
                error_string = error_string + str(file) + '\n'
    if error_string != '':
        cprint(Fore.RED + 'The following file were not converted '
               + ' (nan means no path was found):\n'
               + error_string
               + Fore.RESET)


def convert_clinical_data(bids_dir, path_to_csv):
    from os.path import exists
    # clinical specifications in BIDS
    from os.path import join, split, realpath
    from clinica.iotools.converters.aibl_to_bids.aibl_utils import create_participants_df_AIBL, \
        create_sessions_dict_AIBL

    clinical_spec_path = join(split(realpath(__file__))[0], '../../data/clinical_specifications.xlsx')
    if not exists(clinical_spec_path):
        raise FileNotFoundError(clinical_spec_path + ' file not found ! This is an internal file of Clinica.')

    create_participants_df_AIBL(bids_dir, clinical_spec_path, path_to_csv, delete_non_bids_info=True)
    create_sessions_dict_AIBL(bids_dir, path_to_csv, clinical_spec_path)
