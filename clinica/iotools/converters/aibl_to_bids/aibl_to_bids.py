# coding=utf-8
"""
Covert the AIBL dataset into the BIDS specification.


@author: Simona Bottani

"""

from iotools.converters.aibl_to_bids.aibl_pet_to_bids import av45_paths_to_bids, flute_paths_to_bids, pib_paths_to_bids
from iotools.converters.aibl_to_bids.aibl_t1_to_bids import t1_paths_to_bids
from iotools.converters.aibl_to_bids.clinical_data import create_participants_df_AIBL, create_sessions_dict_AIBL

'''
__author__ = "Simona Bottani"
__copyright__ = "Copyright 2011, The Aramis Lab Team"
__credits__ = ["Simona Bottani"]
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Development"

'''
def convert_all_dataset(path_to_dataset,path_to_csv,bids_dir):
    #Conversion of the entire dataset in BIDS
    
    av45_paths_to_bids(path_to_dataset, path_to_csv, bids_dir)
    flute_paths_to_bids(path_to_dataset, path_to_csv, bids_dir)
    pib_paths_to_bids(path_to_dataset, path_to_csv, bids_dir)
    t1_paths_to_bids(path_to_dataset, path_to_csv, bids_dir)


def convert_clinical_data(bids_dir,clinical_spec_path,path_to_csv):
    #clinical specifications in BIDS
    create_participants_df_AIBL(bids_dir, clinical_spec_path, path_to_csv, delete_non_bids_info=True)
    create_sessions_dict_AIBL(bids_dir, path_to_csv, clinical_spec_path)

if __name__ == '__main__':
    #run the function to convert the entire AIBL dataset

    import sys
    path_to_dataset=sys.argv[1]
    path_to_csv=sys.argv[2]
    bids_dir=sys.argv[3]
    clinical_spec_path=sys.argv[4]
    #path_to_dataset='/Users/simona.bottani/AIBL/AIBL_10zip/AIBL'
    #path_to_csv='/Users/simona.bottani/AIBL/Data_extract_3.2.5'
    #bids_dir='/Users/simona.bottani/AIBL/AIBL_BIDS'
    #clinical_spec_path='/Users/simona.bottani/clinica-aramis/iotools/data/clinical_specifications.xls'
    #  bids_dir='/Volumes/aramis-projects/CLINICA/CLINICA_datasets/BIDS/AIBL_BIDS'

    convert_all_dataset(path_to_dataset,path_to_csv,bids_dir)
    convert_clinical_data(bids_dir, clinical_spec_path, path_to_csv)
