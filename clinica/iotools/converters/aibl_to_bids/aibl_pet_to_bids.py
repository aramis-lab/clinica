# coding=utf-8
"""
Covert the AIBL PET images into the BIDS specification.
There are three pet modalities: av45, pib, flute. All of them are converted in BIDS


@author: Simona Bottani

"""


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



def av45_paths_to_bids(path_to_dataset, path_to_csv, bids_dir):
    '''
        This methods converts all the PET images-av45 in BIDS

        :param path_to_dataset: path_to_dataset
        :param path_to_csv: path_to_csv with clinical data
        :param bids_dir: path to save the AIBL-PET-dataset converted in a BIDS format
        :return: all the images are converted in a BIDS format and saved in the bids_dir        
    '''
    import os
    import pandas
    from os import path
    from numpy import nan
    from iotools.converters.aibl_to_bids.converter_utils import dicom_to_nii, viscode_to_session,find_path_to_pet_modality
    # it reads the dataframe where subject_ID, session_ID and path are saved
    csv_file = pandas.read_csv(os.path.join(path_to_csv, 'aibl_av45meta_28-Apr-2015.csv'))
    images = find_path_to_pet_modality(path_to_dataset, csv_file)
    images.to_csv(os.path.join(bids_dir, 'pet_av45_paths_aibl.tsv'), sep='\t')
    count = 0
    total = images.shape[0]
    for row in images.iterrows():
        # it iterates for each row of the dataframe which contains the T1_paths
        image = row[1]
        subject = image.Subjects_ID
        session = image.Session_ID
        image_path = image.Path_to_pet
        count += 1
        if image_path is nan:
            print 'No path specified for ' + subject + ' in session ' + session
            continue
        print 'Processing subject ' + str(subject) + ' - session ' + session + ', ' + str(count) + ' / ' + str(total)
        session = viscode_to_session(session)
        output_path = path.join(bids_dir, 'sub-AIBL' + subject, 'ses-' + session, 'pet')
        output_filename = 'sub-AIBL' + subject + '_ses-' + session + '_task-rest_acq-av45_pet'
        # image are saved following BIDS specifications
        if path.exists(path.join(output_path,output_filename + '.nii.gz')):
            pass
        else:
            output_image=dicom_to_nii(subject,output_path,output_filename, image_path,dcm2niix='dcm2niix',dcm2nii='dcm2nii',mri_convert='mri_convert')


def pib_paths_to_bids(path_to_dataset, path_to_csv, bids_dir, dcm2niix="dcm2niix", dcm2nii="dcm2nii"):
    '''
        This methods converts all the PET images-pib in BIDS

        :param path_to_dataset: path_to_dataset
        :param path_to_csv: path_to_csv with clinical data
        :param bids_dir: path to save the AIBL-PET-dataset converted in a BIDS format
        :return: all the images are converted in a BIDS format and saved in the bids_dir
    '''
    import os
    import pandas
    from os import path
    from numpy import nan
    from iotools.converters.aibl_to_bids.converter_utils import dicom_to_nii, viscode_to_session, \
        find_path_to_pet_modality
    # it reads the dataframe where subject_ID, session_ID and path are saved
    csv_file = pandas.read_csv(os.path.join(path_to_csv, 'aibl_pibmeta_28-Apr-2015.csv'))
    images = find_path_to_pet_modality(path_to_dataset, csv_file)
    images.to_csv(os.path.join(bids_dir, 'pet_pib_paths_aibl.tsv'), sep='\t')
    count = 0
    total = images.shape[0]
    for row in images.iterrows():
        # it iterates for each row of the dataframe which contains the T1_paths
        image = row[1]
        subject = image.Subjects_ID
        session = image.Session_ID
        image_path = image.Path_to_pet
        count += 1
        if image_path is nan:
            print 'No path specified for ' + subject + ' in session ' + session
            continue
        print 'Processing subject ' + str(subject) + ' - session ' + session + ', ' + str(count) + ' / ' + str(total)
        session = viscode_to_session(session)
        output_path = path.join(bids_dir, 'sub-AIBL' + subject, 'ses-' + session, 'pet')
        output_filename = 'sub-AIBL' + subject + '_ses-' + session + '_task-rest_acq-pib_pet'
        # image are saved following BIDS specifications
        if path.exists(path.join(output_path,output_filename + '.nii.gz')):
            pass
        else:
            output_image = dicom_to_nii(subject, output_path, output_filename, image_path, dcm2niix='dcm2niix',
                                    dcm2nii='dcm2nii', mri_convert='mri_convert')


def flute_paths_to_bids(path_to_dataset, path_to_csv, bids_dir, dcm2niix="dcm2niix", dcm2nii="dcm2nii"):
    '''
        This methods converts all the PET images-flute in BIDS

        :param path_to_dataset: path_to_dataset
        :param path_to_csv: path_to_csv with clinical data
        :param bids_dir: path to save the AIBL-PET-dataset converted in a BIDS format
        :return: all the images are converted in a BIDS format and saved in the bids_dir
    '''
    import os
    import pandas
    from os import path
    from numpy import nan
    from iotools.converters.aibl_to_bids.converter_utils import dicom_to_nii, viscode_to_session, \
        find_path_to_pet_modality
    # it reads the dataframe where subject_ID, session_ID and path are saved

    csv_file = pandas.read_csv(os.path.join(path_to_csv, 'aibl_flutemeta_28-Apr-2015.csv'))
    images = find_path_to_pet_modality(path_to_dataset, csv_file)
    images.to_csv(os.path.join(bids_dir, 'pet_flute_paths_aibl.tsv'), sep='\t')
    count = 0
    total = images.shape[0]
    for row in images.iterrows():
        # it iterates for each row of the dataframe which contains the T1_paths
        image = row[1]
        subject = image.Subjects_ID
        session = image.Session_ID
        image_path = image.Path_to_pet
        count += 1
        if image_path is nan:
            print 'No path specified for ' + subject + ' in session ' + session
            continue
        print 'Processing subject ' + str(subject) + ' - session ' + session + ', ' + str(count) + ' / ' + str(total)

        session = viscode_to_session(session)
        output_path = path.join(bids_dir, 'sub-AIBL' + subject, 'ses-' + session, 'pet')
        output_filename = 'sub-AIBL' + subject + '_ses-' + session + '_task-rest_acq-flute_pet'
        # image are saved following BIDS specifications
        if path.exists(path.join(output_path,output_filename + '.nii.gz')):
            pass
        else:
            output_image = dicom_to_nii(subject, output_path, output_filename, image_path, dcm2niix='dcm2niix',
                                    dcm2nii='dcm2nii', mri_convert='mri_convert')

