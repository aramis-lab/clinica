# coding=utf-8
"""
Covert the AIBL T1 images into the BIDS specification.


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

def t1_paths_to_bids(path_to_dataset,path_to_csv, bids_dir):
    ''' 
        This method converts all the T1 images found in the AIBL dataset downloaded in BIDS

        :param path_to_dataset: path_to_dataset
        :param path_to_csv: path to the csv file containing clinical data
        :param bids_dir: path to save the AIBL-T1-dataset converted in a BIDS format
        :return: all the images are converted in a BIDS format and saved in the bids_dir
    '''
    from os import path
    from numpy import nan
    from iotools.converters.aibl_to_bids.converter_utils import dicom_to_nii, viscode_to_session, find_path_to_T1
    #it reads the dataframe where subject_ID, session_ID and path are saved
    images=find_path_to_T1(path_to_dataset,path_to_csv)
    images.to_csv(path.join(bids_dir, 'T1_MRI_paths.tsv'), sep='\t', index=False)
    count = 0
    total = images.shape[0]

    for row in images.iterrows():
        #it iterates for each row of the dataframe which contains the T1_paths
        image = row[1]
        subject = image.Subjects_ID
        session=image.Session_ID
        image_path=image.Path_to_T1
        count += 1

        if image_path is nan:
            print 'No path specified for ' + subject + ' in session ' + session
            continue
        print 'Processing subject ' + str(subject) + ' - session ' + session + ', ' + str(count) + ' / ' + str(total)

        session=viscode_to_session(session)
        #creation of the path
        output_path = path.join(bids_dir, 'sub-AIBL' + subject, 'ses-' + session, 'anat')
        output_filename = 'sub-AIBL' + subject + '_ses-' + session + '_T1w'
        #image are saved following BIDS specifications
        if path.exists(path.join(output_path,output_filename + '.nii.gz')):
            pass
        else:
            output_image = dicom_to_nii(subject, output_path, output_filename, image_path, dcm2niix='dcm2niix',
                                    dcm2nii='dcm2nii', mri_convert='mri_convert')

