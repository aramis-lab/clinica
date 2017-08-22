# coding=utf-8
"""
Utils to convert AIBL dataset in BIDS


@author: Simona Bottani

"""

'''

__author__ = "Simona Bottani"
__copyright__ = "Copyright 2011, The Aramis Lab Team"
__credits__ = ["Simona Bottani", "Sabrina Fontanella"]
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Development"

'''


def listdir_nohidden(path):
    '''
       
       This method lists all the subdirectories of path except the hidden folders'
       
        :param path: path whose subdirectories are needed
        :return: list of all the subdirectories of path
    '''

    import os

    result = []
    for x in os.listdir(path):
        #check if there is the folder .DS_STORE, if it's present, we do not save them
        if x.startswith("."):
            continue
        else:
            result.append(x)
    return result

def find_T1_folder(subdirectory,path_to_T1_1):
    '''
        
        This method checks if the subdirectory contains a T1 image, and it returns the h

        :param subdirectory: name of the folder
        :return: previous path to arrive to the T1 image
    '''

    import os

    #there are if conditions which checks if the subfolder contain a T1 image
    if subdirectory=='MPRAGE_ADNI_confirmed':
        return os.path.join(path_to_T1_1, 'MPRAGE_ADNI_confirmed')
    elif subdirectory=='MPRAGE':
        return os.path.join(path_to_T1_1,'MPRAGE')
    elif subdirectory=='MPRAGE_ADNI_confirmed_RPT':
        return os.path.join(path_to_T1_1, 'MPRAGE_ADNI_confirmed_RPT')
    elif subdirectory=='MPRAGE_ADNI_confirmed_REPEATX2':
        return os.path.join(path_to_T1_1,'MPRAGE_ADNI_confirmed_REPEATX2')
    elif subdirectory=='MPRAGE_ADNI_confirmed_repeat':
        return os.path.join(path_to_T1_1,'MPRAGE_ADNI_confirmed_repeat')
    elif subdirectory=='MPRAGE_ADNI_confirmed_REPEAT':
        return os.path.join(path_to_T1_1,'MPRAGE_ADNI_confirmed_REPEAT')
    elif subdirectory=='MPRAGE_ADNI_conf_REPEAT':
        return os.path.join(path_to_T1_1,'MPRAGE_ADNI_conf_REPEAT')
    else:
        return 'NaN' #there are no more folders which could contain T1 images

def find_T1_folder_nodata(subdirectory,path_to_T1_1):
    '''
       
       This method checks if the subdirectory contains a T1 image, and it returns the path. This method differs from the 
       find_T1_folder since for these folders the exame_date is not present in the clinical excel file and we will not check if the exame_date corresponds
       to the date stored in the path to the image, but they will be converted anyway 
       
        :param subdirectory: name of the folder
        :return: previous path to arrive to the T1 image
    '''
    import os

    if subdirectory=='MPRAGESAGISOp2ND':
        return os.path.join(path_to_T1_1, 'MPRAGESAGISOp2ND')
    elif subdirectory=='MPRAGE_SAG_ISO_p2_ND':
        return os.path.join(path_to_T1_1, 'MPRAGE_SAG_ISO_p2_ND')
    elif subdirectory=='MPRAGE_SAG_ISO_p2':
        return os.path.join(path_to_T1_1,'MPRAGE_SAG_ISO_p2')
    else:
        return 'NaN' # there are no more folders which could contain T1 images

def find_correspondance_index(i,csv_file):
    '''
        
        This method gives as output the index of the csv file analysed which correspond to the 'i' subject

        :param i: subject_ID
        :param csv_file: csv file where all the information are listed
        :return: index
    '''

    index=[]
    for x in csv_file.RID:
        if i == str(x):
            index=csv_file.RID[csv_file.RID==x].index.tolist()
    return index

def find_correspondance_date(index,csv_file):
    '''
        
        The method returns the dates reported in the csv_file for the i-subject 
        
        :param index: index corresponding to the subject analysed
        :param csv_file: csv file where all the information are listed
        :return date
    '''

    csv_date=[]
    csv_date = csv_file.EXAMDATE[index]
    return csv_date

def match_data(exame_date,i,csv_file):
    '''
        
        This method returns the session_ID. It controls if the dates corresponding to the image (from the name
        of the subdirectory) correspond to one of the dates listed from the csv_file for the subject analysed. The session_ID 
        is the corresponding session for that patient in that date.
        It returns -4 if there are no information.
        
        :param exame_date: date where the image has been taken, it is saved from the name of the corresponding subdirector
        :param i: subject_ID
        :param csv_file: csv file where all the information are listed
        :return session_id of the patient
    '''

    import re

    session_ID=[]
    index = find_correspondance_index(i,csv_file)
    csv_date=find_correspondance_date(index,csv_file)
    for xx in index:
        if str(csv_date[xx])!='-4':
            #check is the date is not '-4'
            m = re.search('([0-9].*)-(.*)-(.*)_(.*)_(.*)_(.*)', str(exame_date)) #string from image directory
            p = re.search('(.*)/(.*)/(.*)', str(csv_date[xx])) #string from the date of the csv_file
            if (p.group(1) == m.group(2)) & (p.group(2) == m.group(3)) & (p.group(3) == m.group(1)):
                session_ID=csv_file.VISCODE[xx]
    if session_ID==[]:
        session_ID='-4'
    return session_ID

def list_of_paths():
    ''' 
        
        It lists all the folders which not contain PET images
    '''
    path_not_for_pet = ['.DS_Store', 'localizer',
              'Space_3D_T2_FLAIR_sag_p2',
              'AXIAL_FLAIR',
              'MPRAGE_ADNI_confirmed_REPEATX2',
              'Axial_PD-T2_TSE',
              'Axial_PD-T2_TSE_repeat',
              'MPRAGE_SAG_ISO_p2_ND',
              'Axial_PD-T2_TSE_confirmed',
              'MPRAGESAGISOp2ND',
              'MPRAGE_ADNI_confirmed',
              'MPRAGE_ADNI_confirmed_repeat',
              'MPRAGE_SAG_ISO_p2',
              'MPRAGE',
              'MPRAGE_ADNI_confirmed_REPEAT',
              'Axial_PD-T2_TSE_confirmed_repeat',
              'MPRAGE_ADNI_conf_REPEAT',
              'Space_3D_T2_FLAIR_sag_p2_REPEAT',
              'MPRAGE_ADNI_confirmed_RPT',
              'Brain_256_1.6_zoom_4_x_4_iter',
              'Space_3D_T2_FLAIR_sag_REPEAT',
              'Axial_PD-T2_TSE_RPTconfirmed',
              'Axial_PD-T2_TSE_RPT_confirmed',
              'Axial_PD-T2_TSE_confirmed_REPEAT',
              'flair_t2_spc_irprep_ns_sag_p2_1mm_iso',
              'localiser'
              ]
    return path_not_for_pet

def check_subdirectories_pet(subdirectories,sub, no_pet):
    ''' 
        
        It returns the correct subdirectories for the PET images, they should belong to the list where there all the possible
        names of the PET images 
        
        :param subdirectories: 
        :param sub: all the possible subdirectories which need to be checked
        :param no pet: list of names of folders which not contain PET images 
        :return subdirectory which is containing a PET image which needs to be converted
    '''

    for j in xrange(len(sub)):
        if (sub[j] not in no_pet) & (sub[j] != '.DS_Store'):
            subdirectories.append(sub[j])
    subdirectories=list(set(subdirectories))
    return subdirectories


def compress_nii(file_path):
    '''
        
        Compress .nii file
        
        :param file_path:
        :return: compressed image
    '''

    from os import path
    import os
    import gzip

    f_in = open(file_path)
    f_out = gzip.open(path.join(file_path + '.gz'), 'wb')
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()

    # Remove the original file in order to have only the compressed image
    os.remove(file_path)

def dicom_to_nii(subject,output_path,output_filename, image_path,dcm2niix='dcm2niix',dcm2nii='dcm2nii',mri_convert='mri_convert'):
    '''
        
        From dicom to nifti converts the dicom images in a nifti files using dicom2nii or mri_convert
    
        :param subject: 
        :param output_path: where nifti image is stored
        :param output_filename: name of the nifti image 
        :param image_path: where dicom files are stored 
        :return: Image in a nifti format 
    '''
    import os
    import gzip

    try:
        os.makedirs(output_path)
    except OSError:
        if not os.path.isdir(output_path):
            raise

    # if image.Is_Dicom:
    command = dcm2niix + ' -b n -z n -o ' + output_path + ' -f ' + output_filename + ' ' + image_path
    os.system(command)
    nifti_file = os.path.join(output_path, output_filename + '.nii')

    # Check if conversion worked (output file exists?)
    if not os.path.isfile(nifti_file):
        command = dcm2nii + ' -a n -d n -e n -i y -g n -p n -m n -r n -x n -o ' + output_path + ' ' + image_path
        os.system(command)
        nifti_file = os.path.join(output_path, subject.replace('_', '') + '.nii')

    if nifti_file==os.path.join(output_path,'DE-IDENTIFIED.nii') or nifti_file==os.path.join(output_path,subject+'.nii'):
        #if the conversion dcm2nii has not worked, freesurfer utils mri_convert is used
        dicom_image=listdir_nohidden(image_path)
        dicom_image=os.path.join(image_path,dicom_image[1])
        #it requires the installation of Freesurfer
        command=mri_convert +' ' +dicom_image +' ' +os.path.join(output_path,output_filename+'.nii')
        os.system(command)

        nifti_file_1 = os.path.join(output_path, output_filename + '.nii')
        f_in = open(nifti_file_1)
        f_out = gzip.open(os.path.join(nifti_file_1 + '.gz'), 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()

        # Remove the original file
        os.remove(nifti_file_1)
        if os.path.exists(os.path.join(output_path,'DE-IDENTIFIED.nii')):
            os.remove(os.path.join(output_path,'DE-IDENTIFIED.nii'))

        output_image=os.path.join(output_path, output_filename + '.nii')
    else:
        output_image = compress_nii(nifti_file)
    #it returns the image .nii.gz
    return output_image

def viscode_to_session(viscode):
    '''
    
        Replace the session label 'bl' with 'M00' or capitalize the session name passed
        as input.

        :param viscode: session name
        :return: M00 if is the baseline session or the original session name capitalized
    '''
    if viscode == 'bl':
        return 'M00'
    else:
        return viscode.capitalize()


def find_path_to_pet_modality(path_to_dataset, csv_file):
    '''
        
        This method creates a Dataframe which contains all the paths to the PET image of a modality (for example AV45 or PIB)
    
        :param path_to_dataset: path to AIBL dataset  
        :param csv_file: file which correspond to the modality
        :return: A dataframe which contains the path for PET images for a single modality and subject_ID and session_ID are 
        reported for each path
    '''

    import os
    import pandas

    no_pet = list_of_paths()
    subjects_ID = listdir_nohidden(path_to_dataset)
    # selection of the subjects_ID from the folder downloaded
    #this subject must be discarded since it is only a sample and not a patient
    if '0151083' in subjects_ID:
        del subjects_ID[subjects_ID.index('0151083')]
    sub_ID = []
    ses_ID = []
    path_pet = []
    for i in subjects_ID:
        # Iteration through all the subjects_ID
        if int(i) in list(csv_file.RID):
            # check if the subject is present in the csv_file for the modality selected
            subdirectories = []
            path_to_pet_1 = os.path.join(path_to_dataset, str(i))
            #subdirectory_all = os.listdir(path_to_pet_1)
            subdirectory_all=listdir_nohidden(path_to_pet_1)
            subdirectories = check_subdirectories_pet(subdirectories, subdirectory_all, no_pet)
            # selection only of the folders which contain PET image
            for j in xrange(len(subdirectories)):
                path_to_pet_2 = os.path.join(path_to_pet_1, subdirectories[j])
                exame_date = listdir_nohidden(path_to_pet_2)
                # exame date of the image which is going to be converted
                for x in xrange(len(exame_date)):
                    # selection of the session_ID matching the data in the csv_file with the one of the image
                    session_ID = match_data(exame_date[x], i, csv_file)
                    if session_ID != '-4':
                        path_to_pet_3 = os.path.join(path_to_pet_2, str(exame_date[x]))
                        #For the RID 1607 there are two PET images of the flute modality, and we select the first
                        if i=='1607':
                            if subdirectories[j]=='Flute_256_1.6_Zoom_plain_4_x_4_Iter':
                                image_ID=['I442930']
                            else:
                                image_ID = listdir_nohidden(path_to_pet_3)
                        else:
                            image_ID = listdir_nohidden(path_to_pet_3)
                        for y in xrange(len(image_ID)):
                            # final path to find the image we want to convert
                            path_to_pet = os.path.join(path_to_pet_3, image_ID[y])#
                            sub_ID.append(i)
                            ses_ID.append(session_ID)
                            path_pet.append(path_to_pet)


    data = pandas.DataFrame({'Subjects_ID': sub_ID,
                             'Session_ID': ses_ID,
                             'Path_to_pet': path_pet})
    # data=final dataframe

    return data

def find_path_to_T1_ADNI(file_mri,subjects_ID,path_to_dataset):
    '''

        This method creates a Dataframe which contains all the paths to the T1 images which are ADNI compliant (as
        explained in the AIBL website). This images differ from the others T1 of the dataset since in the cvs_file is 
        reported the exame date.

        :param file_mri: in the clinical data there are two files which describe the  parameters of the T1 images (MRI 1.5 T 
        and MRI 3T)
        :param subjects_ID: subjects_id in the dataset dowloaded
        :param path_to_dataset: path to AIBL dataset
        :return: A dataframe which contains the path for T1 images and subject_ID and session_ID are reported for each path
    '''
    import os

    sub_ID = []
    ses_ID = []
    path_T1 = []

    for i in subjects_ID:
        for jj in file_mri:
            #it checks all the file_mri
            if int(i) in list(jj.RID):
                # check if the information of the subject are present in the csv_file
                path_to_T1_1 = os.path.join(path_to_dataset, str(i))
                #subdirectories = os.listdir(path_to_T1_1)
                subdirectories = listdir_nohidden(path_to_T1_1)
                for j in xrange(len(subdirectories)):
                    #check if the subdirectory can contain a T1 image
                    path_to_T1_2=find_T1_folder(subdirectories[j],path_to_T1_1)
                    if path_to_T1_2!='NaN':
                        exame_date = listdir_nohidden(path_to_T1_2) # this is the string I need to compare with the csv
                        for x in xrange(len(exame_date)):
                                #check if the corresponding session_ID can be found in the csv_file
                            session_ID = match_data(exame_date[x], i, jj)
                            if session_ID != '-4':
                                path_to_T1_3 = os.path.join(path_to_T1_2, str(exame_date[x]))
                                image_ID = listdir_nohidden(path_to_T1_3)
                                for y in xrange(len(image_ID)):
                                    # compute the final path
                                    path_to_T1 = os.path.join(path_to_T1_3, image_ID[y])
                                    sub_ID.append(i)
                                    ses_ID.append(session_ID)
                                    path_T1.append(path_to_T1)

    return [sub_ID,ses_ID,path_T1]


def find_path_to_T1_SAG(path_to_dataset,subjects_ID,sub_ID,ses_ID,path_T1):
    '''
        
        This method creates a Dataframe which contains all the paths to the T1 images which are not ADNI compliant, 
        they contain the word "SAG" in their name

        :param path_to_dataset: path to AIBL dataset
        :param subjects_ID: subjects_id in the dataset dowloaded
        :param sub_ID: the previous list (from T1_ADNI) where new subjects ID will be appended
        :param ses_ID: the previous list (from T1_ADNI) where new session ID will be appended
        :param path_T1:the previous list (from T1_ADNI) where new paths will be appended
        :return: it completes the list of all the T1 paths including all the images where we didn't find the exame-data but we can fix it with a further analysis
    '''
    import os

    for i in subjects_ID:
        subdirectory_for_subject = []
        path_to_T1_1 = os.path.join(path_to_dataset, str(i))
        #subdirectories = os.listdir(path_to_T1_1)
        subdirectories=listdir_nohidden(path_to_T1_1)
        for j in xrange(len(subdirectories)):
            #we convert only the images which are in this list and we take only one of them for subject
            if subdirectories[j] in ['MPRAGESAGISOp2ND', 'MPRAGE_SAG_ISO_p2_ND', 'MPRAGE_SAG_ISO_p2']:
                subdirectory_for_subject.append(subdirectories[j])
        if not subdirectory_for_subject:
            pass
        else:
            path_to_T1_2 = os.path.join(path_to_T1_1, subdirectory_for_subject[0])

            exame_date = listdir_nohidden(path_to_T1_2)
            if i in [342, 557]:
                session_ID = 'M54'
            else:
                session_ID = 'M00'
            if (i in sub_ID and session_ID != ses_ID[sub_ID.index(i)]) or (i not in sub_ID):
                #if for a subject in the same session we have both this image and the "ADNI" compliant we are converting the second one since the exame-date is more precise
                path_to_T1_3 = os.path.join(path_to_T1_2, str(exame_date[0]))
                image_ID = listdir_nohidden(path_to_T1_3)
                path_to_T1 = os.path.join(path_to_T1_3, image_ID[0])
                #we append the result to the list
                sub_ID.append(i)
                ses_ID.append(session_ID)
                path_T1.append(path_to_T1)

    return [sub_ID, ses_ID, path_T1]

def find_path_to_T1(path_to_dataset,path_to_csv):
    '''
        This method creates a DataFrame for the T1 images, where for each of them the subject ID, the session ID
        and the path to the image are reported
        
        :param path_to_dataset:  path to AIBL dataset 
        :param path_to_csv: path to the csv files downloaded 
        :return: pandas dataframe which contains all the paths for the T1 images, and the correisponding subject_ID and session_ID 
    '''
    import os
    import pandas

    #two csv_files contain information regarding the T1w MRI images
    mri_meta=pandas.read_csv(os.path.join(path_to_csv,"aibl_mrimeta_28-Apr-2015.csv"))
    mri_3meta=pandas.read_csv(os.path.join(path_to_csv,"aibl_mri3meta_28-Apr-2015.csv"))
    file_mri=[mri_meta,mri_3meta]
    subjects_ID = listdir_nohidden(path_to_dataset)
    #list of all the folders which correspond to the subject_ID
    #all the subjects downloaded are taken into account for the conversion, except this sample
    if '0151083' in subjects_ID:
        del subjects_ID[subjects_ID.index('0151083')]
    [sub_ID,ses_ID,path_T1]=find_path_to_T1_ADNI(file_mri,subjects_ID,path_to_dataset)
    [sub_ID, ses_ID, path_T1]=find_path_to_T1_SAG(path_to_dataset,  subjects_ID, sub_ID, ses_ID, path_T1)

    data = pandas.DataFrame({'Subjects_ID': sub_ID,
                    'Session_ID': ses_ID,
                    'Path_to_T1':path_T1})
    #data= final dataframe
    return data



