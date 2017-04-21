def check_two_dcm_folder(dicom_path, bids_folder, image_uid):
    '''
    Check if a folder contains more than one DICOM and if yes, copy the DICOM related to image id passed as parameter into
    a temporary folder called tmp_dicom_folder.


    :param dicom_path: path to the DICOM folder
    :param bids_folder: path to the BIDS folder where the dataset will be stored
    :param image_uid: image id of the fmri
    :return: the path to the original DICOM folder or the path to a temporary folder called tmp_dicom_folder where only
     the DICOM to convert is copied
    '''
    from glob import glob
    from os import path
    from shutil import copy
    import shutil
    import os

    temp_folder_name = 'tmp_dcm_folder'
    dest_path = path.join(bids_folder, temp_folder_name)

    # Check if there is more than one xml file inside the folder
    xml_list = glob(path.join(dicom_path,'*.xml*'))

    if len(xml_list) > 1:
        # Remove the precedent tmp_dcm_folder if is existing
        if os.path.exists(dest_path):
            shutil.rmtree(dest_path)
        os.mkdir(dest_path)
        dmc_to_conv = glob(path.join(dicom_path,'*'+image_uid+'.dcm*'))
        for d in dmc_to_conv:
            copy(d, dest_path)
        return dest_path
    else:
        return dicom_path




