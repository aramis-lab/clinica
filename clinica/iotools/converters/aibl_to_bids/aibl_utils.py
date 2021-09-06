# coding: utf8

"""Utils to convert AIBL dataset in BIDS."""


def listdir_nohidden(path):
    """List all the subdirectories of path except the hidden folders.

    Args:
        path (str): path whose subdirectories are needed

    Returns:
        List(str): list of all the subdirectories of path
    """
    from os import listdir

    return [result for result in listdir(path) if not result.startswith(".")]


def find_T1_folder(subdirectory, path_to_T1_1):
    """Check if the subdirectory contains a T1 image, and it return the h.

    Args:
        subdirectory (str): subdirectory: name of the folder
        path_to_T1_1 (str): [description]

    Returns:
        str: previous path to arrive to the T1 image
    """
    import os

    path_to_convert = {
        "MPRAGE_ADNI_confirmed",
        "MPRAGE",
        "MPRAGE_ADNI_confirmed_RPT",
        "MPRAGE_ADNI_confirmed_REPEATX2",
        "MPRAGE_ADNI_confirmed_repeat",
        "MPRAGE_ADNI_confirmed_REPEAT",
        "MPRAGE_ADNI_conf_REPEAT",
    }
    for j in path_to_convert:
        path = []
        # if conditions which checks if the subfolder contain a T1 image
        if j == subdirectory:
            path = os.path.join(path_to_T1_1, subdirectory)
            return path
    if path == []:
        return "NaN"  # there are no more folders which could contain T1 images


def find_T1_folder_nodata(subdirectory, path_to_T1_1):
    """

    This method checks if the subdirectory contains a T1 image, and it
    returns the path. This method differs from the find_T1_folder since for
    these folders the exame_date is not present in the clinical excel file
    and we will not check if the exame_date corresponds to the date stored
    in the path to the image, but they will be converted anyway

     :param subdirectory: name of the folder

     :return: previous path to arrive to the T1 image
    """
    import os

    path_to_convert = {"MPRAGESAGISOp2ND", "MPRAGE_SAG_ISO_p2_ND", "MPRAGE_SAG_ISO_p2"}
    for j in path_to_convert:
        path = []
        # if conditions which checks if the subfolder contain a T1 image
        if j == subdirectory:
            path = os.path.join(path_to_T1_1, subdirectory)
            return path
    if path == []:
        return "NaN"  # there are no more folders which could contain T1 images


def find_correspondance_index(i, csv_file):
    """Give as output the index of the CSV file analysed which correspond to the 'i' subject.

    :param i: subject_ID
    :param csv_file: CSV file where all the information are listed

    :return: index
    """
    index = []
    for x in csv_file.RID:
        if i == str(x):
            index = csv_file.RID[csv_file.RID == x].index.tolist()
            return index


def find_correspondance_date(index, csv_file):
    """Return the dates reported in the csv_file for the i-subject.

    :param index: index corresponding to the subject analysed
    :param csv_file: csv file where all the information are listed
    :return date
    """
    return csv_file.EXAMDATE[index]


def match_data(exame_date, i, csv_file):
    """

    This method returns the session_ID. It controls if the dates
    corresponding to the image (from the name of the subdirectory)
    correspond to one of the dates listed from the csv_file for the subject
    analysed. The session_ID is the corresponding session for that patient
    in that date.  It returns -4 if there are no information.

    :param exame_date: date where the image has been taken, it is saved
    from the name of the corresponding subdirector
    :param i: subject_ID
    :param csv_file: csv file where all the information are listed

    :return session_id of the patient
    """
    import re

    session_ID = []
    index = find_correspondance_index(i, csv_file)
    csv_date = find_correspondance_date(index, csv_file)
    for xx in index:
        if str(csv_date[xx]) != "-4":
            # check is the date is not '-4'
            m = re.search(
                "([0-9].*)-(.*)-(.*)_(.*)_(.*)_(.*)", str(exame_date)
            )  # string from image directory
            p = re.search(
                "(.*)/(.*)/(.*)", str(csv_date[xx])
            )  # string from the date of the csv_file
            if (
                (p.group(1) == m.group(2))
                & (p.group(2) == m.group(3))
                & (p.group(3) == m.group(1))
            ):
                session_ID = csv_file.VISCODE[xx]
    if session_ID == []:
        session_ID = "-4"
    return session_ID


def list_of_paths():
    """List all the folders which not contain PET images."""
    return [
        ".DS_Store",
        "localizer",
        "Space_3D_T2_FLAIR_sag_p2",
        "AXIAL_FLAIR",
        "MPRAGE_ADNI_confirmed_REPEATX2",
        "Axial_PD-T2_TSE",
        "Axial_PD-T2_TSE_repeat",
        "MPRAGE_SAG_ISO_p2_ND",
        "Axial_PD-T2_TSE_confirmed",
        "MPRAGESAGISOp2ND",
        "MPRAGE_ADNI_confirmed",
        "MPRAGE_ADNI_confirmed_repeat",
        "MPRAGE_SAG_ISO_p2",
        "MPRAGE",
        "MPRAGE_ADNI_confirmed_REPEAT",
        "Axial_PD-T2_TSE_confirmed_repeat",
        "MPRAGE_ADNI_conf_REPEAT",
        "Space_3D_T2_FLAIR_sag_p2_REPEAT",
        "MPRAGE_ADNI_confirmed_RPT",
        "Brain_256_1.6_zoom_4_x_4_iter",
        "Space_3D_T2_FLAIR_sag_REPEAT",
        "Axial_PD-T2_TSE_RPTconfirmed",
        "Axial_PD-T2_TSE_RPT_confirmed",
        "Axial_PD-T2_TSE_confirmed_REPEAT",
        "flair_t2_spc_irprep_ns_sag_p2_1mm_iso",
        "localiser",
    ]


def check_subdirectories_pet(subdirectories, sub, no_pet):
    """

    It returns the correct subdirectories for the PET images, they should
    belong to the list where there all the possible names of the PET images.

    :param subdirectories:
    :param sub: all the possible subdirectories which need to be checked
    :param no pet: list of names of folders which not contain PET images

    :return subdirectory which is containing a PET image which needs to be
    converted
    """

    for j in range(len(sub)):
        if (sub[j] not in no_pet) & (sub[j] != ".DS_Store"):
            subdirectories.append(sub[j])
    subdirectories = list(set(subdirectories))
    return subdirectories


def dicom_to_nii(subject, output_path, output_filename, image_path):
    """Convert the dicom images in a nifti files using dicom2nii or mri_convert.

    :param subject:
    :param output_path: where nifti image is stored
    :param output_filename: name of the nifti image
    :param image_path: where dicom files are stored

    :return: Image in a nifti format
    """
    import os
    import shutil
    import subprocess
    from os.path import exists

    from clinica.iotools.bids_utils import run_dcm2niix
    from clinica.utils.stream import cprint

    try:
        os.makedirs(output_path)
    except OSError:
        if not os.path.isdir(output_path):
            raise

    # if image.Is_Dicom:
    command = f"dcm2niix -b n -z y -o {output_path} -f {output_filename} {image_path}"
    run_dcm2niix(command)

    nifti_file = os.path.join(output_path, output_filename + ".nii.gz")

    if not exists(nifti_file):
        cprint(nifti_file + " should have been created but this did not happen")
    return nifti_file


def viscode_to_session(viscode):
    """Replace the session label 'bl' with 'M00' or capitalize the session name passed as input.

    :param viscode: session name

    :return: M00 if is the baseline session or the original session name
    capitalized
    """
    if viscode == "bl":
        return "M00"
    else:
        return viscode.capitalize()


def find_path_to_pet_modality(path_to_dataset, csv_file):
    """Create a Dataframe which contains all the paths to the PET image of a modality (for example AV45 or PIB).

    :param path_to_dataset: path to AIBL dataset
    :param csv_file: file which correspond to the modality

    :return: A dataframe which contains the path for PET images for a
    single modality and subject_ID and session_ID are reported for each
    path
    """
    import os

    import pandas

    # TODO
    # exclude_subjects = get_exclude_subject(file.txt)

    no_pet = list_of_paths()
    subjects_ID = listdir_nohidden(path_to_dataset)
    # selection of the subjects_ID from the folder downloaded
    # this subject must be discarded since it is only a sample and not a patient
    if "0151083" in subjects_ID:
        del subjects_ID[subjects_ID.index("0151083")]
    sub_ID = []
    ses_ID = []
    path_pet = []

    # Iteration through all the subjects_ID
    def is_int(x):
        for i in x:
            if int(i) in list(csv_file.RID):
                yield i

    #    def append_path(image_ID):

    for i in is_int(subjects_ID):
        # check if the subject is present in the csv_file for the modality selected
        subdirectories = []
        path_to_pet_1 = os.path.join(path_to_dataset, str(i))
        # subdirectory_all = os.listdir(path_to_pet_1)
        subdirectory_all = listdir_nohidden(path_to_pet_1)
        subdirectories = check_subdirectories_pet(
            subdirectories, subdirectory_all, no_pet
        )
        # selection only of the folders which contain PET image
        for j in range(len(subdirectories)):
            path_to_pet_2 = os.path.join(path_to_pet_1, subdirectories[j])
            exame_date = listdir_nohidden(path_to_pet_2)
            # exame date of the image which is going to be converted
            for x in range(len(exame_date)):
                # selection of the session_ID matching the data in the csv_file with the one of the image
                session_ID = match_data(exame_date[x], i, csv_file)
                if session_ID != "-4":
                    path_to_pet_3 = os.path.join(path_to_pet_2, str(exame_date[x]))
                    # For the RID 1607 there are two PET images of the flute modality, and we select the first
                    if i == "1607":
                        if subdirectories[j] == "Flute_256_1.6_Zoom_plain_4_x_4_Iter":
                            image_ID = ["I442930"]
                        else:
                            image_ID = listdir_nohidden(path_to_pet_3)
                    else:
                        image_ID = listdir_nohidden(path_to_pet_3)

                    for y in range(len(image_ID)):
                        # final path to find the image we want to convert
                        path_to_pet = os.path.join(path_to_pet_3, image_ID[y])  #
                        sub_ID.append(i)
                        ses_ID.append(session_ID)
                        path_pet.append(path_to_pet)

    data = pandas.DataFrame(
        {"Subjects_ID": sub_ID, "Session_ID": ses_ID, "Path_to_pet": path_pet}
    )
    # data=final dataframe

    return data


def find_path_to_T1_ADNI(file_mri, subjects_ID, path_to_dataset):
    """

    This method creates a Dataframe which contains all the paths to the T1
    images which are ADNI compliant (as explained in the AIBL website).
    This images differ from the others T1 of the dataset since in the
    cvs_file is reported the exame date.

    :param file_mri: in the clinical data there are two files which
    describe the  parameters of the T1 images (MRI 1.5 T and MRI 3T)
    :param subjects_ID: subjects_id in the dataset dowloaded
    :param path_to_dataset: path to AIBL dataset

    :return: A dataframe which contains the path for T1 images and
    subject_ID and session_ID are reported for each path
    """
    import os

    sub_ID = []
    ses_ID = []
    path_T1 = []

    for i in subjects_ID:
        for jj in file_mri:
            # it checks all the file_mri
            if int(i) in list(jj.RID):
                # check if the information of the subject are present in the csv_file
                path_to_T1_1 = os.path.join(path_to_dataset, str(i))
                # subdirectories = os.listdir(path_to_T1_1)
                subdirectories = listdir_nohidden(path_to_T1_1)
                for j in range(len(subdirectories)):
                    # check if the subdirectory can contain a T1 image
                    path_to_T1_2 = find_T1_folder(subdirectories[j], path_to_T1_1)
                    if path_to_T1_2 != "NaN":
                        exame_date = listdir_nohidden(
                            path_to_T1_2
                        )  # this is the string I need to compare with the csv
                        for x in range(len(exame_date)):
                            # check if the corresponding session_ID can be found in the csv_file
                            session_ID = match_data(exame_date[x], i, jj)
                            if session_ID != "-4":
                                path_to_T1_3 = os.path.join(
                                    path_to_T1_2, str(exame_date[x])
                                )
                                image_ID = listdir_nohidden(path_to_T1_3)
                                for y in range(len(image_ID)):
                                    # compute the final path
                                    path_to_T1 = os.path.join(path_to_T1_3, image_ID[y])
                                    sub_ID.append(i)
                                    ses_ID.append(session_ID)
                                    path_T1.append(path_to_T1)

    return [sub_ID, ses_ID, path_T1]


def find_path_to_T1_SAG(path_to_dataset, subjects_ID, sub_ID, ses_ID, path_T1):
    """

    This method creates a Dataframe which contains all the paths to the T1
    images which are not ADNI compliant, they contain the word "SAG" in
    their name

    :param path_to_dataset: path to AIBL dataset
    :param subjects_ID: subjects_id in the dataset dowloaded
    :param sub_ID: the previous list (from T1_ADNI) where new subjects ID
    will be appended
    :param ses_ID: the previous list (from T1_ADNI) where new session ID
    will be appended
    :param path_T1:the previous list (from T1_ADNI) where new paths will be
    appended

    :return: it completes the list of all the T1 paths including all the
    images where we didn't find the exame-data but we can fix it with a
    further analysis
    """
    import os

    for i in subjects_ID:
        subdirectory_for_subject = []
        path_to_T1_1 = os.path.join(path_to_dataset, str(i))
        # subdirectories = os.listdir(path_to_T1_1)
        subdirectories = listdir_nohidden(path_to_T1_1)
        for j in range(len(subdirectories)):
            # we convert only the images which are in this list and we take only one of them for subject
            if subdirectories[j] in [
                "MPRAGESAGISOp2ND",
                "MPRAGE_SAG_ISO_p2_ND",
                "MPRAGE_SAG_ISO_p2",
            ]:
                subdirectory_for_subject.append(subdirectories[j])
        if not subdirectory_for_subject:
            pass
        else:
            path_to_T1_2 = os.path.join(path_to_T1_1, subdirectory_for_subject[0])

            exame_date = listdir_nohidden(path_to_T1_2)
            if i in [342, 557]:
                session_ID = "M54"
            else:
                session_ID = "M00"
            if (i in sub_ID and session_ID != ses_ID[sub_ID.index(i)]) or (
                i not in sub_ID
            ):
                # if for a subject in the same session we have both this image and the "ADNI" compliant we are converting the second one since the exame-date is more precise
                path_to_T1_3 = os.path.join(path_to_T1_2, str(exame_date[0]))
                image_ID = listdir_nohidden(path_to_T1_3)
                path_to_T1 = os.path.join(path_to_T1_3, image_ID[0])
                # we append the result to the list
                sub_ID.append(i)
                ses_ID.append(session_ID)
                path_T1.append(path_to_T1)

    return [sub_ID, ses_ID, path_T1]


def find_path_to_T1(path_to_dataset, path_to_csv):
    """
    This method creates a DataFrame for the T1 images, where for each of
    them the subject ID, the session ID and the path to the image are
    reported

    :param path_to_dataset:  path to AIBL dataset
    :param path_to_csv: path to the csv files downloaded
    :return: pandas dataframe which contains all the paths for the T1
    images, and the correisponding subject_ID and session_ID
    """
    import glob
    import os

    import pandas

    # two csv_files contain information regarding the T1w MRI images
    mri_meta = pandas.read_csv(
        glob.glob(os.path.join(path_to_csv, "aibl_mrimeta_*.csv"))[0]
    )
    mri_3meta = pandas.read_csv(
        glob.glob(os.path.join(path_to_csv, "aibl_mri3meta_*.csv"))[0]
    )
    file_mri = [mri_meta, mri_3meta]
    subjects_ID = listdir_nohidden(path_to_dataset)
    # list of all the folders which correspond to the subject_ID
    # all the subjects downloaded are taken into account for the conversion, except this sample
    if "0151083" in subjects_ID:
        del subjects_ID[subjects_ID.index("0151083")]
    [sub_ID, ses_ID, path_T1] = find_path_to_T1_ADNI(
        file_mri, subjects_ID, path_to_dataset
    )
    [sub_ID, ses_ID, path_T1] = find_path_to_T1_SAG(
        path_to_dataset, subjects_ID, sub_ID, ses_ID, path_T1
    )

    data = pandas.DataFrame(
        {"Subjects_ID": sub_ID, "Session_ID": ses_ID, "Path_to_T1": path_T1}
    )
    # data= final dataframe
    return data


# Covert the AIBL PET images into the BIDS specification.
# There are three pet modalities: av45, pib, flute. All of them are converted
# in BIDS


def paths_to_bids(path_to_dataset, path_to_csv, bids_dir, modality, overwrite=False):
    """Convert all the T1 images found in the AIBL dataset downloaded in BIDS.

    Args:
        path_to_dataset: path_to_dataset
        path_to_csv: path to the csv file containing clinical data
        bids_dir: path to save the AIBL-T1-dataset converted in a BIDS format
        modality: string 't1', 'av45', 'flute' or 'pib'
        overwrite: if True previous existing outputs will be erased

    Returns:
        list of all the images that are potentially converted in a BIDS format and saved in the bids_dir.
        This does not guarantee existence.
    """
    import glob
    from multiprocessing import Value, cpu_count
    from multiprocessing.dummy import Pool
    from os import makedirs, remove
    from os.path import exists, join

    import pandas as pds
    from numpy import nan

    from clinica.iotools.bids_utils import json_from_dcm
    from clinica.iotools.utils.data_handling import center_nifti_origin
    from clinica.utils.stream import cprint

    if modality.lower() not in ["t1", "av45", "flute", "pib"]:
        # This should never be reached
        raise RuntimeError(f"{modality.lower()} is not supported for conversion")

    counter = None

    makedirs(join(bids_dir, "conversion_info"), exist_ok=True)

    def init(args):
        """Store the counter for later use."""
        global counter
        counter = args

    def create_file(image):
        global counter
        subject = image.Subjects_ID
        session = image.Session_ID
        name_of_path = {
            "t1": "Path_to_T1",
            "av45": "Path_to_pet",
            "flute": "Path_to_pet",
            "pib": "Path_to_pet",
        }
        # depending on the DataFrame, there is different way of accessing
        # the image object
        image_path = image[name_of_path[modality]]
        with counter.get_lock():
            counter.value += 1
            if image_path is nan:
                cprint(f"No path specified for {subject} in session {session}")
            else:
                cprint(
                    f"'[{modality.upper()}] Processing subject {subject} - {session} {session}, {counter.value}/{total}"
                )

        if image_path is nan:
            return nan
        session = viscode_to_session(session)
        # creation of the path
        if modality == "t1":
            output_path = join(bids_dir, f"sub-AIBL{subject}", f"ses-{session}", "anat")
            output_filename = f"sub-AIBL{subject}_ses-{session}_T1w"
        elif modality in ["flute", "pib", "av45"]:
            output_path = join(bids_dir, f"sub-AIBL{subject}", f"ses-{session}", "pet")
            output_filename = (
                f"sub-AIBL{subject}_ses-{session}_task-rest_acq-{modality}_pet"
            )
        # image is saved following BIDS specifications

        if exists(join(output_path, output_filename + ".nii.gz")) and not overwrite:
            cprint(f"Subject {str(subject)} - session {session} already processed.")
            output_image = join(output_path, output_filename + ".nii.gz")
        else:
            if exists(join(output_path, output_filename + ".nii.gz")):
                remove(join(output_path, output_filename + ".nii.gz"))
            output_image = dicom_to_nii(
                subject, output_path, output_filename, image_path
            )
            json_from_dcm(image_path, join(output_path, output_filename + ".json"))

        # Center all images
        center_nifti_origin(output_image, output_image)

        return output_image

    # it reads the DataFrame where subject_ID, session_ID and path are saved
    if modality == "t1":
        images = find_path_to_T1(path_to_dataset, path_to_csv)
    else:
        path_to_csv_pet_modality = glob.glob(
            join(path_to_csv, "aibl_" + modality + "meta_*.csv")
        )[0]
        if not exists(path_to_csv_pet_modality):
            raise FileNotFoundError(
                path_to_csv_pet_modality + " file not found in clinical data folder"
            )
        # Latest version of Flutemetamol CSV file (aibl_flutemeta_01-Jun-2018.csv)
        # has an extra column for some rows. However, each CSV file (regarding PET tracers)
        # contains the same columns. The usecols fixes this issue.
        df_pet = pds.read_csv(
            path_to_csv_pet_modality, sep=",|;", usecols=list(range(0, 36))
        )
        images = find_path_to_pet_modality(path_to_dataset, df_pet)
    images.to_csv(
        join(bids_dir, "conversion_info", modality + "_paths.tsv"),
        index=False,
        sep="\t",
        encoding="utf-8",
    )

    counter = Value("i", 0)
    total = images.shape[0]
    # Reshape inputs to give it as a list to the workers
    images_list = []
    for i in range(total):
        images_list.append(images.iloc[i])

    # intializer are used with the counter variable to keep track of how many
    # files have been processed
    poolrunner = Pool(cpu_count(), initializer=init, initargs=(counter,))
    output_file_treated = poolrunner.map(create_file, images_list)
    del counter
    return output_file_treated


# -- Methods for the clinical data --


def create_participants_df_AIBL(
    input_path, clinical_spec_path, clinical_data_dir, delete_non_bids_info=True
):
    """Create a participants file for the AIBL dataset where information regarding the patients are reported.

    Args:
        input_path: path to the input directory
        clinical_spec_path: path to the clinical file
        clinical_data_dir: directory to the clinical data files
        delete_non_bids_info: if True delete all the rows of the subjects
        that are not available in the BIDS dataset

    Returns:
        a pandas DataFrame that contains the participants data
    """
    import glob
    import os
    import re
    from os import path

    import numpy as np
    import pandas as pd

    fields_bids = ["participant_id"]
    fields_dataset = []
    prev_location = ""
    index_to_drop = []

    location_name = "AIBL location"
    clinical_spec_path = clinical_spec_path + "_participant.tsv"

    if not os.path.exists(clinical_spec_path):
        raise FileNotFoundError(clinical_spec_path + " not found in clinical data.")
    participants_specs = pd.read_csv(clinical_spec_path, sep="\t")
    participant_fields_db = participants_specs["AIBL"]
    field_location = participants_specs[location_name]
    participant_fields_bids = participants_specs["BIDS CLINICA"]

    # Extract the list of the available fields for the dataset (and the corresponding BIDS version)
    for i in range(0, len(participant_fields_db)):
        if not pd.isnull(participant_fields_db[i]):
            fields_bids.append(participant_fields_bids[i])
            fields_dataset.append(participant_fields_db[i])

    # Init the dataframe that will be saved in the file participant.tsv
    participant_df = pd.DataFrame(columns=fields_bids)

    csv_files = []
    for i in range(0, len(participant_fields_db)):
        # If a field not empty is found
        if not pd.isnull(participant_fields_db[i]):
            # Extract the file location of the field and read the value from the file
            tmp = field_location[i].split("/")
            location = tmp[0]
            # If a sheet is available
            if len(tmp) > 1:
                sheet = tmp[1]
            else:
                sheet = ""
            # Check if the file to open for a certain field it's the same of the previous field
            if location == prev_location and sheet == prev_sheet:
                pass
            else:
                file_ext = os.path.splitext(location)[1]
                file_to_read_path = path.join(clinical_data_dir, location)

                if file_ext == ".xlsx":
                    file_to_read = pd.read_excel(
                        glob.glob(file_to_read_path)[0], sheet_name=sheet
                    )
                elif file_ext == ".csv":
                    file_to_read = pd.read_csv(glob.glob(file_to_read_path)[0])
                prev_location = location
                prev_sheet = sheet

            field_col_values = []
            # For each field in fields_dataset extract all the column values
            for j in range(0, len(file_to_read)):
                # Convert the alternative_id_1 to string if is an integer/float
                if participant_fields_bids[i] == "alternative_id_1" and (
                    file_to_read[participant_fields_db[i]].dtype == np.float64
                    or file_to_read[participant_fields_db[i]].dtype == np.int64
                ):
                    if not pd.isnull(file_to_read.at[j, participant_fields_db[i]]):
                        # value_to_append = str(file_to_read.get_value(j, participant_fields_db[i])).rstrip('.0')
                        value_to_append = str(
                            file_to_read.at[j, participant_fields_db[i]]
                        )

                    else:
                        value_to_append = "n/a"
                else:
                    value_to_append = file_to_read.at[j, participant_fields_db[i]]
                field_col_values.append(value_to_append)
            # Add the extracted column to the participant_df
            participant_df[participant_fields_bids[i]] = pd.Series(field_col_values)

    # Adding participant_id column with BIDS ids
    for i in range(0, len(participant_df)):
        value = participant_df["alternative_id_1"][i]
        participant_df["participant_id"][i] = "sub-AIBL" + str(value)
        participant_df["date_of_birth"][i] = re.search(
            "/([0-9].*)", str(participant_df["date_of_birth"][i])
        ).group(1)
        if participant_df["sex"][i] == 1:
            participant_df["sex"][i] = "M"
        else:
            participant_df["sex"][i] = "F"

    participant_df.replace("-4", "n/a")

    # Delete all the rows of the subjects that are not available in the BIDS dataset
    if delete_non_bids_info:
        participant_df = participant_df.drop(index_to_drop)

    participant_df.to_csv(
        os.path.join(input_path, "participants.tsv"),
        sep="\t",
        index=False,
        encoding="utf8",
    )

    return participant_df


def create_sessions_dict_AIBL(input_path, clinical_data_dir, clinical_spec_path):
    """Extract the information regarding the sessions and store them in a dictionary.

    :param input_path: path to the input folder
    :param clinical_spec_path: path to the clinical file
    :param clinical_data_dir: directory to the clinical data files
    :return: A dataframe saved in a tsv file which contains information for each session
    """
    import glob
    from os import path

    import pandas as pd

    # Load data
    location = "AIBL location"
    sessions = pd.read_csv(clinical_spec_path + "_sessions.tsv", sep="\t")
    sessions_fields = sessions["AIBL"]
    field_location = sessions[location]
    sessions_fields_bids = sessions["BIDS CLINICA"]
    fields_dataset = []
    fields_bids = []
    sessions_dict = {}

    for i in range(0, len(sessions_fields)):
        if not pd.isnull(sessions_fields[i]):
            fields_bids.append(sessions_fields_bids[i])
            fields_dataset.append(sessions_fields[i])

    sessions_df = pd.DataFrame(columns=fields_bids)

    files_to_read = []
    sessions_fields_to_read = []
    for i in range(0, len(sessions_fields)):
        # If the i-th field is available
        if not pd.isnull(sessions_fields[i]):
            # Load the file
            tmp = field_location[i]
            location = tmp[0]
            file_to_read_path = path.join(clinical_data_dir, tmp)
            files_to_read.append(glob.glob(file_to_read_path)[0])
            sessions_fields_to_read.append(sessions_fields[i])

    rid = pd.read_csv(files_to_read[0], dtype={"text": str}, low_memory=False).RID
    rid = list(set(rid))
    for r in rid:
        dict = []
        for i in files_to_read:
            file_to_read = pd.read_csv(i, dtype={"text": str})
            if len(file_to_read.columns) == 1:
                file_to_read = pd.read_csv(i, sep=";", low_memory=False)

            # information are written following the BIDS specifications
            viscode = file_to_read.loc[(file_to_read["RID"] == r), "VISCODE"]
            for j in sessions_fields_to_read:
                if j in list(file_to_read.columns.values) and j == "MMSCORE":
                    MMSCORE = file_to_read.loc[(file_to_read["RID"] == r), j]
                    MMSCORE[MMSCORE == -4] = "n/a"
                elif j in list(file_to_read.columns.values) and j == "CDGLOBAL":
                    CDGLOBAL = file_to_read.loc[(file_to_read["RID"] == r), j]
                    CDGLOBAL[CDGLOBAL == -4] = "n/a"
                elif j in list(file_to_read.columns.values) and j == "DXCURREN":
                    DXCURREN = file_to_read.loc[(file_to_read["RID"] == r), j]
                    DXCURREN[DXCURREN == -4] = "n/a"
                    DXCURREN[DXCURREN == 1] = "CN"
                    DXCURREN[DXCURREN == 2] = "MCI"
                    DXCURREN[DXCURREN == 3] = "AD"
                elif j in list(file_to_read.columns.values) and j == "EXAMDATE":
                    EXAMDATE = file_to_read.loc[(file_to_read["RID"] == r), j]
                elif j in list(file_to_read.columns.values) and j == "PTDOB":
                    PTDOB = file_to_read.loc[(file_to_read["RID"] == r), j]

        examdates = get_examdates(
            r, EXAMDATE.to_list(), viscode.to_list(), clinical_data_dir
        )
        age = get_ages(PTDOB.values[0], examdates)

        viscode[viscode == "bl"] = "M00"
        viscode = viscode.str.upper()

        dict = pd.DataFrame(
            {
                "session_id": "ses-" + viscode,
                "age": age,
                "MMS": MMSCORE,
                "cdr_global": CDGLOBAL,
                "diagnosis": DXCURREN,
                "examination_date": examdates,
            }
        )

        cols = dict.columns.tolist()
        dict = dict[cols[-1:] + cols[:-1]]

        bids_paths = path.join(input_path, "sub-AIBL" + str(r))
        if path.exists(bids_paths):
            dict.to_csv(
                path.join(
                    input_path,
                    "sub-AIBL" + str(r),
                    "sub-AIBL" + str(r) + "_sessions.tsv",
                ),
                sep="\t",
                index=False,
                encoding="utf8",
            )


def create_scans_dict_AIBL(input_path, clinical_data_dir, clinical_spec_path):
    """Create scans.tsv files for AIBL.

    Args:
        input_path: path to the input folder
        clinical_spec_path: path to the clinical file
        clinical_data_dir: directory to the clinical data files
    """
    import glob
    from os import path

    import pandas as pd

    import clinica.iotools.bids_utils as bids

    # Load data
    location = "AIBL location"
    scans = pd.read_csv(clinical_spec_path + "_scans.tsv", sep="\t")
    scans_fields = scans["AIBL"]
    field_location = scans[location]
    scans_fields_bids = scans["BIDS CLINICA"]
    fields_dataset = []
    fields_bids = []

    # Keep only fields for which there are AIBL fields
    for i in range(0, len(scans_fields)):
        if not pd.isnull(scans_fields[i]):
            fields_bids.append(scans_fields_bids[i])
            fields_dataset.append(scans_fields[i])

    files_to_read = []
    sessions_fields_to_read = []
    for i in range(0, len(scans_fields)):
        # If the i-th field is available
        if not pd.isnull(scans_fields[i]):
            # Load the file
            tmp = field_location[i]
            file_to_read_path = path.join(clinical_data_dir, tmp)
            files_to_read.append(glob.glob(file_to_read_path)[0])
            sessions_fields_to_read.append(scans_fields[i])

    bids_ids = [
        path.basename(sub_path)
        for sub_path in glob.glob(path.join(input_path, "sub-AIBL*"))
    ]

    # This dictionary should be automatically computed from the dataset
    #
    ses_dict = bids.get_sessions_map_AIBL(bids_ids, input_path)

    scans_dict = bids.create_scans_dict(
        clinical_data_dir,
        "AIBL",
        clinical_spec_path,
        bids_ids,
        "RID",
        "VISCODE",
        ses_dict,
    )
    bids.write_scans_tsv(input_path, bids_ids, scans_dict)


def get_examdates(rid, examdates, viscodes, clinical_data_dir):
    import glob
    from datetime import datetime
    from os import path

    import pandas as pd
    from dateutil.relativedelta import relativedelta

    res_examdates = []
    csv_list = [
        glob.glob(path.join(clinical_data_dir, "aibl_mri3meta_*.csv"))[0],
        glob.glob(path.join(clinical_data_dir, "aibl_mrimeta_*.csv"))[0],
        glob.glob(path.join(clinical_data_dir, "aibl_cdr_*.csv"))[0],
        glob.glob(path.join(clinical_data_dir, "aibl_flutemeta_*.csv"))[0],
        glob.glob(path.join(clinical_data_dir, "aibl_mmse_*.csv"))[0],
        glob.glob(path.join(clinical_data_dir, "aibl_pibmeta_*.csv"))[0],
    ]

    for e in range(len(examdates)):
        exam = examdates[e]

        if exam != "-4":
            res_examdates.append(exam)
            continue

        # If EXAMDATE does not exist (-4) we try to obtain it from another .csv file
        for csv_file in csv_list:
            if "aibl_flutemeta" in csv_file:
                csv_data = pd.read_csv(
                    csv_file, low_memory=False, usecols=list(range(0, 36))
                )
            else:
                csv_data = pd.read_csv(csv_file, low_memory=False)
            exam_date = csv_data[
                (csv_data.RID == rid) & (csv_data.VISCODE == viscodes[e])
            ]
            if not exam_date.empty and exam_date.iloc[0].EXAMDATE != "-4":
                exam = exam_date.iloc[0].EXAMDATE
                break

        # If EXAMDATE still does not exist (-4) we add the session months to baseline date
        if exam == "-4":
            bl_index = viscodes.index("bl")
            if bl_index > -1:
                bl_date = examdates[bl_index]
                bl_examdate = datetime.strptime(bl_date, "%m/%d/%Y")
                if viscodes[e] != "bl":
                    months = int(viscodes[e][1:])
                    examdate = bl_examdate + relativedelta(months=+months)
                    exam = examdate.strftime("%m/%d/%Y")

        if exam == "-4":
            print(f"No EXAMDATE for subject %{rid}, at session {viscodes[e]}")

        res_examdates.append(exam)

    return res_examdates


def get_ages(pt_dob, examdates):
    """Calculate age as time passed by since DOB to EXAMDATE.

    :param pt_dob: string - Date of birth of patient ("/%Y" format)
    :param examdates: list - Exam dates ("%m/%d/%Y" format)
    :return: list - Age at each exam date
    """
    from datetime import datetime

    age = []
    dob = datetime.strptime(pt_dob, "/%Y")

    for exam in examdates:
        examdate = datetime.strptime(exam, "%m/%d/%Y")
        delta = examdate - dob
        age.append(round(delta.days / 365.25, 1))

    return age
