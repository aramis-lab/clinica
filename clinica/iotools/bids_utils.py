# coding: utf-8

"""Methods used by BIDS converters."""


# -- Methods for the clinical data --
# @ToDo:test this function
def create_participants_df(
    study_name,
    clinical_spec_path,
    clinical_data_dir,
    bids_ids,
    delete_non_bids_info=True,
):
    """Create the file participants.tsv.

    Args:
        study_name: name of the study (Ex. ADNI)
        clinical_spec_path: path to the clinical file
        clinical_data_dir: path to the directory where the clinical data are stored
        bids_ids: list of bids ids
        delete_non_bids_info: if True delete all the rows of the subjects that
        are not available in the BIDS dataset

    Returns: a pandas dataframe that contains the participants data
    """
    import os
    from os import path

    import numpy as np
    import pandas as pd

    from clinica.utils.stream import cprint

    fields_bids = ["participant_id"]
    prev_location = ""
    index_to_drop = []
    subjects_to_drop = []
    location_name = study_name + " location"

    # Load the data from the clincal specification file
    participants_specs = pd.read_csv(clinical_spec_path + "_participant.tsv", sep="\t")
    participant_fields_db = participants_specs[study_name]
    field_location = participants_specs[location_name]
    participant_fields_bids = participants_specs["BIDS CLINICA"]

    # Extract the list of the available BIDS fields for the dataset
    for i in range(0, len(participant_fields_db)):
        if not pd.isnull(participant_fields_db[i]):
            fields_bids.append(participant_fields_bids[i])

    # Init the dataframe that will be saved in the file participants.tsv
    participant_df = pd.DataFrame(columns=fields_bids)

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
            # Check if the file to open for a certain field is the same of the previous field
            if location == prev_location and sheet == prev_sheet:
                pass
            else:
                file_ext = os.path.splitext(location)[1]
                file_to_read_path = path.join(clinical_data_dir, location)

                if file_ext == ".xlsx":
                    file_to_read = pd.read_excel(file_to_read_path, sheet_name=sheet)
                elif file_ext == ".csv":
                    file_to_read = pd.read_csv(file_to_read_path)
                prev_location = location
                prev_sheet = sheet

            field_col_values = []
            # For each field in fields_dataset extract all the column values
            for j in range(0, len(file_to_read)):
                # Convert the alternative_id_1 to string if is an integer/float
                value_to_read = file_to_read[participant_fields_db[i]]
                if participant_fields_bids[i] == "alternative_id_1" and (
                    value_to_read.dtype == np.float64 or value_to_read.dtype == np.int64
                ):
                    if not pd.isnull(file_to_read.at[j, participant_fields_db[i]]):
                        value_to_append = str(
                            file_to_read.at[j, participant_fields_db[i]]
                        ).rstrip(".0")
                    else:
                        value_to_append = np.NaN
                else:
                    value_to_append = file_to_read.at[j, participant_fields_db[i]]
                field_col_values.append(value_to_append)
            # Add the extracted column to the participant_df
            participant_df[participant_fields_bids[i]] = pd.Series(field_col_values)

    if study_name == "ADNI" or study_name == "AIBL":
        # ADNImerge contains one row for each visits so there are duplicates
        participant_df = participant_df.drop_duplicates(
            subset=["alternative_id_1"], keep="first"
        )
    elif study_name == "OASIS":
        # OASIS provides several MRI for the same session
        participant_df = participant_df[
            ~participant_df.alternative_id_1.str.endswith("_MR2")
        ]
    participant_df.reset_index(inplace=True, drop=True)

    # Adding participant_id column with BIDS ids
    for i in range(0, len(participant_df)):
        if study_name == "OASIS":
            value = (participant_df["alternative_id_1"][i].split("_"))[1]
        elif study_name == "OASIS3":
            value = participant_df["alternative_id_1"][i].replace("OAS3", "")
        else:
            value = remove_space_and_symbols(participant_df["alternative_id_1"][i])

        bids_id = [s for s in bids_ids if value in s]

        if len(bids_id) == 0:
            index_to_drop.append(i)
            subjects_to_drop.append(value)
        else:
            participant_df.at[i, "participant_id"] = bids_id[0]

    if len(subjects_to_drop) > 0:
        cprint(
            msg=(
                "The following subjects of dataset directory were not found in your BIDS folder :\n"
                + ", ".join(subjects_to_drop)
            ),
            lvl="info",
        )
    # Delete all the rows of the subjects that are not available in the BIDS dataset
    if delete_non_bids_info:
        participant_df = participant_df.drop(index_to_drop)

    participant_df = participant_df.fillna("n/a")

    return participant_df


def create_sessions_dict(
    clinical_data_dir,
    study_name,
    clinical_spec_path,
    bids_ids,
    name_column_ids,
    subj_to_remove=[],
):
    """Extract the information regarding the sessions and store them in a dictionary (session M00 only).

    Args:
        clinical_data_dir: path to the input folder
        study_name: name of the study (Ex: ADNI)
        clinical_spec_path:  path to the clinical file
        bids_ids: list of bids ids
        name_column_ids: name of the column where the subject ids are stored
        subj_to_remove: subjects to remove
    """
    import os
    from os import path

    import numpy as np
    import pandas as pd

    # Load data
    location = study_name + " location"
    sessions = pd.read_csv(clinical_spec_path + "_sessions.tsv", sep="\t")
    sessions_fields = sessions[study_name]
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

    for i in range(0, len(sessions_fields)):
        # If the i-th field is available
        if not pd.isnull(sessions_fields[i]):
            # Load the file
            tmp = field_location[i].split("/")
            location = tmp[0]
            if len(tmp) > 1:
                sheet = tmp[1]
            else:
                sheet = ""

            file_to_read_path = path.join(clinical_data_dir, location)
            file_ext = os.path.splitext(location)[1]

            if file_ext == ".xlsx":
                file_to_read = pd.read_excel(file_to_read_path, sheet_name=sheet)
            elif file_ext == ".csv":
                file_to_read = pd.read_csv(file_to_read_path)

            for r in range(0, len(file_to_read.values)):
                row = file_to_read.iloc[r]
                # Extracts the subject ids columns from the dataframe
                subj_id = row[name_column_ids]
                if hasattr(subj_id, "dtype"):
                    if subj_id.dtype == np.int64:
                        subj_id = str(subj_id)
                # Removes all the - from
                subj_id_alpha = remove_space_and_symbols(subj_id)
                if study_name == "OASIS":
                    subj_id_alpha = str(subj_id[0:3] + "IS" + subj_id[3] + subj_id[5:9])
                if study_name == "OASIS3":
                    subj_id_alpha = str(subj_id[0:3] + "IS" + subj_id[3:])

                # Extract the corresponding BIDS id and create the output file if doesn't exist
                subj_bids = [s for s in bids_ids if subj_id_alpha in s]
                if len(subj_bids) == 0:
                    # If the subject is not an excluded one
                    if subj_id not in subj_to_remove:
                        print(
                            f"{sessions_fields[i]} for {subj_id} not found in the BIDS converted."
                        )
                else:
                    subj_bids = subj_bids[0]
                    sessions_df[sessions_fields_bids[i]] = row[sessions_fields[i]]

                    subj_dir = path.join(
                        path.dirname(path.dirname(clinical_data_dir)),
                        "out",
                        "bids",
                        subj_bids,
                    )
                    session_names = get_bids_subjs_list(subj_dir)
                    for s in session_names:
                        s_name = s.replace("ses-", "")
                        if subj_bids not in sessions_dict:
                            sessions_dict.update({subj_bids: {}})
                        if s_name not in sessions_dict[subj_bids].keys():
                            sessions_dict[subj_bids].update({s_name: {"session_id": s}})
                        (sessions_dict[subj_bids][s_name]).update(
                            {sessions_fields_bids[i]: row[sessions_fields[i]]}
                        )

    return sessions_dict


def create_scans_dict(
    clinical_data_dir,
    study_name,
    clinic_specs_path,
    bids_ids,
    name_column_ids,
    name_column_ses,
    ses_dict,
):
    """[summary].

    Args:
        clinical_data_dir: path to the directory where the clinical data are stored
        study_name: name of the study (Ex ADNI)
        clinic_specs_path: path to the clinical specification file
        bids_ids: list of bids ids
        name_column_ids: name of the column where the subject id is contained
        name_column_ses: name of the column where the viscode of the session is contained
        ses_dict: links the session id to the viscode of the session.

    Returns: a pandas DataFrame that contains the scans information for all sessions of all participants.
    """
    import datetime
    import glob
    from os import path

    import pandas as pd

    scans_dict = {}
    prev_file = ""
    prev_sheet = ""

    if study_name not in get_supported_dataset():
        raise Exception(
            "Dataset not supported. Supported datasets are:", get_supported_dataset()
        )

    # Init the dictionary with the subject ids
    for bids_id in bids_ids:
        scans_dict[bids_id] = dict()
        for session_id in {"ses-" + key for key in ses_dict[bids_id].keys()}:
            scans_dict[bids_id][session_id] = {
                "T1/DWI/fMRI/FMAP": {},
                "PIB": {},
                "AV45": {},
                "FLUTE": {},
                "FDG": {},
            }

    scans_specs = pd.read_csv(clinic_specs_path + "_scans.tsv", sep="\t")
    fields_dataset = []
    fields_location = []
    fields_bids = []
    fields_mod = []

    # Extract the fields available and the corresponding bids name, location and type
    for i in range(0, len(scans_specs[study_name])):
        field = scans_specs[study_name][i]
        if not pd.isnull(field):
            fields_dataset.append(field)
            fields_bids.append(scans_specs["BIDS CLINICA"][i])
            fields_location.append(scans_specs[study_name + " location"][i])
            fields_mod.append(scans_specs["Modalities related"][i])

    # For each field available extract the original name, extract from the file all the values and fill a data structure
    for i in range(0, len(fields_dataset)):
        # Location is composed by file/sheet
        location = fields_location[i].split("/")
        file_name = location[0]
        if len(location) > 1:
            sheet = location[1]
        else:
            sheet = ""
        # Check if the file to read is already opened
        if file_name == prev_file and sheet == prev_sheet:
            pass
        else:
            file_to_read_path = path.join(clinical_data_dir, file_name)
            file_ext = path.splitext(file_name)[1]
            if file_ext == ".xlsx":
                file_to_read = pd.read_excel(
                    glob.glob(file_to_read_path)[0], sheet_name=sheet
                )
            elif file_ext == ".csv":
                file_to_read = pd.read_csv(
                    glob.glob(file_to_read_path)[0], sep=None, engine="python"
                )
            prev_file = file_name
            prev_sheet = sheet

        for bids_id in bids_ids:
            original_id = bids_id.replace("sub-" + study_name, "")
            for session_name in {"ses-" + key for key in ses_dict[bids_id].keys()}:
                row_to_extract = file_to_read[
                    (file_to_read[name_column_ids] == int(original_id))
                    & (file_to_read[name_column_ses] == ses_dict[session_name])
                ].index.tolist()
                if len(row_to_extract) > 0:
                    row_to_extract = row_to_extract[0]
                    # Fill the dictionary with all the information
                    value = file_to_read.iloc[row_to_extract][fields_dataset[i]]

                    if study_name == "AIBL":  # Deal with special format in AIBL
                        if value == "-4":
                            value = "n/a"
                        elif fields_bids[i] == "acq_time":
                            date_obj = datetime.datetime.strptime(value, "%m/%d/%Y")
                            value = date_obj.strftime("%Y-%m-%dT%H:%M:%S")

                    scans_dict[bids_id][session_name][fields_mod[i]][
                        fields_bids[i]
                    ] = value
                else:
                    print(f"Scans information for {bids_id} {session_name} not found.")
                    scans_dict[bids_id][session_name][fields_mod[i]][
                        fields_bids[i]
                    ] = "n/a"

    return scans_dict


def write_modality_agnostic_files(study_name, bids_dir):
    """
    Write the files README, dataset_description.json, .bidsignore and .bids-validator-config.json
    at the root of the BIDS directory.

    Args:
        study_name: name of the study (Ex ADNI)
        bids_dir: path to the bids directory
    """
    import json
    from os import path

    import clinica

    dataset_dict = {"Name": study_name, "BIDSVersion": "1.4.1", "DatasetType": "raw"}
    dataset_json = json.dumps(dataset_dict, skipkeys=True, indent=4)
    f = open(path.join(bids_dir, "dataset_description.json"), "w")
    f.write(dataset_json)
    f.close()

    file = open(path.join(bids_dir, "README"), "w")
    file.write(
        f"This BIDS directory was generated with Clinica v{clinica.__version__}.\n"
        f"More information on http://www.clinica.run\n"
    )
    file.close()

    validator_dict = {
        "ignore": [
            # Possibly dcm2nii(x) errors
            "NIFTI_UNIT",
            "INCONSISTENT_PARAMETERS",
            # fMRI-specific errors
            "SLICE_TIMING_NOT_DEFINED",
            "NIFTI_PIXDIM4",
            "BOLD_NOT_4D",
            "REPETITION_TIME_MUST_DEFINE",
            "TASK_NAME_MUST_DEFINE",
            # Won't fix errors
            "MISSING_SESSION",  # Allows subjects to have different sessions
            "INCONSISTENT_SUBJECTS",  # Allows subjects to have different modalities
            "SCANS_FILENAME_NOT_MATCH_DATASET",  # Necessary until PET is added to BIDS standard
            "CUSTOM_COLUMN_WITHOUT_DESCRIPTION",  # We won't create these JSON files as clinical description
            # is already done in TSV files of clinica.
            "NO_AUTHORS",  # Optional field in dataset_description.json
        ],
        "warn": [],
        "error": [],
        "ignoredFiles": [],
    }
    validator_json = json.dumps(validator_dict, skipkeys=True, indent=4)
    f = open(path.join(bids_dir, ".bids-validator-config.json"), "w")
    f.write(validator_json)
    f.close()

    f = open(path.join(bids_dir, ".bidsignore"), "w")
    f.write(
        "pet/\nconversion_info/\n"
    )  # pet/ is necessary until PET is added to BIDS standard
    f.close()


def write_sessions_tsv(bids_dir, sessions_dict):
    """Create <participant_id>_sessions.tsv files.

    Write the content of the function create scans dict in several TSV files
    following the BIDS specification.

    Args:
        bids_dir: path to the bids directory
        sessions_dict: output of the function create_sessions_dict
    """
    import os
    from glob import glob
    from os import path

    import pandas as pd

    bids_paths = glob(path.join(bids_dir, "sub-*"))

    for sp in bids_paths:
        bids_id = sp.split(os.sep)[-1]

        if bids_id in sessions_dict:
            session_df = pd.DataFrame.from_dict(sessions_dict[bids_id], orient="index")
            cols = session_df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            session_df = session_df[cols]
        else:
            print(f"No session data available for {sp}")
            session_df = pd.DataFrame(columns=["session_id"])
            session_df["session_id"] = pd.Series("M00")

        session_df = session_df.set_index("session_id").fillna("n/a")
        session_df.to_csv(
            path.join(sp, bids_id + "_sessions.tsv"),
            sep="\t",
            encoding="utf8",
        )


def write_scans_tsv(bids_dir, bids_ids, scans_dict):
    """Write the scans dict into TSV files.

    Args:
        bids_dir:  path to the BIDS directory
        bids_ids: list of bids ids
        scans_dict:  the output of the function create_scans_dict

    """
    import os
    from glob import glob
    from os import path

    import pandas as pd

    for bids_id in bids_ids:
        sessions_paths = glob(path.join(bids_dir, bids_id, "ses-*"))
        for session_path in sessions_paths:
            session_name = session_path.split(os.sep)[-1]
            scans_df = pd.DataFrame()
            tsv_name = bids_id + "_" + session_name + "_scans.tsv"
            # If the file already exists, remove it
            if os.path.exists(path.join(bids_dir, bids_id, session_name, tsv_name)):
                os.remove(path.join(bids_dir, bids_id, session_name, tsv_name))

            mod_available = glob(path.join(bids_dir, bids_id, session_name, "*"))
            for mod in mod_available:
                mod_name = os.path.basename(mod)
                # Grab scan files excluding their sidecar JSON.
                files = [
                    file
                    for file in glob(path.join(mod, "*"))
                    if "json" not in os.path.splitext(file)[1]
                ]
                for file in files:
                    file_name = os.path.basename(file)
                    if mod_name == "anat" or mod_name == "dwi" or mod_name == "func":
                        f_type = "T1/DWI/fMRI/FMAP"
                    elif mod_name == "pet":
                        description_dict = {
                            carac.split("-")[0]: carac.split("-")[1]
                            for carac in file_name.split("_")
                            if "-" in carac
                        }
                        f_type = description_dict["trc"].upper()
                    elif mod_name == "swi":
                        pass
                    else:
                        continue

                    row_to_append = pd.DataFrame(
                        scans_dict[bids_id][session_name][f_type], index=[0]
                    )
                    # Insert the column filename as first value
                    row_to_append.insert(0, "filename", path.join(mod_name, file_name))
                    scans_df = scans_df.append(row_to_append)

            scans_df = scans_df.set_index("filename").fillna("n/a")
            scans_df.to_csv(
                path.join(bids_dir, bids_id, session_name, tsv_name),
                sep="\t",
                encoding="utf8",
            )


# -- Other methods --
def contain_dicom(folder_path):
    """Check if a folder contains DICOM images.

    Args:
        folder_path: path to the folder

    Returns:
        True if DICOM files are found inside the folder, False otherwise
    """
    from glob import glob
    from os import path

    dcm_files = glob(path.join(folder_path, "*.dcm"))
    if len(dcm_files) > 0:
        return True

    return False












def get_supported_dataset():
    """Return the list of supported datasets."""
    return ["ADNI", "CLINAD", "PREVDEMALS", "INSIGHT", "OASIS", "AIBL"]


def get_bids_subjs_list(bids_path):
    """Given a BIDS compliant dataset, return the list of all the subjects available."""
    import os
    from os import path

    return [d for d in os.listdir(bids_path) if os.path.isdir(path.join(bids_path, d))]


def get_bids_subjs_paths(bids_path):
    """Given a BIDS compliant dataset, returns the list of all paths to the subjects folders."""
    import os
    from os import path

    return [
        path.join(bids_path, d)
        for d in os.listdir(bids_path)
        if os.path.isdir(path.join(bids_path, d))
    ]


def compute_new_subjects(original_ids, bids_ids):
    """Check for new subject to convert.

    This function checks for news subjects to convert to the BIDS version i.e. subjects
    contained in the unorganised version that are not available in the BIDS version.

    Args:
        original_ids: list of all the ids of the unorganized folder.
        bids_ids: list of all the BIDS ids contained inside the BIDS converted version of the dataset

    Returns:
        a list containing the original_ids of the subjects that are not available in the BIDS converted version
    """
    to_return = []
    original_ids = remove_space_and_symbols(original_ids)

    for s in original_ids:
        if not any(s in id for id in bids_ids):
            to_return.append(s)

    return to_return


def remove_space_and_symbols(data):
    """Remove spaces and  - _ from a list (or a single) of strings.

    Args:
        data: list of strings or a single string to clean

    Returns:
        data: list of strings or a string without space and symbols _ and -
    """
    import re

    if type(data) is list:
        for i in range(0, len(data)):
            data[i] = re.sub("[-_ ]", "", data[i])
    else:
        data = re.sub("[-_ ]", "", data)

    return data


def get_ext(file_path):
    import os

    root, ext = os.path.splitext(file_path)
    if ext in ".gz":
        file_ext = os.path.splitext(root)[1] + ext
    else:
        file_ext = ext
    return file_ext


def compress_nii(file_path):
    """Compress nii files.

    Args:
        file_path (str): path to the file to convert
    """
    import gzip
    import shutil
    from os import remove

    with open(file_path, "rb") as f_in:
        with gzip.open(file_path + ".gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    remove(file_path)


def json_from_dcm(dcm_dir, json_path):
    """
    Writes descriptive JSON file from DICOM header

    Args:
        dcm_dir (str): Path to the DICOM directory
        json_path (str): Path to the output JSON file
    """
    import json
    from glob import glob
    from os import path

    from pydicom import dcmread
    from pydicom.tag import Tag

    from clinica.utils.stream import cprint

    fields_dict = {
        "DeviceSerialNumber": Tag(("0018", "1000")),
        "Manufacturer": Tag(("0008", "0070")),
        "ManufacturersModelName": Tag(("0008", "1090")),
        "SoftwareVersions": Tag(("0018", "1020")),
        "BodyPart": Tag(("0018", "0015")),
        "Units": Tag(("0054", "1001")),
        # Institution
        "InstitutionName": Tag(("0008", "0080")),
        "InstitutionAddress": Tag(("0008", "0081")),
        "InstitutionalDepartmentName": Tag(("0008", "1040")),
        # MRI
        "MagneticFieldStrength": Tag(("0018", "0087")),
        # PET
        "InjectedRadioactivity": Tag(("0018", "1074")),
        "MolarActivity": Tag(("0018", "1077")),
        "InjectionStart": Tag(("0018", "1042")),
        "FrameDuration": Tag(("0018", "1242")),
    }

    try:
        dcm_path = glob(path.join(dcm_dir, "*.dcm"))[0]
        ds = dcmread(dcm_path)
        json_dict = dict()
        for key, tag in fields_dict.items():
            if tag in ds.keys():
                json_dict[key] = ds.get(tag).value

        json = json.dumps(json_dict, skipkeys=True, indent=4)
        with open(json_path, "w") as f:
            f.write(json)
    except IndexError:
        cprint(msg=f"No DICOM found at {dcm_dir}", lvl="warning")
