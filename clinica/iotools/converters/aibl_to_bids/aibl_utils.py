"""Utils to convert AIBL dataset in BIDS."""
import warnings
from dataclasses import dataclass
from functools import partial
from os import PathLike
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import pandas as pd


def _listdir_nohidden(path: PathLike) -> List[str]:
    """List all the subdirectories of path except the hidden folders.

    Parameters
    ----------
    path: str
        Path to list subdirectories from.

    Returns
    -------
    list of str:
        Subdirectories found within path.
    """
    from pathlib import Path

    return [
        str(p.name)
        for p in Path(path).iterdir()
        if p.is_dir() and not p.name.startswith(".")
    ]


def _find_t1_in_paths(
    subdirectory: str,
    path_to_t1_images: str,
    paths_to_convert: Sequence[str],
) -> Optional[str]:
    """Find the directory containing T1 images.

    Parameters
    ----------
    subdirectory : str
        Name of the folder.

    path_to_t1_images : str
        Path to T1 images.

    paths_to_convert : Sequence[str]
        Sequence of paths to convert.

    Return
    ------
    str or None:
        Previous path to arrive to the T1 image.
        Return None if no path found.
    """
    import os

    for folder in paths_to_convert:
        if folder == subdirectory:
            return os.path.join(path_to_t1_images, subdirectory)
    return None


_find_t1_folder = partial(
    _find_t1_in_paths,
    paths_to_convert=[
        "MPRAGE_ADNI_confirmed",
        "MPRAGE",
        "MPRAGE_ADNI_confirmed_RPT",
        "MPRAGE_ADNI_confirmed_REPEATX2",
        "MPRAGE_ADNI_confirmed_repeat",
        "MPRAGE_ADNI_confirmed_REPEAT",
        "MPRAGE_ADNI_conf_REPEAT",
    ],
)


_find_t1_folder_nodata = partial(
    _find_t1_in_paths,
    paths_to_convert=[
        "MPRAGESAGISOp2ND",
        "MPRAGE_SAG_ISO_p2_ND",
        "MPRAGE_SAG_ISO_p2",
    ],
)


def _find_correspondence_index(subject_id: str, df: pd.DataFrame) -> List[str]:
    """Returns the index of the CSV file analysed for a given subject.

    Parameters
    ----------
    subject_id : str
        The subject's identifier.

    df : pd.DataFrame
        The dataframe in which to look for subject data.

    Return
    ------
    list of str :
        List of indexes for the subject

    Raises
    ------
    ValueError :
        If no correspondence was found.
    """
    for rid in df.RID:
        if subject_id == str(rid):
            return df.RID[df.RID == rid].index.tolist()
    raise ValueError(
        f"Unable to find correspondence between the subject ID {subject_id} "
        "and RID in provided dataframe."
    )


def _find_correspondence_date(index: list, df: pd.DataFrame) -> pd.Series:
    """Return the dates reported in the csv_file for a given index.

    Parameters
    ----------
    index : list
        List of index.
    df : pd.DataFrame
        DataFrame where all the information are listed

    Return
    ------
    pd.Series :
        Exam dates at index.
    """
    return df.EXAMDATE[index]


def match_data(exam_date: str, subject_id: str, csv_file: pd.DataFrame) -> str:
    """

    This method returns the session_ID. It controls if the dates
    corresponding to the image (from the name of the subdirectory)
    correspond to one of the dates listed from the csv_file for the subject
    analysed. The session_ID is the corresponding session for that patient
    in that date.  It returns -4 if there are no information.

    :param exam_date: date where the image has been taken, it is saved
    from the name of the corresponding subdirector
    :type exam_date: str
    :param subject_id: Subject's identifier
    :type subject_id: str
    :param csv_file: csv file where all the information are listed
    :type csv_file: pd.DataFrame
    :return session_id of the patient
    :rtype: str
    """
    import re

    session_id = None
    indices = _find_correspondence_index(subject_id, csv_file)
    csv_date = _find_correspondence_date(indices, csv_file)
    for index in indices:
        if str(csv_date[index]) == "-4":
            continue
        m = re.search(
            "([0-9].*)-(.*)-(.*)_(.*)_(.*)_(.*)", exam_date
        )  # string from image directory
        p = re.search(
            "(.*)/(.*)/(.*)", str(csv_date[index])
        )  # string from the date of the csv_file
        if (
            (p.group(1) == m.group(2))
            & (p.group(2) == m.group(3))
            & (p.group(3) == m.group(1))
        ):
            session_id = csv_file.VISCODE[index]
    return session_id or "-4"


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
    :param no_pet: list of names of folders which not contain PET images

    :return subdirectory which is containing a PET image which needs to be
    converted
    """

    for j in range(len(sub)):
        if (sub[j] not in no_pet) & (sub[j] != ".DS_Store"):
            subdirectories.append(sub[j])
    subdirectories = list(set(subdirectories))
    return subdirectories


def _dicom_to_nii(subject, output_path, output_filename, image_path):
    """Convert the DICOM images to NIfTI files using dcm2niix.

    :param subject:
    :param output_path: where NIfTI image is stored
    :param output_filename: name of the NIfTI image
    :param image_path: where DICOM files are stored

    :return: Image in NIfTI format
    """
    import os
    from os.path import exists

    from clinica.iotools.bids_utils import run_dcm2niix
    from clinica.utils.stream import cprint

    try:
        os.makedirs(output_path)
    except OSError:
        if not os.path.isdir(output_path):
            raise

    # if image.Is_Dicom:
    run_dcm2niix(
        input_dir=image_path,
        output_dir=output_path,
        output_fmt=output_filename,
        compress=True,
        bids_sidecar=False,
    )

    nifti_file = os.path.join(output_path, output_filename + ".nii.gz")

    if not exists(nifti_file):
        cprint(nifti_file + " should have been created but this did not happen")
    return nifti_file


def _find_path_to_pet_modality(
    path_to_dataset: Path, csv_file: pd.DataFrame
) -> pd.DataFrame:
    """Create a Dataframe which contains all the paths to the PET image of a modality (for example AV45 or PIB).

    :param path_to_dataset: path to AIBL dataset
    :param csv_file: file which correspond to the modality

    :return: A dataframe which contains the path for PET images for a
    single modality and subject_ID and session_ID are reported for each
    path
    """
    no_pet = list_of_paths()
    subjects_id = _listdir_nohidden(path_to_dataset)

    if "0151083" in subjects_id:
        del subjects_id[subjects_id.index("0151083")]

    scans = ScanCollection([], "pet")

    for subject_id in subjects_id:
        if int(subject_id) in list(csv_file.RID):
            # check if the subject is present in the csv_file for the modality selected
            subdirectories = []
            path_to_pet = path_to_dataset / str(subject_id)
            subdirectory_all = _listdir_nohidden(path_to_pet)
            subdirectories = check_subdirectories_pet(
                subdirectories, subdirectory_all, no_pet
            )
            # selection only of the folders which contain PET image
            for folder in subdirectories:
                path_to_pet = path_to_pet / folder
                exam_dates = _listdir_nohidden(path_to_pet)
                # exam date of the image which is going to be converted
                for exam_date in exam_dates:
                    # selection of the session_ID matching the data in the csv_file with the one of the image
                    session_id = match_data(str(exam_date), subject_id, csv_file)
                    if session_id != "-4":
                        path_to_pet = path_to_pet / str(exam_date)
                        # For the RID 1607 there are two PET images of the flute modality, and we select the first
                        if (
                            subject_id == "1607"
                            and folder == "Flute_256_1.6_Zoom_plain_4_x_4_Iter"
                        ):
                            image_ids = ["I442930"]
                        else:
                            image_ids = _listdir_nohidden(path_to_pet)
                        for image_id in image_ids:
                            scans.add_scan(
                                subject_id, session_id, path_to_pet / image_id
                            )
    return scans.to_df()


@dataclass
class Scan:
    """Class for keeping track of the relationship between a subject ID,
    a session ID, and a path to a modality scan.
    """

    subject_id: str
    session_id: str
    path: Path


@dataclass
class ScanCollection:
    """
    Maintain coherence in the collection : it is impossible to add
    a scan for an already existing subject session pair.
    """

    scans: List[Scan]
    modality: str

    def _is_new_scan_valid(self, new_scan: Scan) -> bool:
        for scan in self.scans:
            if (scan.subject_id, scan.session_id) == (
                new_scan.subject_id,
                new_scan.session_id,
            ):
                warnings.warn(
                    f"Already a scan for {self.modality} in the collection for subject {scan.subject_id} "
                    f"and session {scan.session_id} : {scan.path}. "
                    f"Cannot add another scan {new_scan.path} for this subject session pair."
                )
                return False
        return True

    def add(self, scan: Scan):
        if self._is_new_scan_valid(scan):
            self.scans.append(scan)

    def add_scan(self, subject_id: str, session_id: str, path: Path):
        """Add a scan to the collection."""
        self.add(Scan(subject_id, session_id, path))

    def to_df(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "Subjects_ID": [scan.subject_id for scan in self.scans],
                "Session_ID": [scan.session_id for scan in self.scans],
                f"Path_to_{self.modality}": [str(scan.path) for scan in self.scans],
            }
        )


def _find_path_to_t1_adni(
    mri_files: Tuple[pd.DataFrame, pd.DataFrame],
    subject_ids: List[str],
    path_to_dataset: Path,
) -> ScanCollection:
    """Creates a Dataframe which contains all the paths to the T1
    images which are ADNI compliant (as explained in the AIBL website).
    These images differ from the others T1 of the dataset since the exam
    date is reported in the CSV file.

    Parameters
    ----------
    mri_files : tuple of 2 DataFrames
        List of dataframes extracted from the clinical data.
        Indeed, there are two files which describe the parameters of
        the T1 images : MRI 1.5 T and MRI 3T.

    subjects_ids : list of str
        List of subjects IDs in the downloaded dataset.

    path_to_dataset : Path
        Path to AIBL dataset.

    Return
    ------
    ScanCollection :
        The collected T1 scan paths.
    """
    scans = ScanCollection([], "T1")

    for subject_id in subject_ids:
        for mri_file in mri_files:
            if int(subject_id) in list(mri_file.RID):
                # check if the information of the subject are present in the csv_file
                path_to_t1 = path_to_dataset / str(subject_id)
                for folder in _listdir_nohidden(path_to_t1):
                    if path_to_t1 := _find_t1_folder(folder, path_to_t1):
                        for exam_date in _listdir_nohidden(path_to_t1):
                            # check if the corresponding session_ID can be found in the csv_file
                            session_id = match_data(exam_date, subject_id, mri_file)
                            if session_id != "-4":
                                path_to_t1 = path_to_t1 / str(exam_date)
                                for image_id in _listdir_nohidden(path_to_t1):
                                    scans.add_scan(
                                        subject_id, session_id, path_to_t1 / image_id
                                    )
    return scans


def _find_path_to_t1_sag(
    subjects_id: List[str],
    path_to_dataset: Path,
    t1_scans: ScanCollection,
) -> ScanCollection:
    """Creates a ScanCollection which contains all the paths to the T1 images which are not
    ADNI compliant, they contain the word "SAG" in their name.

    Parameters
    ----------
    subjects_id : list of str
        Subjects ID in the downloaded dataset.

    path_to_dataset : Path
        Path to AIBL dataset.

    t1_scans : ScanCollection
        The previously computed collection of T1 scan paths.
        The function will potentially add new scan paths to this collection.

    Return
    ------
    ScanCollection :
        The collected T1 scan paths.

    Notes
    -----
    This function completes the list of all the T1 paths including all the images where
    we didn't find the exam data, but we can fix it with further analysis.
    """
    for subject_id in subjects_id:
        path_to_t1 = path_to_dataset / str(subject_id)
        directories_to_analyze = (
            "MPRAGESAGISOp2ND",
            "MPRAGE_SAG_ISO_p2_ND",
            "MPRAGE_SAG_ISO_p2",
        )
        subdirectory_for_subject = [
            d for d in _listdir_nohidden(path_to_t1) if d in directories_to_analyze
        ]
        if not subdirectory_for_subject:
            continue
        path_to_t1 = path_to_t1 / subdirectory_for_subject[0]
        exam_date = _listdir_nohidden(path_to_t1)
        session_id = "M54" if subject_id in (342, 557) else "M00"
        # if for a subject in the same session we have both this image
        # and the "ADNI" compliant we are converting the second one
        # since the exam date is more precise.
        path_to_t1 = path_to_t1 / str(exam_date[0])
        image_id = _listdir_nohidden(path_to_t1)
        path_to_t1 = path_to_t1 / image_id[0]
        t1_scans.add_scan(subject_id, session_id, path_to_t1)

    return t1_scans


def _find_path_to_t1(path_to_dataset: Path, path_to_csv: str) -> pd.DataFrame:
    """Creates a DataFrame for the T1 images.
    For each image, the subject ID, the session ID and the path to the image are
    reported.

    Parameters
    ----------
    path_to_dataset : str
        Path to AIBL dataset.

    path_to_csv : str
        Path to the folder with the CSV files downloaded.

    Return
    ------
    pd.DataFrame :
        Pandas dataframe which contains all the paths for the T1
        images, and the corresponding subject ID, and session ID.
    """
    import glob
    import os

    import pandas as pd

    # two csv_files contain information regarding the T1w MRI images
    mri_meta_df = pd.read_csv(
        glob.glob(os.path.join(path_to_csv, "aibl_mrimeta_*.csv"))[0]
    )
    mri_3meta_df = pd.read_csv(
        glob.glob(os.path.join(path_to_csv, "aibl_mri3meta_*.csv"))[0]
    )
    mri_dfs = (mri_meta_df, mri_3meta_df)
    subjects_id = _listdir_nohidden(path_to_dataset)
    # list of all the folders which correspond to the subject ID
    # all the subjects downloaded are taken into account for the conversion, except this sample
    if "0151083" in subjects_id:
        del subjects_id[subjects_id.index("0151083")]
    t1_scans = _find_path_to_t1_adni(mri_dfs, subjects_id, path_to_dataset)
    t1_scans = _find_path_to_t1_sag(subjects_id, path_to_dataset, t1_scans)
    return t1_scans.to_df()


# Convert the AIBL PET images into the BIDS specification.
# There are three pet modalities: av45, pib, flute. All of them are converted
# in BIDS
def _create_file(image, modality, bids_dir, overwrite):
    from os import remove
    from os.path import exists, join

    from numpy import nan

    from clinica.iotools.bids_utils import json_from_dcm
    from clinica.iotools.converter_utils import viscode_to_session
    from clinica.iotools.utils.data_handling import center_nifti_origin
    from clinica.utils.pet import Tracer
    from clinica.utils.stream import cprint

    subject = image.Subjects_ID
    participant_id = f"sub-AIBL{image.Subjects_ID}"
    session_id = image.Session_ID
    name_of_path = {
        "t1": "Path_to_T1",
        "av45": "Path_to_pet",
        "flute": "Path_to_pet",
        "pib": "Path_to_pet",
    }

    image_path = image[name_of_path[modality]]

    if image_path == nan:
        cprint(
            msg=(
                f"[{modality.upper()}] No path specified for {subject} "
                f"in session {session_id}"
            ),
            lvl="info",
        )
        return nan
    else:
        cprint(
            msg=(
                f"[{modality.upper()}] Processing subject {subject} "
                f"in session {session_id}"
            ),
            lvl="info",
        )

    session_id = viscode_to_session(session_id)

    # creation of the path
    if modality == "t1":
        output_path = join(bids_dir, f"{participant_id}", f"{session_id}", "anat")
        output_filename = f"{participant_id}_{session_id}_T1w"
    elif modality in ["flute", "pib", "av45"]:
        tracer = {"flute": Tracer.FMM, "pib": Tracer.PIB, "av45": Tracer.AV45}[modality]
        output_path = join(bids_dir, f"{participant_id}", f"{session_id}", "pet")
        output_filename = f"{participant_id}_{session_id}_trc-{tracer}_pet"
    else:
        return None

    # image is saved following BIDS specifications
    if exists(join(output_path, output_filename + ".nii.gz")) and not overwrite:
        cprint(f"Subject {str(subject)} - session {session_id} already processed.")
        output_image = join(output_path, output_filename + ".nii.gz")
    else:
        if exists(join(output_path, output_filename + ".nii.gz")):
            remove(join(output_path, output_filename + ".nii.gz"))
        output_image = _dicom_to_nii(subject, output_path, output_filename, image_path)
        json_from_dcm(image_path, join(output_path, output_filename + ".json"))

    # Center all images
    center_nifti_origin(output_image, output_image)

    return output_image


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
    from functools import partial
    from multiprocessing import Pool, cpu_count
    from os import makedirs
    from os.path import exists, join
    from pathlib import Path

    import pandas as pds

    path_to_dataset = Path(path_to_dataset)

    if modality.lower() not in ["t1", "av45", "flute", "pib"]:
        # This should never be reached
        raise RuntimeError(f"{modality.lower()} is not supported for conversion")

    makedirs(join(bids_dir, "conversion_info"), exist_ok=True)

    # it reads the DataFrame where subject_ID, session_ID and path are saved
    if modality == "t1":
        images = _find_path_to_t1(path_to_dataset, path_to_csv)
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
            path_to_csv_pet_modality,
            sep=",|;",
            usecols=list(range(0, 36)),
            engine="python",
        )
        images = _find_path_to_pet_modality(path_to_dataset, df_pet)

    images.to_csv(
        join(bids_dir, "conversion_info", modality + "_paths.tsv"),
        index=False,
        sep="\t",
        encoding="utf-8",
    )

    images_list = list([data for _, data in images.iterrows()])

    with Pool(processes=max(cpu_count() - 1, 1)) as pool:
        create_file_ = partial(
            _create_file,
            modality=modality,
            bids_dir=bids_dir,
            overwrite=overwrite,
        )
        output_file_treated = pool.map(create_file_, images_list)

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
    from os import path

    import numpy as np
    import pandas as pd

    fields_bids = ["participant_id"]
    fields_dataset = []
    prev_location = ""
    prev_sheet = ""
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

    # Compute BIDS-compatible participant ID.
    participant_df["participant_id"] = "sub-AIBL" + participant_df["alternative_id_1"]

    # Keep year-of-birth only.
    participant_df["date_of_birth"] = participant_df["date_of_birth"].str.extract(
        r"/(\d{4}).*"
    )

    # Normalize sex value.
    participant_df["sex"] = participant_df["sex"].map({1: "M", 2: "F"}).fillna("n/a")

    # Normalize known NA values.
    participant_df.replace(-4, "n/a", inplace=True)

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

    for i in range(0, len(sessions_fields)):
        if not pd.isnull(sessions_fields[i]):
            fields_bids.append(sessions_fields_bids[i])
            fields_dataset.append(sessions_fields[i])

    files_to_read = []
    sessions_fields_to_read = []
    for i in range(0, len(sessions_fields)):
        # If the i-th field is available
        if not pd.isnull(sessions_fields[i]):
            # Load the file
            tmp = field_location[i]
            file_to_read_path = path.join(clinical_data_dir, tmp)
            files_to_read.append(glob.glob(file_to_read_path)[0])
            sessions_fields_to_read.append(sessions_fields[i])

    rid = pd.read_csv(files_to_read[0], dtype={"text": str}, low_memory=False).RID
    rid = list(set(rid))
    for r in rid:
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

        viscode[viscode == "bl"] = "M000"
        viscode = viscode.str.upper()

        sessions = pd.DataFrame(
            {
                "months": viscode.str[1:],
                "age": age,
                "MMS": MMSCORE,
                "cdr_global": CDGLOBAL,
                "diagnosis": DXCURREN,
                "examination_date": examdates,
            }
        )
        sessions = sessions.assign(
            session_id=lambda df: df.months.apply(lambda x: f"ses-M{int(x):03d}")
        )
        cols = sessions.columns.tolist()
        sessions = sessions[cols[-1:] + cols[:-1]]

        bids_paths = path.join(input_path, "sub-AIBL" + str(r))
        if path.exists(bids_paths):
            sessions.to_csv(
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
    ses_dict = {
        bids_id: {"M000": "bl", "M018": "m18", "M036": "m36", "M054": "m54"}
        for bids_id in bids_ids
    }
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
