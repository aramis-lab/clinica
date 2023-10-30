import warnings
from dataclasses import dataclass
from enum import Enum
from functools import partial
from os import PathLike
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import pandas as pd

from clinica.utils.pet import Tracer


class Modality(str, Enum):
    T1 = "t1"
    AV45 = "av45"
    FLUTE = "flute"
    PIB = "pib"

    @property
    def name_of_path(self) -> str:
        return "Path_to_T1" if self == "t1" else "Path_to_pet"

    @property
    def tracer(self) -> Tracer:
        if self == "t1":
            raise ValueError("T1 modality does not have a tracer.")
        if self == "flute":
            return Tracer.FMM
        if self == "pib":
            return Tracer.PIB
        if self == "av45":
            return Tracer.AV45

    @property
    def bids_folder(self) -> str:
        return "anat" if self == "t1" else "pet"

    @property
    def suffix(self) -> str:
        return "T1w" if self == "t1" else f"trc-{self.tracer}_pet"


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
    modality: Modality

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
                self.modality.name_of_path: [str(scan.path) for scan in self.scans],
            }
        )


def paths_to_bids(
    path_to_dataset: Path,
    path_to_csv: Path,
    bids_dir: Path,
    modality: Modality,
    overwrite: bool = False,
    n_procs: Optional[int] = 1,
) -> List[str]:
    """Convert all the T1 images found in the AIBL dataset downloaded in BIDS.

    Parameters
    ----------
    path_to_dataset : Path
        Path to the AIBL dataset.

    path_to_csv : Path
        Path to the csv file containing clinical data.

    bids_dir : Path
        Path to save the AIBL-T1-dataset converted in a BIDS format.

    modality : Modality
        Modality to convert.

    overwrite : bool
        If True previous existing outputs will be erased.

    n_procs : int, optional
        The requested number of processes.
        If specified, it should be between 1 and the number of available CPUs.
        Default=1.

    Return
    ------
    list of str :
        list of all the images that are potentially converted in a BIDS format
        and saved in the bids_dir. This does not guarantee existence.
    """
    from functools import partial
    from multiprocessing import Pool

    (bids_dir / "conversion_info").mkdir(parents=True, exist_ok=True)

    if modality == Modality.T1:
        images = _find_path_to_t1(path_to_dataset, path_to_csv)
    else:
        images = _find_path_to_pet_modality(path_to_dataset, path_to_csv, modality)

    images.to_csv(
        bids_dir / "conversion_info" / f"{modality}_paths.tsv",
        index=False,
        sep="\t",
        encoding="utf-8",
    )

    images_list = list([data for _, data in images.iterrows()])

    with Pool(processes=n_procs) as pool:
        create_file_ = partial(
            _create_file,
            modality=modality,
            bids_dir=bids_dir,
            overwrite=overwrite,
        )
        output_file_treated = pool.map(create_file_, images_list)

    return output_file_treated


def _find_path_to_t1(path_to_dataset: Path, path_to_csv: Path) -> pd.DataFrame:
    """Creates a DataFrame for the T1 images.
    For each image, the subject ID, the session ID and the path to the image are
    reported.

    Parameters
    ----------
    path_to_dataset : Path
        Path to AIBL dataset.

    path_to_csv : Path
        Path to the folder with the CSV files downloaded.

    Return
    ------
    pd.DataFrame :
        Pandas dataframe which contains all the paths for the T1
        images, and the corresponding subject ID, and session ID.
    """
    import pandas as pd

    # two csv_files contain information regarding the T1w MRI images
    mri_dfs = tuple(
        [
            pd.read_csv(_get_first_file_matching_pattern(path_to_csv, pattern))
            for pattern in ("aibl_mrimeta_*.csv", "aibl_mri3meta_*.csv")
        ]
    )
    subjects_id = _listdir_nohidden(path_to_dataset)
    # all the subjects downloaded are taken into account for the conversion, except this sample
    if "0151083" in subjects_id:
        del subjects_id[subjects_id.index("0151083")]
    t1_scans = _find_path_to_t1_adni(mri_dfs, subjects_id, path_to_dataset)
    t1_scans = _find_path_to_t1_sag(subjects_id, path_to_dataset, t1_scans)
    return t1_scans.to_df()


def _get_first_file_matching_pattern(folder: Path, pattern: str) -> Path:
    """Return the first file in given folder matching the provided pattern.
    Raise a ValueError if no such file found.
    """
    if pattern == "":
        raise ValueError("Pattern is not valid.")
    try:
        return sorted(list(folder.glob(pattern)))[0]
    except IndexError:
        raise ValueError(f"No file matching pattern {folder}/{pattern}.")


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

    subject_ids : list of str
        List of subjects IDs in the downloaded dataset.

    path_to_dataset : Path
        Path to AIBL dataset.

    Return
    ------
    ScanCollection :
        The collected T1 scan paths.
    """
    scans = ScanCollection([], Modality.T1)

    for subject_id in subject_ids:
        for mri_file in mri_files:
            if int(subject_id) in list(mri_file.RID):
                path_to_subject = path_to_dataset / str(subject_id)
                for folder in _listdir_nohidden(path_to_subject):
                    if path_to_t1_folder := _find_t1_folder(folder, path_to_subject):
                        for exam_date in _listdir_nohidden(path_to_t1_folder):
                            session_id = _find_patient_session_id(
                                exam_date, subject_id, mri_file
                            )
                            if session_id != "-4":
                                path_to_t1 = path_to_t1_folder / str(exam_date)
                                for image_id in _listdir_nohidden(path_to_t1):
                                    scans.add_scan(
                                        subject_id, session_id, path_to_t1 / image_id
                                    )
    return scans


def _find_t1_in_paths(
    subdirectory: str,
    path_to_t1_images: Path,
    paths_to_convert: Sequence[str],
) -> Optional[Path]:
    """Find the directory containing T1 images.

    Parameters
    ----------
    subdirectory : str
        Name of the folder.

    path_to_t1_images : Path
        Path to T1 images.

    paths_to_convert : Sequence[str]
        Sequence of paths to convert.

    Return
    ------
    str or None:
        Previous path to arrive to the T1 image.
        Return None if no path found.
    """
    for folder in paths_to_convert:
        if folder == subdirectory:
            return path_to_t1_images / subdirectory
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


def _find_patient_session_id(
    exam_date: str, subject_id: str, metadata_df: pd.DataFrame
) -> str:
    """Returns the session_ID for the requested subject.

    It controls if the dates corresponding to the image (from the name of the subdirectory)
    correspond to one of the dates listed from the csv_file for the subject analysed.
    The session_ID is the corresponding session for that patient in that date.

    It returns -4 if there are no information.

    Parameters
    ----------
    exam_date : str
        Date where the image has been taken, it is saved from the name
        of the corresponding subdirectory.

    subject_id : str
        The subject's identifier.

    metadata_df : pd.DataFrame
        DataFrame obtained from the CSV file where all the information are listed.

    Return
    ------
    str :
        Session ID of the patient.
    """
    import re

    session_id = None
    indices = _find_correspondence_index(subject_id, metadata_df)
    csv_date = metadata_df.EXAMDATE[indices]

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
            session_id = metadata_df.VISCODE[index]

    return session_id or "-4"


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


def _find_path_to_pet_modality(
    path_to_dataset: Path,
    path_to_csv: Path,
    modality: Modality,
) -> pd.DataFrame:
    """Create a Dataframe which contains all the paths to the PET
    image of a modality (for example AV45 or PIB).

    Parameters
    ----------
    path_to_dataset : Path
        Path to AIBL dataset.

    path_to_csv : Path
        Path to CSV file.

    modality : Modality
        The PET modality for which to find paths.

    Return
    ------
    pd.DataFrame :
        A dataframe which contains the path for PET images for a single
        modality and subject_ID and session_ID are reported for each path.
    """
    df = _load_pet_csv(path_to_csv, modality)
    subjects_id = _listdir_nohidden(path_to_dataset)
    scans = ScanCollection([], modality)

    if "0151083" in subjects_id:
        del subjects_id[subjects_id.index("0151083")]

    for subject_id in subjects_id:
        if int(subject_id) in list(df.RID):
            path_to_subject = path_to_dataset / str(subject_id)
            subject_subdirectories_pet = _check_subdirectories_pet(
                _listdir_nohidden(path_to_subject),
                _list_folder_without_pet(),
            )
            for folder in subject_subdirectories_pet:
                path_to_pet_folder = path_to_subject / folder
                exam_dates = _listdir_nohidden(path_to_pet_folder)
                # exam date of the image which is going to be converted
                for exam_date in exam_dates:
                    # selection of the session_ID matching the data in the csv_file with the one of the image
                    session_id = _find_patient_session_id(
                        str(exam_date), subject_id, df
                    )
                    if session_id != "-4":
                        path_to_pet = path_to_pet_folder / str(exam_date)
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


def _load_pet_csv(path_to_csv: Path, modality: Modality) -> pd.DataFrame:
    import pandas as pd

    path_to_csv_pet_modality = _get_first_file_matching_pattern(
        path_to_csv, f"aibl_{modality}meta_*.csv"
    )

    if not path_to_csv_pet_modality.exists():
        raise FileNotFoundError(
            f"{path_to_csv_pet_modality} file not found in clinical data folder."
        )
    # Latest version of Flutemetamol CSV file (aibl_flutemeta_01-Jun-2018.csv)
    # has an extra column for some rows. However, each CSV file (regarding PET tracers)
    # contains the same columns. The usecols fixes this issue.
    return pd.read_csv(
        path_to_csv_pet_modality,
        sep=",|;",
        usecols=list(range(0, 36)),
        engine="python",
    )


def _check_subdirectories_pet(
    subdirectories: List[str], no_pet: List[str]
) -> List[str]:
    """Returns the correct subdirectories for the PET images, they should
    belong to the list where there all the possible names of the PET images.

    Parameters
    ----------
    subdirectories : list of str
        The possible subdirectories which need to be checked.

    no_pet : list of str
        List of folder names which do not contain PET images.

    Return
    ------
    list of str :
        List of subdirectories which contain a PET image which
        needs to be converted.
    """
    return list(
        set([s for s in subdirectories if (s not in no_pet) & (s != ".DS_Store")])
    )


def _list_folder_without_pet() -> List[str]:
    """Return a list of all the folders which do not contain PET images."""
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


def _create_file(
    image: pd.Series, modality: Modality, bids_dir: Path, overwrite: bool
) -> Optional[str]:
    """Convert the AIBL PET images into the BIDS specification.

    There are three pet modalities: av45, pib, flute. All of them are converted in BIDS.

    Parameters
    ----------
    image : pd.Series
    modality : Modality
    bids_dir : Path
    overwrite : bool

    Return
    ------
    str or None :
        Path to file
    """
    from numpy import nan

    from clinica.iotools.bids_utils import json_from_dcm
    from clinica.iotools.converter_utils import viscode_to_session
    from clinica.iotools.utils.data_handling import center_nifti_origin
    from clinica.utils.stream import cprint

    participant_id = f"sub-AIBL{image.Subjects_ID}"
    session_id = image.Session_ID
    image_path = image[modality.name_of_path]

    if image_path == nan:
        cprint(
            msg=(
                f"[{modality.upper()}] No path specified for {image.Subjects_ID} "
                f"in session {session_id}"
            ),
            lvl="info",
        )
        return nan
    else:
        cprint(
            msg=(
                f"[{modality.upper()}] Processing subject {image.Subjects_ID} "
                f"in session {session_id}"
            ),
            lvl="info",
        )
    session_id = viscode_to_session(session_id)
    output_path = (
        bids_dir
        / f"{participant_id}"
        / f"{session_id}"
        / modality.bids_folder
        / f"{participant_id}_{session_id}_{modality.suffix}"
    )
    # image is saved following BIDS specifications
    if output_path.with_suffix(".nii.gz").exists() and not overwrite:
        cprint(f"Subject {image.Subjects_ID} - session {session_id} already processed.")
        output_image = output_path.with_suffix(".nii.gz")
    else:
        if output_path.with_suffix(".nii.gz").exists():
            output_path.with_suffix(".nii.gz").unlink()
        output_image = _dicom_to_nii(output_path, image_path)
        json_from_dcm(image_path, str(output_path.with_suffix(".json")))

    center_nifti_origin(output_image, output_image)

    return str(output_image)


def _dicom_to_nii(output_path: Path, image_path: str) -> Path:
    """Convert the DICOM images to NIfTI files using dcm2niix.

    Parameters
    ----------
    output_path : Path
     Path to the NIfTI image without the extension.

    image_path : str
        Path to where the DICOM files are stored.

    Return
    ------
    Path :
        Path to the image in NIfTI format.
    """
    from clinica.iotools.bids_utils import run_dcm2niix
    from clinica.utils.stream import cprint

    try:
        output_path.parent.mkdir(parents=True)
    except OSError:
        if not output_path.parent.is_dir():
            raise

    # if image.Is_Dicom:
    run_dcm2niix(
        input_dir=image_path,
        output_dir=str(output_path.parent),
        output_fmt=output_path.name,
        compress=True,
        bids_sidecar=False,
    )

    nifti_file = output_path.with_suffix(".nii.gz")

    if not nifti_file.exists():
        cprint(f"{nifti_file} should have been created but this did not happen")
    return nifti_file
