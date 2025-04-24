from pathlib import Path
from typing import Iterable, Optional, Union

from clinica.converters.abstract_converter import Converter
from clinica.utils.filemanip import UserProvidedPath

from ._utils import ADNIModality

__all__ = ["convert"]


def convert(
    path_to_dataset: UserProvidedPath,
    bids_dir: UserProvidedPath,
    path_to_clinical: UserProvidedPath,
    clinical_data_only: bool = False,
    subjects: Optional[UserProvidedPath] = None,
    modalities: Optional[Iterable[Union[str, ADNIModality]]] = None,
    xml_path: Optional[UserProvidedPath] = None,
    mod_to_update: bool = False,
    n_procs: Optional[int] = 1,
):
    from .._utils import validate_input_path

    path_to_dataset = validate_input_path(path_to_dataset)
    bids_dir = validate_input_path(bids_dir, check_exist=False)
    path_to_clinical = validate_input_path(path_to_clinical)
    if xml_path:
        xml_path = validate_input_path(xml_path)
    if subjects:
        subjects = validate_input_path(subjects)
    modalities = _validate_input_modalities(modalities)

    adni_to_bids = AdniToBids()
    adni_to_bids.check_adni_dependencies()

    if not clinical_data_only:
        adni_to_bids.convert_images(
            source_dir=path_to_dataset,
            clinical_dir=path_to_clinical,
            dest_dir=bids_dir,
            modalities=modalities,
            subjects=subjects,
            mod_to_update=mod_to_update,
            n_procs=n_procs,
        )
    adni_to_bids.convert_clinical_data(
        clinical_data_dir=path_to_clinical,
        out_path=bids_dir,
        clinical_data_only=clinical_data_only,
        subjects=subjects,
        xml_path=xml_path,
    )


class AdniToBids(Converter):
    @classmethod
    def check_adni_dependencies(cls) -> None:
        """Check the dependencies of ADNI converter."""
        from clinica.utils.check_dependency import ThirdPartySoftware, check_software

        check_software(ThirdPartySoftware.DCM2NIIX)

    @staticmethod
    def _create_modality_agnostic_files(bids_dir: Path):
        from clinica.converters.study_models import StudyName
        from clinica.utils.stream import cprint

        from .._utils import write_modality_agnostic_files

        cprint("Creating modality agnostic files...", lvl="info")
        readme_data = {
            "link": "http://adni.loni.usc.edu",
            "desc": (
                "ADNI is a global research effort that actively supports the investigation and development of "
                "treatments that slow or stop the progression of Alzheimer's disease (AD).This multisite, longitudinal "
                "study assesses clinical, imaging, genetic and biospecimen biomarkers through the process of normal "
                "aging to mild cognitive impairment (MCI) and AD dementia.With established, standardized methods for "
                "imaging and biomarker collection and analysis, ADNI facilitates a way for scientists to conduct "
                "cohesive research and share compatible data with other researchers around the world."
            ),
        }
        write_modality_agnostic_files(
            study_name=StudyName.ADNI,
            readme_data=readme_data,
            bids_dir=bids_dir,
        )

    def convert_clinical_data(
        self,
        clinical_data_dir: Path,
        out_path: Path,
        clinical_data_only: bool = False,
        subjects: Optional[Path] = None,
        xml_path: Optional[Path] = None,
    ):
        """Convert the clinical data of ADNI specified into the file clinical_specifications_adni.xlsx.

        Args:
            clinical_data_dir:  path to the clinical data directory
            out_path: path to the BIDS directory
            clinical_data_only: process clinical data only
            subjects: restrict processing to this manifest of subjects
            xml_path: path to the XML metadata files
        """
        from clinica.converters.study_models import StudyName
        from clinica.dataset import (
            get_paths_to_subjects_in_bids_dataset,
            get_subjects_from_bids_dataset,
        )
        from clinica.utils.stream import cprint

        from .._utils import create_participants_df
        from ._json import create_json_metadata
        from ._utils import (
            correct_diagnosis_sc_adni3,
            create_adni_scans_files,
            create_adni_sessions_dict,
        )

        clinical_specifications_folder = Path(__file__).parents[0] / "specifications"
        if not out_path.exists():
            msg = f"BIDS folder {out_path} not found."
            cprint(msg, lvl="error")
            raise FileNotFoundError(msg)
        conversion_path = out_path / "conversion_info"
        self._create_modality_agnostic_files(out_path)

        if clinical_data_only:
            bids_ids, bids_subjects_paths = _get_bids_subjects_info(
                clinical_data_dir=clinical_data_dir,
                out_path=out_path,
                subjects=subjects,
            )
        else:
            bids_ids = get_subjects_from_bids_dataset(out_path)
            bids_subjects_paths = get_paths_to_subjects_in_bids_dataset(out_path)

        # -- Creation of participant.tsv --
        cprint("Creating participants.tsv...", lvl="info")
        participants_df = create_participants_df(
            StudyName.ADNI,
            clinical_specifications_folder,
            clinical_data_dir,
            bids_ids,
        )
        # Replace the original values with the standard defined by the AramisTeam
        participants_df["sex"] = participants_df["sex"].replace("Male", "M")
        participants_df["sex"] = participants_df["sex"].replace("Female", "F")

        # Correction of diagnosis_sc for ADNI3 participants
        participants_df = correct_diagnosis_sc_adni3(clinical_data_dir, participants_df)
        participants_df.to_csv(
            out_path / "participants.tsv",
            sep="\t",
            index=False,
            encoding="utf-8",
        )
        # -- Creation of sessions.tsv --
        cprint("Creating sessions files...", lvl="info")
        create_adni_sessions_dict(
            bids_ids,
            clinical_specifications_folder,
            clinical_data_dir,
            bids_subjects_paths,
        )
        # -- Creation of scans files --
        if conversion_path.exists():
            cprint("Creating scans files...", lvl="info")
            create_adni_scans_files(conversion_path, bids_subjects_paths)

        if xml_path is not None:
            if xml_path.exists():
                create_json_metadata(bids_subjects_paths, bids_ids, xml_path)
            else:
                cprint(
                    msg=(
                        f"Clinica was unable to find {xml_path}, "
                        "skipping xml metadata extraction."
                    ),
                    lvl="warning",
                )

    def convert_images(
        self,
        source_dir: Path,
        clinical_dir: Path,
        dest_dir: Path,
        modalities: Iterable[ADNIModality],
        subjects: Optional[Path] = None,
        mod_to_update: bool = False,
        n_procs: Optional[int] = 1,
    ):
        """Convert the images of ADNI.

        Args:
            source_dir: path to the ADNI directory
            clinical_dir: path to the clinical data directory
            dest_dir: path to the BIDS directory
            subjects: Path to list of subjects to process
            modalities: modalities to convert (T1, PET_FDG, PET_AMYLOID, PET_TAU, DWI, FLAIR, fMRI)
            mod_to_update: if given pre-existing images in the BIDS directory will be erased and extracted again.
        """
        from clinica.utils.stream import cprint

        from ._utils import get_subjects_list
        from .modality_converters import modality_converter_factory

        modalities = modalities or ADNIModality
        subjects = get_subjects_list(source_dir, clinical_dir, subjects)

        # Create the output folder if is not already existing
        (dest_dir / "conversion_info").mkdir(parents=True, exist_ok=True)
        version_number = len([f for f in (dest_dir / "conversion_info").iterdir()])
        conversion_dir = dest_dir / "conversion_info" / f"v{version_number}"
        conversion_dir.mkdir(parents=True, exist_ok=True)
        cprint(f"Destination folder = {dest_dir}", lvl="debug")

        for modality in modalities:
            for converter in modality_converter_factory(modality):
                converter(
                    source_dir=source_dir,
                    csv_dir=clinical_dir,
                    destination_dir=dest_dir,
                    conversion_dir=conversion_dir,
                    subjects=subjects,
                    mod_to_update=mod_to_update,
                    n_procs=n_procs,
                )


def _validate_input_modalities(
    modalities: Optional[Iterable[str]] = None,
) -> tuple[ADNIModality, ...]:
    if modalities:
        return tuple((ADNIModality(modality) for modality in modalities))
    return tuple(ADNIModality)


def _get_bids_subjects_info(
    clinical_data_dir: Path,
    out_path: Path,
    subjects: Optional[Path] = None,
) -> tuple[list[str], list[Path]]:
    from clinica.converters.study_models import StudyName, bids_id_factory

    from .._utils import load_clinical_csv

    # Read optional list of participants.
    if subjects:
        subjects = [subject for subject in subjects.read_text().split("\n") if subject]
    # Load all participants from ADNIMERGE.
    adni_merge = load_clinical_csv(clinical_data_dir, "ADNIMERGE")
    participants = adni_merge["PTID"].unique()
    # Filter participants if requested.
    participants = sorted(participants & subjects if subjects else participants)
    # Compute their corresponding BIDS IDs and paths.
    bids_ids = [
        bids_id_factory(StudyName.ADNI).from_original_study_id(p) for p in participants
    ]
    bids_paths = [out_path / bids_id for bids_id in bids_ids]

    return bids_ids, bids_paths
