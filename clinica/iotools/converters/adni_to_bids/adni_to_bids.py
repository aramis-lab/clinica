from pathlib import Path
from typing import Optional

from clinica.iotools.abstract_converter import Converter

__all__ = ["convert"]


def convert(
    path_to_dataset: Path,
    bids_dir: Path,
    path_to_clinical: Path,
    clinical_data_only: bool = False,
    subjects: Optional[Path] = None,
    modalities: Optional[list[str]] = None,
    xml_path: Optional[Path] = None,
    force_new_extraction: bool = False,
    n_procs: Optional[int] = 1,
):
    adni_to_bids = AdniToBids()
    adni_to_bids.check_adni_dependencies()
    if not path_to_dataset.exists():
        raise FileNotFoundError(
            f"The provided dataset path {path_to_dataset} does not exist."
        )
    if not path_to_clinical.exists():
        raise FileNotFoundError(
            f"The provided path to clinical data {path_to_clinical} does not exist."
        )
    if not clinical_data_only:
        adni_to_bids.convert_images(
            source_dir=path_to_dataset,
            clinical_dir=path_to_clinical,
            dest_dir=bids_dir,
            subjects=subjects,
            modalities=modalities,
            force_new_extraction=force_new_extraction,
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
    def get_modalities_supported(cls) -> list[str]:
        """Return a list of modalities supported.

        Returns: a list containing the modalities supported by the converter
        (T1, PET_FDG, PET_AMYLOID, PET_TAU, DWI, FLAIR, fMRI, FMAP)
        """
        return [
            "T1",
            "PET_FDG",
            "PET_AMYLOID",
            "PET_TAU",
            "DWI",
            "FLAIR",
            "fMRI",
            "FMAP",
        ]

    @classmethod
    def check_adni_dependencies(cls) -> None:
        """Check the dependencies of ADNI converter."""
        from clinica.utils.check_dependency import ThirdPartySoftware, check_software

        check_software(ThirdPartySoftware.DCM2NIIX)

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
        from clinica.iotools.bids_utils import (
            StudyName,
            create_participants_df,
            get_bids_subjs_list,
            get_bids_subjs_paths,
            write_modality_agnostic_files,
        )
        from clinica.iotools.converters.adni_to_bids.adni_utils import (
            correct_diagnosis_sc_adni3,
            create_adni_scans_files,
            create_adni_sessions_dict,
        )
        from clinica.utils.stream import cprint

        from .adni_json import create_json_metadata

        clinical_specifications_folder = Path(__file__).parents[0] / "specifications"
        if not out_path.exists():
            raise IOError("BIDS folder not found.")
        conversion_path = out_path / "conversion_info"

        if clinical_data_only:
            bids_ids, bids_subjects_paths = _get_bids_subjs_info(
                clinical_data_dir=clinical_data_dir,
                out_path=out_path,
                subjects=subjects,
            )
        else:
            bids_ids = get_bids_subjs_list(out_path)
            bids_subjects_paths = get_bids_subjs_paths(out_path)

        # -- Creation of modality agnostic files --
        cprint("Creating modality agnostic files...")
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
            bids_dir=out_path,
        )
        # -- Creation of participant.tsv --
        cprint("Creating participants.tsv...")
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
        cprint("Creating sessions files...")
        create_adni_sessions_dict(
            bids_ids,
            clinical_specifications_folder,
            clinical_data_dir,
            bids_subjects_paths,
        )
        # -- Creation of scans files --
        if conversion_path.exists():
            cprint("Creating scans files...")
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
        subjects: Optional[Path] = None,
        modalities: Optional[list[str]] = None,
        force_new_extraction: bool = False,
        n_procs: Optional[int] = 1,
    ):
        """Convert the images of ADNI.

        Args:
            source_dir: path to the ADNI directory
            clinical_dir: path to the clinical data directory
            dest_dir: path to the BIDS directory
            subjects: Path to list of subjects to process
            modalities: modalities to convert (T1, PET_FDG, PET_AMYLOID, PET_TAU, DWI, FLAIR, fMRI)
            force_new_extraction: if given pre-existing images in the BIDS directory will be erased and extracted again.
        """
        from copy import copy

        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_av45_fbb_pet as adni_av45_fbb
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_dwi as adni_dwi
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet as adni_fdg
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_flair as adni_flair
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap as adni_fmap
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmri as adni_fmri
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_pib_pet as adni_pib
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_t1 as adni_t1
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_tau_pet as adni_tau
        from clinica.iotools.converters.adni_to_bids.adni_utils import get_subjects_list
        from clinica.utils.stream import cprint

        modalities = modalities or self.get_modalities_supported()
        adni_merge = load_clinical_csv(clinical_dir, "ADNIMERGE")
        subjects = get_subjects_list(source_dir, clinical_dir, subjs_list_path)

        # Create the output folder if is not already existing
        (dest_dir / "conversion_info").mkdir(parents=True, exist_ok=True)
        version_number = len([f for f in (dest_dir / "conversion_info").iterdir()])
        conversion_dir = dest_dir / "conversion_info" / f"v{version_number}"
        conversion_dir.mkdir(parents=True, exist_ok=True)
        cprint(f"Destination folder = {dest_dir}", lvl="debug")

        converters = {
            "T1": [adni_t1.convert_adni_t1],
            "PET_FDG": [
                adni_fdg.convert_adni_fdg_pet,
                adni_fdg.convert_adni_fdg_pet_uniform,
            ],
            "PET_AMYLOID": [
                adni_pib.convert_adni_pib_pet,
                adni_av45_fbb.convert_adni_av45_fbb_pet,
            ],
            "PET_TAU": [adni_tau.convert_adni_tau_pet],
            "DWI": [adni_dwi.convert_adni_dwi],
            "FLAIR": [adni_flair.convert_adni_flair],
            "fMRI": [adni_fmri.convert_adni_fmri],
            "FMAP": [adni_fmap.convert_adni_fmap],
        }

        for modality in modalities:
            if modality not in converters:
                raise Exception(f"{modality} is not a valid input modality")
            for converter in converters[modality]:
                converter(
                    source_dir=source_dir,
                    csv_dir=clinical_dir,
                    destination_dir=dest_dir,
                    conversion_dir=conversion_dir,
                    subjects=subjects,
                    mod_to_update=force_new_extraction,
                    n_procs=n_procs,
                )


def _get_bids_subjs_info(
    clinical_data_dir: Path,
    out_path: Path,
    subjects: Optional[Path] = None,
) -> tuple[list[str], list[Path]]:
    from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv

    # Read optional list of participants.
    if subjects:
        subjects = [subject for subject in subjects.read_text().split("\n") if subject]
    # Load all participants from ADNIMERGE.
    adni_merge = load_clinical_csv(clinical_data_dir, "ADNIMERGE")
    participants = adni_merge["PTID"].unique()
    # Filter participants if requested.
    participants = sorted(participants & subjects if subjects else participants)
    # Compute their corresponding BIDS IDs and paths.
    bids_ids = [f"sub-ADNI{p.replace('_', '')}" for p in participants]
    bids_paths = [out_path / bids_id for bids_id in bids_ids]

    return bids_ids, bids_paths
