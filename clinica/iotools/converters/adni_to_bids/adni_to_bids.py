from os import PathLike
from pathlib import Path
from typing import List, Optional

from clinica.iotools.abstract_converter import Converter


def get_bids_subjs_info(
    clinical_data_dir: str,
    out_path: str,
    subjects_list_path: Optional[str] = None,
):
    from os import path

    from pandas import read_csv

    # Read optional list of participants.
    subjects_list = (
        set([line.rstrip("\n") for line in open(subjects_list_path)])
        if subjects_list_path
        else None
    )

    # Load all participants from ADNIMERGE.
    adni_merge_path = path.join(clinical_data_dir, "ADNIMERGE.csv")
    participants = set(
        read_csv(adni_merge_path, sep=",", usecols=["PTID"], squeeze=True).unique()
    )

    # Filter participants if requested.
    participants = sorted(
        participants & subjects_list if subjects_list else participants
    )

    # Compute their corresponding BIDS IDs and paths.
    bids_ids = [f"sub-ADNI{p.replace('_', '')}" for p in participants]
    bids_paths = [path.join(out_path, bids_id) for bids_id in bids_ids]

    return bids_ids, bids_paths


class AdniToBidsConverter(Converter):
    study_name: str = "ADNI"
    # link: str = "http://adni.loni.usc.edu"
    # description: str = (
    #    "ADNI is a global research effort that actively supports the investigation "
    #    "and development of treatments that slow or stop the progression of Alzheimer's "
    #    "disease (AD).This multisite, longitudinal study assesses clinical, imaging, "
    #    "genetic and biospecimen biomarkers through the process of normal aging to mild "
    #    "cognitive impairment (MCI) and AD dementia.With established, standardized methods "
    #    "for imaging and biomarker collection and analysis, ADNI facilitates a way for "
    #    "scientists to conduct cohesive research and share compatible data with other "
    #    "researchers around the world."
    # )

    def __init__(
        self,
        source_dataset: PathLike,
        destination_dataset: PathLike,
        clinical_data_directory: PathLike,
        xml_data_directory: Optional[PathLike] = None,
        clinical_data_only: bool = False,
        force_new_extraction: bool = False,
    ):
        from clinica.utils.exceptions import ClinicaParserError

        super().__init__(
            source_dataset,
            destination_dataset,
            clinical_data_directory,
            clinical_data_only,
        )
        self.force_new_extraction = force_new_extraction
        self.xml_data_directory = (
            Path(xml_data_directory) if xml_data_directory else None
        )
        if self.clinical_data_only and self.force_new_extraction:
            raise ClinicaParserError(
                "[ADNI2BIDS] Arguments `clinical_data_only` and `force_new_extraction` "
                "are mutually exclusive."
            )

    @classmethod
    def get_modalities_supported(cls) -> List[str]:
        """Return a list of modalities supported.

        Returns: a list containing the modalities supported by the converter
        (T1, PET_FDG, PET_AMYLOID, PET_TAU, DWI, FLAIR, fMRI)
        """
        return ["T1", "PET_FDG", "PET_AMYLOID", "PET_TAU", "DWI", "FLAIR", "fMRI"]

    def convert_clinical_data(
        self, subjects_list_path: Optional[PathLike] = None
    ) -> None:
        """Convert the clinical data of ADNI specified into the file
        clinical_specifications_adni.xlsx.

        Parameters
        ----------
        subjects_list_path : str, optional
            If specified, restrict the processing to the subjects specified
            in the corresponding file.
            If not specified, all subjects will be handled.
        """
        import os
        from os import path

        import clinica.iotools.bids_utils as bids
        import clinica.iotools.converters.adni_to_bids.adni_utils as adni_utils
        from clinica.utils.stream import cprint

        from .adni_json import create_json_metadata

        clinic_specs_path = path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            "data",
            "clinical_specifications_adni",
        )

        conversion_path = self.destination_dataset / "conversion_info"

        if self.clinical_data_only:
            bids_ids, bids_subjs_paths = get_bids_subjs_info(
                clinical_data_dir=self.clinical_data_directory,
                out_path=self.destination_dataset,
                subjects_list_path=subjects_list_path,
            )
        else:
            bids_ids = bids.get_bids_subjs_list(self.destination_dataset)
            bids_subjs_paths = bids.get_bids_subjs_paths(self.destination_dataset)

        # -- Creation of participant.tsv --
        cprint("Creating participants.tsv...")
        participants_df = bids.create_participants_df(
            self.study_name,
            clinic_specs_path,
            self.clinical_data_directory,
            bids_ids,
        )

        # Replace the original values with the standard defined by the AramisTeam
        participants_df["sex"] = participants_df["sex"].replace("Male", "M")
        participants_df["sex"] = participants_df["sex"].replace("Female", "F")

        # Correction of diagnosis_sc for ADNI3 participants
        participants_df = adni_utils.correct_diagnosis_sc_adni3(
            self.clinical_data_directory, participants_df
        )

        participants_df.to_csv(
            self.destination_dataset / "participants.tsv",
            sep="\t",
            index=False,
            encoding="utf-8",
        )

        # -- Creation of sessions.tsv --
        cprint("Creating sessions files...")
        adni_utils.create_adni_sessions_dict(
            bids_ids, clinic_specs_path, self.clinical_data_directory, bids_subjs_paths
        )

        # -- Creation of scans files --
        if conversion_path.exists():
            cprint("Creating scans files...")
            adni_utils.create_adni_scans_files(conversion_path, bids_subjs_paths)

        if self.xml_data_directory:
            if self.xml_data_directory.exists():
                create_json_metadata(
                    bids_subjs_paths, bids_ids, self.xml_data_directory
                )
            else:
                cprint(
                    msg=(
                        f"Clinica was unable to find {self.xml_data_directory}, "
                        "skipping xml metadata extraction."
                    ),
                    lvl="warning",
                )
        super().convert_clinical_data(subjects_list_path)

    def convert_images(
        self,
        subjects_list_path: Optional[PathLike] = None,
        modalities: Optional[List[str]] = None,
    ):
        """Convert the images of ADNI.

        Parameters
        ----------
        subjects_list_path : str, optional
            list of subjects to process

        modalities : list of str, optional
            List of modalities to convert among:
            T1, PET_FDG, PET_AMYLOID, PET_TAU, DWI, FLAIR, fMRI
            If not specified, all modalities will be handled.
        """
        from copy import copy

        import pandas as pd

        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_av45_fbb_pet as adni_av45_fbb
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_dwi as adni_dwi
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet as adni_fdg
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_flair as adni_flair
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmri as adni_fmri
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_pib_pet as adni_pib
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_t1 as adni_t1
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_tau_pet as adni_tau
        from clinica.utils.stream import cprint

        modalities = modalities or self.get_modalities_supported()

        adni_merge = pd.read_csv(self.clinical_data_directory / "ADNIMERGE.csv")

        # Load a file with subjects list or compute all the subjects
        if subjects_list_path:
            cprint("Loading a subjects lists provided by the user...")
            subjs_list = [line.rstrip("\n") for line in open(subjects_list_path)]
            subjs_list_copy = copy(subjs_list)

            # Check that there are no errors in subjs_list given by the user
            for subj in subjs_list_copy:
                adnimerge_subj = adni_merge[adni_merge.PTID == subj]
                if len(adnimerge_subj) == 0:
                    cprint(
                        msg=f"Subject with PTID {subj} does not exist. Please check your subjects list.",
                        lvl="warning",
                    )
                    subjs_list.remove(subj)
            del subjs_list_copy

        else:
            cprint("Using all the subjects contained into the ADNIMERGE.csv file...")
            subjs_list = list(adni_merge["PTID"].unique())

        conversion_info = self.destination_dataset / "conversion_info"
        conversion_info.mkdir(exist_ok=True)
        version_number = len(list(conversion_info.iterdir()))
        conversion_dir = conversion_info / f"v{version_number}"
        conversion_dir.mkdir()
        cprint(self.destination_dataset, lvl="debug")

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
        }

        for modality in modalities:
            if modality not in converters:
                raise Exception(f"{modality} is not a valid input modality")
            for converter in converters[modality]:
                converter(
                    source_dir=self.source_dataset,
                    csv_dir=self.clinical_data_directory,
                    destination_dir=self.destination_dataset,
                    conversion_dir=conversion_dir,
                    subjects=subjs_list,
                    mod_to_update=self.force_new_extraction,
                )
