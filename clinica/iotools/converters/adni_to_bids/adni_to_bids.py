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


class AdniToBids(Converter):
    @classmethod
    def get_modalities_supported(cls) -> List[str]:
        """Return a list of modalities supported.

        Returns: a list containing the modalities supported by the converter
        (T1, PET_FDG, PET_AMYLOID, PET_TAU, DWI, FLAIR, fMRI)
        """
        return ["T1", "PET_FDG", "PET_AMYLOID", "PET_TAU", "DWI", "FLAIR", "fMRI"]

    @classmethod
    def check_adni_dependencies(cls) -> None:
        """Check the dependencies of ADNI converter."""
        from clinica.utils.check_dependency import check_dcm2niix

        check_dcm2niix()

    def convert_clinical_data(
        self,
        clinical_data_dir: str,
        out_path: str,
        clinical_data_only: bool = False,
        subjects_list_path: Optional[str] = None,
        xml_path: Optional[str] = None,
    ):
        """Convert the clinical data of ADNI specified into the file clinical_specifications_adni.xlsx.

        Args:
            clinical_data_dir:  path to the clinical data directory
            out_path: path to the BIDS directory
            clinical_data_only: process clinical data only
            subjects_list_path: restrict processing to this manifest of subjects
            xml_path: path to the XML metadata files
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
        try:
            os.path.exists(out_path)
        except IOError:
            print("BIDS folder not found.")
            raise

        conversion_path = path.join(out_path, "conversion_info")

        if clinical_data_only:
            bids_ids, bids_subjs_paths = get_bids_subjs_info(
                clinical_data_dir=clinical_data_dir,
                out_path=out_path,
                subjects_list_path=subjects_list_path,
            )
        else:
            bids_ids = bids.get_bids_subjs_list(out_path)
            bids_subjs_paths = bids.get_bids_subjs_paths(out_path)

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
        bids.write_modality_agnostic_files(
            study_name="ADNI", readme_data=readme_data, bids_dir=out_path
        )

        # -- Creation of participant.tsv --
        cprint("Creating participants.tsv...")
        participants_df = bids.create_participants_df(
            "ADNI", clinic_specs_path, clinical_data_dir, bids_ids
        )

        # Replace the original values with the standard defined by the AramisTeam
        participants_df["sex"] = participants_df["sex"].replace("Male", "M")
        participants_df["sex"] = participants_df["sex"].replace("Female", "F")

        # Correction of diagnosis_sc for ADNI3 participants
        participants_df = adni_utils.correct_diagnosis_sc_adni3(
            clinical_data_dir, participants_df
        )

        participants_df.to_csv(
            path.join(out_path, "participants.tsv"),
            sep="\t",
            index=False,
            encoding="utf-8",
        )

        # -- Creation of sessions.tsv --
        cprint("Creating sessions files...")
        adni_utils.create_adni_sessions_dict(
            bids_ids, clinic_specs_path, clinical_data_dir, bids_subjs_paths
        )

        # -- Creation of scans files --
        if os.path.exists(conversion_path):
            cprint("Creating scans files...")
            adni_utils.create_adni_scans_files(conversion_path, bids_subjs_paths)

        if xml_path is not None:
            if os.path.exists(xml_path):
                create_json_metadata(bids_subjs_paths, bids_ids, xml_path)
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
        source_dir,
        clinical_dir,
        dest_dir,
        subjs_list_path=None,
        modalities=None,
        force_new_extraction=False,
    ):
        """Convert the images of ADNI.

        Args:
            source_dir: path to the ADNI directory
            clinical_dir: path to the clinical data directory
            dest_dir: path to the BIDS directory
            subjs_list_path: list of subjects to process
            modalities: modalities to convert (T1, PET_FDG, PET_AMYLOID, PET_TAU, DWI, FLAIR, fMRI)
            force_new_extraction: if given pre-existing images in the BIDS directory will be erased and extracted again.
        """
        import os
        from copy import copy
        from os import path

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

        adni_merge_path = path.join(clinical_dir, "ADNIMERGE.csv")
        adni_merge = pd.read_csv(adni_merge_path)

        # Load a file with subjects list or compute all the subjects
        if subjs_list_path is not None:
            cprint("Loading a subjects lists provided by the user...")
            subjs_list = [line.rstrip("\n") for line in open(subjs_list_path)]
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

        # Create the output folder if is not already existing
        os.makedirs(dest_dir, exist_ok=True)
        os.makedirs(path.join(dest_dir, "conversion_info"), exist_ok=True)
        version_number = len(os.listdir(path.join(dest_dir, "conversion_info")))
        conversion_dir = path.join(dest_dir, "conversion_info", f"v{version_number}")
        os.makedirs(conversion_dir)
        cprint(dest_dir, lvl="debug")

        converters = {
            "T1": [adni_t1.convert_adni_t1],
            "PET_FDG": [adni_fdg.convert_adni_fdg_pet],
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
                    source_dir,
                    clinical_dir,
                    dest_dir,
                    conversion_dir,
                    subjs_list,
                    force_new_extraction,
                )
