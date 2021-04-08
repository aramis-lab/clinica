# coding: utf-8

from clinica.iotools.abstract_converter import Converter


class AdniToBids(Converter):
    def get_modalities_supported(self):
        """Return a list of modalities supported.

        Returns: a list containing the modalities supported by the converter
        (T1, PET_FDG, PET_AMYLOID, PET_TAU, DWI, FLAIR, fMRI)
        """
        return ["T1", "PET_FDG", "PET_AMYLOID", "PET_TAU", "DWI", "FLAIR", "fMRI"]

    def check_adni_dependencies(self):
        """Check the dependencies of ADNI converter."""
        from clinica.utils.check_dependency import check_dcm2nii, check_dcm2niix

        check_dcm2nii()
        check_dcm2niix()

    def convert_clinical_data(self, clinical_data_dir, out_path):
        """Convert the clinical data of ADNI specified into the file clinical_specifications_adni.xlsx.

        Args:
            clinical_data_dir:  path to the clinical data directory
            out_path: path to the BIDS directory
        """
        import os
        from os import path

        import pandas as pd

        import clinica.iotools.bids_utils as bids
        import clinica.iotools.converters.adni_to_bids.adni_utils as adni_utils
        from clinica.utils.stream import cprint

        clinic_specs_path = path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            "data",
            "clinical_specifications_adni.xlsx",
        )
        try:
            os.path.exists(out_path)
        except IOError:
            print("BIDS folder not found.")
            raise

        bids_ids = bids.get_bids_subjs_list(out_path)
        bids_subjs_paths = bids.get_bids_subjs_paths(out_path)
        conversion_path = path.join(out_path, "conversion_info")

        if not bids_ids:
            adni_merge_path = path.join(clinical_data_dir, "ADNIMERGE.csv")
            adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=",")
            bids_ids = [
                "sub-ADNI" + subj.replace("_", "")
                for subj in list(adni_merge.PTID.unique())
            ]
            bids_subjs_paths = [path.join(out_path, subj) for subj in bids_ids]

        # -- Creation of modality agnostic files --
        cprint("Creating modality agnostic files...")
        bids.write_modality_agnostic_files("ADNI", out_path)

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
        cprint("Creating scans files...")
        adni_utils.create_adni_scans_files(
            conversion_path,
            bids_subjs_paths,
        )

    def convert_images(
        self,
        source_dir,
        clinical_dir,
        dest_dir,
        subjs_list_path=None,
        modalities=["T1", "PET_FDG", "PET_AMYLOID", "PET_TAU", "DWI", "FLAIR", "fMRI"],
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
        from colorama import Fore

        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_av45_fbb_pet as adni_av45_fbb
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_dwi as adni_dwi
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet as adni_fdg
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_flair as adni_flair
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmri as adni_fmri
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_pib_pet as adni_pib
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_t1 as adni_t1
        import clinica.iotools.converters.adni_to_bids.adni_modalities.adni_tau_pet as adni_tau
        from clinica.utils.stream import cprint

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
                        f"{Fore.RED}Subject with PTID {subj} does not exist. Please check your subjects list.{Fore.RESET}"
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
        cprint(dest_dir)

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
