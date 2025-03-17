"""Convert OASIS dataset (https://sites.wustl.edu/oasisbrains/) to BIDS."""

from pathlib import Path
from typing import Optional

import nibabel as nb
import numpy as np

from clinica.converters.abstract_converter import Converter
from clinica.utils.filemanip import UserProvidedPath

__all__ = ["convert"]


def convert(
    path_to_dataset: UserProvidedPath,
    bids_dir: UserProvidedPath,
    path_to_clinical: UserProvidedPath,
    subjects: Optional[UserProvidedPath] = None,
    n_procs: Optional[int] = 1,
    **kwargs,
):
    from clinica.converters.factory import get_converter_name
    from clinica.converters.study_models import StudyName
    from clinica.utils.stream import cprint

    from .._utils import validate_input_path

    path_to_dataset = validate_input_path(path_to_dataset)
    bids_dir = validate_input_path(bids_dir, check_exist=False)
    path_to_clinical = validate_input_path(path_to_clinical)
    if subjects:
        """
        cprint(
            (
                f"Subject filtering is not yet implemented in {get_converter_name(StudyName.OASIS)} converter. "
                "All subjects available will be converted."
            ),
            lvl="warning",
        )
        """
        subjects = validate_input_path(subjects)

    OasisToBids().convert(
        path_to_dataset,
        bids_dir,
        path_to_clinical,
        subjects=subjects,
        n_procs=n_procs,
    )


class OasisToBids(Converter):
    def convert(
        self,
        source_dir: Path,
        destination_dir: Path,
        clinical_data_dir: Path,
        subjects: Optional[UserProvidedPath] = None,
        n_procs: Optional[int] = 1,
    ):
        self._create_modality_agnostic_files(destination_dir)
        self.convert_images(
            source_dir, destination_dir, subjects=subjects, n_procs=n_procs
        )
        self.convert_clinical_data(clinical_data_dir, destination_dir)

    def convert_clinical_data(self, clinical_data_dir: Path, bids_dir: Path):
        """Convert the clinical data defined inside the clinical_specifications.xlx into BIDS.

        Args:
            clinical_data_dir: path to the folder with the original clinical data
            bids_dir: path to the BIDS directory
        """
        from clinica.dataset import get_subjects_from_bids_dataset
        from clinica.utils.stream import cprint

        cprint("Converting clinical data...", lvl="info")
        bids_ids = get_subjects_from_bids_dataset(bids_dir)
        self._create_participants_tsv(clinical_data_dir, bids_dir, bids_ids)
        self._create_sessions_tsv(clinical_data_dir, bids_dir, bids_ids)
        self._create_scans_tsv(bids_dir)

    def _create_participants_tsv(
        self,
        clinical_data_dir: Path,
        bids_dir: Path,
        bids_ids: list[str],
    ):
        from clinica.converters._utils import create_participants_df
        from clinica.converters.study_models import StudyName

        participants_df = create_participants_df(
            study_name=StudyName.OASIS,
            clinical_specifications_folder=Path(__file__).parents[1] / "specifications",
            clinical_data_dir=clinical_data_dir,
            bids_ids=bids_ids,
        )
        # Replace the values of the diagnosis_bl column
        participants_df["diagnosis_bl"].replace([0.0, np.nan], "CN", inplace=True)
        participants_df["diagnosis_bl"].replace(
            [0.5, 1.0, 1.5, 2.0], "AD", inplace=True
        )
        # Following line has no sense
        # participants_df['diagnosis_bl'].replace(participants_df['diagnosis_bl']>0.0, 'AD', inplace=True)
        participants_df = participants_df.fillna("n/a")
        participants_df.to_csv(
            bids_dir / "participants.tsv",
            sep="\t",
            index=False,
            encoding="utf-8",
        )

    @staticmethod
    def _create_sessions_tsv(
        clinical_data_dir: Path,
        bids_dir: Path,
        bids_ids: list[str],
    ) -> None:
        from ._utils import create_sessions_df, write_sessions_tsv

        sessions_df = create_sessions_df(
            clinical_data_dir=clinical_data_dir,
            clinical_specifications_folder=Path(__file__).parents[1] / "specifications",
            bids_ids=bids_ids,
        )

        write_sessions_tsv(bids_dir, sessions_df)

    @staticmethod
    def _create_scans_tsv(bids_dir: Path) -> None:
        from ._utils import write_scans_tsv

        write_scans_tsv(bids_dir)

    @staticmethod
    def _create_modality_agnostic_files(bids_dir: Path):
        from clinica.converters._utils import write_modality_agnostic_files
        from clinica.converters.study_models import StudyName

        if not bids_dir.exists():
            bids_dir.mkdir(parents=True)
        readme_data = {
            "link": "https://sites.wustl.edu/oasisbrains/#access",
            "desc": (
                "This set consists of a cross-sectional collection of 416 subjects aged 18 to 96. For each subject, 3 "
                "or 4 individual T1-weighted MRI scans obtained in single scan sessions are included. The subjects are "
                "all right-handed and include both men and women. 100 of the included subjects over the age of 60 have "
                "been clinically diagnosed with very mild to moderate Alzheimer’s disease (AD). Additionally, a "
                "reliability data set is included containing 20 nondemented subjects imaged on a subsequent visit "
                "within 90 days of their initial session."
            ),
        }
        write_modality_agnostic_files(
            study_name=StudyName.OASIS,
            readme_data=readme_data,
            bids_dir=bids_dir,
        )

    @staticmethod
    def convert_single_subject(subj_folder: Path, dest_dir: Path):
        from clinica.cmdline import setup_clinica_logging
        from clinica.converters.study_models import StudyName, bids_id_factory
        from clinica.utils.stream import cprint

        from ._utils import get_first_image, get_image_with_good_orientation

        # This function is executed in a multiprocessing context
        # such that we need to re-configure the clinica logger in the child processes.
        # Note that logging messages could easily be lost (for example when logging
        # to a file from two different processes). A better solution would be to
        # implement a logging process consuming logging messages from a multiprocessing.Queue...
        setup_clinica_logging("INFO")

        t1_folder = subj_folder / "PROCESSED" / "MPRAGE" / "SUBJ_111"
        cprint(f"Converting {subj_folder.name}", lvl="info")
        participant_id = bids_id_factory(StudyName.OASIS).from_original_study_id(
            subj_folder.name
        )
        bids_subj_folder = dest_dir / participant_id
        if not bids_subj_folder.is_dir():
            bids_subj_folder.mkdir(parents=True)
        session_folder = bids_subj_folder / "ses-M000"
        if not session_folder.is_dir():
            (session_folder / "anat").mkdir(parents=True, exist_ok=True)

        # In order do convert the Analyze format to Nifti the path to the .img file is required
        nb.save(
            get_image_with_good_orientation(get_first_image(t1_folder)),
            session_folder / "anat" / f"{participant_id}_ses-M000_T1w.nii.gz",
        )

    def convert_images(
        self,
        source_dir: Path,
        dest_dir: Path,
        subjects: Optional[Path] = None,
        n_procs: Optional[int] = 1,
    ):
        """Convert T1w images to BIDS.

        Parameters
        ----------
        source_dir: path to the OASIS dataset

        dest_dir: path to the BIDS directory

        subjects: path to list of subjects to process

        n_procs : int, optional
            The requested number of processes.
            If specified, it should be between 1 and the number of available CPUs.
            Default=1.

        Notes
        -----
        Previous version of this method used mri_convert from FreeSurfer to convert
        Analyze data from OASIS-1. To remove this strong dependency, NiBabel is used instead.
        """
        from functools import partial
        from multiprocessing import Pool

        from clinica.utils.stream import cprint

        from ._utils import get_subjects_list

        if not dest_dir.exists():
            dest_dir.mkdir(parents=True)

        subjects_folders = get_subjects_list(source_dir, subjects)

        func = partial(self.convert_single_subject, dest_dir=dest_dir)
        # If n_procs==1 do not rely on a Process Pool to enable classical debugging
        if n_procs == 1:
            for folder in subjects_folders:
                func(folder)
        else:
            with Pool(processes=n_procs) as pool:
                pool.map(func, subjects_folders)
