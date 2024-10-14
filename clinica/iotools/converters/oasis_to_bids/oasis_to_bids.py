"""Convert OASIS dataset (https://sites.wustl.edu/oasisbrains/) to BIDS."""

from pathlib import Path
from typing import Optional

import nibabel as nb
import numpy as np

from clinica.iotools.abstract_converter import Converter
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
    from clinica.iotools.bids_utils import StudyName
    from clinica.iotools.converters.factory import get_converter_name
    from clinica.utils.stream import cprint

    from ..utils import validate_input_path

    path_to_dataset = validate_input_path(path_to_dataset)
    bids_dir = validate_input_path(bids_dir, check_exist=False)
    path_to_clinical = validate_input_path(path_to_clinical)
    if subjects:
        cprint(
            (
                f"Subject filtering is not yet implemented in {get_converter_name(StudyName.OASIS)} converter. "
                "All subjects available will be converted."
            ),
            lvl="warning",
        )
    OasisToBids().convert(
        path_to_dataset,
        bids_dir,
        path_to_clinical,
        n_procs=n_procs,
    )


class OasisToBids(Converter):
    def convert(
        self,
        source_dir: Path,
        destination_dir: Path,
        clinical_data_dir: Path,
        n_procs: Optional[int] = 1,
    ):
        self.convert_images(source_dir, destination_dir, n_procs=n_procs)
        self.convert_clinical_data(clinical_data_dir, destination_dir)

    def convert_clinical_data(self, clinical_data_dir: Path, bids_dir: Path):
        """Convert the clinical data defined inside the clinical_specifications.xlx into BIDS.

        Args:
            clinical_data_dir: path to the folder with the original clinical data
            bids_dir: path to the BIDS directory
        """
        from clinica.iotools.bids_utils import get_bids_subjs_list
        from clinica.utils.stream import cprint

        cprint("Converting clinical data...", lvl="info")
        bids_ids = get_bids_subjs_list(bids_dir)
        self._create_participants_tsv(clinical_data_dir, bids_dir, bids_ids)
        self._create_sessions_tsv(clinical_data_dir, bids_dir, bids_ids)
        self._create_scans_tsv(bids_dir)
        self._create_modality_agnostic_files(bids_dir)

    def _create_participants_tsv(
        self,
        clinical_data_dir: Path,
        bids_dir: Path,
        bids_ids: list[str],
    ):
        from clinica.iotools.bids_utils import StudyName, create_participants_df

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

    def _create_sessions_tsv(
        self,
        clinical_data_dir: Path,
        bids_dir: Path,
        bids_ids: list[str],
    ) -> None:
        from .oasis_to_bids_utils import create_sessions_dict, write_sessions_tsv

        sessions_dict = create_sessions_dict(
            clinical_data_dir=clinical_data_dir,
            clinical_specifications_folder=Path(__file__).parents[1] / "specifications",
            bids_ids=bids_ids,
        )

        write_sessions_tsv(bids_dir, sessions_dict)

    def _create_scans_tsv(
        self,
        bids_dir: Path,
    ) -> None:
        from clinica.iotools.converters.oasis_to_bids.oasis_to_bids_utils import (
            write_scans,
        )

        write_scans(bids_dir)

    def _create_modality_agnostic_files(self, bids_dir: Path):
        from clinica.iotools.bids_utils import StudyName, write_modality_agnostic_files

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
        from clinica.iotools.bids_utils import StudyName, bids_id_factory
        from clinica.utils.stream import cprint

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
        try:
            img_file_path = next(t1_folder.glob("*.img"))
        except StopIteration:
            msg = f"No file ending in .img found in {t1_folder}."
            cprint(msg, lvl="error")
            raise FileNotFoundError(msg)
        nb.save(
            _get_image_with_good_orientation(img_file_path),
            session_folder / "anat" / f"{participant_id}_ses-M000_T1w.nii.gz",
        )

    def convert_images(
        self, source_dir: Path, dest_dir: Path, n_procs: Optional[int] = 1
    ):
        """Convert T1w images to BIDS.

        Parameters
        ----------
        source_dir: path to the OASIS dataset

        dest_dir: path to the BIDS directory

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

        if not dest_dir.exists():
            dest_dir.mkdir(parents=True)

        subjects_folders = [
            path
            for path in source_dir.rglob("OAS1_*")
            if path.is_dir() and path.name.endswith("_MR1")
        ]
        func = partial(self.convert_single_subject, dest_dir=dest_dir)
        # If n_procs==1 do not rely on a Process Pool to enable classical debugging
        if n_procs == 1:
            for folder in subjects_folders:
                func(folder)
        else:
            with Pool(processes=n_procs) as pool:
                pool.map(func, subjects_folders)


def _get_image_with_good_orientation(image_path: Path) -> nb.Nifti1Image:
    # First, convert to Nifti so that we can extract the s_form with NiBabel
    # (NiBabel creates an 'Spm2AnalyzeImage' object that does not contain 'get_sform' method
    img_with_wrong_orientation_analyze = nb.load(image_path)

    # OASIS-1 images have the same header but sform is incorrect
    # To solve this issue, we use header from images converted with FreeSurfer
    # to generate a 'clean hard-coded' header
    # affine:
    # [[   0.    0.   -1.   80.]
    #  [   1.    0.    0. -128.]
    #  [   0.    1.    0. -128.]
    #  [   0.    0.    0.    1.]]
    # fmt: off
    affine = np.array(
        [
            0, 0, -1, 80,
            1, 0, 0, -128,
            0, 1, 0, -128,
            0, 0, 0, 1
        ]
    ).reshape(4, 4)
    # fmt: on
    s_form = affine.astype(np.int16)

    hdr = nb.Nifti1Header()
    hdr.set_data_shape((256, 256, 160))
    hdr.set_data_dtype(np.int16)
    hdr["bitpix"] = 16
    hdr.set_sform(s_form, code="scanner")
    hdr.set_qform(s_form, code="scanner")
    hdr["extents"] = 16384
    hdr["xyzt_units"] = 10

    img_with_good_orientation_nifti = nb.Nifti1Image(
        np.round(img_with_wrong_orientation_analyze.get_fdata(dtype="float32")).astype(
            np.int16
        ),
        s_form,
        header=hdr,
    )
    # Header correction to obtain dim0 = 3
    return nb.funcs.four_to_three(img_with_good_orientation_nifti)[0]
