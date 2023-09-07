"""Convert OASIS dataset (http://www.oasis-brains.org/) to BIDS."""
from os import PathLike
from pathlib import Path
from typing import List, Optional

from clinica.iotools.abstract_converter import Converter


class OasisToBidsConverter(Converter):
    study_name: str = "OASIS"
    link: str = "https://www.oasis-brains.org/#access"
    description: str = (
        "This set consists of a cross-sectional collection of 416 subjects aged 18 to 96. "
        "For each subject, 3 or 4 individual T1-weighted MRI scans obtained in single scan "
        "sessions are included. The subjects are all right-handed and include both men and "
        "women. 100 of the included subjects over the age of 60 have been clinically diagnosed "
        "with very mild to moderate Alzheimer’s disease (AD). Additionally, a reliability data "
        "set is included containing 20 nondemented subjects imaged on a subsequent visit within "
        "90 days of their initial session."
    )

    def convert_clinical_data(
        self, subjects_list_path: Optional[PathLike] = None
    ) -> None:
        """Convert the clinical data defined inside the clinical_specifications.xlx into BIDS."""
        import os
        from os import path

        import numpy as np

        import clinica.iotools.bids_utils as bids
        from clinica.utils.stream import cprint

        cprint("Converting clinical data...")
        bids_ids = bids.get_bids_subjs_list(self.destination_dataset)

        iotools_folder = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        clinic_specs_path = path.join(iotools_folder, "data", "clinical_specifications")

        # --Create participants.tsv--
        participants_df = bids.create_participants_df(
            self.study_name,
            clinic_specs_path,
            self.clinical_data_directory,
            bids_ids,
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
            self.destination_dataset / "participants.tsv",
            sep="\t",
            index=False,
            encoding="utf-8",
        )

        # --Create sessions files--
        sessions_dict = bids.create_sessions_dict_OASIS(
            self.clinical_data_directory,
            self.destination_dataset,
            self.study_name,
            clinic_specs_path,
            bids_ids,
            "ID",
        )
        for y in bids_ids:
            if sessions_dict[y]["M000"]["diagnosis"] > 0:
                sessions_dict[y]["M000"]["diagnosis"] = "AD"
            else:
                sessions_dict[y]["M000"]["diagnosis"] = "CN"

        bids.write_sessions_tsv(self.destination_dataset, sessions_dict)

        # --Create scans files--
        # Note: We have no scans information for OASIS
        scans_dict = bids.create_scans_dict(
            self.clinical_data_directory,
            self.study_name,
            clinic_specs_path,
            bids_ids,
            "ID",
            "",
            sessions_dict,
        )
        bids.write_scans_tsv(self.destination_dataset, bids_ids, scans_dict)
        super().convert_clinical_data(subjects_list_path)

    @staticmethod
    def convert_single_subject(subject_folder: Path, destination_dataset: Path):
        import nibabel as nb
        import numpy as np

        t1_folder = subject_folder / "PROCESSED" / "MPRAGE" / "SUBJ_111"
        subject_id = subject_folder.name
        print("Converting ", subject_id)
        numerical_id = (subject_id.split("_"))[1]
        participant_id = "sub-OASIS1" + str(numerical_id)
        bids_subject_folder = destination_dataset / participant_id
        if not bids_subject_folder.is_dir():
            bids_subject_folder.mkdir(parents=True, exist_ok=True)

        session_folder = bids_subject_folder / "ses-M000"
        if not session_folder.is_dir():
            session_folder.mkdir(parents=True, exist_ok=True)
            (session_folder / "anat").mkdir()

        # In order do convert the Analyze format to Nifti the path to the .img file is required
        img_file_path = str(next(t1_folder.glob("*.img")))
        output_path = session_folder / "anat" / f"{participant_id}_ses-M000_T1w.nii.gz"

        # First, convert to Nifti so that we can extract the s_form with NiBabel
        # (NiBabel creates an 'Spm2AnalyzeImage' object that does not contain 'get_sform' method
        img_with_wrong_orientation_analyze = nb.load(img_file_path)

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
            np.round(
                img_with_wrong_orientation_analyze.get_fdata(dtype="float32")
            ).astype(np.int16),
            s_form,
            header=hdr,
        )
        # Header correction to obtain dim0 = 3
        img_with_good_dimension = nb.funcs.four_to_three(
            img_with_good_orientation_nifti
        )[0]

        nb.save(img_with_good_dimension, output_path)

    def convert_images(
        self,
        subjects_list_path: Optional[PathLike] = None,
        modalities: Optional[List[str]] = None,
    ) -> None:
        """Convert T1w images to BIDS.

        Notes
        -----
        Previous version of this method used mri_convert from FreeSurfer to convert
        Analyze data from OASIS-1. To remove this strong dependency, NiBabel is used instead.
        """
        from functools import partial
        from multiprocessing import Pool, cpu_count

        subjects_folders = [
            path
            for path in self.source_dataset.rglob("OAS1_*")
            if path.is_dir() and path.name.endswith("_MR1")
        ]

        with Pool(processes=max(cpu_count() - 1, 1)) as pool:
            func = partial(
                self.convert_single_subject,
                destination_dataset=self.destination_dataset,
            )
            pool.map(func, subjects_folders)
