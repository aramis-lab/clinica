"""Convert OASIS dataset (http://www.oasis-brains.org/) to BIDS."""

from clinica.iotools.abstract_converter import Converter


class OasisToBids(Converter):
    def convert_clinical_data(self, clinical_data_dir, bids_dir):
        """Convert the clinical data defined inside the clinical_specifications.xlx into BIDS.

        Args:
            clinical_data_dir: path to the folder with the original clinical data
            bids_dir: path to the BIDS directory
        """
        import os
        from os import path

        import numpy as np

        import clinica.iotools.bids_utils as bids
        from clinica.utils.stream import cprint

        cprint("Converting clinical data...")
        bids_ids = bids.get_bids_subjs_list(bids_dir)

        iotools_folder = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        clinic_specs_path = path.join(iotools_folder, "data", "clinical_specifications")

        # --Create participants.tsv--
        participants_df = bids.create_participants_df(
            "OASIS", clinic_specs_path, clinical_data_dir, bids_ids
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
            path.join(bids_dir, "participants.tsv"),
            sep="\t",
            index=False,
            encoding="utf-8",
        )

        # --Create sessions files--
        sessions_dict = bids.create_sessions_dict_OASIS(
            clinical_data_dir, bids_dir, "OASIS", clinic_specs_path, bids_ids, "ID"
        )
        for y in bids_ids:
            if sessions_dict[y]["M00"]["diagnosis"] > 0:
                sessions_dict[y]["M00"]["diagnosis"] = "AD"
            else:
                sessions_dict[y]["M00"]["diagnosis"] = "CN"

        bids.write_sessions_tsv(bids_dir, sessions_dict)

        # --Create scans files--
        # Note: We have no scans information for OASIS
        scans_dict = bids.create_scans_dict(
            clinical_data_dir,
            "OASIS",
            clinic_specs_path,
            bids_ids,
            "ID",
            "",
            sessions_dict,
        )
        bids.write_scans_tsv(bids_dir, bids_ids, scans_dict)

        # -- Creation of modality agnostic files --
        readme_data = {
            "link": "https://www.oasis-brains.org/#access",
            "desc": (
                "This set consists of a cross-sectional collection of 416 subjects aged 18 to 96. For each subject, 3 "
                "or 4 individual T1-weighted MRI scans obtained in single scan sessions are included. The subjects are "
                "all right-handed and include both men and women. 100 of the included subjects over the age of 60 have "
                "been clinically diagnosed with very mild to moderate Alzheimer’s disease (AD). Additionally, a "
                "reliability data set is included containing 20 nondemented subjects imaged on a subsequent visit "
                "within 90 days of their initial session."
            ),
        }
        bids.write_modality_agnostic_files(
            study_name="OASIS-1", readme_data=readme_data, bids_dir=bids_dir
        )

    @staticmethod
    def convert_single_subject(subj_folder, dest_dir):
        import os
        from os import path
        from pathlib import Path

        import nibabel as nb
        import numpy as np

        t1_folder = path.join(subj_folder, "PROCESSED", "MPRAGE", "SUBJ_111")
        subj_id = path.basename(subj_folder)
        print("Converting ", subj_id)
        numerical_id = (subj_id.split("_"))[1]
        bids_id = "sub-OASIS1" + str(numerical_id)
        bids_subj_folder = path.join(dest_dir, bids_id)
        if not path.isdir(bids_subj_folder):
            os.mkdir(bids_subj_folder)

        session_folder = path.join(bids_subj_folder, "ses-M00")
        if not os.path.isdir(session_folder):
            os.mkdir(path.join(session_folder))
            os.mkdir(path.join(session_folder, "anat"))

        # In order do convert the Analyze format to Nifti the path to the .img file is required
        img_file_path = str(next(Path(t1_folder).glob("*.img")))
        output_path = path.join(session_folder, "anat", bids_id + "_ses-M00_T1w.nii.gz")

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

    def convert_images(self, source_dir, dest_dir):
        """Convert T1w images to BIDS.

        Args:
            source_dir: path to the OASIS dataset
            dest_dir: path to the BIDS directory

        Note:
            Previous version of this method used mri_convert from FreeSurfer to convert
            Analyze data from OASIS-1. To remove this strong dependency, NiBabel is used instead.
        """
        import os
        from functools import partial
        from multiprocessing import Pool, cpu_count
        from pathlib import Path

        if not os.path.isdir(dest_dir):
            os.mkdir(dest_dir)

        subjs_folders = [
            path
            for path in Path(source_dir).rglob("OAS1_*")
            if path.is_dir() and path.name.endswith("_MR1")
        ]

        with Pool(processes=max(cpu_count() - 1, 1)) as pool:
            func = partial(self.convert_single_subject, dest_dir=dest_dir)
            pool.map(func, subjs_folders)
