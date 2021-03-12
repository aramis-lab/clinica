# coding: utf8

"""
Convert OASIS dataset (http://www.oasis-brains.org/) to BIDS.
"""

from clinica.iotools.abstract_converter import Converter


class Oasis3ToBids(Converter):
    def convert_clinical_data(self, clinical_data_dir, bids_dir):
        """
        Convert the clinical data defined inside the clinical_specifications.xlx into BIDS

        Args:
            clinical_data_dir: path to the folder with the original clinical data
            bids_dir: path to the BIDS directory
        """
        from os import path
        import os
        import numpy as np
        import clinica.iotools.bids_utils as bids
        from clinica.utils.stream import cprint
        cprint('Converting clinical data...')
        bids_ids = bids.get_bids_subjs_list(bids_dir)

        iotools_folder = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        clinic_specs_path = path.join(iotools_folder, 'data',
                                      'clinical_specifications.xlsx')

        # --Create participants.tsv--
        participants_df = bids.create_participants_df('OASIS', clinic_specs_path, clinical_data_dir, bids_ids)

        # Replace the values of the diagnosis_bl column
        participants_df['diagnosis_bl'].replace([0.0, np.nan], 'CN', inplace=True)
        participants_df['diagnosis_bl'].replace([0.5, 1.0, 1.5, 2.0], 'AD', inplace=True)
        # Following line has no sense
        # participants_df['diagnosis_bl'].replace(participants_df['diagnosis_bl']>0.0, 'AD', inplace=True)
        participants_df.to_csv(path.join(bids_dir, 'participants.tsv'), sep='\t', index=False, encoding='utf-8')

        # --Create sessions files--
        sessions_dict = bids.create_sessions_dict(clinical_data_dir, 'OASIS', clinic_specs_path, bids_ids, 'ID')
        for y in bids_ids:
            if sessions_dict[y]['M00']['diagnosis'] > 0:
                sessions_dict[y]['M00']['diagnosis'] = 'AD'
            else:
                sessions_dict[y]['M00']['diagnosis'] = 'CN'

        bids.write_sessions_tsv(bids_dir, sessions_dict)

        # --Create scans files--
        # Note: We have no scans information for OASIS
        scans_dict = bids.create_scans_dict(clinical_data_dir, 'OASIS', clinic_specs_path, bids_ids, 'ID')
        bids.write_scans_tsv(bids_dir, bids_ids, scans_dict)

    def convert_images(self, source_dir, dest_dir):
        """
        Convert T1 images to BIDS

        Args:
            source_dir: path to the OASIS dataset
            dest_dir: path to the BIDS directory

        Note:
            Previous version of this method used mri_convert from FreeSurfer to convert Analyze data from OASIS-1.
            To remove this strong dependency, NiBabel is used instead.
        """
        from os import path
        from glob import glob
        import os
        from multiprocessing.dummy import Pool
        from multiprocessing import cpu_count

        def convert_single_subject(subj_folder):
            import os
            import re
            import shutil
            import nibabel as nb
            import numpy as np

            def copy_file(source, target):
                try:
                    shutil.copy(source, target)
                except IOError as e:
                    print("Unable to copy file. %s" % e)

            t1_folder = path.join(subj_folder, 'anat1')
            subj_id = os.path.basename(subj_folder)
            print('Converting ', subj_id)
            old_sess_fold = subj_id.split("_")
            numerical_id = re.sub(r'^OAS3', '', old_sess_fold[0])
            bids_id = 'sub-OASIS3' + str(numerical_id)
            bids_subj_folder = path.join(dest_dir, bids_id)
            if not os.path.isdir(bids_subj_folder):
                os.mkdir(bids_subj_folder)

            #session_name = "".join(old_sess_fold[1:])
            session_name = old_sess_fold[2]
            session_folder = path.join(bids_subj_folder, 'ses-' + session_name)
            if not os.path.isdir(session_folder):
                os.mkdir(path.join(session_folder))
                os.mkdir(path.join(session_folder, 'anat'))

            # In order do convert the Analyze format to Nifti the path to the .img file is required
            img_file_path = glob(path.join(t1_folder, '**/*.nii.gz'), recursive=True)[0]
            #TODO: if not nifit ( len(img_file_path) = 0), then convert
            #img_file_path = glob(path.join(t1_folder, '**/*.img'), recursive=True)[0]

            base_file_name = os.path.basename(img_file_path).split(".")[0]
            json_file_path = glob(path.join(t1_folder, '**/' + base_file_name + ".json"), recursive=True)[0]

            output_path = path.join(session_folder, 'anat')
            copy_file(json_file_path, output_path)
            # First, convert to Nifti so that we can extract the s_form with NiBabel
            # (NiBabel creates an 'Spm2AnalyzeImage' object that does not contain 'get_sform' method
            # img_with_wrong_orientation_analyze = nb.load(img_file_path)

            # OASIS-1 images have the same header but sform is incorrect
            # To solve this issue, we use header from images converted with FreeSurfer
            # to generate a 'clean hard-coded' header
            # affine:
            # [[   0.    0.   -1.   80.]
            #  [   1.    0.    0. -128.]
            #  [   0.    1.    0. -128.]
            #  [   0.    0.    0.    1.]]
            # affine = np.array([0, 0, -1, 80, 1, 0, 0, -128, 0, 1, 0, -128, 0, 0, 0, 1]).reshape(4, 4)
            # s_form = affine.astype(np.int16)
            #
            # hdr = nb.Nifti1Header()
            # hdr.set_data_shape((256, 256, 160))
            # hdr.set_data_dtype(np.int16)
            # hdr['bitpix'] = 16
            # hdr.set_sform(s_form, code='scanner')
            # hdr.set_qform(s_form, code='scanner')
            # hdr['extents'] = 16384
            # hdr['xyzt_units'] = 10
            #
            # img_with_good_orientation_nifti = nb.Nifti1Image(
            #     np.round(img_with_wrong_orientation_analyze.get_data()).astype(np.int16),
            #     s_form,
            #     header=hdr
            # )
            # nb.save(img_with_good_orientation_nifti, output_path)

        if not os.path.isdir(dest_dir):
            os.mkdir(dest_dir)

        subjs_folders = glob(path.join(source_dir, 'OAS3*'))
        # subjs_folders = [subj_folder for subj_folder in subjs_folders if subj_folder.endswith('_MR1')]
        poolrunner = Pool(max(os.cpu_count()-1, 1))
        poolrunner.map(convert_single_subject, subjs_folders)
