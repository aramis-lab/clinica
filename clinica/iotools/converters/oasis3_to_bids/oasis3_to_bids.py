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
        participants_df = bids.create_participants_df('OASIS3', clinic_specs_path, clinical_data_dir, bids_ids)

        # Replace the values of the diagnosis_bl column
        participants_df['diagnosis_bl'].replace([0.0, np.nan], 'CN', inplace=True)
        participants_df['diagnosis_bl'].replace([0.5, 1.0, 1.5, 2.0], 'AD', inplace=True)
        # Following line has no sense
        # participants_df['diagnosis_bl'].replace(participants_df['diagnosis_bl']>0.0, 'AD', inplace=True)
        participants_df.to_csv(path.join(bids_dir, 'participants.tsv'), sep='\t', index=False, encoding='utf-8')

        # --Create sessions files--
        sessions_dict = bids.create_sessions_dict(clinical_data_dir, 'OASIS3', clinic_specs_path, bids_ids, 'ID')
        for y in bids_ids:
            for z in sessions_dict[y].keys():
                if sessions_dict[y][z]['diagnosis'] > 0:
                    sessions_dict[y][z]['diagnosis'] = 'AD'
                else:
                    sessions_dict[y][z]['diagnosis'] = 'CN'

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

        def convert_single_subject(subj_folder):
            import os
            import re
            import shutil

            def copy_file(source, target):
                try:
                    new_filename = re.sub("OAS3", "OASIS3", path.basename(source))
                    new_filename = re.sub("_echo-", "_run-", new_filename)
                    if re.search(".*_pet.(nii.gz|json)", new_filename) is not None:
                        new_filename = re.sub("task-rest_", "", new_filename)
                        new_filename = re.sub("acq", "trc", new_filename)
                        target = path.join(path.dirname(target), "pet")
                    elif re.search(".*_(minIP|swi|dwi).(nii.gz|json)", new_filename) is not None:
                        target = path.join(path.dirname(target), "swi")
                    if not os.path.isdir(target):
                        os.mkdir(path.join(target))
                    shutil.copy(source, target)
                    os.rename(path.join(target, path.basename(source)), path.join(target, new_filename))
                except IOError as e:
                    print("Unable to copy file. %s" % e)

            subj_id = os.path.basename(subj_folder)
            print('Converting ', subj_id)
            old_sess_fold = subj_id.split("_")
            numerical_id = re.sub(r'^OAS3', '', old_sess_fold[0])
            bids_id = 'sub-OASIS3' + str(numerical_id)
            bids_subj_folder = path.join(dest_dir, bids_id)
            if not os.path.isdir(bids_subj_folder):
                os.mkdir(bids_subj_folder)

            session_name = old_sess_fold[2]
            session_folder = path.join(bids_subj_folder, 'ses-' + session_name)
            if not os.path.isdir(session_folder):
                os.mkdir(path.join(session_folder))

            # In order do convert the Analyze format to Nifti the path to the .img file is required
            nifti_folders = glob(path.join(subj_folder, '**/*.nii.gz'), recursive=True)
            #TODO: if not nifit ( len(img_file_path) = 0), then convert

            for nifti_folder in nifti_folders:
                #  Get the json files
                base_file_name = os.path.basename(nifti_folder).split(".")[0]
                json_file_paths = glob(path.join(subj_folder, '**/' + base_file_name + ".json"), recursive=True)

                # Looking for an existing func or anat folder
                dir_name = path.dirname(nifti_folder)
                regex = re.compile("(?P<folder>anat|func)[1-9]")
                if regex.search(dir_name) is not None:
                    output_folder = regex.search(dir_name).group("folder")
                elif regex.search(path.abspath(dir_name)) is not None:
                    output_folder = regex.search(path.abspath(dir_name)).group("folder")
                else:
                    output_folder = "anat"

                # Moove nifti and json files to the right folder
                output_path = path.join(session_folder, output_folder)
                copy_file(nifti_folder, output_path)
                if len(json_file_paths) > 0:
                    copy_file(json_file_paths[0], output_path)

        if not os.path.isdir(dest_dir):
            os.mkdir(dest_dir)

        subjs_folders = glob(path.join(source_dir, 'OAS3*'))
        for subj_folder in subjs_folders:
            convert_single_subject(subj_folder)
