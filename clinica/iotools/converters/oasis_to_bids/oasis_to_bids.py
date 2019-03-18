# coding: utf8

"""
Convert OASIS dataset (http://www.oasis-brains.org/) to BIDS.
"""

from clinica.iotools.abstract_converter import Converter

__author__ = "Sabrina Fontanella"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Completed"


class OasisToBids(Converter):
    def convert_clinical_data(self, clinical_data_dir, bids_dir):
        """
        Convert the clinical data defined inside the clinical_specifications.xlx into BIDS

        Args:
            clinical_data_dir: path to the folder with the original clinical data
            bids_dir: path to the bids directory

        """
        import clinica.iotools.bids_utils as bids
        from os import path
        import os
        import numpy as np
        print('Converting clinical data...')
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

        """
        from os import path
        from glob import glob
        import os
        from multiprocessing.dummy import Pool
        from multiprocessing import cpu_count

        def convert_single_subject(subj_folder):
            import os
            import subprocess

            t1_folder = path.join(subj_folder, 'PROCESSED', 'MPRAGE',
                                  'SUBJ_111')
            subj_id = os.path.basename(subj_folder)
            print('Converting ', subj_id)
            numerical_id = (subj_id.split("_"))[1]
            bids_id = 'sub-OASIS1' + str(numerical_id)
            bids_subj_folder = path.join(dest_dir, bids_id)
            if not os.path.isdir(bids_subj_folder):
                os.mkdir(bids_subj_folder)

            session_folder = path.join(bids_subj_folder, 'ses-M00')
            if not os.path.isdir(session_folder):
                os.mkdir(path.join(session_folder))
                os.mkdir(path.join(session_folder, 'anat'))

            # In order do convert the Analyze format to NIFTI the path to the .img file is required
            img_file_path = glob(path.join(t1_folder, '*.img'))[0]
            output_path = path.join(session_folder, 'anat',
                                    bids_id + '_ses-M00_T1w.nii.gz')
            subprocess.run('mri_convert' + ' ' + img_file_path + ' ' + output_path,
                           shell=True,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)

        if not os.path.isdir(dest_dir):
            os.mkdir(dest_dir)

        subjs_folders = glob(path.join(source_dir, 'OAS1_*'))
        poolrunner = Pool(cpu_count() - 1)
        poolrunner.map(convert_single_subject, subjs_folders)
