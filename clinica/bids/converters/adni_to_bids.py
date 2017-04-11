
from clinica.bids.abstract_converter import Converter
from clinica.engine.cmdparser import CmdParser
import logging

__author__ = "Jorge Samper and Sabrina Fontanella"
__copyright__ = "Copyright 2017, The Aramis Lab Team"
__credits__ = [""]
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Sabrina Fontanella"
__email__ = "sabrina.fontanella@icm-institute.org"
__status__ = "Development"


class ADNI_TO_BIDS(Converter, CmdParser):

    def define_name(self):
        self._name = 'adni-to-bids'

    def define_options(self):
        self._args.add_argument("dataset_directory",
                               help='Path of the unorganized ADNI directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("--clinical-only", type=bool, default=False,
                                dest='clinical_only',
                                help='(Optional) Given an already existing BIDS output folder, convert only the clinical data.')

    def convert_clinical_data(self, input_path, out_path):
        """

        Convert clinical data of ADNI dataset.

        :param src:
        :param out_path:
        :return:

        {'sub-ADNI009S1354': { 'bl': {u'diagnosis': 'Dementia', 'session_id': 'ses-bl'},
                              'm06': {u'diagnosis': 'Dementia', 'session_id': 'ses-m06'}},...
         ....
         }
        """
        from os import path
        import os
        import logging
        import clinica.bids.bids_utils as bids
        import pandas as pd
        from glob import glob
        from os.path import normpath

        logging.basicConfig(filename=path.join(out_path, 'conversion_clinical.log'),
                            format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG,
                            datefmt='%m/%d/%Y %I:%M')

        clinic_specs_path = path.join(os.path.dirname(os.path.dirname(__file__)), 'data',
                                      'clinical_specifications_adni.xlsx')


        try:
            os.path.exists(out_path)
        except IOError:
            print 'BIDS folder not found.'
            raise

        bids_ids = bids.get_bids_subjs_list(out_path)
        bids_subjs_paths = bids.get_bids_subjs_paths(out_path)

        # -- Creation of participant.tsv --
        participants_df = bids.create_participants_df(input_path, out_path, 'ADNI', 'clinical_specifications_adni.xlsx', bids_ids)
        participants_df.to_csv(path.join(out_path, 'participant.tsv'), sep='\t', index=False)

        # -- Creation of sessions.tsv --
        logging.info("--Creation of sessions files. --")
        print("\nCreation of sessions files...")
        # Load data
        sessions = pd.read_excel(clinic_specs_path, sheetname='sessions.tsv')
        sessions_fields = sessions['ADNI']
        field_location = sessions['ADNI location']
        sessions_fields_bids = sessions['BIDS CLINICA']
        file_to_read = pd.read_csv(path.join(input_path, 'clinicalData', 'ADNIMERGE.csv'))
        fields_dataset = []
        fields_bids = []
        sessions_dict = {}
        subj_to_remove = []

        for i in range(0, len(sessions_fields)):
            if not pd.isnull(sessions_fields[i]):
                fields_bids.append(sessions_fields_bids[i])
                fields_dataset.append(sessions_fields[i])

        sessions_df = pd.DataFrame(columns=fields_bids)

        # For each line of ADNIMERGE check if the subjects is available in the BIDS dataset and extract the information
        for r in range(0, len(file_to_read.values)):
            row = file_to_read.iloc[r]
            id_ref = 'PTID'
            subj_id = row[id_ref.decode('utf-8')]
            # Removes all the - from
            subj_id_alpha = bids.remove_space_and_symbols(subj_id)
            subj_bids = [s for s in bids_ids if subj_id_alpha in s]
            if len(subj_bids) == 0:
                pass
            else:
                subj_bids = subj_bids[0]
                for i in range(0, len(sessions_fields)):
                    # If the i-th field is available
                    if not pd.isnull(sessions_fields[i]):
                        sessions_df[sessions_fields_bids[i]] = row[sessions_fields[i]]
                        visit_id = row['VISCODE']
                        # If the dictionary already contain the subject add or update information regarding a specific session,
                        # otherwise create the enty
                        if sessions_dict.has_key(subj_bids):
                            sess_available = sessions_dict[subj_bids].keys()
                            if visit_id in sess_available:
                                sessions_dict[subj_bids][visit_id].update({sessions_fields_bids[i]: row[sessions_fields[i]]})
                            else:
                                sessions_dict[subj_bids].update({visit_id: {'session_id': 'ses-'+visit_id,
                                                                            sessions_fields_bids[i]: row[sessions_fields[i]]}})
                        else:
                            sessions_dict.update({subj_bids: {visit_id: {'session_id': 'ses-'+visit_id,
                                                                         sessions_fields_bids[i]: row[sessions_fields[i]]}}})

        # Write the information contained inside the dictionary to the proper session file
        for sp in bids_subjs_paths:
            bids_id = sp.split(os.sep)[-1]
            sessions_df = pd.DataFrame(columns=fields_bids)
            if sessions_dict.has_key(bids_id):
                sess_aval = sessions_dict[bids_id].keys()
                for ses in sess_aval:
                    sessions_df = sessions_df.append(pd.DataFrame(sessions_dict[bids_id][ses], index=['i', ]))

                sessions_df.to_csv(path.join(sp, bids_id + '_sessions.tsv'), sep='\t', index=False)

        # -- Creation of scans files --
        print 'Creation of scans files...'
        scans_dict = {}

        for bids_id in bids_ids:
            scans_dict.update({bids_id: {'T1/DWI/fMRI': {}, 'FDG': {}}})

        scans_specs = pd.read_excel(clinic_specs_path, sheetname='scans.tsv')
        scans_fields_db = scans_specs['ADNI']
        scans_fields_bids = scans_specs['BIDS CLINICA']
        scans_fields_mod = scans_specs['Modalities related']
        fields_bids = ['filename']

        for i in range(0, len(scans_fields_db)):
            if not pd.isnull(scans_fields_db[i]):
                fields_bids.append(scans_fields_bids[i])

        scans_df = pd.DataFrame(columns=(fields_bids))

        for bids_subj_path in bids_subjs_paths:
            # Create the file
            bids_id = os.path.basename(normpath(bids_subj_path))

            sessions_paths = glob(path.join(bids_subj_path, 'ses-*'))
            for session_path in sessions_paths:
                session_name = session_path.split(os.sep)[-1]
                tsv_name = bids_id + '_' + session_name + "_scans.tsv"

                # If the file already exists, remove it
                if os.path.exists(path.join(session_path, tsv_name)):
                    os.remove(path.join(session_path, tsv_name))

                scans_tsv = open(path.join(session_path, tsv_name), 'a')
                scans_df.to_csv(scans_tsv, sep='\t', index=False)

                # Extract modalities available for each subject
                mod_available = glob(path.join(session_path, '*'))
                for mod in mod_available:
                    mod_name = os.path.basename(mod)
                    files = glob(path.join(mod, '*'))
                    for file in files:
                        file_name = os.path.basename(file)
                        if mod == "anat" or mod == "dwi" or mod == "func":
                            type_mod = 'T1/DWI/fMRI'
                        else:
                            type_mod = 'FDG'

                        scans_df['filename'] = pd.Series(path.join(mod_name, file_name))
                        scans_df.to_csv(scans_tsv, header=False, sep='\t', index=False)

                scans_df = pd.DataFrame(columns=(fields_bids))

        print '-- Scans files created for each subject. --'

    def replace_sequence_chars(self, sequence_name):
        sequence = sequence_name.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_').replace(
                        '(',
                        '_').replace(
                        ')', '_').replace(':', '_')

        return sequence

    def convert_from_dicom(self, input_path, output_path, bids_name, mod_type):
        """

        :param t1_path:
        :param output_path:
        :param bids_name:
        :return:
        """

        import os
        import clinica.bids.bids_utils as bids
        from os import path
        from glob import glob

        if not os.path.exists(output_path):
            os.mkdir(output_path)
        os.system('dcm2niix -b n -z y -o ' + output_path + ' -f ' + bids_name + bids.get_bids_suff(mod_type) + ' ' + input_path)

        # If dcm2niix didn't work use dcm2nii
        if not os.path.exists(path.join(output_path, bids_name + bids.get_bids_suff(mod_type) + '.nii.gz')):
            print 'Conversion with dcm2niix failed, trying with dcm2nii'
            #os.system('dcm2nii -a n -d n -e n -i y -g n -p n -m n -r n -x n -o ' + output_path + ' ' + image_path')

        # If the conversion failed with both tools
        if not os.path.exists(path.join(output_path, bids_name + bids.get_bids_suff(mod_type) + '.nii.gz')):
            print 'Conversion of the dicom failed for ', input_path

    def fill_zeros(self, s, length):
        return ('0' * (length - len(str(s)))) + str(s)

    def days_between(self, d1, d2):
        from datetime import datetime
        d1 = datetime.strptime(d1, "%Y-%m-%d")
        d2 = datetime.strptime(d2, "%Y-%m-%d")
        return abs((d2 - d1).days)

    def convert_images(self, source_dir, clinical_dir, dest_dir, subjs_list_path='', mod_to_add='', mod_to_update='', add_subjs=False):
        """
        The function first computes the paths of the right image to be converted and
        :param source_dir:
        :param dest_dir:
        :return:
        """

        import clinica.bids.bids_utils as bids
        import os
        from os import path
        import pandas as pd
        import shutil
        from glob import glob
        import numpy as np

        print "*******************************"
        print "ADNI to BIDS converter"
        print "*******************************"

        subjs_list = []
        adni_merge_path = path.join(clinical_dir, 'ADNIMERGE.csv')
        adni_merge = pd.read_csv(adni_merge_path)

        # Options to use
        convert_func = (mod_to_add == '' and mod_to_update == '') or mod_to_add == 'func' or mod_to_update == 'func'
        print convert_func


        # Load a file with subjects list or compute all the subjects
        if subjs_list_path != '':
            subjs_list = [line.rstrip('\n') for line in open(subjs_list_path)]
        else:
            print 'Using all the subjects contained into ADNI merge...'
            subjs_list = adni_merge['PTID'].unique()

        print 'Subjects found:', len(subjs_list)

        # Extract the list of sessions from VISCODE column of adnimerge file, remove duplicate and convert to a list
        sess_list = adni_merge['VISCODE'].drop_duplicates().tolist()

        bids_ids = []
        alpha_ids = []

        if (mod_to_add == '' and mod_to_update == '') or (mod_to_add != '' and os.path.exists(dest_dir) == False):
            print 'Creating the output folder'
            os.mkdir(dest_dir)
            os.mkdir(path.join(dest_dir, 'conversion_info'))
        #
        # # Compute anat paths
        # if (mod_to_add == '' or mod_to_add == 'anat') and (mod_to_update == '' or mod_to_update == 'anat'):
        #     t1_paths = self.compute_t1_paths(source_dir, clinical_dir, dest_dir, subjs_list)
        #
        # # Compute pet paths
        # if (mod_to_add == '' or mod_to_add == 'pet') and (mod_to_update == '' or mod_to_update == 'pet'):
        #     print 'Calculating paths for PET FDG...'
        #     pet_fdg_paths = self.compute_fdg_pet_paths(source_dir, clinical_dir, dest_dir, subjs_list)
        #     print 'Done!'
        #     pet_av45_paths = self.compute_av45_pet_paths(source_dir, clinical_dir, dest_dir, subjs_list)

        # Compute func paths
        if (mod_to_add == '' or mod_to_add == 'func') and (mod_to_update == '' or mod_to_update == 'func'):
            print 'Calculating paths for FMRI...'
            print subjs_list
            fmri_paths = self.compute_fmri_path(source_dir, clinical_dir, dest_dir, subjs_list)
            #fmri_paths = pd.read_csv(path.join(dest_dir,'paths','fmri_paths.tsv'), sep='\t')
            print 'Done!'


        # Create subjects folders
        for subj in subjs_list:
            alpha_id = bids.remove_space_and_symbols(subj)
            bids_id = 'sub-ADNI' + alpha_id
            alpha_ids.append(alpha_id)
            bids_ids.append(bids_id)
            if (mod_to_add == '' and mod_to_update == '') or (
                    mod_to_add != '' and not os.path.exists(path.join(dest_dir, bids_id))):
                os.mkdir(path.join(dest_dir, bids_id))


        # For each subject, extract the info from the paths files and convert the modalities
        for i in range(0, len(subjs_list)):
            print 'Converting ', subjs_list[i]

            if convert_func:
                # Extract the list of sessions available
                sess_list = mri_info = fmri_paths[(fmri_paths['Subject_ID'] == subjs_list[i])]['VISCODE'].values

                # For each session available, create the folder is doesn't exist and convert them
                for ses in sess_list:
                    bids_ses_id = 'ses-' + ses
                    bids_file_name = bids_ids[i] + '_ses-' + ses
                    ses_path = path.join(dest_dir, bids_ids[i], bids_ses_id)

                    if not os.path.exists(ses_path):
                        os.mkdir(ses_path)

                    fmri_info = fmri_paths[(fmri_paths['Subject_ID'] == subjs_list[i]) & (fmri_paths['VISCODE'] == ses)]
                    if not fmri_info['Path'].empty:
                        if type(fmri_info['Path'].values[0]) != float:
                            if not os.path.exists(path.join(ses_path, 'func')):
                                os.mkdir(path.join(ses_path, 'func'))
                                fmri_path = fmri_info['Path'].values[0]
                                bids.convert_fmri(fmri_path, path.join(ses_path, 'func'), bids_file_name)


        #
        # for i in range(0, len(subjs_list)):
        #     print 'Converting ', subjs_list[i]
        #     for ses in sess_list:
        #         bids_file_name = bids_ids[i] + '_ses-' + ses
        #         bids_ses_id = 'ses-' + ses
        #         ses_path = path.join(dest_dir, bids_ids[i], bids_ses_id)
        #         # if (mod_to_add == '' and mod_to_update == '') or (mod_to_add != '' and not os.path.exists(ses_path)):
        #         #     print 'Creating ses folder'
        #         #     os.mkdir(ses_path)
        #         if not os.path.exists(ses_path):
        #             os.mkdir(ses_path)
        #
        #         if mod_to_add != '' and mod_to_add != 'anat':
        #             pass
        #         elif (mod_to_add == '' and mod_to_update == '') or mod_to_add == 'anat' or mod_to_update == 'anat':
        #             # Convert T1
        #             t1_info = t1_paths[(t1_paths['Subject_ID'] == subjs_list[i]) & (t1_paths['VISCODE'] == ses)]
        #             if len(t1_info) == 0:
        #                 print 'No T1 found for the visit ' + ses
        #             else:
        #                 t1_path = (t1_info['Path'].values[0]).replace(' ', '\ ')
        #                 # Check if the pet folder already exist
        #                 if os.path.isdir(path.join(ses_path, 'anat')):
        #                     if mod_to_add != '':
        #                         raise IOError('anat modality found. For updating the dataset use the flag -updated_mod')
        #
        #                     print 'Removing the old anat folder...'
        #                     shutil.rmtree(path.join(ses_path, 'anat'))
        #                     print 'Removed!'
        #
        #                 os.mkdir(path.join(ses_path, 'anat'))
        #                 bids.convert_T1(t1_path, path.join(ses_path, 'anat'), bids_file_name)
        #                 t1_bids_path = glob(path.join(path.join(ses_path, 'anat'),
        #                                               bids_file_name + bids.get_bids_suff('T1') + '.nii.gz'))[0]
        #                 # Correct the t1 image and overwrite the previous one
        #                 self.center_nifti_origin(t1_bids_path, t1_bids_path)
        #         # Convert pet (fdg and av45)
        #         if mod_to_add != '' and mod_to_add != 'pet':
        #             pass
        #         elif (mod_to_add == '' and mod_to_update == '') or mod_to_add == 'pet' or mod_to_update == 'pet':
        #             pet_fdg_info = pet_fdg_paths
        #             # Check if the pet folder already exist
        #             if os.path.isdir(path.join(ses_path, 'pet')):
        #                 if mod_to_add != '':
        #                     raise IOError('PET modality found. For updating the dataset use the flag -updated_mod')
        #
        #                 print 'Removing the old PET folder...'
        #                 shutil.rmtree(path.join(ses_path, 'pet'))
        #
        #             pet_fdg_info = pet_fdg_paths[
        #                 (pet_fdg_paths['Subject_ID'] == subjs_list[i]) & (pet_fdg_paths['VISCODE'] == ses)]
        #             print pet_fdg_info
        #             if len(pet_fdg_info['Path']) != 0:
        #                 original = pet_fdg_info['Original'].values
        #                 if len(original) > 0:
        #                     if pet_fdg_info['Original'].values[0]:
        #                         os.mkdir(path.join(ses_path, 'pet'))
        #                         pet_path = pet_fdg_info['Path'].values[0]
        #
        #                         self.center_nifti_origin(pet_path,
        #                                                  path.join(path.join(ses_path, 'pet'),
        #                                                            bids_file_name + '_task-rest_acq-fdg_pet.nii.gz'))
        #                     else:
        #                         pet_path = pet_fdg_info['Path'].values[0].replace(' ', '\ ')
        #                         self.convert_from_dicom(pet_path, path.join(ses_path, 'pet'),
        #                                                 bids_file_name + '_task-rest_acq-fdg', 'pet')
        #                 else:
        #                     print 'Original not found for ', subjs_list[i]
        #             else:
        #                 print 'No pet_fdg path found for subject: ' + subj + ' visit: ' + ses
        #
        #                 # Convert pet_av45
        #             pet_av45_info = pet_av45_paths
        #
        #             # Check if the pet folder already exist
        #             if os.path.isdir(path.join(ses_path, 'pet')):
        #                 if mod_to_add != '':
        #                     raise IOError('PET modality found. For updating the dataset use the flag -updated_mod')
        #
        #             pet_av45_info = pet_av45_paths[
        #                 (pet_av45_paths['Subject_ID'] == subjs_list[i]) & (pet_av45_paths['VISCODE'] == ses)]
        #             original = pet_av45_info['Original'].values
        #             if len(pet_av45_info['Path']) != 0:
        #                 if len(original) > 0:
        #                     if pet_av45_info['Original'].values[0] == True:
        #                         if not os.path.exists(path.join(ses_path, 'pet')):
        #                             os.mkdir(path.join(ses_path, 'pet'))
        #                         pet_path = pet_av45_info['Path'].values[0]
        #
        #                         self.center_nifti_origin(pet_path,
        #                                                  path.join(path.join(ses_path, 'pet'),
        #                                                            bids_file_name + '_task-rest_acq-av45_pet.nii.gz'))
        #                     else:
        #                         pet_path = pet_av45_info['Path'].values[0].replace(' ', '\ ')
        #                         self.convert_from_dicom(pet_path, path.join(ses_path, 'pet'),
        #                                                 bids_file_name + '_task-rest_acq-av45', 'pet')
        #                 else:
        #                     print 'Original not found for ', subjs_list[i]
        #             else:
        #                 print 'No pet av45 path found for subject: ' + subj + ' visit: ' + ses
        #
        #         # Convert func
        #         if mod_to_add != '' and mod_to_add != 'func':
        #             pass
        #         elif (mod_to_add == '' and mod_to_update == '') or mod_to_add == 'func' or mod_to_update == 'func':
        #             if os.path.isdir(path.join(ses_path, 'func')):
        #                 if mod_to_add != '':
        #                     raise IOError('Func modality found. For updating the dataset use the flag -updated_mod')
        #
        #                 print 'Removing the old func folder...'
        #                 shutil.rmtree(path.join(ses_path, 'func'))
        #             fmri_info = fmri_paths[(fmri_paths['Subject_ID'] == subjs_list[i]) & (fmri_paths['VISCODE'] == ses)]
        #             if not fmri_info['Path'].empty:
        #                 if type(fmri_info['Path'].values[0]) != float:
        #                     if not os.path.exists(path.join(ses_path, 'func')):
        #                         os.mkdir(path.join(ses_path, 'func'))
        #                     fmri_path = fmri_info['Path'].values[0]
        #                     bids.convert_fmri(fmri_path, path.join(ses_path, 'func'), bids_file_name)

    def center_nifti_origin(self, input_image, output_image):
        import nibabel as nib
        import numpy as np

        img = nib.load(input_image)
        canonical_img = nib.as_closest_canonical(img)
        hd = canonical_img.header
        if hd['quatern_b'] != 0 or hd['quatern_c'] != 0 or hd['quatern_d'] != 0:
            print 'Warning: Not all values in quatern are equal to zero'
        qform = np.zeros((4, 4))
        for i in range(1, 4):
            qform[i - 1, i - 1] = hd['pixdim'][i]
            qform[i - 1, 3] = -1.0 * hd['pixdim'][i] * hd['dim'][i] / 2.0
        new_img = nib.Nifti1Image(canonical_img.get_data(caching='unchanged'), qform)
        nib.save(new_img, output_image)

    def adni1_image(self, subject_id, timepoint, visit_str, mprage_meta_subj, ida_meta_subj):

        # Get the preferred scan (image series that has been Scaled)
        filtered_mprage = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Processed')
                                           & (mprage_meta_subj.Visit == visit_str)
                                           & (mprage_meta_subj.Sequence.map(lambda x: x.endswith('Scaled')))]

        # If there is not a preferred image we use ADNI2 processing (get the original) preferring 1.5T images
        if filtered_mprage.shape[0] < 1:
            mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']
            return self.adni2_image(subject_id, timepoint, visit_str, mprage_meta_subj_orig, preferred_field_strength=1.5)

        filtered_mprage_mag = filtered_mprage
        if len(filtered_mprage.MagStrength.unique()) > 1:
            filtered_mprage_mag = filtered_mprage[filtered_mprage.MagStrength == 1.5]  # Select 1.5T images

        scan = filtered_mprage_mag.iloc[0]
        series_id = scan.SeriesID

        filtered_scan = ida_meta_subj[ida_meta_subj.LONIUID == series_id]

        if filtered_scan.shape[0] < 1:
            # If no IDA_META for 1.5T try for 3T
            filtered_mprage_mag = filtered_mprage[filtered_mprage.MagStrength == 3.0]
            scan = filtered_mprage_mag.iloc[0]
            series_id = scan.SeriesID

            filtered_scan = ida_meta_subj[ida_meta_subj.LONIUID == series_id]

            if filtered_scan.shape[0] < 1:
                print 'NO IDA Meta: ' + subject_id + ' for visit ' + timepoint + ' - ' + visit_str
                return None

        original = True
        ida_scan = filtered_scan.iloc[0]
        if ida_scan.Scanner.find('Philips') > -1:

            scan = (mprage_meta_subj[
                        (mprage_meta_subj['Orig/Proc'] == 'Original') & (mprage_meta_subj.SeriesID == series_id)]).iloc[
                0]
            sequence = scan.Sequence

        else:  # scan already selected above
            sequence = scan.Sequence[:scan.Sequence.find('N3') - 2]
            original = False

        sequence = sequence.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_').replace('(',
                                                                                                            '_').replace(
            ')', '_').replace(':', '_')

        return {'Subject_ID': subject_id,
                'VISCODE': timepoint,
                'Visit': visit_str,
                'Sequence': sequence,
                'Scan_Date': scan.ScanDate,
                'Study_ID': str(scan.StudyID),
                'Series_ID': str(scan.SeriesID),
                'Field_Strength': scan.MagStrength,
                'Original': original}

    def adni2_image(self, subject_id, timepoint, visit_str, mprage_meta_subj_orig, preferred_field_strength=3.0):

        cond_mprage = ((mprage_meta_subj_orig.Visit == visit_str) & mprage_meta_subj_orig.Sequence.map(lambda x: ((
                                                                                                                  x.lower().find(
                                                                                                                      'mprage') > -1) | (
                                                                                                                  x.lower().find(
                                                                                                                      'mp-rage') > -1) | (
                                                                                                                  x.lower().find(
                                                                                                                      'mp rage') > -1)) & (
                                                                                                                 x.find(
                                                                                                                     '2') < 0) & (
                                                                                                                 x.lower().find(
                                                                                                                     'repeat') < 0)))

        cond_spgr = ((mprage_meta_subj_orig.Visit == visit_str) & mprage_meta_subj_orig.Sequence.map(
            lambda x: (x.lower().find('spgr') > -1) & (x.lower().find('acc') < 0) & (x.lower().find('repeat') < 0)))

        filtered_scan = mprage_meta_subj_orig[cond_mprage | cond_spgr]
        if filtered_scan.shape[0] < 1:

            # TODO Improve this code. Don't make a double verification for the whole condition.
            # Invert order of filtering: less to more restrictive, check for the repeated as for the MagStrength
            cond_mprage_rep = ((mprage_meta_subj_orig.Visit == visit_str) & mprage_meta_subj_orig.Sequence.map(
                lambda x: ((x.lower().find('mprage') > -1) | (x.lower().find('mp-rage') > -1) | (
                x.lower().find('mp rage') > -1)) & (x.find('2') < 0)))

            cond_spgr_rep = ((mprage_meta_subj_orig.Visit == visit_str) & mprage_meta_subj_orig.Sequence.map(
                lambda x: (x.lower().find('spgr') > -1) & (x.lower().find('acc') < 0)))

            filtered_scan = mprage_meta_subj_orig[cond_mprage_rep | cond_spgr_rep]
            if filtered_scan.shape[0] < 1:
                print 'NO MPRAGE Meta2: ' + subject_id + ' for visit ' + timepoint + ' - ' + visit_str
                return None

        if len(filtered_scan.MagStrength.unique()) > 1:
            filtered_scan = filtered_scan[
                filtered_scan.MagStrength == preferred_field_strength]  # Select preferred_field_strength images

        scan = filtered_scan.iloc[0]

        sequence = scan.Sequence.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_').replace('(',
                                                                                                                 '_').replace(
            ')', '_').replace(':', '_')

        return {'Subject_ID': subject_id,
                'VISCODE': timepoint,
                'Visit': visit_str,
                'Sequence': sequence,
                'Scan_Date': scan.ScanDate,
                'Study_ID': str(scan.StudyID),
                'Series_ID': str(scan.SeriesID),
                'Field_Strength': scan.MagStrength,
                'Original': True}

    def adnigo_image(self,subject_id, timepoint, visit_str, mprage_meta_subj, ida_meta_subj, original_phase):

        if original_phase == 'ADNI1':
            filtered_mprage = mprage_meta_subj[(mprage_meta_subj['Orig/Proc'] == 'Processed')
                                               & (mprage_meta_subj.MagStrength == 1.5)
                                               & (mprage_meta_subj.Visit == visit_str)
                                               & (mprage_meta_subj.Sequence.map(lambda x: x.endswith('Scaled')))]
            if filtered_mprage.shape[0] > 0:
                return self.adni1_image(subject_id, timepoint, visit_str, mprage_meta_subj, ida_meta_subj)

        mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']
        return self.adni2_image(subject_id, timepoint, visit_str, mprage_meta_subj_orig)

    def select_image_qc(self, id_list, mri_qc_subj):
        import numpy as np

        if len(id_list) == 0:
            return None

        selected_image = None
        image_ids = ['I' + str(imageuid) for imageuid in id_list]
        int_ids = [int(imageuid) for imageuid in id_list]
        images_qc = mri_qc_subj[mri_qc_subj.loni_image.isin(image_ids)]

        if images_qc.shape[0] < 1:
            return min(int_ids)

        if np.sum(images_qc.series_selected) == 1:
            selected_image = images_qc[images_qc.series_selected == 1].iloc[0].loni_image[1:]
        else:
            images_not_rejected = images_qc[images_qc.series_quality < 4]

            if images_not_rejected.shape[0] < 1:

                # There are no images that passed the qc
                # so we'll try to see if there are other images without qc,
                # otherwise return None
                qc_ids = set([int(qc_id[1:]) for qc_id in images_qc.loni_image.unique()])
                no_qc_ids = list(set(int_ids) - qc_ids)

                if len(no_qc_ids) == 0:
                    return None
                else:
                    return min(no_qc_ids)

            series_quality = [q if q > 0 else 4 for q in list(images_not_rejected.series_quality)]
            best_q = np.amin(series_quality)

            images_best_qc = images_not_rejected[images_not_rejected.series_quality == best_q]
            if images_best_qc.shape[0] == 1:
                selected_image = images_best_qc.iloc[0].loni_image[1:]
            else:
                selected_image = min(int_ids)

        return int(selected_image)

    def dti_image(self, subject_id, timepoint, visit_str, ida_meta_scans, mri_qc_subj, enhanced):

        sel_image = self.select_image_qc(list(ida_meta_scans.IMAGEUID), mri_qc_subj)
        if sel_image is None:
            return None

        sel_scan = ida_meta_scans[ida_meta_scans.IMAGEUID == sel_image].iloc[0]

        sequence = sel_scan.Sequence
        sequence = sequence.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_').replace('(', '_').replace(')', '_').replace(':', '_')

        image_dict = {'Subject_ID': subject_id,
                      'VISCODE': timepoint,
                      'Visit': visit_str,
                      'Sequence':  sequence,
                      'Scan_Date': sel_scan['Scan Date'],
                      'Study_ID': str(int(sel_scan.LONISID)),
                      'Series_ID': str(int(sel_scan.LONIUID)),
                      'Image_ID': str(int(sel_scan.IMAGEUID)),
                      'Field_Strength': sel_scan.MagStrength,
                      'Scanner': sel_scan.Scanner,
                      'Enhanced': enhanced}

        return image_dict

    def visits_to_timepoints_t1(self, subject, mprage_meta_subj_orig, adnimerge_subj):
        from datetime import datetime

        mprage_meta_subj_orig = mprage_meta_subj_orig[mprage_meta_subj_orig['Visit'] != 'ADNI Baseline']

        visits = dict()

        unique_visits = list(mprage_meta_subj_orig.Visit.unique())

        pending_timepoints = []

        # We try to obtain the corresponding image Visit for a given VISCODE
        for adni_row in adnimerge_subj.iterrows():  # (adnimerge_subj[adnimerge_subj.FLDSTRENG.map(lambda x: x is not '')]).iterrows():
            visit = adni_row[1]
            if visit.ORIGPROT == 'ADNI2':
                if visit.VISCODE == 'bl':
                    preferred_visit_name = 'ADNI2 Screening MRI-New Pt'
                elif visit.VISCODE == 'm03':
                    preferred_visit_name = 'ADNI2 Month 3 MRI-New Pt'
                elif visit.VISCODE == 'm06':
                    preferred_visit_name = 'ADNI2 Month 6-New Pt'
                else:
                    year = str(int(visit.VISCODE[1:]) / 12)
                    preferred_visit_name = 'ADNI2 Year ' + year + ' Visit'
            else:
                if visit.VISCODE == 'bl':
                    if visit.ORIGPROT == 'ADNI1':
                        preferred_visit_name = 'ADNI Screening'
                    else:  # ADNIGO
                        preferred_visit_name = 'ADNIGO Screening MRI'
                elif visit.VISCODE == 'm03':  # Only for ADNIGO Month 3
                    preferred_visit_name = 'ADNIGO Month 3 MRI'
                else:
                    month = int(visit.VISCODE[1:])
                    if month < 54:
                        preferred_visit_name = 'ADNI1/GO Month ' + str(month)
                    else:
                        preferred_visit_name = 'ADNIGO Month ' + str(month)

            if preferred_visit_name in unique_visits:
                key_preferred_visit = (visit.VISCODE, visit.COLPROT, visit.ORIGPROT)
                if key_preferred_visit not in visits.keys():
                    visits[key_preferred_visit] = preferred_visit_name
                elif visits[key_preferred_visit] != preferred_visit_name:
                    print 'Multiple visits for one timepoint!'
                    print subject
                    print key_preferred_visit
                    print visits[key_preferred_visit]
                    print visit
                unique_visits.remove(preferred_visit_name)
                continue

            pending_timepoints.append(visit)

        # Then for images.Visit non matching the expected labels we find the closest date in visits list
        for visit in unique_visits:
            image = (mprage_meta_subj_orig[mprage_meta_subj_orig.Visit == visit]).iloc[0]
            min_db = 100000
            min_db2 = 0
            min_visit = None
            min_visit2 = None

            for timepoint in pending_timepoints:
                db = self.days_between(image.ScanDate, timepoint.EXAMDATE)
                if db < min_db:
                    min_db2 = min_db
                    min_visit2 = min_visit

                    min_db = db
                    min_visit = timepoint

            if min_visit is None:
                print 'No corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.Visit
                print image
                continue

            if min_visit2 is not None and min_db > 90:
                print 'More than 60 days for corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.Visit + ' on ' + image.ScanDate
                print 'Timepoint 1: ' + min_visit.VISCODE + ' - ' + min_visit.ORIGPROT + ' on ' + min_visit.EXAMDATE + ' (Distance: ' + str(
                    min_db) + ' days)'
                print 'Timepoint 2: ' + min_visit2.VISCODE + ' - ' + min_visit2.ORIGPROT + ' on ' + min_visit2.EXAMDATE + ' (Distance: ' + str(
                    min_db2) + ' days)'

                # If image is too close to the date between two visits we prefer the earlier visit
                if (datetime.strptime(min_visit.EXAMDATE, "%Y-%m-%d")
                        > datetime.strptime(image.ScanDate, "%Y-%m-%d")
                        > datetime.strptime(min_visit2.EXAMDATE, "%Y-%m-%d")):
                    dif = self.days_between(min_visit.EXAMDATE, min_visit2.EXAMDATE)
                    if abs((dif / 2.0) - min_db) < 30:
                        min_visit = min_visit2

                print 'We prefer ' + min_visit.VISCODE

            key_min_visit = (min_visit.VISCODE, min_visit.COLPROT, min_visit.ORIGPROT)
            if key_min_visit not in visits.keys():
                visits[key_min_visit] = image.Visit
            elif visits[key_min_visit] != image.Visit:
                print 'Multiple visits for one timepoint!'
                print subject
                print key_min_visit
                print visits[key_min_visit]
                print image.Visit

        return visits

    def visits_to_timepoints_dti(self, subject, ida_meta_subj, adnimerge_subj):
        from datetime import datetime

        visits = dict()
        unique_visits = list(ida_meta_subj.Visit.unique())
        pending_timepoints = []

        # We try to obtain the corresponding image Visit for a given VISCODE
        for adni_row in adnimerge_subj.iterrows():
            visit = adni_row[1]
            if visit.ORIGPROT == 'ADNI2':
                if visit.VISCODE == 'bl':
                    preferred_visit_name = 'ADNI2 Screening MRI-New Pt'
                elif visit.VISCODE == 'm03':
                    preferred_visit_name = 'ADNI2 Month 3 MRI-New Pt'
                elif visit.VISCODE == 'm06':
                    preferred_visit_name = 'ADNI2 Month 6-New Pt'
                else:
                    year = str(int(visit.VISCODE[1:]) / 12)
                    preferred_visit_name = 'ADNI2 Year ' + year + ' Visit'
            else:
                if visit.VISCODE == 'bl':
                    if visit.ORIGPROT == 'ADNI1':
                        preferred_visit_name = 'ADNI Screening'
                    else:  # ADNIGO
                        preferred_visit_name = 'ADNIGO Screening MRI'
                elif visit.VISCODE == 'm03':  # Only for ADNIGO Month 3
                    preferred_visit_name = 'ADNIGO Month 3 MRI'
                else:
                    month = int(visit.VISCODE[1:])
                    if month < 54:
                        preferred_visit_name = 'ADNI1/GO Month ' + str(month)
                    else:
                        preferred_visit_name = 'ADNIGO Month ' + str(month)

            if preferred_visit_name in unique_visits:
                key_preferred_visit = (visit.VISCODE, visit.COLPROT, visit.ORIGPROT)
                if key_preferred_visit not in visits.keys():
                    visits[key_preferred_visit] = preferred_visit_name
                elif visits[key_preferred_visit] != preferred_visit_name:
                    print 'Multiple visits for one timepoint!'
                    print subject
                    print key_preferred_visit
                    print visits[key_preferred_visit]
                    print visit
                unique_visits.remove(preferred_visit_name)
                continue

            pending_timepoints.append(visit)

        # Then for images.Visit non matching the expected labels we find the closest date in visits list
        for visit in unique_visits:
            image = (ida_meta_subj[ida_meta_subj.Visit == visit]).iloc[0]
            min_db = 100000
            min_db2 = 0
            min_visit = None
            min_visit2 = None

            for timepoint in pending_timepoints:
                db = self.days_between(image['Scan Date'], timepoint.EXAMDATE)
                if db < min_db:
                    min_db2 = min_db
                    min_visit2 = min_visit

                    min_db = db
                    min_visit = timepoint

            if min_visit is None:
                print 'No corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.Visit
                print image
                continue

            if min_visit2 is not None and min_db > 90:
                print 'More than 60 days for corresponding timepoint in ADNIMERGE for subject ' + subject + ' in visit ' + image.Visit + ' on ' + image.ScanDate
                print 'Timepoint 1: ' + min_visit.VISCODE + ' - ' + min_visit.ORIGPROT + ' on ' + min_visit.EXAMDATE + ' (Distance: ' + str(
                    min_db) + ' days)'
                print 'Timepoint 2: ' + min_visit2.VISCODE + ' - ' + min_visit2.ORIGPROT + ' on ' + min_visit2.EXAMDATE + ' (Distance: ' + str(
                    min_db2) + ' days)'

                # If image is too close to the date between two visits we prefer the earlier visit
                if (datetime.strptime(min_visit.EXAMDATE, "%Y-%m-%d")
                        > datetime.strptime(image.ScanDate, "%Y-%m-%d")
                        > datetime.strptime(min_visit2.EXAMDATE, "%Y-%m-%d")):
                    dif = self.days_between(min_visit.EXAMDATE, min_visit2.EXAMDATE)
                    if abs((dif / 2.0) - min_db) < 30:
                        min_visit = min_visit2

                print 'We prefer ' + min_visit.VISCODE

            key_min_visit = (min_visit.VISCODE, min_visit.COLPROT, min_visit.ORIGPROT)
            if key_min_visit not in visits.keys():
                visits[key_min_visit] = image.Visit
            elif visits[key_min_visit] != image.Visit:
                print 'Multiple visits for one timepoint!'
                print subject
                print key_min_visit
                print visits[key_min_visit]
                print image.Visit

        return visits

    def compute_t1_paths(self, source_dir, clinical_dir, dest_dir, subjs_list):
        """

        :param source_dir:
        :param clinical_dir:
        :param dest_dir:
        :param subjs_list:
        :return:
        """

        import pandas as pd
        from os import path, walk

        t1_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                     'Study_ID', 'Field_Strength', 'Series_ID', 'Original']

        t1_df = pd.DataFrame(columns=t1_col_df)
        adni_merge_path = path.join(clinical_dir, 'ADNIMERGE.csv')
        adni_screening_path = path.join(clinical_dir, 'ADNI_ScreeningList_8_22_12.csv')
        ida_meta_path = path.join(clinical_dir, 'IDA_MR_METADATA_Listing.csv')
        mprage_meta_path = path.join(clinical_dir, 'MPRAGEMETA.csv')

        adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
        ida_meta = pd.io.parsers.read_csv(ida_meta_path, sep=',')
        mprage_meta = pd.io.parsers.read_csv(mprage_meta_path, sep=',')

        for subj in subjs_list:
            adnimerge_subj = adni_merge[adni_merge.PTID == subj]
            # Sort the values by examination date
            adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')

            mprage_meta_subj = mprage_meta[mprage_meta.SubjectID == subj]
            mprage_meta_subj = mprage_meta_subj.sort_values('ScanDate')

            ida_meta_subj = ida_meta[ida_meta.Subject == subj]

            mprage_meta_subj_orig = mprage_meta_subj[mprage_meta_subj['Orig/Proc'] == 'Original']
            visits = self.visits_to_timepoints_t1(subj, mprage_meta_subj, adnimerge_subj)

            keys = visits.keys()
            keys.sort()
            for visit_info in visits.keys():
                if visit_info[1] == 'ADNI1':
                    image_dict = self.adni1_image(subj, visit_info[0], visits[visit_info], mprage_meta_subj,
                                                  ida_meta_subj)
                elif visit_info[1] == 'ADNIGO':
                    image_dict = self.adnigo_image(subj, visit_info[0], visits[visit_info], mprage_meta_subj,
                                                   ida_meta_subj,
                                                   visit_info[2])
                else:  # ADNI2
                    image_dict = self.adni2_image(subj, visit_info[0], visits[visit_info], mprage_meta_subj_orig)

                if image_dict is None:
                    image_dict = {'Subject_ID': subj,
                                  'VISCODE': visit_info[0],
                                  'Visit': visits[visit_info],
                                  'Sequence': '',
                                  'Scan_Date': '',
                                  'Study_ID': '',
                                  'Series_ID': '',
                                  'Field_Strength': '',
                                  'Original': True}

                row_to_append = pd.DataFrame(image_dict, index=['i', ])
                t1_df = t1_df.append(row_to_append, ignore_index=True)

        images = t1_df
        is_dicom = []
        nifti_paths = []
        count = 0

        for row in images.iterrows():

            image = row[1]
            seq_path = path.join(source_dir, str(image.Subject_ID), image.Sequence)

            count += 1
            # print 'Processing Subject ' + str(image.Subject_ID) + ' - Session ' + image.VISCODE + ', ' + str(
            #     count) + ' / ' + str(total)

            series_path = ''
            for (dirpath, dirnames, filenames) in walk(seq_path):
                found = False
                for d in dirnames:
                    if d == 'S' + str(image.Series_ID):
                        series_path = path.join(dirpath, d)
                        found = True
                        break
                if found:
                    break

            nifti_path = series_path
            dicom = True

            for (dirpath, dirnames, filenames) in walk(series_path):
                for f in filenames:
                    if f.endswith(".nii"):
                        dicom = False
                        nifti_path = path.join(dirpath, f)
                        break

            is_dicom.append(dicom)
            nifti_paths.append(nifti_path)

        images.loc[:, 'Is_Dicom'] = pd.Series(is_dicom, index=images.index)
        images.loc[:, 'Path'] = pd.Series(nifti_paths, index=images.index)

        # Drop all the lines that have the Path section empty
        images = images.drop(images[images.Path == ''].index)

        # Store the paths inside a file called conversion_info inside the input directory
        images.to_csv(path.join(dest_dir, 'conversion_info', 't1_paths.tsv'), sep='\t', index=False)

        return images

    def compute_fdg_pet_paths(self, source_dir, clinical_dir, dest_dir, subjs_list):
        import pandas as pd
        import os
        from os import walk, path
        from numpy import nan, argsort

        pet_fdg_col = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID',
                       'Series_ID', 'Image_ID', 'Original']

        pet_fdg_df = pd.DataFrame(columns=pet_fdg_col)
        petqc_path = path.join(clinical_dir, 'PETQC.csv')
        pet_meta_list_path = path.join(clinical_dir, 'PET_META_LIST.csv')
        petqc = pd.io.parsers.read_csv(petqc_path, sep=',')
        pet_meta_list = pd.io.parsers.read_csv(pet_meta_list_path, sep=',')

        for subj in subjs_list:
            pet_qc_subj = petqc[(petqc.PASS == 1) & (petqc.RID == int(subj[-4:]))]
            subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subj]
            if subject_pet_meta.shape[0] < 1:
                # print 'NO Screening: Subject - ' + subject + ' for visit ' + qc_visit.VISCODE2
                continue
            for visit in list(pet_qc_subj.VISCODE2.unique()):
                pet_qc_visit = pet_qc_subj[pet_qc_subj.VISCODE2 == visit]
                if pet_qc_visit.shape[0] > 1:
                    normal_images = []
                    normal_meta = []
                    for row in pet_qc_visit.iterrows():
                        image = row[1]
                        pet_meta_image = \
                            subject_pet_meta[(subject_pet_meta['Image ID'] == int(image.LONIUID[1:]))].iloc[0]
                        if pet_meta_image.Sequence.lower().find('early') < 0:
                            normal_images.append(image)
                            normal_meta.append(pet_meta_image)
                    if len(normal_images) == 0:
                        print 'No regular FDG-PET image: Subject - ' + subj + ' for visit ' + visit
                        continue
                    if len(normal_images) == 1:
                        qc_visit = normal_images[0]
                    else:
                        qc_visit = None
                        index = argsort([x['Series ID'] for x in normal_meta])
                        for i in index[::-1]:
                            coreg_avg = subject_pet_meta[(subject_pet_meta['Sequence'] == 'Co-registered, Averaged')
                                                         & (
                                                             subject_pet_meta['Series ID'] == normal_meta[i][
                                                                 'Series ID'])]
                            if coreg_avg.shape[0] > 0:
                                qc_visit = normal_images[i]
                                break
                        if qc_visit is None:
                            qc_visit = normal_images[index[len(index) - 1]]
                else:
                    qc_visit = pet_qc_visit.iloc[0]
                int_image_id = int(qc_visit.LONIUID[1:])
                original_pet_meta = subject_pet_meta[
                    (subject_pet_meta['Orig/Proc'] == 'Original') & (subject_pet_meta['Image ID'] == int_image_id)]
                if original_pet_meta.shape[0] < 1:
                    original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                         & (subject_pet_meta.Sequence.map(
                        lambda x: (x.lower().find('fdg') > -1)))
                                                         & (subject_pet_meta['Scan Date'] == qc_visit.EXAMDATE)]
                    if original_pet_meta.shape[0] < 1:
                        print 'NO Screening: Subject - ' + subj + ' for visit ' + qc_visit.VISCODE2
                        continue
                original_image = original_pet_meta.iloc[0]
                averaged_pet_meta = subject_pet_meta[(subject_pet_meta['Sequence'] == 'Co-registered, Averaged') & (
                    subject_pet_meta['Series ID'] == original_image['Series ID'])]
                if averaged_pet_meta.shape[0] < 1:
                    sel_image = original_image
                    original = True
                else:
                    sel_image = averaged_pet_meta.iloc[0]
                    original = False
                visit = sel_image.Visit
                sequence = sel_image.Sequence
                sequence = sequence.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_').replace('(',
                                                                                                                    '_').replace(
                    ')', '_').replace(':', '_')
                date = sel_image['Scan Date']
                study_id = sel_image['Study ID']
                series_id = sel_image['Series ID']
                image_id = sel_image['Image ID']

                row_to_append = pd.DataFrame(
                    [[subj, qc_visit.VISCODE2, str(visit), sequence, date, str(study_id), str(series_id), str(image_id),
                      original]],
                    columns=pet_fdg_col)
                pet_fdg_df = pet_fdg_df.append(row_to_append, ignore_index=True)

        images = pet_fdg_df
        count = 0
        total = images.shape[0]
        image_folders = []
        for row in images.iterrows():
            image = row[1]
            seq_path = path.join(source_dir, str(image.Subject_ID), image.Sequence)
            count += 1
            print 'Processing Subject ' + str(image.Subject_ID) + ' - Session ' + image.VISCODE + ', ' + str(
                count) + ' / ' + str(total)
            image_path = ''
            for (dirpath, dirnames, filenames) in walk(seq_path):
                found = False
                for d in dirnames:
                    if d == 'I' + str(image.Image_ID):
                        image_path = path.join(dirpath, d)
                        found = True
                        break
                if found:
                    break
            if image.Original:
                for (dirpath, dirnames, filenames) in walk(image_path):
                    for f in filenames:
                        if f.endswith(".nii"):
                            image_path = path.join(dirpath, f)
                            break
            image_folders.append(image_path)
            if image_path == '':
                print 'Not found ' + str(image.Subject_ID)
        images.loc[:, 'Path'] = pd.Series(image_folders, index=images.index)

        fdg_csv_path = path.join(dest_dir, 'conversion_info')
        if not os.path.exists(fdg_csv_path):
            os.mkdir(fdg_csv_path)
        images.to_csv(path.join(fdg_csv_path, 'fdg_pet_paths.tsv'), sep='\t', index=False)

        return images

    def compute_av45_pet_paths(self, source_dir, clinical_dir, dest_dir, subjs_list):
        import pandas as pd
        from numpy import nan
        from os import walk, path
        import os

        pet_av45_col = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                        'Study_ID', 'Series_ID', 'Image_ID', 'Original']

        pet_av45_df = pd.DataFrame(columns=pet_av45_col)

        pet_visits_corr = {'ADNI Baseline': 'bl', 'ADNI2 Baseline-New Pt': 'bl'}
        av45qc_path = path.join(clinical_dir, 'AV45QC.csv')
        av45qc = pd.io.parsers.read_csv(av45qc_path, sep=',')
        baseline = av45qc[(av45qc.VISCODE2 == 'bl') & (av45qc.PASS == 1) & av45qc.RID.isin(
            [int(s[-4:]) for s in subjs_list])]

        pet_meta_list_path = path.join(clinical_dir, 'PET_META_LIST.csv')
        pet_meta_list = pd.io.parsers.read_csv(pet_meta_list_path, sep=',')
        pet_meta_list = pet_meta_list[pet_meta_list.Visit.map(lambda x: x.find('Baseline') > -1)]

        for subject in subjs_list:
            pet_qc_subj = av45qc[(av45qc.PASS == 1) & (av45qc.RID == int(subject[-4:]))]
            subject_pet_meta = pet_meta_list[pet_meta_list['Subject'] == subject]

            if subject_pet_meta.shape[0] < 1:
                # print 'NO Screening: Subject - ' + subject + ' for visit ' + qc_visit.VISCODE2
                continue

            for visit in list(pet_qc_subj.VISCODE2.unique()):
                pet_qc_visit = pet_qc_subj[pet_qc_subj.VISCODE2 == visit]
                if pet_qc_visit.shape[0] > 1:
                    normal_images = []
                    normal_meta = []
                    for row in pet_qc_visit.iterrows():
                        image = row[1]
                        pet_meta_image = \
                            subject_pet_meta[(subject_pet_meta['Image ID'] == int(image.LONIUID[1:]))].iloc[0]
                        if pet_meta_image.Sequence.lower().find('early') < 0:
                            normal_images.append(image)
                            normal_meta.append(pet_meta_image)

                    if len(normal_images) == 0:
                        print 'No regular AV45-PET image: Subject - ' + subject + ' for visit ' + visit
                        continue
                    if len(normal_images) == 1:
                        qc_visit = normal_images[0]
                    else:
                        qc_visit = None
                        index = nan.argsort([x['Series ID'] for x in normal_meta])
                        for i in index[::-1]:
                            coreg_avg = subject_pet_meta[
                                (subject_pet_meta['Sequence'] == 'AV45 Co-registered, Averaged')
                                & (subject_pet_meta['Series ID'] == normal_meta[i]['Series ID'])]
                            if coreg_avg.shape[0] > 0:
                                qc_visit = normal_images[i]
                                break
                        if qc_visit is None:
                            qc_visit = normal_images[index[len(index) - 1]]
                else:
                    qc_visit = pet_qc_visit.iloc[0]
                int_image_id = int(qc_visit.LONIUID[1:])
                original_pet_meta = subject_pet_meta[
                    (subject_pet_meta['Orig/Proc'] == 'Original') & (subject_pet_meta['Image ID'] == int_image_id)]
                if original_pet_meta.shape[0] < 1:
                    original_pet_meta = subject_pet_meta[(subject_pet_meta['Orig/Proc'] == 'Original')
                                                         & (subject_pet_meta.Sequence.map(
                        lambda x: (x.lower().find('av45') > -1)))
                                                         & (subject_pet_meta['Scan Date'] == qc_visit.EXAMDATE)]
                    if original_pet_meta.shape[0] < 1:
                        print 'NO Screening: Subject - ' + subject + ' for visit ' + qc_visit.VISCODE2
                        continue
                original_image = original_pet_meta.iloc[0]
                averaged_pet_meta = subject_pet_meta[
                    (subject_pet_meta['Sequence'] == 'AV45 Co-registered, Averaged') & (
                        subject_pet_meta['Series ID'] == original_image['Series ID'])]
                if averaged_pet_meta.shape[0] < 1:
                    sel_image = original_image
                    original = True
                else:
                    sel_image = averaged_pet_meta.iloc[0]
                    original = False
                visit = sel_image.Visit
                sequence = sel_image.Sequence
                sequence = sequence.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_').replace('(',
                                                                                                                    '_').replace(
                    ')', '_').replace(':', '_')
                viscode = qc_visit['VISCODE2']
                date = sel_image['Scan Date']
                study_id = sel_image['Study ID']
                series_id = sel_image['Series ID']
                image_id = sel_image['Image ID']

                row_to_append = pd.DataFrame(
                    [[subject, viscode, str(visit), sequence, date, str(study_id), str(series_id), str(image_id),
                      original]],
                    columns=pet_av45_col)
                pet_av45_df = pet_av45_df.append(row_to_append, ignore_index=True)

        subjects = pet_av45_df
        count = 0
        total = subjects.shape[0]

        image_folders = []
        for row in subjects.iterrows():

            subject = row[1]
            seq_path = path.join(source_dir, str(subject.Subject_ID), subject.Sequence)

            count += 1
            print 'Processing subject ' + str(subject.Subject_ID) + ', ' + str(count) + ' / ' + str(total)

            image_path = ''
            for (dirpath, dirnames, filenames) in walk(seq_path):
                found = False
                for d in dirnames:
                    if d == 'I' + str(subject.Image_ID):
                        image_path = path.join(dirpath, d)
                        found = True
                        break
                if found:
                    break

            if subject.Original:
                for (dirpath, dirnames, filenames) in walk(image_path):
                    for f in filenames:
                        if f.endswith(".nii"):
                            image_path = path.join(dirpath, f)
                            break

            image_folders.append(image_path)

            if image_path == '':
                print 'Not found ' + str(subject.Subject_ID)

        subjects.loc[:, 'Path'] = pd.Series(image_folders, index=subjects.index)

        av45_csv_path = path.join(dest_dir, 'conversion_info')
        if not os.path.exists(av45_csv_path):
            os.mkdir(av45_csv_path)
        subjects.to_csv(path.join(av45_csv_path, 'av45_pet_paths.tsv'), sep='\t', index=False)

        return subjects

    def compute_fmri_path(self, source_dir, clinical_dir, dest_dir, subjs_list):
        '''
        Compute the path for fmri images following the criteria described in the document Bids for AramisLab.

        :param source_dir: path to the ADNI image folder
        :param clinical_dir: path to the directory with all the clinical data od ADNI
        :param dest_dir: path to the output_folder
        :param subjs_list: subjects list
        :return: pandas Dataframe containing the path for each fmri
        '''
        import os
        from os import path
        from os import walk
        import pandas as pd

        fmri_col = ['Subject_ID', 'VISCODE', 'Visit', 'IMAGEUID', 'Sequence', 'Scan Date', 'LONIUID', 'Scanner',
                    'MagStregth', 'Path']

        fmri_df = pd.DataFrame(columns=fmri_col)

        # Load the requested clinical data
        mayo_mri_fmri_path = path.join(clinical_dir, 'MAYOADIRL_MRI_FMRI_09_15_16.csv')
        mayo_mri_imageqc_path = path.join(clinical_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')
        ida_mr_metadata_path = path.join(clinical_dir, 'IDA_MR_Metadata_Listing.csv')

        mayo_mri_fmri = pd.io.parsers.read_csv(mayo_mri_fmri_path, sep=',')
        ida_mr_metadata = pd.io.parsers.read_csv(ida_mr_metadata_path, sep=',')
        mayo_mri_imageqc = pd.io.parsers.read_csv(mayo_mri_imageqc_path, sep=',')

        for subj in subjs_list:
            print subj
            fmri_subjs_info = mayo_mri_fmri[(mayo_mri_fmri.RID == int(subj[-4:]))]
            # Extract visits available
            visits_list = fmri_subjs_info['VISCODE2'].tolist()
            # Removing duplicates
            visits_list = list(set(visits_list))

            if len(visits_list) != 0:
                for viscode in visits_list:
                    scan_date = ''
                    sequence = ''
                    loni_uid = ''
                    visit = ''
                    mag_strenght = ''
                    image_path = ''

                    fmri_subj = fmri_subjs_info[fmri_subjs_info['VISCODE2'] == viscode]

                    if not fmri_subj.empty:
                        # If there are multiple scans for the same session, check what is the one selected for the usage (field 'series_selected') or
                        # choose the one with the best quality
                        if len(fmri_subj) > 1:
                            fmri_imageuid = fmri_subj['IMAGEUID'].tolist()
                            loni_uid = ['I' + str(imageuid) for imageuid in fmri_imageuid]
                            images_qc = mayo_mri_imageqc[mayo_mri_imageqc.loni_image.isin(loni_uid)]
                            series_selected_values = images_qc['series_selected'].tolist()
                            sum_series_selected = sum(series_selected_values)
                            if sum_series_selected == 1:
                                imageuid_to_select = images_qc[images_qc['series_selected'] > 0]['loni_image'].iloc[
                                    0].replace('I', '')
                            else:
                                imageuid_to_select = self.select_image_qc(fmri_imageuid, images_qc)

                            fmri_subj = fmri_subj[fmri_subj['IMAGEUID'] == int(imageuid_to_select)].iloc[0]
                        else:
                            fmri_subj = fmri_subj.iloc[0]

                        fmri_imageuid = fmri_subj['IMAGEUID']

                        # Discard scans made with non Philips scanner and with a bad quality
                        fmri_metadata = ida_mr_metadata[ida_mr_metadata['IMAGEUID'] == fmri_imageuid]
                        if not fmri_metadata.empty:
                            fmri_metadata = fmri_metadata.iloc[0]

                            if not 'Philips' in fmri_metadata['Scanner']:
                                print 'No Philips scanner for ', subj, 'visit', viscode, '. Skipped.'
                                continue

                            elif 4 in mayo_mri_imageqc[mayo_mri_imageqc['loni_image'] == 'I' + str(fmri_imageuid)][
                                'series_quality'].values:
                                print 'Bad scan quality for ', subj, 'visit', viscode, '. Skipped.'
                                continue

                            scan_date = fmri_subj.SCANDATE
                            sequence = self.replace_sequence_chars(fmri_subj.SERDESC)
                            scanner = fmri_metadata['Scanner']
                            loni_uid = fmri_metadata['LONIUID']
                            visit = fmri_metadata['Visit']
                            mag_strenght = fmri_metadata['MagStrength']

                            # Calculate the path
                            seq_path = path.join(source_dir, str(subj), sequence)
                            for (dirpath, dirnames, filenames) in walk(seq_path):
                                found = False
                                for d in dirnames:
                                    if d == 'S' + str(loni_uid):
                                        image_path = path.join(dirpath, d)
                                        # Check if the path exists
                                        if not os.path.isdir(image_path):
                                            print 'Path not existing for subject:', subj, 'visit:', visit
                                        found = True
                                        break
                                if found:
                                    break

                            # The session scmri correspond to the baseline
                            if viscode == 'scmri':
                                viscode = 'bl'
                        else:
                            logging.info('Missing visit, sequence, scan date and loniuid for ', subj, 'visit', visit)
                        row_to_append = pd.DataFrame(
                            [[subj, str(viscode), visit, str(fmri_imageuid), sequence, scan_date, str(loni_uid),
                              scanner, mag_strenght, image_path]], columns=fmri_col)

                        fmri_df = fmri_df.append(row_to_append, ignore_index=True)
                    else:
                        logging.info('Missing fMRI for ', subj, 'visit', visit)

        fmri_df.to_csv(path.join(dest_dir, 'conversion_info', 'fmri_paths.tsv'), sep='\t', index=False)
        return fmri_df

    def compute_dti_paths(self, csv_dir, adni_dir, subjs_list):

        import pandas as pd
        from os import path, walk

        dti_col_df = ['Subject_ID', 'VISCODE', 'Visit', 'Sequence', 'Scan_Date',
                      'Study_ID', 'Series_ID', 'Image_ID', 'Field_Strength', 'Scanner', 'Enhanced']

        dti_df = pd.DataFrame(columns=dti_col_df)
        # subjs_list_path = path.join(source_dir, 'clinicalData', 'subjects_list.xlsx')

        adni_merge_path = path.join(csv_dir, 'ADNIMERGE.csv')
        ida_meta_path = path.join(csv_dir, 'IDA_MR_METADATA_Listing.csv')
        mri_qc_path = path.join(csv_dir, 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')

        # adni_merge_path = path.join(source_dir, 'clinicalData', 'ADNIMERGE.csv')
        # ida_meta_path = path.join(source_dir, 'clinicalData', 'IDA_MR_METADATA_Listing.csv')
        # mri_qc_path = path.join(source_dir, 'clinicalData', 'MAYOADIRL_MRI_IMAGEQC_12_08_15.csv')

        # subjs_list_excel = pd.read_excel(subjs_list_path)
        # subjs_list = subjs_list_excel['PTID']

        adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
        ida_meta = pd.io.parsers.read_csv(ida_meta_path, sep=',')
        ida_meta = ida_meta[ida_meta.Sequence.map(lambda x: x.lower().find('dti') > -1)]

        mri_qc = pd.io.parsers.read_csv(mri_qc_path, sep=',')
        mri_qc = mri_qc[mri_qc.series_type == 'DTI']

        for subj in subjs_list:

            adnimerge_subj = adni_merge[adni_merge.PTID == subj]
            # Sort the values by examination date
            adnimerge_subj = adnimerge_subj.sort_values('EXAMDATE')

            ida_meta_subj = ida_meta[ida_meta.Subject == subj]
            ida_meta_subj = ida_meta_subj.sort_values('Scan Date')

            # print subj
            # print type(subj)

            mri_qc_subj = mri_qc[mri_qc.RID == int(subj[-4:])]

            visits = self.visits_to_timepoints_dti(subj, ida_meta_subj, adnimerge_subj)

            keys = visits.keys()
            keys.sort()

            for visit_info in visits.keys():

                # visit_info = (VISCODE, COLPROT, ORIGPROT)

                visit_str = visits[visit_info]
                visit_ida_meta = ida_meta_subj[ida_meta_subj.Visit == visit_str]

                axial_ida_meta = visit_ida_meta[visit_ida_meta.Sequence.map(lambda x: x.lower().find('enhanced') < 0)]
                enhanced_ida_meta = visit_ida_meta[visit_ida_meta.Sequence.map(lambda x: x.lower().find('enhanced') > -1)]

                axial = self.dti_image(subj, visit_info[0], visits[visit_info], axial_ida_meta, mri_qc_subj, False)
                enhanced = self.dti_image(subj, visit_info[0], visits[visit_info], enhanced_ida_meta, mri_qc_subj, True)

                if axial is not None:
                    row_to_append = pd.DataFrame(axial, index=['i', ])
                    dti_df = dti_df.append(row_to_append, ignore_index=True)

                if enhanced is not None:
                    row_to_append = pd.DataFrame(enhanced, index=['i', ])
                    dti_df = dti_df.append(row_to_append, ignore_index=True)

        images = dti_df
        is_dicom = []
        nifti_paths = []
        count = 0

        for row in images.iterrows():

            image = row[1]
            seq_path = path.join(adni_dir, str(image.Subject_ID), image.Sequence)

            count += 1
            # print 'Processing Subject ' + str(image.Subject_ID) + ' - Session ' + image.VISCODE + ', ' + str(
            #     count) + ' / ' + str(total)

            series_path = ''
            s = 'S' + str(image.Series_ID)

            for (dirpath, dirnames, filenames) in walk(seq_path):
                found = False
                for d in dirnames:
                    if d == s:
                        series_path = path.join(dirpath, d)
                        found = True
                        break
                if found:
                    break

            nifti_path = series_path
            dicom = True

            # for (dirpath, dirnames, filenames) in walk(series_path):
            #     for f in filenames:
            #         if f.endswith(".nii"):
            #             dicom = False
            #             nifti_path = path.join(dirpath, f)
            #             break

            is_dicom.append(dicom)
            nifti_paths.append(nifti_path)

        images.loc[:, 'Is_Dicom'] = pd.Series(is_dicom, index=images.index)
        images.loc[:, 'Path'] = pd.Series(nifti_paths, index=images.index)

        # Drop all the lines that have the Path section empty
        # images = images.drop(images[images.Path == ''].index)
        # Store the paths inside a file called t1_paths inside the input directory
        images.to_csv(path.join('/Users/jorge.samper/Workspace', 'dti_paths.tsv'), sep='\t', index=False)

        return images

    def run_pipeline(self, args):
        if args.modality is True:
            self.convert_clinical_data(args.dataset_directory, args.bids_directory)



    #TODO To integrate
    def viscode_to_session(viscode):
        if viscode == 'bl':
            return 'M00'
        else:
            return viscode.capitalize()

    def dti_to_bids(self, dti_paths, bids_dir):

        #TODO Remove imports from here?
        from clinica.bids.bids_utils import remove_space_and_symbols, dcm_to_nii
        import pandas as pd
        from os import path, makedirs
        from numpy import nan

        dti_images = pd.io.parsers.read_csv(dti_paths, sep='\t')

        for row in dti_images.iterrows():
            image = row[1]
            if image.Path is nan:
                continue

            subject = 'sub-' + remove_space_and_symbols(image.Subject_ID)
            session = 'ses-' + viscode_to_session(image.VISCODE)

            output_path = path.join(bids_dir, subject, session, 'dwi')
            #TODO Define the standard notation for acquisition name: axial and axialEnhanced? Or enhancedAxial?
            bids_name = subject + '_' + session + '_acq-' + ('axialEnhanced' if image.Enhanced else 'axial') + '_dwi'

            try:
                makedirs(output_path)
            except OSError:
                if not path.isdir(output_path):
                    raise
            dcm_to_nii(image.Path, output_path, bids_name)







