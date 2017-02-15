from clinica.bids.abstract_converter import Converter
from clinica.engine.cmdparser import CmdParser

__author__ = "Jorge Samper and Sabrina Fontanella"
__copyright__ = "Copyright 2017, The Aramis Lab Team"
__credits__ = [""]
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Sabrina Fontanella"
__email__ = "sabrina.fontanella@icm-institute.org"
__status__ = "Development"


class ADNI_TO_BIDS(Converter, CmdParser) :

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
        bids.create_participants(input_path, out_path, 'ADNI', 'clinical_specifications_adni.xlsx', bids_ids)

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

    def convert_images(self, source_dir, dest_dir, mod_to_add = '', mod_to_update = ''):
        """
        :param source_dir:
        :param dest_dir:
        :return:
        """

        import clinica.bids.bids_utils as bids
        import os
        from os import path
        import pandas as pd
        import shutil

        # Compute the path to the t1 images
        if (mod_to_add == '' or mod_to_add == 't1') and ( mod_to_update == '' or mod_to_update == 't1' ):
            t1_paths = self.compute_t1_paths(source_dir)


        subjs_list_path = path.join(source_dir, 'clinicalData', 'subjects_list_from_adnimerge.xlsx')
        subjs_list_excel = pd.read_excel(subjs_list_path)
        subjs_list = subjs_list_excel['PTID']
        adni_merge_path = path.join(source_dir, 'clinicalData', 'ADNIMERGE.csv')
        adni_merge = pd.read_csv(adni_merge_path)
        # Extract the list of sessions from VISCODE column of adnimerge file, remove duplicate and convert to a list
        session_list = adni_merge['VISCODE'].drop_duplicates().tolist()
        print session_list
        sess_list = ['bl']
        bids_ids = []
        alpha_ids = []

        pet_fdg_paths = self.compute_fdg_pet_paths(source_dir, subjs_list)

        if mod_to_add == '' and mod_to_update=='':
            os.mkdir(dest_dir)

        for subj in subjs_list:
            alpha_id = bids.remove_space_and_symbols(subj)
            bids_id = 'sub-ADNI' + alpha_id
            alpha_ids.append(alpha_id)
            bids_ids.append(bids_id)
            if mod_to_add == '' and mod_to_update=='':
                os.mkdir(path.join(dest_dir, bids_id))

        for i in range(0, len(subjs_list)):
            for ses in sess_list:
                bids_file_name = bids_ids[i] + '_ses-' + ses
                bids_ses_id = 'ses-' + ses
                ses_path = path.join(dest_dir, bids_ids[i], bids_ses_id)
                if mod_to_add == '' and mod_to_update=='':
                    os.mkdir(ses_path)

                if mod_to_add != '' and mod_to_add !='t1':
                    pass
                elif (mod_to_add=='' and mod_to_update=='') or mod_to_add == 't1' or mod_to_update == 't1':
                    # Convert T1
                    t1_info = t1_paths[(t1_paths['subj_id'] == subjs_list[i]) & (t1_paths['session'] == ses)]
                    if t1_info['dicom'].values[0] == False:
                        t1_path = t1_info['path'].values[0]
                        if len(t1_path) == 0:
                            print 'No path'
                        else:
                            os.mkdir(path.join(ses_path, 'anat'))
                            bids.convert_T1(t1_path, path.join(ses_path, 'anat'), bids_file_name)
                    else:
                        # Convert the image using dcm2nii
                        print 'Dicom found, needs to be converted'
                        os.mkdir(path.join(ses_path, 'anat'))
                        t1_path = (t1_info['path'].values[0]).replace(' ', '\ ')
                        if len(t1_path)!=0:
                            self.convert_from_dicom(t1_path, path.join(ses_path, 'anat'), bids_file_name, 'T1')
                        else:
                            print 'No path for dicom'

                if mod_to_add != '' and mod_to_add !='pet_fdg' and mod_to_update!='pet_fdg':
                    pass
                elif mod_to_add=='' or mod_to_add == 'pet_fdg' or mod_to_update == 'pet_fdg':
                    pet_fdg_info = pet_fdg_paths
                    visit_session = {'ADNI Baseline': 'bl', 'ADNI2 Baseline-New Pt': 'bl'}


                    images = pet_fdg_info
                    count = 0
                    total = images.shape[0]

                    # Check if the pet folder already exist
                    if os.path.isdir(path.join(ses_path, 'pet')):
                        if mod_to_add != '':
                            raise IOError('PET modality found. For updating the dataset use the flag -updated_mod')

                        print 'Removing the old PET folder...'
                        shutil.rmtree(path.join(ses_path, 'pet'))


                    pet_fdg_info = pet_fdg_paths[(pet_fdg_paths['Subject_ID'] == subjs_list[i])]
                    original = pet_fdg_info['Original'].values
                    if len(original) > 0:
                        if pet_fdg_info['Original'].values[0] == True:
                            print 'is original'
                        else:
                            pet_path = pet_fdg_info['Path'].values[0].replace(' ', '\ ')
                            self.convert_from_dicom(pet_path,  path.join(ses_path, 'pet'), bids_file_name, 'pet')
                    else:
                        print 'Original not found for ', subjs_list[i]

    def compute_t1_paths(self, source_dir):
        """
        Select the T1 to use for each subject.

        Adapted from a script of Jorge Samper
        :param source_dir:
        :return: Dictionary that has the following structure
        { subj_id1 : { session1 : { modality: file_name, ...},
                       session2 : {...}}
          subj_id2 : ....}
        """

        import pandas as pd
        from os import path, walk

        t1_col_df = ['subj_id', 'session','filename', 'date', 'series_id']
        sessions_list = ['bl']

        t1_df = pd.DataFrame(columns = t1_col_df)
        row_to_append = pd.DataFrame(columns = t1_col_df)
        subjs_list_path = path.join(source_dir, 'clinicalData', 'subjects_list_from_adnimerge.xlsx')
        adni_merge_path = path.join(source_dir, 'clinicalData', 'ADNIMERGE.csv')
        adni_screening_path = path.join(source_dir, 'clinicalData', 'ADNI_ScreeningList_8_22_12.csv')
        ida_meta_path = path.join(source_dir, 'clinicalData', 'IDA_MR_METADATA_Listing.csv')
        mprage_meta_path = path.join(source_dir, 'clinicalData', 'MPRAGEMETA.csv')

        subjs_list_excel = pd.read_excel(subjs_list_path)
        subjs_list = subjs_list_excel['PTID']

        adni_merge = pd.io.parsers.read_csv(adni_merge_path, sep=',')
        adni1_screening = pd.io.parsers.read_csv(adni_screening_path, sep=',')
        ida_meta = pd.io.parsers.read_csv(ida_meta_path, sep=',')
        mprage_meta = pd.io.parsers.read_csv(mprage_meta_path, sep=',')


        for ses in sessions_list:
            adni_subj = adni_merge[(adni_merge.VISCODE == ses) & adni_merge.PTID.isin(subjs_list)]

            for row in adni_subj.iterrows():
                subject = row[1]
                subj_id = subject.PTID
                original = True

                if subject.ORIGPROT == 'ADNI1':

                    # If the subject came from ADNI1, 1.5T image will be selected
                    filtered_screening = adni1_screening[adni1_screening.PTID == subj_id]
                    if filtered_screening.shape[0] < 1:
                        print 'NO Screening: ' + subj_id
                        continue

                    sel_image = filtered_screening.iloc[0]
                    series_id = sel_image['Series.ID']

                    filtered_scan = ida_meta[ida_meta.LONIUID == series_id]
                    if filtered_scan.shape[0] < 1:
                        print 'NO IDA Meta: ' + subj_id
                        continue
                    ida_scan = filtered_scan.iloc[0]

                    if ida_scan.Scanner.find('Philips') > -1:

                        scan = \
                            (mprage_meta[(mprage_meta['Orig/Proc'] == 'Original') & (mprage_meta.SeriesID == series_id)]).iloc[
                                0]
                        sequence = scan.Sequence

                    else:
                        scan = \
                            (mprage_meta[(mprage_meta['Orig/Proc'] == 'Processed') & (mprage_meta.SeriesID == series_id)]).iloc[
                                0]
                        sequence = scan.Sequence[:scan.Sequence.find('N3') - 2]
                        original = False
                else:
                    visits = ('ADNIGO Screening MRI', 'ADNI2 Screening MRI-New Pt')

                    cond_mprage = ((mprage_meta.SubjectID == subj_id) & mprage_meta.Visit.isin(visits) & (
                    mprage_meta['Orig/Proc'] == 'Original')
                                   & (mprage_meta.Sequence.map(lambda x: x.find('MPRAGE') > -1) & mprage_meta.Sequence.map(
                        lambda x: x.find('2') < 0)))

                    cond_spgr = ((mprage_meta.SubjectID == subj_id) & (mprage_meta.Visit.isin(visits)) & (
                    mprage_meta['Orig/Proc'] == 'Original')
                                 & (mprage_meta.Sequence.map(lambda x: x.find('SPGR') > -1) & mprage_meta.Sequence.map(
                        lambda x: x.lower().find('acc') < 0) & mprage_meta.Sequence.map(lambda x: x.find('REPEAT') < 0)))

                    filtered_scan = mprage_meta[cond_mprage | cond_spgr]
                    if filtered_scan.shape[0] < 1:
                        print 'NO MPRAGE Meta: ' + subj_id
                        continue

                    scan = filtered_scan.iloc[0]
                    sequence = scan.Sequence
                    series_id = scan.SeriesID

                date = scan.ScanDate
                sequence = sequence.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_')
                row_to_append = pd.DataFrame([[subj_id, ses, sequence, date, str(series_id)]], columns= t1_col_df)
                t1_df = t1_df.append(row_to_append, ignore_index=True)

            subjects = t1_df
            is_dicom = []
            nifti_paths = []

            for row in subjects.iterrows():
                    subject = row[1]

                    seq_path = path.join(source_dir, subject['subj_id'], subject['filename'])

                    series_path = ''
                    for (dirpath, dirnames, filenames) in walk(seq_path):
                        found = False
                        for d in dirnames:

                            if d == 'S' + str(subject['series_id']):
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

                    if nifti_path == '':
                        print 'Not found ' + str(subject['subj_id'])

            subjects.loc[:, 'dicom'] = pd.Series(is_dicom, index=subjects.index)
            subjects.loc[:, 'path'] = pd.Series(nifti_paths, index=subjects.index)

        return subjects

    def compute_fdg_pet_paths(self, source_dir, subjs_list):
        import pandas as pd
        import csv
        from os import walk, path

        pet_fdg_col = ['Subject_ID', 'Visit', 'Sequence', 'Scan_Date', 'Study_ID',
                                                             'Series_ID', 'Image_ID', 'Original']
        pet_fdg_df = pd.DataFrame(columns=pet_fdg_col)
        petqc_path = path.join(source_dir, 'clinicalData', 'PETQC.csv')
        pet_meta_list_path = path.join(source_dir,'clinicalData', 'PET_META_LIST.csv')
        out_images_filename = '/Users/sabrina.fontanella/clinica_io/phase1_pet.csv'
        output_paths_list = '/Users/sabrina.fontanella/clinica_io/phase2_pet.csv'

        try:
            petqc = pd.io.parsers.read_csv(petqc_path, sep=',')
        except IOError:
            print ('File PETQC.csv not found in folder clinicalData of ADNI.')

        baseline = petqc[
            (petqc.VISCODE2 == 'bl') & (petqc.PASS == 1) & petqc.RID.isin([int(s[-4:]) for s in subjs_list])]

        try:
            pet_meta_list = pd.io.parsers.read_csv(pet_meta_list_path, sep=',')
        except IOError:
            print ('File PET_META_LIST.csv not found in folder clinicalData of ADNI.')

        pet_meta_list = pet_meta_list[pet_meta_list.Visit.map(lambda x: x.find('Baseline') > -1)]


        for row in baseline.iterrows():
            subject = row[1]
            zeros_rid = self.fill_zeros(subject.RID, 4)
            int_image_id = int(subject.LONIUID[1:])
            filtered_pet_meta = pet_meta_list[pet_meta_list['Subject'].map(lambda x: x.endswith(zeros_rid))]

            if filtered_pet_meta.shape[0] < 1:
                print 'NO Screening: RID - ' + subject.RID
                continue

            original_pet_meta = filtered_pet_meta[
                (filtered_pet_meta['Orig/Proc'] == 'Original') & (filtered_pet_meta['Image ID'] == int_image_id)]
            if original_pet_meta.shape[0] < 1:
                print 'NO Screening: RID - ' + subject.RID
                continue

            original_image = original_pet_meta.iloc[0]

            averaged_pet_meta = filtered_pet_meta[(filtered_pet_meta['Sequence'] == 'Co-registered, Averaged') & (
            filtered_pet_meta['Series ID'] == original_image['Series ID'])]
            if averaged_pet_meta.shape[0] < 1:
                sel_image = original_image
                original = True
            else:
                sel_image = averaged_pet_meta.iloc[0]
                original = False

            subj_id = sel_image.Subject
            visit = sel_image.Visit
            sequence = sel_image.Sequence
            sequence = sequence.replace(' ', '_').replace('/', '_').replace(';', '_').replace('*', '_').replace('(',
                                                                                                                    '_').replace(
                ')', '_').replace(':', '_')
            date = sel_image['Scan Date']
            study_id = sel_image['Study ID']
            series_id = sel_image['Series ID']
            image_id = sel_image['Image ID']


            row_to_append = pd.DataFrame([[subj_id, str(visit), sequence, date, study_id, str(series_id), str(image_id), original]], columns=pet_fdg_col)
            #print row_to_append
            pet_fdg_df = pet_fdg_df.append(row_to_append, ignore_index=True)


        pet_fdg_df.to_csv(out_images_filename, '\t', index=False)
        subjects = pet_fdg_df
        count = 0
        total = subjects.shape[0]

        image_folders = []
        for row in subjects.iterrows():

            subject = row[1]
            seq_path = path.join(source_dir, str(subject.Subject_ID), subject.Sequence)
            print seq_path

            count += 1
            print 'Processing subject ' + str(subject.Subject_ID) + ', ' + str(count) + ' / ' + str(total)

            image_path = ''
            for (dirpath, dirnames, filenames) in walk(seq_path):
                found = False
                for d in dirnames:
                    # print 'd =', d
                    # print 'image ID', subject.Image_ID
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
                #subjects.to_csv(output_paths_list, sep='\t', index=False)

        subjects.to_csv(output_paths_list, sep='\t', index=False)
        return subjects

    def run_pipeline(self, args):
        if args.modality is True:
            self.convert_clinical_data(args.dataset_directory, args.bids_directory)








