from clinica.bids.abstract_converter import Converter
from clinica.engine.cmdparser import CmdParser

__author__ = "Sabrina Fontanella"
__copyright__ = "Copyright 2017, The Aramis Lab Team"
__credits__ = ["Jorge Samper"]
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
        bids_subjs_paths = bids.get_bids_subjs_pats(out_path)

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
                sessions_df = pd.DataFrame(sessions_dict[bids_id]['bl'], index=['i', ])
                for ses in sess_aval:
                    sessions_df = sessions_df.append(pd.DataFrame(sessions_dict[bids_id][ses], index=['i', ]))

                sessions_df.to_csv(path.join(sp, bids_id + '_sessions.tsv'), sep='\t', index=False)


    def convert_t1_from_dicom(self, t1_path, output_path, bids_name):
        """

        :param t1_path:
        :param output_path:
        :param bids_name:
        :return:
        """

        import os
        import clinica.bids.bids_utils as bids

        if not os.path.exists(output_path):
            os.mkdir(output_path)
        os.system('dcm2niix -b n -z y -o ' + output_path + ' -f ' + bids_name + bids.get_bids_suff('T1') + ' ' + t1_path )


    def convert_images(self, source_dir, dest_dir):
        """
        :param source_dir:
        :param dest_dir:
        :return:
        """

        import clinica.bids.bids_utils as bids
        import os
        from os import path
        import pandas as pd

        t1_paths = self.compute_t1_paths(source_dir)
        subjs_list_path = path.join(source_dir, 'clinicalData', 'subjects_list_from_adnimerge.xlsx')
        subjs_list_excel = pd.read_excel(subjs_list_path)
        subjs_list = subjs_list_excel['PTID']
        adni_merge_path = path.join(source_dir, 'clinicalData', 'ADNIMERGE.csv')
        sess_list = ['bl']
        bids_ids = []
        alpha_ids = []

        os.mkdir(dest_dir)

        for subj in subjs_list:
            alpha_id = bids.remove_space_and_symbols(subj)
            bids_id = 'sub-ADNI' + alpha_id
            alpha_ids.append(alpha_id)
            bids_ids.append(bids_id)
            os.mkdir(path.join(dest_dir, bids_id))

        for i in range(0, len(subjs_list)):
            for ses in sess_list:
                bids_file_name = bids_ids[i] + '_ses-' + ses
                bids_ses_id = 'ses-' + ses
                ses_path = path.join(dest_dir, bids_ids[i], bids_ses_id)
                os.mkdir(ses_path)

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
                        self.convert_t1_from_dicom(t1_path, path.join(ses_path, 'anat'), bids_file_name)
                    else:
                        print 'No path for dicom'

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
        from os import path

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
        print subjects

    def run_pipeline(self, args):
        if args.modality is True:
            self.convert_clinical_data(args.dataset_directory, args.bids_directory)








