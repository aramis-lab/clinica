from os import path, walk, remove, makedirs
from glob import glob
import os
import logging
from clinica.bids.bids_utils import  remove_space_and_symbols
import clinica.bids.bids_utils as bids
import pandas as pd
import csv
from clinica.bids.converter_utils import  remove_space_and_symbols
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

    #def __init__(self):
    #    super(ADNI_TO_BIDS, self).__init__()

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

    def convert_clinical_data(self, src, dest_dir):pass

    def convert_images(self, source_dir, dest_dir):
        t1_paths = compute_t1_paths(source_dir)
        # print t1_table
        subjs_list_path = path.join(source_dir, 'clinicalData', 'subjects_list_from_adnimerge.xlsx')
        subjs_list_excel = pd.read_excel(subjs_list_path)
        subjs_list = subjs_list_excel['PTID']
        adni_merge_path = path.join(source_dir, 'clinicalData', 'ADNIMERGE.csv')
        sess_list = ['bl']
        bids_ids = []
        alpha_ids = []

        for subj in subjs_list:
            alpha_id = remove_space_and_symbols(subj)
            bids_id = 'sub-' + alpha_id
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
                        bids.convert_T1(t1_path, path.join(ses_path, 'anat'), bids_file_name)
                else:
                    # Convert the image using dcm2nii
                    print 'Dicom found, needs to be converted'



    def compute_t1_paths(source_dir):
        """
        Select the T1 to use for each subject.

        Adapted from a script of Jorge Samper
        :param source_dir:
        :return: Dictionary that has the following structure
        { subj_id1 : {session1 : {modality: file_name, ...}, subj_id2 : ....}
        """
        t1_dict = {}
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

    def run_pipeline(self, args):
        if args.modality is True:
            self.convert_clinical_data(args.dataset_directory, args.bids_directory)

        #self.convert_images(args.xx)






